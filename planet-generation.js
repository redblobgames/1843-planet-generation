/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 *
 * Adapting mapgen4 code for a sphere. Quick & dirty, for procjam2018
 */
'use strict';

const SEED = 123;

const SimplexNoise = require('simplex-noise');
const FlatQueue = require('flatqueue');
const colormap = require('./colormap');
const {vec3, mat4} = require('gl-matrix');
const {makeRandInt, makeRandFloat} = require('@redblobgames/prng');
const SphereMesh = require('./sphere-mesh');

const regl = require('regl')({
    canvas: "#output",
    extensions: ['OES_element_index_uint', 'OES_standard_derivatives']
});

const u_colormap = regl.texture({
    width: colormap.width,
    height: colormap.height,
    data: colormap.data,
    wrapS: 'clamp',
    wrapT: 'clamp'
});


/* UI parameters */
let N = 10000;
let P = 20;
let jitter = 0.75;
let rotation = -1;
let drawMode = 'centroid';
let draw_plateVectors = false;
let draw_plateBoundaries = false;

window.setN = newN => { N = newN; generateMesh(); };
window.setP = newP => { P = newP; generateMap(); };
window.setJitter = newJitter => { jitter = newJitter; generateMesh(); };
window.setRotation = newRotation => { rotation = newRotation; draw(); };
window.setDrawMode = newMode => { drawMode = newMode; draw(); };
window.setDrawPlateVectors = flag => { draw_plateVectors = flag; draw(); };
window.setDrawPlateBoundaries = flag => { draw_plateBoundaries = flag; draw(); };

const renderPoints = regl({
    frag: `
precision mediump float;
void main() {
   gl_FragColor = vec4(0, 0, 0, 1);
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
uniform float u_pointsize;
attribute vec3 a_xyz;
void main() {
  gl_Position = u_projection * vec4(a_xyz, 1);
  gl_PointSize = gl_Position.z > 0.0? 0.0 : u_pointsize;
}
`,

    depth: {
        enable: false,
    },
    
    uniforms: {
        u_projection: regl.prop('u_projection'),
        u_pointsize: regl.prop('u_pointsize'),
    },

    primitive: 'points',
    count: regl.prop('count'),
    attributes: {
        a_xyz: regl.prop('a_xyz'),
    },
});


const renderLines = regl({
    frag: `
precision mediump float;
uniform vec4 u_multiply_rgba, u_add_rgba;
varying vec4 v_rgba;
void main() {
   gl_FragColor = v_rgba * u_multiply_rgba + u_add_rgba;
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
attribute vec3 a_xyz;
attribute vec4 a_rgba;
varying vec4 v_rgba;
void main() {
  vec4 pos = u_projection * vec4(a_xyz, 1);
  v_rgba = (-2.0 * pos.z) * a_rgba;
  gl_Position = pos;
}
`,

    depth: {
        enable: false,
    },
    
    uniforms: {
        u_projection: regl.prop('u_projection'),
        u_multiply_rgba: regl.prop('u_multiply_rgba'),
        u_add_rgba: regl.prop('u_add_rgba'),
    },

    blend: {
        enable: true,
        func: {src: 'one', dst: 'one minus src alpha'},
        equation: {
            rgb: 'add',
            alpha: 'add'
        },
        color: [0, 0, 0, 0],
    },
    primitive: 'lines',
    count: regl.prop('count'),
    attributes: {
        a_xyz: regl.prop('a_xyz'),
        a_rgba: regl.prop('a_rgba'),
    },
});


const renderTriangles = regl({
    frag: `
precision mediump float;
uniform sampler2D u_colormap;
varying vec2 v_tm;
void main() {
   float e = v_tm.x > 0.0? 0.5 * (v_tm.x * v_tm.x + 1.0) : 0.5 * (v_tm.x + 1.0);
   gl_FragColor = texture2D(u_colormap, vec2(e, v_tm.y));
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
attribute vec3 a_xyz;
attribute vec2 a_tm;
varying vec2 v_tm;
void main() {
  v_tm = a_tm;
  gl_Position = u_projection * vec4(a_xyz, 1);
}
`,

    uniforms: {
        u_colormap: u_colormap,
        u_projection: regl.prop('u_projection'),
    },

    count: regl.prop('count'),
    attributes: {
        a_xyz: regl.prop('a_xyz'),
        a_tm: regl.prop('a_tm'),
    },
});


const renderIndexedTriangles = regl({
    frag: `
#extension GL_OES_standard_derivatives : enable

precision mediump float;

uniform sampler2D u_colormap;
uniform vec2 u_light_angle;
uniform float u_inverse_texture_size, u_slope, u_flat, u_c, u_d, u_outline_strength;

varying vec2 v_tm;
void main() {
   float e = v_tm.x > 0.0? 0.5 * (v_tm.x * v_tm.x + 1.0) : 0.5 * (v_tm.x + 1.0);
   float dedx = dFdx(v_tm.x);
   float dedy = dFdy(v_tm.x);
   vec3 slope_vector = normalize(vec3(dedy, dedx, u_d * 2.0 * u_inverse_texture_size));
   vec3 light_vector = normalize(vec3(u_light_angle, mix(u_slope, u_flat, slope_vector.z)));
   float light = u_c + max(0.0, dot(light_vector, slope_vector));
   float outline = 1.0 + u_outline_strength * max(dedx,dedy);
   gl_FragColor = vec4(texture2D(u_colormap, vec2(e, v_tm.y)).rgb * light / outline, 1);
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
attribute vec3 a_xyz;
attribute vec2 a_tm;
varying vec2 v_tm;
void main() {
  v_tm = a_tm;
  gl_Position = u_projection * vec4(a_xyz, 1);
}
`,

    uniforms: {
        u_colormap: u_colormap,
        u_projection: regl.prop('u_projection'),
        u_light_angle: [Math.cos(Math.PI/3), Math.sin(Math.PI/3)],
        u_inverse_texture_size: 1.0 / 2048,
        u_d: 60,
        u_c: 0.15,
        u_slope: 6,
        u_flat: 2.5,
        u_outline_strength: 5,
    },

    elements: regl.prop('elements'),
    attributes: {
        a_xyz: regl.prop('a_xyz'),
        a_tm: regl.prop('a_tm'),
    },
});

/**********************************************************************
 * Geometry
 */

let _randomNoise = new SimplexNoise(makeRandFloat(SEED));
const persistence = 2/3;
const amplitudes = Array.from({length: 5}, (_, octave) => Math.pow(persistence, octave));

function fbm_noise(nx, ny, nz) {
    let sum = 0, sumOfAmplitudes = 0;
    for (let octave = 0; octave < amplitudes.length; octave++) {
        let frequency = 1 << octave;
        sum += amplitudes[octave] * _randomNoise.noise3D(nx * frequency, ny * frequency, nz * frequency);
        sumOfAmplitudes += amplitudes[octave];
    }
    return sum / sumOfAmplitudes;
}

/* Calculate the centroid and push it onto an array */
function pushCentroidOfTriangle(out, ax, ay, az, bx, by, bz, cx, cy, cz) {
    // TODO: renormalize to radius 1
    out.push((ax+bx+cx)/3, (ay+by+cy)/3, (az+bz+cz)/3);
}


function generateTriangleCenters(mesh, {r_xyz}) {
    let {numTriangles} = mesh;
    let t_xyz = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = mesh.s_begin_r(3*t),
            b = mesh.s_begin_r(3*t+1),
            c = mesh.s_begin_r(3*t+2);
        pushCentroidOfTriangle(t_xyz,
                 r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
                 r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
                 r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]);
    }
    return t_xyz;
}

function generateVoronoiGeometry(mesh, {r_xyz, t_xyz}, r_color_fn) {
    const {numSides} = mesh;
    let xyz = [], tm = [];

    for (let s = 0; s < numSides; s++) {
        let inner_t = mesh.s_inner_t(s),
            outer_t = mesh.s_outer_t(s),
            begin_r = mesh.s_begin_r(s);
        let rgb = r_color_fn(begin_r);
        xyz.push(t_xyz[3*inner_t], t_xyz[3*inner_t+1], t_xyz[3*inner_t+2],
                      t_xyz[3*outer_t], t_xyz[3*outer_t+1], t_xyz[3*outer_t+2],
                      r_xyz[3*begin_r], r_xyz[3*begin_r+1], r_xyz[3*begin_r+2]);
        tm.push(rgb, rgb, rgb);
    }
    return {xyz, tm};
}

class QuadGeometry {
    constructor () {
        /* xyz = position in 3-space;
           tm = temperature, moisture
           I = indices for indexed drawing mode */
    }

    setMesh({numSides, numRegions, numTriangles}) {
        this.I = new Int32Array(3 * numSides);
        this.xyz = new Float32Array(3 * (numRegions + numTriangles));
        this.tm = new Float32Array(2 * (numRegions + numTriangles));
    }

    setMap(mesh, {r_xyz, t_xyz, r_color_fn, s_flow, r_elevation, t_elevation, r_moisture, t_moisture}) {
        const V = 0.95;
        const {numSides, numRegions, numTriangles} = mesh;
        const {xyz, tm, I} = this;

        xyz.set(r_xyz);
        xyz.set(t_xyz, r_xyz.length);
        // TODO: multiply all the r, t points by the elevation, taking V into account

        let p = 0;
        for (let r = 0; r < numRegions; r++) {
            tm[p++] = r_elevation[r];
            tm[p++] = r_moisture[r];
        }
        for (let t = 0; t < numTriangles; t++) {
            tm[p++] = t_elevation[t];
            tm[p++] = t_moisture[t];
        }

        let i = 0, count_valley = 0, count_ridge = 0;
        let {_halfedges, _triangles} = mesh;
        for (let s = 0; s < numSides; s++) {
            let opposite_s = mesh.s_opposite_s(s),
                r1 = mesh.s_begin_r(s),
                r2 = mesh.s_begin_r(opposite_s),
                t1 = mesh.s_inner_t(s),
                t2 = mesh.s_inner_t(opposite_s);
            
            // Each quadrilateral is turned into two triangles, so each
            // half-edge gets turned into one. There are two ways to fold
            // a quadrilateral. This is usually a nuisance but in this
            // case it's a feature. See the explanation here
            // https://www.redblobgames.com/x/1725-procedural-elevation/#rendering
            let coast = r_elevation[r1] < 0.0 || r_elevation[r2] < 0.0;
            if (coast || s_flow[s] > 0 || s_flow[opposite_s] > 0) {
                // It's a coastal or river edge, forming a valley
                I[i++] = r1; I[i++] = numRegions+t2; I[i++] = numRegions+t1;
                count_valley++;
            } else {
                // It's a ridge
                I[i++] = r1; I[i++] = r2; I[i++] = numRegions+t1;
                count_ridge++;
            }
        }

        console.log('ridge=', count_ridge, ', valley=', count_valley);
    }
}

/**********************************************************************
 * Plates
 */

function pickRandomRegions(mesh, N, randInt) {
    let {numRegions} = mesh;
    let chosen_r = new Set();
    while (chosen_r.size < N && chosen_r.size < numRegions) {
        chosen_r.add(randInt(numRegions));
    }
    return chosen_r;
}


function generatePlates(mesh, r_xyz) {
    let r_plate = new Int32Array(mesh.numRegions);
    r_plate.fill(-1);
    let plate_r = pickRandomRegions(mesh, Math.min(P, N), makeRandInt(SEED));
    let queue = Array.from(plate_r);
    for (let r of queue) { r_plate[r] = r; }
    let out_r = [];
    const randInt = makeRandInt(SEED);
    for (let queue_out = 0; queue_out < mesh.numRegions; queue_out++) {
        let pos = queue_out + randInt(queue.length - queue_out);
        let current_r = queue[pos];
        queue[pos] = queue[queue_out];
        mesh.r_circulate_r(out_r, current_r);
        for (let neighbor_r of out_r) {
            if (r_plate[neighbor_r] === -1) {
                r_plate[neighbor_r] = r_plate[current_r];
                queue.push(neighbor_r);
            }
        }
    }

    // Assign a random movement vector for each plate
    let plate_vec = [];
    for (let center_r of plate_r) {
        let neighbor_r = mesh.r_circulate_r([], center_r)[0];
        let p0 = r_xyz.slice(3 * center_r, 3 * center_r + 3),
            p1 = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
        plate_vec[center_r] = vec3.normalize([], vec3.subtract([], p1, p0));
    }

    return {plate_r, r_plate, plate_vec};
}


/* Distance from any point in seeds_r to all other points, but 
 * don't go past any point in stop_r */
function assignDistanceField(mesh, seeds_r, stop_r) {
    const randInt = makeRandInt(SEED);
    let {numRegions} = mesh;
    let r_distance = new Float32Array(numRegions);
    r_distance.fill(Infinity);
    
    let queue = [];
    for (let r of seeds_r) {
        queue.push(r);
        r_distance[r] = 0;
    }
    
    let out_r = [];
    for (let queue_out = 0; queue_out < mesh.numRegions; queue_out++) {
        let pos = queue_out + randInt(queue.length - queue_out);
        let current_r = queue[pos];
        queue[pos] = queue[queue_out];
        mesh.r_circulate_r(out_r, current_r);
        for (let neighbor_r of out_r) {
            if (r_distance[neighbor_r] === Infinity && !stop_r.has(neighbor_r)) {
                r_distance[neighbor_r] = r_distance[current_r] + 1;
                queue.push(neighbor_r);
            }
        }
    }
    return r_distance;
    // TODO: possible enhancement: keep track of which seed is closest
    // to this point, so that we can assign variable mountain/ocean
    // elevation to each seed instead of them always being +1/-1
}


/* Calculate the collision measure, which is the amount
 * that any neighbor's plate vector is pushing against 
 * the current plate vector. */
const COLLISION_THRESHOLD = 0.75;
function findCollisions(mesh, r_xyz, plate_is_ocean, r_plate, plate_vec) {
    const deltaTime = 1e-2; // simulate movement
    let {numRegions} = mesh;
    let mountain_r = new Set(),
        coastline_r = new Set(),
        ocean_r = new Set();
    let r_out = [];
    /* For each region, I want to know how much it's being compressed
       into an adjacent region. The "compression" is the change in
       distance as the two regions move. I'm looking for the adjacent
       region from a different plate that pushes most into this one*/
    for (let current_r = 0; current_r < numRegions; current_r++) {
        let bestCompression = Infinity, best_r = -1;
        mesh.r_circulate_r(r_out, current_r);
        for (let neighbor_r of r_out) {
            if (r_plate[current_r] !== r_plate[neighbor_r]) {
                /* sometimes I regret storing xyz in a compact array... */
                let current_pos = r_xyz.slice(3 * current_r, 3 * current_r + 3),
                    neighbor_pos = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
                /* simulate movement for deltaTime seconds */
                let distanceBefore = vec3.distance(current_pos, neighbor_pos),
                    distanceAfter = vec3.distance(vec3.add([], current_pos, vec3.scale([], plate_vec[r_plate[current_r]], deltaTime)),
                                                  vec3.add([], neighbor_pos, vec3.scale([], plate_vec[r_plate[neighbor_r]], deltaTime)));
                /* how much closer did these regions get to each other? */
                let compression = distanceBefore - distanceAfter;
                /* keep track of the adjacent region that gets closest */
                if (compression < bestCompression) {
                    best_r = neighbor_r;
                    bestCompression = compression;
                }
            }
        }
        if (best_r !== -1) {
            /* at this point, bestCompression tells us how much closer
               we are getting to the region that's pushing into us the most */
            let collided = bestCompression > COLLISION_THRESHOLD * deltaTime;
            if (plate_is_ocean.has(current_r) && plate_is_ocean.has(best_r)) {
                (collided? coastline_r : ocean_r).add(current_r);
            } else if (!plate_is_ocean.has(current_r) && !plate_is_ocean.has(best_r)) {
                if (collided) mountain_r.add(current_r);
            } else {
                (collided? mountain_r : coastline_r).add(current_r);
            }
        }
    }
    return {mountain_r, coastline_r, ocean_r};
}


function assignRegionElevation(mesh, {r_xyz, plate_is_ocean, r_plate, plate_vec, /* out */ r_elevation}) {
    const epsilon = 1e-3;
    let {numRegions} = mesh;

    let {mountain_r, coastline_r, ocean_r} = findCollisions(
        mesh, r_xyz, plate_is_ocean, r_plate, plate_vec);

    for (let r = 0; r < numRegions; r++) {
        if (r_plate[r] === r) {
            (plate_is_ocean.has(r)? ocean_r : coastline_r).add(r);
        }
    }

    let stop_r = new Set();
    for (let r of mountain_r) { stop_r.add(r); }
    for (let r of coastline_r) { stop_r.add(r); }
    for (let r of ocean_r) { stop_r.add(r); }

    console.log('seeds mountain/coastline/ocean:', mountain_r.size, coastline_r.size, ocean_r.size, 'plate_is_ocean', plate_is_ocean.size,'/', P);
    let r_distance_a = assignDistanceField(mesh, mountain_r, ocean_r);
    let r_distance_b = assignDistanceField(mesh, ocean_r, coastline_r);
    let r_distance_c = assignDistanceField(mesh, coastline_r, stop_r);

    for (let r = 0; r < numRegions; r++) {
        let a = r_distance_a[r] + epsilon,
            b = r_distance_b[r] + epsilon,
            c = r_distance_c[r] + epsilon;
        if (a === Infinity && b === Infinity) {
            r_elevation[r] = 0.1;
        } else {
            r_elevation[r] = (1/a - 1/b) / (1/a + 1/b + 1/c);
        }
        r_elevation[r] += 0.1 * fbm_noise(r_xyz[3*r], r_xyz[3*r+1], r_xyz[3*r+2]);
    }
}



/**********************************************************************
 * Rivers - from mapgen4
 */

function assignTriangleValues(mesh, {r_elevation, r_moisture, /* out */ t_elevation, t_moisture}) {
    const {numTriangles} = mesh;
    for (let t = 0; t < numTriangles; t++) {
        let s0 = 3*t;
        let r1 = mesh.s_begin_r(s0),
            r2 = mesh.s_begin_r(s0+1),
            r3 = mesh.s_begin_r(s0+2);
        t_elevation[t] = 1/3 * (r_elevation[r1] + r_elevation[r2] + r_elevation[r3]);
        t_moisture[t] = 1/3 * (r_moisture[r1] + r_moisture[r2] + r_moisture[r3]);
    }
}


let _queue = new FlatQueue();
function assignDownflow(mesh, {t_elevation, /* out */ t_downflow_s, /* out */ order_t}) {
    /* Use a priority queue, starting with the ocean triangles and
     * moving upwards using elevation as the priority, to visit all
     * the land triangles */
    let {numTriangles} = mesh,
        queue_in = 0;
    t_downflow_s.fill(-999);
    /* Part 1: ocean triangles get downslope assigned to the lowest neighbor */
    for (let t = 0; t < numTriangles; t++) {
        if (t_elevation[t] < 0) {
            let best_s = -1, best_e = t_elevation[t];
            for (let j = 0; j < 3; j++) {
                let s = 3 * t + j,
                    e = t_elevation[mesh.s_outer_t(s)];
                if (e < best_e) {
                    best_e = e;
                    best_s = s;
                }
            }
            order_t[queue_in++] = t;
            t_downflow_s[t] = best_s;
            _queue.push(t, t_elevation[t]);
        }
    }
    /* Part 2: land triangles get visited in elevation priority */
    for (let queue_out = 0; queue_out < numTriangles; queue_out++) {
        let current_t = _queue.pop();
        for (let j = 0; j < 3; j++) {
            let s = 3 * current_t + j;
            let neighbor_t = mesh.s_outer_t(s); // uphill from current_t
            if (t_downflow_s[neighbor_t] === -999 && t_elevation[neighbor_t] >= 0.0) {
                t_downflow_s[neighbor_t] = mesh.s_opposite_s(s);
                order_t[queue_in++] = neighbor_t;
                _queue.push(neighbor_t, t_elevation[neighbor_t]);
            }
        }
    }
}


function assignFlow(mesh, {order_t, t_elevation, t_moisture, t_downflow_s, /* out */ t_flow, /* out */ s_flow}) {
    let {numTriangles, _halfedges} = mesh;
    s_flow.fill(0);
    for (let t = 0; t < numTriangles; t++) {
        if (t_elevation[t] >= 0.0) {
            t_flow[t] = 0.5 * t_moisture[t] * t_moisture[t];
        } else {
            t_flow[t] = 0;
        }
    }
    for (let i = order_t.length-1; i >= 0; i--) {
        let tributary_t = order_t[i];
        let flow_s = t_downflow_s[tributary_t];
        let trunk_t = (_halfedges[flow_s] / 3) | 0;
        if (flow_s >= 0) {
            t_flow[trunk_t] += t_flow[tributary_t];
            s_flow[flow_s] += t_flow[tributary_t]; // TODO: isn't s_flow[flow_s] === t_flow[?]
            if (t_elevation[trunk_t] > t_elevation[tributary_t]) {
                t_elevation[trunk_t] = t_elevation[tributary_t];
            }
        }
    }
}




/**********************************************************************
 * Main
 */

// ugh globals, sorry
var mesh, map = {};
var quadGeometry = new QuadGeometry();

function generateMesh() {
    let result = SphereMesh.makeSphere(N, jitter, makeRandFloat(SEED));
    mesh = result.mesh;
    quadGeometry.setMesh(mesh);
    
    map.r_xyz = result.r_xyz;
    map.t_xyz = generateTriangleCenters(mesh, map);
    map.r_elevation = new Float32Array(mesh.numRegions);
    map.t_elevation = new Float32Array(mesh.numTriangles);
    map.r_moisture = new Float32Array(mesh.numRegions);
    map.t_moisture = new Float32Array(mesh.numTriangles);
    map.t_downflow_s = new Int32Array(mesh.numTriangles);
    map.order_t = new Int32Array(mesh.numTriangles);
    map.t_flow = new Float32Array(mesh.numTriangles);
    map.s_flow = new Float32Array(mesh.numSides);

    generateMap();
}

function generateMap() {
    let result = generatePlates(mesh, map.r_xyz);
    map.plate_r = result.plate_r;
    map.r_plate = result.r_plate;
    map.plate_vec = result.plate_vec;
    map.plate_is_ocean = new Set();
    for (let r of map.plate_r) {
        if (makeRandInt(r)(10) < 5) {
            map.plate_is_ocean.add(r);
            // TODO: either make tiny plates non-ocean, or make sure tiny plates don't create seeds for rivers
        }
    }
    assignRegionElevation(mesh, map);
    // TODO: assign region moisture in a better way!
    for (let r = 0; r < mesh.numRegions; r++) {
        map.r_moisture[r] = (map.r_plate[r] % 10) / 10.0;
    }
    assignTriangleValues(mesh, map);
    assignDownflow(mesh, map);
    assignFlow(mesh, map);

    quadGeometry.setMap(mesh, map);
    draw();
}


function drawPlateVectors(u_projection, mesh, {r_xyz, r_plate, plate_vec}) {
    let line_xyz = [], line_rgba = [];
    
    for (let r = 0; r < mesh.numRegions; r++) {
        line_xyz.push(r_xyz.slice(3 * r, 3 * r + 3));
        line_rgba.push([1, 1, 1, 1]);
        line_xyz.push(vec3.add([], r_xyz.slice(3 * r, 3 * r + 3),
                               vec3.scale([], plate_vec[r_plate[r]], 2 / Math.sqrt(N))));
        line_rgba.push([1, 0, 0, 0]);
    }

    renderLines({
        u_projection,
        u_multiply_rgba: [1, 1, 1, 1],
        u_add_rgba: [0, 0, 0, 0],
        a_xyz: line_xyz,
        a_rgba: line_rgba,
        count: line_xyz.length,
    });
}

function drawPlateBoundaries(u_projection, mesh, {t_xyz, r_plate}) {
    let line_xyz = [], line_rgba = [];
    for (let s = 0; s < mesh.numSides; s++) {
        let begin_r = mesh.s_begin_r(s),
            end_r = mesh.s_end_r(s);
        if (r_plate[begin_r] !== r_plate[end_r]) {
            let inner_t = mesh.s_inner_t(s),
                outer_t = mesh.s_outer_t(s);
            line_xyz.push(t_xyz.slice(3 * inner_t, 3 * inner_t + 3),
                          t_xyz.slice(3 * outer_t, 3 * outer_t + 3));
            line_rgba.push([1, 1, 1, 1], [1, 1, 1, 1]);
        }
    }
    renderLines({
        u_projection,
        u_multiply_rgba: [1, 1, 1, 1],
        u_add_rgba: [0, 0, 0, 0],
        a_xyz: line_xyz,
        a_rgba: line_rgba,
        count: line_xyz.length,
    });
}

function drawRivers(u_projection, mesh, {t_xyz, s_flow}) {
    let line_xyz = [], line_rgba = [];

    for (let s = 0; s < mesh.numSides; s++) {
        if (s_flow[s] > 1) {
            let flow = 0.1 * Math.sqrt(s_flow[s]);
            let inner_t = mesh.s_inner_t(s),
                outer_t = mesh.s_outer_t(s);
            line_xyz.push(t_xyz.slice(3 * inner_t, 3 * inner_t + 3),
                          t_xyz.slice(3 * outer_t, 3 * outer_t + 3));
            if (flow > 1) flow = 1;
            let rgba_premultiplied = [0.2 * flow, 0.5 * flow, 0.7 * flow, flow];
            line_rgba.push(rgba_premultiplied, rgba_premultiplied);
        }
    }
    renderLines({
        u_projection,
        u_multiply_rgba: [1, 1, 1, 1],
        u_add_rgba: [0, 0, 0, 0],
        a_xyz: line_xyz,
        a_rgba: line_rgba,
        count: line_xyz.length,
    });
}

function draw() {
    let u_pointsize = 0.1 + 100 / Math.sqrt(N);
    let u_projection = mat4.create();
    mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
    mat4.rotate(u_projection, u_projection, -rotation, [0.1, 1, 0]);
    mat4.rotate(u_projection, u_projection, -Math.PI/2+0.2, [1, 0, 0]);

    function r_color_fn(r) {
        let m = map.r_moisture[r];
        let e = map.r_elevation[r];
        return [e, m];
    }

    if (drawMode === 'centroid') {
        let triangleGeometry = generateVoronoiGeometry(mesh, map, r_color_fn);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.xyz,
            a_tm: triangleGeometry.tm,
            count: triangleGeometry.xyz.length / 3,
        });
    } else if (drawMode === 'quads') {
        renderIndexedTriangles({
            u_projection,
            a_xyz: quadGeometry.xyz,
            a_tm: quadGeometry.tm,
            elements: quadGeometry.I,
        });
    }

    drawRivers(u_projection, mesh, map);
    
    if (draw_plateVectors) {
        drawPlateVectors(u_projection, mesh, map);
    }
    if (draw_plateBoundaries) {
        drawPlateBoundaries(u_projection, mesh, map);
    }
    
    // renderPoints({
    //     u_projection,
    //     u_pointsize,
    //     a_xyz: map.r_xyz,
    //     count: mesh.numRegions,
    // });
}

generateMesh();
