/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

const SEED = 123;

const SimplexNoise = require('simplex-noise');
const colormap = require('../../maps/mapgen4/colormap');
const {vec3, mat4} = require('gl-matrix');
const {makeRandInt, makeRandFloat} = require('@redblobgames/prng');
const SphereMesh = require('./sphere-mesh');

const regl = require('regl')({
    canvas: "#output",
    extensions: ['OES_element_index_uint']
});

const u_colormap = regl.texture({
    width: colormap.width,
    height: colormap.height,
    data: colormap.data,
    wrapS: 'clamp',
    wrapT: 'clamp'
});

const renderPoints = regl({
    frag: `
precision mediump float;
varying vec3 v_rgb;
void main() {
   gl_FragColor = vec4(v_rgb, 1);
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
uniform float u_pointsize;
uniform float u_brightness;
attribute vec3 a_xyz;
varying vec3 v_rgb;
void main() {
  v_rgb = (a_xyz + 1.0) / 2.0 * vec3(1.0, 0.8, 0.9) * u_brightness;
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
        u_brightness: regl.prop('u_brightness'),
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
  v_rgba = pos.z > 0.0? vec4(0, 0, 0, 0) : a_rgba;
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
        func: {src:'one', dst:'one minus src alpha'},
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
varying vec3 v_rgb;
void main() {
   gl_FragColor = texture2D(u_colormap, v_rgb.xy);
   // gl_FragColor = vec4(v_rgb, 1);
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
attribute vec3 a_xyz;
attribute vec3 a_rgb;
varying vec3 v_rgb;
void main() {
  v_rgb = a_rgb;
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
        a_rgb: regl.prop('a_rgb'),
    },
});


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

function randomColor(r, r_xyz) {
    let x = r_xyz[3*r], y = r_xyz[3*r+1], z = r_xyz[3*r+2];
    let e = fbm_noise(x, y, z) - 0.1,
        m = fbm_noise(3*x+5, 3*y+5, z+5) + 2.5 * (0.5 - Math.abs(z));
    if (e > 0) e = e * e + z*z*z*z;
    return [
        0.5 * (1 + e),
        0.5 * (1 + m),
        (r / 100) % 1.0,
    ];
}

function generateDelaunayGeometry(r_xyz, r_color_fn, mesh) {
    const {numTriangles} = mesh;
    let geometry = [], colors = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = mesh.s_begin_r(3*t), b = mesh.s_begin_r(3*t+1), c = mesh.s_begin_r(3*t+2);
        geometry.push(
            r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
            r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
            r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]
        );
        colors.push(r_color_fn(a),
                    r_color_fn(b),
                    r_color_fn(c));
    }
    return {geometry, colors};
}


/* Calculate the centroid and push it onto an array */
function pushCentroidOfTriangle(out, ax, ay, az, bx, by, bz, cx, cy, cz) {
    // TODO: renormalize to radius 1
    out.push((ax+bx+cx)/3, (ay+by+cy)/3, (az+bz+cz)/3);
}


/* Calculate the circumenter and push it onto an array */
function pushCircumcenterOfTriangle(out, ax, ay, az, bx, by, bz, cx, cy, cz) {
    // https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
    let a = [ax, ay, az],
        b = [bx, by, bz],
        c = [cx, cy, cz],
        ac = vec3.subtract([], c, a),
        ab = vec3.subtract([], b, a),
        abXac = vec3.cross([], ab, ac),
        numerator = vec3.add([],
                             vec3.scale([], vec3.cross([], abXac, ab), vec3.squaredLength(ac)),
                             vec3.scale([], vec3.cross([], ac, abXac), vec3.squaredLength(ab))),
        toCircumsphereCenter = vec3.scale([],
                                          numerator,
                                          1 / (2 * vec3.squaredLength(abXac))
                                         ),
        circumcenter = vec3.add([], a, toCircumsphereCenter);
    // TODO: renormalize to radius 1
    out.push(circumcenter[0], circumcenter[1], circumcenter[2]);
}


function generateTriangleCenters(mesh, r_xyz, pushCenter) {
    let {numTriangles} = mesh;
    let t_xyz = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = mesh.s_begin_r(3*t),
            b = mesh.s_begin_r(3*t+1),
            c = mesh.s_begin_r(3*t+2);
        pushCenter(t_xyz,
                 r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
                 r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
                 r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]);
    }
    return t_xyz;
}

function generateVoronoiGeometry(mesh, r_xyz, t_xyz, r_color_fn) {
    const {numSides} = mesh;
    let geometry = [], colors = [];

    for (let s = 0; s < numSides; s++) {
        let inner_t = mesh.s_inner_t(s),
            outer_t = mesh.s_outer_t(s),
            begin_r = mesh.s_begin_r(s);
        let rgb = r_color_fn(begin_r);
        geometry.push(t_xyz[3*inner_t], t_xyz[3*inner_t+1], t_xyz[3*inner_t+2],
                      t_xyz[3*outer_t], t_xyz[3*outer_t+1], t_xyz[3*outer_t+2],
                      r_xyz[3*begin_r], r_xyz[3*begin_r+1], r_xyz[3*begin_r+2]);
        colors.push(rgb, rgb, rgb);
    }
    return {geometry, colors};
}


function pickRandomRegions(mesh, N, randInt) {
    let {numRegions} = mesh;
    let chosen_r = new Set();
    while (chosen_r.size < N && chosen_r.size < numRegions) {
        chosen_r.add(randInt(numRegions));
    }
    return chosen_r;
}


let algorithm = 2;
let N = 10000;
let P = 20;
let jitter = 0.5;
let rotation = -1.5;
let drawMode = 'centroid';

window.setAlgorithm = newAlgorithm => { algorithm = newAlgorithm; generate(); };
window.setN = newN => { N = newN; generate(); };
window.setP = newP => { P = newP; generate(); };
window.setJitter = newJitter => { jitter = newJitter; generate(); };
window.setRotation = newRotation => { rotation = newRotation; draw(); };
window.setDrawMode = newMode => { drawMode = newMode; draw(); };

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
    let epsilon = 1e-2;
    let {numRegions} = mesh;
    let mountain_r = new Set(),
        coastline_r = new Set(),
        ocean_r = new Set();
    let r_out = [];
    for (let current_r = 0; current_r < numRegions; current_r++) {
        let bestCollision = Infinity, best_r = -1;
        mesh.r_circulate_r(r_out, current_r);
        for (let neighbor_r of r_out) {
            if (r_plate[current_r] !== r_plate[neighbor_r]) {
                let current_pos = r_xyz.slice(3 * current_r, 3 * current_r + 3),
                    neighbor_pos = r_xyz.slice(3 * neighbor_r, 3 * neighbor_r + 3);
                let distanceBefore = vec3.distance(current_pos, neighbor_pos),
                    distanceAfter = vec3.distance(vec3.add([], current_pos, vec3.scale([], plate_vec[r_plate[current_r]], epsilon)),
                                                  vec3.add([], neighbor_pos, vec3.scale([], plate_vec[r_plate[neighbor_r]], epsilon)));
                let collision = distanceBefore - distanceAfter;
                if (collision < bestCollision) {
                    best_r = neighbor_r;
                    bestCollision = collision;
                }
            }
        }
        if (best_r !== -1) {
            let collided = bestCollision > COLLISION_THRESHOLD * epsilon;
            if (plate_is_ocean.has(current_r) && plate_is_ocean.has(best_r)) {
                (collided? ocean_r : ocean_r).add(current_r);
            } else if (!plate_is_ocean.has(current_r) && !plate_is_ocean.has(best_r)) {
                if (collided) mountain_r.add(current_r);
            } else {
                (collided? mountain_r : coastline_r).add(current_r);
            }
        }
    }
    return {mountain_r, coastline_r, ocean_r};
}


function assignElevation(mesh, r_xyz, plate_is_ocean, r_plate, plate_vec) {
    const epsilon = 1e-3;
    let {numRegions} = mesh;
    let r_elevation = new Float32Array(numRegions);

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

    console.log(mountain_r.size, coastline_r.size, ocean_r.size);
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
    }
    return r_elevation;
}

var mesh, r_xyz, t_xyz, plate_r, r_plate, plate_vec, plate_is_ocean, r_elevation;

function generate() {
    let result = SphereMesh.makeSphere(algorithm, N, jitter);
    mesh = result.mesh;
    r_xyz = result.r_xyz;
    t_xyz = generateTriangleCenters(
        mesh, r_xyz,
        drawMode === 'centroid'? pushCentroidOfTriangle : pushCircumcenterOfTriangle
    );
                                        
    result = generatePlates(mesh, r_xyz);
    plate_r = result.plate_r;
    r_plate = result.r_plate;
    plate_vec = result.plate_vec;
    plate_is_ocean = new Set();
    for (let r of plate_r) {
        if (makeRandInt(r)(10) < 1) {
            plate_is_ocean.add(r);
        }
    }
    r_elevation = assignElevation(mesh, r_xyz, plate_is_ocean, r_plate, plate_vec);

    draw();
}

function draw() {
    function r_color_fn(r) {
        // return [plate_is_ocean.has(r_plate[r])? 0.25 : 0.6, 0, 0];
        let color = randomColor(r_plate[r], r_xyz);
        let e = r_elevation[r] - 0.2;
        if (e > 0) e = e * e;
        return [0.5 + 0.5 * e + 0.1 * (color[0] - 0.5), color[1], 1];
    }

    let u_pointsize = (drawMode === 'points'? 5 : 1) * (0.1 + 100 / Math.sqrt(N));
    let u_projection = mat4.create();
    mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
    mat4.rotate(u_projection, u_projection, -rotation, [0.1, 1, 0]);
    mat4.rotate(u_projection, u_projection, -Math.PI/2+0.2, [1, 0, 0]);
    
    if (drawMode === 'delaunay') {
        let triangleGeometry = generateDelaunayGeometry(r_xyz, r_color_fn, mesh);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.geometry,
            a_rgb: triangleGeometry.colors,
            count: triangleGeometry.geometry.length / 3,
        });
    } else if (drawMode === 'voronoi' || drawMode === 'centroid') {
        let triangleGeometry = generateVoronoiGeometry(mesh, r_xyz, t_xyz, r_color_fn);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.geometry,
            a_rgb: triangleGeometry.colors,
            count: triangleGeometry.geometry.length / 3,
        });
    }

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

    line_xyz = []; line_rgba = [];
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
                          
    renderPoints({
        u_projection,
        u_pointsize,
        u_brightness: drawMode === 'points'? 1.0 : 0.0,
        a_xyz: r_xyz,
        count: mesh.numRegions,
    });
}

generate();
