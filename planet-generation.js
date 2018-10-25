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


const renderTriangles = regl({
    frag: `
precision mediump float;
uniform sampler2D u_colormap;
varying vec3 v_rgb;
void main() {
   gl_FragColor = texture2D(u_colormap, v_rgb.xy);
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


/* use global variable drawMode to decide which type of triangle center to push */
function pushCenterOfTriangle(out, ax, ay, az, bx, by, bz, cx, cy, cz) {
    ((drawMode === 'centroid')? pushCentroidOfTriangle : pushCircumcenterOfTriangle)(out, ax, ay, az, bx, by, bz, cx, cy, cz);
}

function generateVoronoiGeometry(r_xyz, r_color_fn, mesh) {
    const {numTriangles, numSides} = mesh;
    let centers = [];
    let geometry = [], colors = [];

    for (let t = 0; t < numTriangles; t++) {
        // TODO: flag for circumcenter vs centroid
        let a = mesh.s_begin_r(3*t), b = mesh.s_begin_r(3*t+1), c = mesh.s_begin_r(3*t+2);
        pushCenterOfTriangle(centers,
                 r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
                 r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
                 r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]);
    }

    for (let s = 0; s < numSides; s++) {
        let inner_t = mesh.s_inner_t(s),
            outer_t = mesh.s_outer_t(s),
            begin_r = mesh.s_begin_r(s);
        let rgb = r_color_fn(begin_r);
        geometry.push(centers[3*inner_t], centers[3*inner_t+1], centers[3*inner_t+2],
                      centers[3*outer_t], centers[3*outer_t+1], centers[3*outer_t+2],
                      r_xyz[3*begin_r], r_xyz[3*begin_r+1], r_xyz[3*begin_r+2]);
        colors.push(rgb, rgb, rgb);
    }
    return {geometry, colors, centers};
}


function pickRandomRegions(mesh, N, randInt) {
    let {numRegions} = mesh;
    let chosen_r = new Set();
    while (chosen_r.size < N) {
        chosen_r.add(randInt(numRegions));
    }
    return chosen_r;
}


let algorithm = 2;
let N = 2000;
let jitter = 0.5;
let rotation = -1.5;
let drawMode = 'centroid';

window.setAlgorithm = newAlgorithm => { algorithm = newAlgorithm; draw(); };
window.setN = newN => { N = newN; draw(); };
window.setJitter = newJitter => { jitter = newJitter; draw(); };
window.setRotation = newRotation => { rotation = newRotation; draw(); };
window.setDrawMode = newMode => { drawMode = newMode; draw(); };

function draw() {
    let u_pointsize = (drawMode === 'points'? 5 : 1) * (0.1 + 100 / Math.sqrt(N));
    let u_projection = mat4.create();
    mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
    mat4.rotate(u_projection, u_projection, -rotation, [0.1, 1, 0]);
    mat4.rotate(u_projection, u_projection, -Math.PI/2+0.2, [1, 0, 0]);
    
    let {mesh, r_xyz} = SphereMesh.makeSphere(algorithm, N, jitter);
    function r_color_fn(r) {
        return randomColor(r_plate_r[r], r_xyz);
    }

    if (drawMode === 'delaunay') {
        let triangleGeometry = generateDelaunayGeometry(r_xyz, r_color_fn, mesh);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.geometry,
            a_rgb: triangleGeometry.colors,
            count: triangleGeometry.geometry.length / 3,
        });
    } else if (drawMode === 'voronoi' || drawMode === 'centroid') {
        let triangleGeometry = generateVoronoiGeometry(r_xyz, r_color_fn, mesh);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.geometry,
            a_rgb: triangleGeometry.colors,
            count: triangleGeometry.geometry.length / 3,
        });
    }
    
    renderPoints({
        u_projection,
        u_pointsize,
        u_brightness: drawMode === 'points'? 1.0 : 0.0,
        a_xyz: r_xyz,
        count: mesh.numRegions,
    });
}

draw();
