/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

const Delaunator = require('delaunator');
const {vec3, mat4} = require('gl-matrix');
const regl = require('regl')({
    canvas: "#output",
    extensions: ['OES_element_index_uint']
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
varying vec3 v_rgb;
void main() {
   gl_FragColor = vec4(v_rgb, 1);
}
`,

    vert: `
precision mediump float;
uniform mat4 u_projection;
attribute vec3 a_xyz;
attribute vec3 a_rgb;
varying vec3 v_rgb;
void main() {
  v_rgb = mix((a_xyz + 1.0) / 2.0 * vec3(1.5, 1.0, 1.3), a_rgb, 0.33);
  gl_Position = u_projection * vec4(a_xyz, 1);
}
`,

    uniforms: {
        u_projection: regl.prop('u_projection'),
    },

    count: regl.prop('count'),
    attributes: {
        a_xyz: regl.prop('a_xyz'),
        a_rgb: regl.prop('a_rgb'),
    },
});


let _randomLat = [], _randomLon = [];
function generateFibonacciSphere1(N, jitter) {
    let a_latlong = [];

    // First algorithm from http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
    const s = 3.6/Math.sqrt(N);
    const dz = 2.0/N;
    for (let k = 0, long = 0, z = 1 - dz/2; k !== N; k++, z -= dz) {
        let r = Math.sqrt(1-z*z);
        let latDeg = Math.asin(z) * 180 / Math.PI;
        let lonDeg = long * 180 / Math.PI;
        if (_randomLat[k] === undefined) _randomLat[k] = Math.random() - Math.random();
        if (_randomLon[k] === undefined) _randomLon[k] = Math.random() - Math.random();
        latDeg += jitter * _randomLat[k] * (latDeg - Math.asin(Math.max(-1, z - dz * 2 * Math.PI * r / s)) * 180 / Math.PI);
        lonDeg += jitter * _randomLon[k] * (s/r * 180 / Math.PI);
        a_latlong.push(latDeg, lonDeg % 360.0);
        long += s/r;
    }
    return a_latlong;
}

function generateFibonacciSphere2(N, jitter) {
    let a_latlong = [];

    // Second algorithm from http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
    const s = 3.6/Math.sqrt(N);
    const dlong = Math.PI * (3-Math.sqrt(5));  /* ~2.39996323 */
    const dz = 2.0 / N;
    for (let k = 0, long = 0, z = 1 - dz/2; k !== N; k++, z -= dz) {
        let r = Math.sqrt(1 - z*z);
        let latDeg = Math.asin(z) * 180 / Math.PI;
        let lonDeg = long * 180 / Math.PI;
        if (_randomLat[k] === undefined) _randomLat[k] = Math.random() - Math.random();
        if (_randomLon[k] === undefined) _randomLon[k] = Math.random() - Math.random();
        latDeg += jitter * _randomLat[k] * (latDeg - Math.asin(Math.max(-1, z - dz * 2 * Math.PI * r / s)) * 180 / Math.PI);
        lonDeg += jitter * _randomLon[k] * (s/r * 180 / Math.PI);
        a_latlong.push(latDeg, lonDeg % 360.0);
        long += dlong;
    }
    return a_latlong;
}


/* calculate x,y,z from spherical coordinates lat,lon and then push
 * them onto out array; for one-offs pass [] as the first argument */
function pushCartesianFromSpherical(out, latDeg, lonDeg) {
    let latRad = latDeg / 180.0 * Math.PI,
        lonRad = lonDeg / 180.0 * Math.PI;
    out.push(Math.cos(latRad) * Math.cos(lonRad),
             Math.cos(latRad) * Math.sin(lonRad),
             Math.sin(latRad));
    return out;
}



/** Add south pole back into the mesh.
 *
 * We run the Delaunay Triangulation on all points *except* the south
 * pole, which gets mapped to infinity with the stereographic
 * projection. This function adds the south pole into the
 * triangulation. The Delaunator guide explains how the halfedges have
 * to be connected to make the mesh work.
 * <https://mapbox.github.io/delaunator/>
 *
 * Returns the new {triangles, halfedges} for the triangulation with
 * one additional point added around the convex hull.
 */
function addSouthPoleToMesh(southPoleId, {triangles, halfedges}) {
    // This logic is from <https://github.com/redblobgames/dual-mesh>,
    // where I use it to insert a "ghost" region on the "back" side of
    // the planar map. The same logic works here. In that code I use
    // "s" for edges ("sides"), "r" for regions ("points"), t for triangles
    let numSides = triangles.length;
    function s_next_s(s) { return (s % 3 == 2) ? s-2 : s+1; }

    let numUnpairedSides = 0, firstUnpairedSide = -1;
    let pointIdToSideId = []; // seed to side
    for (let s = 0; s < numSides; s++) {
        if (halfedges[s] === -1) {
            numUnpairedSides++;
            pointIdToSideId[triangles[s]] = s;
            firstUnpairedSide = s;
        }
    }
    
    let newTriangles = new Int32Array(numSides + 3 * numUnpairedSides);
    let newHalfedges = new Int32Array(numSides + 3 * numUnpairedSides);
    newTriangles.set(triangles);
    newHalfedges.set(halfedges);

    for (let i = 0, s = firstUnpairedSide;
         i < numUnpairedSides;
         i++, s = pointIdToSideId[newTriangles[s_next_s(s)]]) {

        // Construct a pair for the unpaired side s
        let newSide = numSides + 3 * i;
        newHalfedges[s] = newSide;
        newHalfedges[newSide] = s;
        newTriangles[newSide] = newTriangles[s_next_s(s)];
        
        // Construct a triangle connecting the new side to the south pole
        newTriangles[newSide + 1] = newTriangles[s];
        newTriangles[newSide + 2] = southPoleId;
        let k = numSides + (3 * i + 4) % (3 * numUnpairedSides);
        newHalfedges[newSide + 2] = k;
        newHalfedges[k] = newSide + 2;
    }

    return {
        triangles: newTriangles,
        halfedges: newHalfedges,
    };
}


/* Pick some pastel colors per index and cache them.
 * Not all operations will benefit from this, but some
 * (like rotation) do, so might as well cache. */
let _randomColor = [];
function randomColor(index) {
    if (!_randomColor[index]) {
        _randomColor[index] = [
            0.5 + Math.random() * 0.5,
            0.6 + Math.random() * 0.3,
            0.5 + Math.random() * 0.5,
        ];
    }
    return _randomColor[index];
}

function generateDelaunayGeometry(r_xyz, delaunay) {
    let {triangles} = delaunay;
    let numTriangles = triangles.length / 3;
    let geometry = [], colors = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = triangles[3*t], b = triangles[3*t+1], c = triangles[3*t+2];
        let rgb = randomColor(a+b+c);
        geometry.push(
            r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
            r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
            r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]
        );
        colors.push(rgb, rgb, rgb);
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

function generateVoronoiGeometry(r_xyz, delaunay) {
    let {triangles, halfedges} = delaunay;
    let numPoints = r_xyz.length / 3,
        numSides = triangles.length,
        numTriangles = numSides / 3;
    let centers = [];
    let geometry = [], colors = [];

    for (let t = 0; t < numTriangles; t++) {
        let a = triangles[3*t], b = triangles[3*t+1], c = triangles[3*t+2];
        // TODO: flag for circumcenter vs centroid
        pushCenterOfTriangle(centers,
                 r_xyz[3*a], r_xyz[3*a+1], r_xyz[3*a+2],
                 r_xyz[3*b], r_xyz[3*b+1], r_xyz[3*b+2],
                 r_xyz[3*c], r_xyz[3*c+1], r_xyz[3*c+2]);
    }

    for (let s = 0; s < numSides; s++) {
        let inner_t = (s / 3) | 0,
            outer_t = (halfedges[s] / 3) | 0,
            begin_r = triangles[s];
        let rgb = randomColor(begin_r);
        geometry.push(centers[3*inner_t], centers[3*inner_t+1], centers[3*inner_t+2],
                      centers[3*outer_t], centers[3*outer_t+1], centers[3*outer_t+2],
                      r_xyz[3*begin_r], r_xyz[3*begin_r+1], r_xyz[3*begin_r+2]);
        colors.push(rgb, rgb, rgb);
    }
    return {geometry, colors, centers};
}


function stereographicProjection(r_xyz) {
    // See <https://en.wikipedia.org/wiki/Stereographic_projection>
    const degToRad = Math.PI / 180;
    let numPoints = r_xyz.length / 3;
    let r_XY = [];
    for (let r = 0; r < numPoints; r++) {
        let x = r_xyz[3*r],
            y = r_xyz[3*r + 1],
            z = r_xyz[3*r + 2];
        let X = x / (1-z),
            Y = y / (1-z);
        r_XY.push(X, Y);
    }
    return r_XY;
}

let algorithm = 1;
function generateFibonacciSphere(N, jitter) {
    return [null, generateFibonacciSphere1, generateFibonacciSphere2][algorithm](N, jitter);
}

let N = 1000;
let jitter = 0.0;
let rotation = -4;
let drawMode = 'points';

window.setAlgorithm = newAlgorithm => { algorithm = newAlgorithm; draw(); };
window.setN = newN => { N = newN; draw(); };
window.setJitter = newJitter => { jitter = newJitter; draw(); };
window.setRotation = newRotation => { rotation = newRotation; draw(); };
window.setDrawMode = newMode => { drawMode = newMode; draw(); };

function draw() {
    let u_pointsize = (drawMode === 'points'? 5 : 1) * (0.1 + 100 / Math.sqrt(N));
    let u_projection = mat4.create();
    mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
    mat4.rotate(u_projection, u_projection, -rotation, [1, 0.5, 0]);
    
    let latlong = generateFibonacciSphere(N, jitter);
    let r_xyz = [];
    for (let r = 0; r < latlong.length/2; r++) {
        pushCartesianFromSpherical(r_xyz, latlong[2*r], latlong[2*r+1]);
    }

    let delaunay = new Delaunator(stereographicProjection(r_xyz));

    /* TODO: rotate an existing point into this spot instead of creating one */
    r_xyz.push(0, 0, 1);
    delaunay = addSouthPoleToMesh(r_xyz.length/3 - 1, delaunay);

    if (drawMode === 'delaunay') {
        let triangleGeometry = generateDelaunayGeometry(r_xyz, delaunay);
        renderTriangles({
            u_projection,
            a_xyz: triangleGeometry.geometry,
            a_rgb: triangleGeometry.colors,
            count: triangleGeometry.geometry.length / 3,
        });
    } else if (drawMode === 'voronoi' || drawMode === 'centroid') {
        let triangleGeometry = generateVoronoiGeometry(r_xyz, delaunay);
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
        count: r_xyz.length / 3,
    });
}

draw();
