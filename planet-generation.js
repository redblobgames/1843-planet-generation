/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

/* global Delaunator, createREGL, vec3, mat4 */

const regl = createREGL({
    canvas: "#output",
    extensions: ['OES_element_index_uint']
});

const renderPoints = regl({
    frag: `
precision mediump float;
varying vec2 v_uv;
void main() {
   gl_FragColor = vec4(0, 0, 0, 1);
}
`,

    vert: `
#define M_PI 3.1415926535897932384626433832795
precision highp float;
uniform mat4 u_projection;
uniform float u_pointsize;
attribute vec2 a_latlong;
varying vec2 v_uv;
void main() {
  v_uv = vec2(a_latlong.r / 180.0 + 0.5, fract((a_latlong.g + 90.0) / 360.0));
  float lat = a_latlong.r / 180.0 * M_PI;
  float lon = a_latlong.g / 180.0 * M_PI;
  gl_Position = u_projection * 
                 vec4(cos(lat) * cos(lon),
                      cos(lat) * sin(lon),
                      sin(lat),
                      1);
  gl_PointSize = gl_Position.z > 0.0? 1.0 : u_pointsize;
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
        a_latlong: regl.prop('a_latlong'),
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
precision highp float;
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


let randomLat = [], randomLon = [];
function generateFibonacciSphere1(N, jitter) {
    let a_latlong = [];

    // First algorithm from http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
    const s = 3.6/Math.sqrt(N);
    const dz = 2.0/N;
    for (let k = 0, long = 0, z = 1 - dz/2; k !== N; k++, z -= dz) {
        let r = Math.sqrt(1-z*z);
        let latDeg = Math.asin(z) * 180 / Math.PI;
        let lonDeg = long * 180 / Math.PI;
        if (randomLat[k] === undefined) randomLat[k] = Math.random() - Math.random();
        if (randomLon[k] === undefined) randomLon[k] = Math.random() - Math.random();
        latDeg += jitter * randomLat[k] * (latDeg - Math.asin(Math.max(-1, z - dz * 2 * Math.PI * r / s)) * 180 / Math.PI);
        lonDeg += jitter * randomLon[k] * (s/r * 180 / Math.PI);
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
        if (randomLat[k] === undefined) randomLat[k] = Math.random() - Math.random();
        if (randomLon[k] === undefined) randomLon[k] = Math.random() - Math.random();
        latDeg += jitter * randomLat[k] * (latDeg - Math.asin(Math.max(-1, z - dz * 2 * Math.PI * r / s)) * 180 / Math.PI);
        lonDeg += jitter * randomLon[k] * (s/r * 180 / Math.PI);
        a_latlong.push(latDeg, lonDeg % 360.0);
        long += dlong;
    }
    return a_latlong;
}


/** Add south pole back into the mesh.
 *
 * We run the Delaunay Triangulation on all points *except* the south pole,
 * which gets mapped to infinity with the stereographic projection. This
 * function adds the south pole into the triangulation.
 *
 * Returns the new {triangles, halfedges} for the triangulation with
 * one additional point added around the convex hull.
 */
function addSouthPoleToMesh(southPoleId, {triangles, halfedges}) {
    // This logic is from <https://github.com/redblobgames/dual-mesh>,
    // where I use it to insert a "ghost" region on the "back" side
    // of the planar map. The same logic works here.
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


function generateDelaunayGeometry(latlong, delaunay) {
    let {triangles} = delaunay;
    let numTriangles = triangles.length / 3;
    let geometry = [], colors = [];
    for (let t = 0; t < numTriangles; t++) {
        let a = triangles[3*t], b = triangles[3*t+1], c = triangles[3*t+2];
        let rgb = randomColor(a+b+c);
        pushCartesianFromSpherical(geometry, latlong[2*a], latlong[2*a+1]);
        pushCartesianFromSpherical(geometry, latlong[2*b], latlong[2*b+1]);
        pushCartesianFromSpherical(geometry, latlong[2*c], latlong[2*c+1]);
        colors.push(rgb, rgb, rgb);
    }
    return {geometry, colors};
}

function stereographicProjection(latlong) {
    const degToRad = Math.PI / 180;
    let numPoints = latlong.length / 2;
    let points = [];
    for (let r = 0; r < numPoints; r++) {
        let lat = degToRad * latlong[2*r],
            lon = degToRad * latlong[2*r + 1];
        // See <https://en.wikipedia.org/wiki/Stereographic_projection>
        let x = Math.cos(lat) * Math.cos(lon),
            y = Math.cos(lat) * Math.sin(lon),
            z = Math.sin(lat);
        let X = x / (1-z),
            Y = y / (1-z);
        points.push(X, Y);
    }
    return points;
}

let algorithm = 2;
function generateFibonacciSphere(N, jitter) {
    return [null, generateFibonacciSphere1, generateFibonacciSphere2][algorithm](N, jitter);
}
      
let N = 1000;
let jitter = 0.0;
let rotation = -1.5;

function setAlgorithm(newAlgorithm) { algorithm = newAlgorithm; draw(); }
function setN(newN) { N = newN; draw(); }
function setJitter(newJitter) { jitter = newJitter; draw(); }
function setRotation(newRotation) { rotation = newRotation; draw(); }

function draw() {
    let u_pointsize = 0.1 + 100 / Math.sqrt(N);
    let u_projection = mat4.create();
    mat4.scale(u_projection, u_projection, [1, 1, 0.5, 1]); // avoid clipping
    mat4.rotate(u_projection, u_projection, -rotation, [1, 0.5, 0]);
    
    let a_latlong = generateFibonacciSphere(N, jitter);
    let delaunay = new Delaunator(stereographicProjection(a_latlong));

    /* TODO: Better approach would be to rotate the sphere so that
       a_latlong[N-1] is at the south pole, then run Delaunator on
       everything except that point, then stitch that point back into
       the mesh. However, for now I'm adding a new point into the mesh
       at the south pole and hoping there's nothing there. */
    a_latlong.push(90, 0);
    delaunay = addSouthPoleToMesh(a_latlong.length/2 - 1, delaunay);
    
    let triangleGeometry = generateDelaunayGeometry(a_latlong, delaunay);

    renderTriangles({
        u_projection,
        a_xyz: triangleGeometry.geometry,
        a_rgb: triangleGeometry.colors,
        count: triangleGeometry.geometry.length / 3,
    });
    renderPoints({
        u_projection,
        u_pointsize,
        a_latlong,
        count: a_latlong.length / 2,
    });
}

draw();
