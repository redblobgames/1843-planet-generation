/*
 * From https://www.redblobgames.com/x/1843-planet-generation/
 * Copyright 2018 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

/* global createREGL, mat4 */

const regl = createREGL({
    canvas: "#output",
    extensions: ['OES_element_index_uint']
});

let u_projection = mat4.fromRotation([], 0.5, [1, 1, 0]);

const render = regl({
    frag: `
precision mediump float;
// uniform sampler2D u_texture;
varying vec2 v_uv;
void main() {
   // gl_FragColor = texture2D(u_texture, v_uv);
   gl_FragColor = vec4(v_uv.r, 0.5, v_uv.g, 1);
}
`,

    vert: `
#define M_PI 3.1415926535897932384626433832795
precision highp float;
uniform mat4 u_projection;
attribute vec2 a_latlong;
varying vec2 v_uv;
void main() {
  v_uv = a_latlong / 360.0 + 0.5;
  float theta = a_latlong.r / 180.0 * M_PI;
  float rho = a_latlong.g / 180.0 * M_PI;
  gl_Position = u_projection * vec4(sin(theta) * cos(rho),
                     sin(theta) * sin(rho),
                     cos(theta),
                     1);
}
`,

    uniforms: {
        u_projection,
    },

    count: regl.prop('count'),
    attributes: {
        a_latlong: regl.prop('a_latlong'),
    },
});


let count = 0;
let a_latlong = [];
for (let lat = -180; lat < 180; lat += 20) {
    for (let long = -180; long < 180; long += 20) {
        a_latlong.push([
            lat, long, lat+20, long, lat, long+20,
            lat+20, long, lat, long+20, lat+20, long+20,
        ]);
        count += 2;
    }
}

render({
    a_latlong,
    count,
});
