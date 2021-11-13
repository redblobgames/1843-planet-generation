#!/bin/sh
esbuild planet-generation.js --bundle --minify --sourcemap --outfile=build/_bundle.js
