// Vendored libraries for the ExpressOM Quarto dashboard engine.
// Bundled with esbuild into tools/xq/xq-vendor.min.js. Rebuild:
//   cd tools/vendor && npm install && npx esbuild vendor-entry.js \
//     --bundle --minify --format=iife -o ../xq/xq-vendor.min.js
// then regenerate the embedded copy: python3 tools/embed_xq_assets.py
import * as Plot from "@observablehq/plot";
import * as Inputs from "@observablehq/inputs";
import * as htl from "htl";
import { scaleLinear, scaleLog, scaleSqrt, scaleOrdinal, scaleSequential } from "d3-scale";
import { axisBottom, axisLeft } from "d3-axis";
import { select } from "d3-selection";
import { extent, bisector, ticks, quantileSorted, max, min, mean } from "d3-array";
import { format } from "d3-format";
import { unzlibSync } from "fflate";
import createScatterplot from "regl-scatterplot";
const d3 = { scaleLinear, scaleLog, scaleSqrt, scaleOrdinal, scaleSequential, axisBottom, axisLeft, select, extent, bisector, ticks, quantileSorted, max, min, mean, format };
globalThis.XQV = { Plot, Inputs, htl, d3, createScatterplot, unzlibSync };
