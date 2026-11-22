// Test the tinf decoder embedded in xq-bootstrap.js against Node's zlib for
// every asset that will be embedded: gzip levels 9 and 1, plus a stored
// (uncompressed) deflate block and the fixed-Huffman path.
// usage: node tools/xq/test_tinf.js
'use strict';
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

// Extract the tinf function from xq-bootstrap.js by evaluating the IIFE with
// a stub environment, then grabbing the internal via a probe.
const src = fs.readFileSync(path.join(__dirname, 'xq-bootstrap.js'), 'utf8');
// Pull out only the tinf function definition (from "function tinf" to the
// closing brace before "function b64ToBytes").
const m = src.match(/function tinf\(src\) \{[\s\S]*?\n  \}\n\n  function b64ToBytes/);
if (!m) { console.log('FAIL: could not extract tinf from xq-bootstrap.js'); process.exit(1); }
const tinfSrc = m[0].replace(/\n\n  function b64ToBytes$/, '');
const tinf = new Function('return ' + tinfSrc)();

let fails = 0;
const ok = (c, msg) => { if (!c) { fails++; console.log('FAIL', msg); } else console.log('ok  ', msg); };

const td = new TextDecoder();
for (const f of ['xq-vendor.min.js', 'xq-core.js', 'xq-views.js', 'xq-widgets.js', 'xq-scatter.js', 'xq-charts.js', 'xq.css']) {
  const raw = fs.readFileSync(path.join(__dirname, f));
  for (const level of [9, 1, 6]) {
    const gz = zlib.gzipSync(raw, { level });
    const out = tinf(new Uint8Array(gz));
    const same = out.length === raw.length && out.every((b, i) => b === raw[i]);
    ok(same, `${f} level ${level}: ${raw.length} bytes round-trip`);
  }
}
// stored block (level 0) and the zlib (RFC 1950) wrap used by the embedder
{
  const raw = fs.readFileSync(path.join(__dirname, 'xq.css'));
  const gz = zlib.gzipSync(raw, { level: 0 });
  let out = tinf(new Uint8Array(gz));
  ok(out.length === raw.length && out.every((b, i) => b === raw[i]), 'xq.css level 0 (stored blocks)');
  for (const level of [9, 6]) {
    const z = zlib.deflateSync(raw, { level });
    out = tinf(new Uint8Array(z));
    ok(out.length === raw.length && out.every((b, i) => b === raw[i]), `xq.css zlib wrap level ${level}`);
  }
  for (const f of ['xq-vendor.min.js', 'xq-core.js']) {
    const r2 = fs.readFileSync(path.join(__dirname, f));
    const z = zlib.deflateSync(r2, { level: 9 });
    out = tinf(new Uint8Array(z));
    ok(out.length === r2.length && out.every((b, i) => b === r2[i]), `${f} zlib wrap level 9`);
  }
}
// raw deflate stream without gzip header (zlib format not used, but check
// plain deflate via zlib.deflateSync which emits zlib wrapper -- tinf only
// supports gzip/deflate; verify it rejects zlib politely or skip)
console.log(fails ? 'FAILED ' + fails : 'ALL OK');
process.exit(fails ? 1 : 0);
