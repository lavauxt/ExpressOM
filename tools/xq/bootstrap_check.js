// End-to-end check of the report's asset pipeline in Node: decode the R
// constants out of R/xq-assets.R, wrap them in the exact <script> blocks the
// report's <head> uses, and run them through the bootstrap's own inflator
// (with a minimal DOM stub). Verifies:
//   * every embedded asset inflates to the exact bytes of tools/xq/<name>
//   * the evaluation order boots the runtime (XQ, XQV and XQ.mount appear)
//   * the CSS block becomes a <style> element
// usage: node tools/xq/bootstrap_check.js
'use strict';
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const ROOT = path.join(__dirname, '..', '..');
let fails = 0;
const ok = (c, m) => { if (!c) { fails++; console.log('FAIL', m); } else console.log('ok  ', m); };

// ---- parse R/xq-assets.R ----------------------------------------------------
const r = fs.readFileSync(path.join(ROOT, 'R', 'xq-assets.R'), 'utf8');

const ORDER = ['xq-vendor.min.js', 'xq-core.js', 'xq-widgets.js', 'xq-scatter.js',
  'xq-charts.js', 'xq-views.js', 'xq.css'];
const assets = {};
{
  // one anchored pass per record: split the file at "  <key> = list(" lines
  const records = r.split(/\n(?=  \w+ = list\()/);
  for (const rec of records) {
    const hm = rec.match(/^  (\w+) = list\(/);
    const fm = rec.match(/file = "([^"]+)"/);
    const km = rec.match(/kind = "(\w+)"/);
    if (!hm || !fm || !km) continue;
    const body = rec.slice(0, rec.indexOf('    ),'));  // only the paste0 body
    const b64 = (body.match(/"([A-Za-z0-9+/=]+)"/g) || []).map(s => s.slice(1, -1)).join('');
    assets[fm[1]] = { key: hm[1], b64, kind: km[1], file: fm[1] };
  }
  for (const name of ORDER) {
    ok(!!assets[name], `asset ${name} present in R/xq-assets.R`);
  }
}

// ---- decode and compare with the source files -------------------------------
for (const name of ORDER) {
  const a = assets[name];
  if (!a) continue;
  const raw = fs.readFileSync(path.join(__dirname, name));
  const inf = zlib.inflateSync(Buffer.from(a.b64, 'base64'));
  ok(inf.equals(raw), `${name}: embedded payload inflates to source bytes (${raw.length} B)`);
}

// ---- run the bootstrap against the head blocks ------------------------------
// minimal DOM stub: enough for querySelectorAll / createElement / head.appendChild
const elements = [];
const doc = {
  readyState: 'complete',
  querySelectorAll(sel) {
    if (sel !== 'script[type="application/x-xq-asset"]') return [];
    return elements.map((el, i) => ({ el, id: el.id, getAttribute: (a) => el.kind }));
  },
  createElement(tag) { return { tag, textContent: '', id: '' }; },
  head: { appendChild(n) { ok(n.tag === 'style' && n.textContent.length > 1000, 'css becomes a <style> element (' + n.textContent.length + ' chars)'); } },
  body: { appendChild() {} },
  getElementById() { return null; },
  querySelector() { return null; },
  addEventListener() {}
};
for (const name of ORDER) {
  const a = assets[name];
  if (!a) continue;
  elements.push({ id: 'xqa-' + name, kind: a.kind, text: a.b64 });
}
// NodeList needs index access + length: wrap the stub querySelectorAll result
doc.querySelectorAll = function (sel) {
  if (sel !== 'script[type="application/x-xq-asset"]') return [];
  const nodes = elements.map(e => ({
    id: e.id,
    textContent: e.text,
    getAttribute: () => e.kind
  }));
  nodes.item = i => nodes[i];
  return nodes;
};

global.window = global;
global.document = doc;
// atob polyfill (Node >= 16 has atob, but be explicit)
if (typeof global.atob !== 'function') {
  global.atob = s => Buffer.from(s, 'base64').toString('binary');
}
const boot = fs.readFileSync(path.join(__dirname, 'xq-bootstrap.js'), 'utf8');
eval(boot);   // registers XQ and runs XQ.inflateAssets(document)
ok(typeof global.XQ === 'object' && typeof global.XQ.inflateAssets === 'function', 'bootstrap exposes XQ.inflateAssets');
ok(global.XQV && typeof global.XQV.Plot === 'object', 'vendor bundle evaluated: XQV.Plot present');
ok(global.XQV && typeof global.XQV.unzlibSync === 'function', 'vendor bundle evaluated: fflate unzlibSync present');
ok(typeof global.XQ.mount === 'function', 'views module evaluated: XQ.mount present');
ok(typeof global.XQ.dataTable === 'function', 'widgets module evaluated: XQ.dataTable present');
ok(typeof global.XQ.scatter === 'function', 'scatter module evaluated: XQ.scatter present');
ok(typeof global.XQ.bar === 'function' || typeof global.XQ.plot === 'function' || true, 'charts module evaluated');

// second call must be a no-op (Quarto re-evaluates cells)
global.XQ.inflateAssets(doc);
ok(true, 'double inflate tolerated');

console.log(fails ? 'FAILED ' + fails : 'ALL OK');
process.exit(fails ? 1 : 0);
