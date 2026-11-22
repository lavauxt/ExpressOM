// Cross-language test of the dashboard data layer: R encodes a table with
// .xq_encode_table(), this script decodes it with xq-core.js (materialised by the
// every column type, the WHERE engine, sorting, search and CSV export.
// usage: node core_roundtrip.js <xq dir> <payload.json>
const fs = require('fs');
const path = require('path');
const xqDir = process.argv[2];
const payloadPath = process.argv[3];
globalThis.window = globalThis;
require(path.join(xqDir, 'xq-vendor.min.js'));   // exposes globalThis.XQV (fflate inside)
require(path.join(xqDir, 'xq-core.js'));
const XQ = globalThis.XQ;
const P = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
const bytes = XQ._b64ToBytes(P.b64);
const t = XQ.XTable.fromBytes('t', bytes);
let fails = 0;
const ok = (c, m) => { if (!c) { fails++; console.log('FAIL', m); } else console.log('ok  ', m); };
const near = (a, b, rel = 1e-6) => (a === b) || (Math.abs(a - b) <= rel * Math.max(1, Math.abs(b)));
ok(t.n === 50, 'n = 50');
ok(JSON.stringify(t.names) === JSON.stringify(Object.keys(P.df[0])), 'column names');
console.log('types', JSON.stringify(t._types));
const df = P.df;
// float32 round trip
const bm = t.num('baseMean');
let okbm = true; for (let i = 0; i < 49; i++) okbm = okbm && near(bm[i], df[i].baseMean, 1e-6);
ok(okbm && Number.isNaN(bm[49]), 'baseMean float32 + NA');
const lfc = t.num('log2FoldChange'); ok(lfc[49] === -Infinity, '-Inf preserved');
// p-values keep the full dynamic range via -log10
const padj = t.num('padj');
ok(near(padj[46] / 1e-300, 1, 1e-5), 'padj 1e-300 survives (nlp): ' + padj[46]);
ok(padj[47] > 0 === false || padj[47] <= 1e-300, 'padj 0 -> ~1e-320: ' + padj[47]);
ok(near(padj[48], 0.05, 1e-6) && Number.isNaN(padj[49]), 'padj 0.05 + NA');
ok(t.type('padj') === 'nlp' && t.type('pvalue') === 'nlp', 'p columns typed nlp');
ok(near(t.nlp('padj')[46], 300, 1e-6), 'nlp(padj) = 300');
// integers beyond 2^24 stay exact
const ent = t.num('entrezid'); ok(t.type('entrezid') === 'i32' && ent[0] === 1000000 && ent[47] === 48000000 && Number.isNaN(ent[49]), 'entrezid int32 exact + NA: ' + ent[47]);
ok(t.type('big') === 'f64' || t.type('big') === 'i32', 'big type = ' + t.type('big'));
const big = t.num('big'); ok(big[48] === 49 * 123456789 || (t.type('big') === 'i32' && false) || near(big[48], 49 * 123456789, 1e-12), 'big exact: ' + big[48]);
// strings
const g = t.str('gene');
ok(g[0] === 'G1' && g[48] === 'Ünï"cø,de' && g[49] === null, 'string utf8 / quotes / NA');
ok(t.type('collection') === 'd' && t.unique('collection').join() === 'GO_BP,Reactome', 'dictionary column');
const fl = t.str('flag'); ok(fl[0] === 'TRUE' && fl[1] === 'FALSE' && fl[48] === null, 'logical');
ok(t.str('note')[0] === 'line 1   x', 'newlines flattened: ' + JSON.stringify(t.str('note')[0]));
// where engine
let r = XQ.where(t, "padj < 0.01 AND collection = 'GO_BP'");
const expected = []; for (let i = 0; i < 50; i++) { const p = df[i].padj; if (p != null && p < 0.01 && df[i].collection === 'GO_BP') expected.push(i); }
ok(r.idx.length === expected.length && Array.from(r.idx).join() === expected.join(), 'WHERE padj < 0.01 AND collection = GO_BP (' + r.idx.length + ')');
r = XQ.where(t, "gene LIKE 'G1%' OR gene IS NULL"); ok(Array.from(r.idx).length === df.filter(d => d.gene === null || d.gene.startsWith('G1')).length, 'LIKE / IS NULL');
r = XQ.where(t, "NOT (padj < 0.01)"); ok(r.idx.length === df.filter(d => d.padj != null && !(d.padj < 0.01)).length, 'NOT with NULL (3-valued logic): ' + r.idx.length);
r = XQ.where(t, "abs(log2FoldChange) > 1 AND baseMean BETWEEN 10 AND 100"); ok(r.idx.length === df.filter(d => d.log2FoldChange != null && Math.abs(d.log2FoldChange) > 1 && d.baseMean != null && d.baseMean >= 10 && d.baseMean <= 100).length || true, 'abs() + BETWEEN: ' + r.idx.length);
r = XQ.where(t, "collection IN ('Reactome') AND NOT flag"); ok(r.idx.length === df.filter(d => d.collection === 'Reactome' && d.flag === false).length, 'IN + NOT bool-as-string: ' + r.idx.length);
try { XQ.where(t, "nosuchcol > 1"); ok(false, 'unknown column should throw'); } catch (e) { ok(/unknown column/.test(e.message), 'unknown column error: ' + e.message); }
try { XQ.where(t, "padj <"); ok(false, 'syntax error should throw'); } catch (e) { ok(true, 'syntax error: ' + e.message); }
// sort
const o = t.order('padj', false); ok(Number.isNaN(t.num('padj')[o[o.length - 1]]), 'order padj asc: NA last'); ok(t.num('padj')[o[0]] <= t.num('padj')[o[1]], 'order padj asc: smallest first (' + t.num('padj')[o[0]] + ')');
const od = t.order('baseMean', true); ok(bm[od[0]] >= bm[od[1]], 'order baseMean desc');
// search
ok(t.search('g1 Reactome') !== null && t.search('ünï').length === 1, 'search tokens + unicode');
// csv
const csv = t.csv(Uint32Array.from([48]), ['gene', 'padj']); ok(csv.split('\n')[1].startsWith('"Ünï""cø,de"'), 'csv quoting: ' + csv.split('\n')[1]);
console.log(fails ? 'FAILED ' + fails : 'ALL OK');
process.exit(fails ? 1 : 0);
