/* xq-views.js -- one function per dashboard card.
 *
 * R writes a small JSON spec per card (<script type="application/json"
 * id="xqj-v_<id>">) and the Quarto cell is just  XQ.mount("<id>", invalidation).
 * Nothing from R is ever interpolated into JavaScript source, and a card only
 * builds when it is first visible. Controls are Observable Inputs; large point
 * clouds are WebGL (regl-scatterplot); the filter boxes run on typed arrays
 * with DuckDB-Wasm loaded on demand.
 */
(function (g) {
  'use strict';
  const XQ = g.XQ;
  const doc = g.document;
  const el = XQ.el, F = XQ.fmt, C = XQ.charts;
  const V = (XQ.views = {});
  const Inputs = function () { return g.XQV.Inputs; };

  // ------------------------------------------------------------------ helpers
  const trunc = C.trunc;
  const grid = function () { const d = el('div', 'xq-grid2'); Array.prototype.forEach.call(arguments, function (n) { if (n) d.append(n); }); return d; };
  const stack = function () { const d = el('div', 'xq-stack'); Array.prototype.forEach.call(arguments, function (n) { if (n) d.append(n); }); return d; };
  const h4 = function (t) { return el('h4', 'xq-h4', t); };
  const fig = function (title, node, o) { return XQ.frame(node, Object.assign({ title: title }, o || {})); };
  const controls = function () { const d = el('div', 'xq-controls'); Array.prototype.forEach.call(arguments, function (n) { if (n) d.append(n); }); return d; };
  const num = function (x, d) { return F.num(x, d == null ? 3 : d); };

  function disposeTree(node) {
    if (!node || !node.querySelectorAll) return;
    const all = [node].concat(Array.from(node.querySelectorAll('*')));
    all.forEach(function (n) {
      try { if (n._xqDispose) n._xqDispose(); } catch (e) { /* noop */ }
    });
  }

  /** The only function the Quarto cells call. */
  XQ.mount = function (id, invalidation) {
    const spec = XQ.json('v_' + id);
    if (!spec) return el('div', 'xq-err', 'This panel has no data (view "' + id + '" is missing).');
    const fn = V[spec.view];
    if (!fn) return el('div', 'xq-err', 'Unknown view type: ' + spec.view);
    const host = XQ.lazy(function () { return fn(spec); }, { minHeight: spec.minHeight || 180, label: spec.label || 'Preparing this panel…' });
    if (invalidation && invalidation.then) invalidation.then(function () { disposeTree(host); });
    // The "How to read this" note of the card travels with its view, so it is
    // always rendered ABOVE the chart (Quarto may reorder sibling blocks).
    const hs = doc.getElementById('xqh-' + id);
    if (!hs) return host;
    const wrap = el('div', 'xq-mounted');
    const slot = el('div', 'xq-help-slot');
    slot.innerHTML = hs.textContent;
    wrap.append(slot, host);
    return wrap;
  };

  function groupIdx(tbl, col) {
    const s = tbl.str(col), m = new Map();
    for (let i = 0; i < s.length; i++) { const k = s[i]; if (!m.has(k)) m.set(k, []); m.get(k).push(i); }
    const out = new Map();
    m.forEach(function (v, k) { out.set(k, Uint32Array.from(v)); });
    return out;
  }
  const firstOf = function (tbl, cands) { return tbl.first(cands); };
  function labelsOf(tbl, col) { const s = col ? tbl.str(col) : null; return function (i) { return s && s[i] != null ? s[i] : '#' + (i + 1); }; }
  function select(options, o) { return Inputs().select(options, o); }
  function onInput(node, fn, ms) { node.addEventListener('input', ms ? XQ.debounce(fn, ms) : fn); }
  const padjChoices = function (cur) {
    const base = [0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.000001];
    if (cur != null && base.indexOf(cur) < 0) { base.push(cur); base.sort(function (a, b) { return b - a; }); }
    return base;
  };

  // --------------------------------------------------------------- tables
  /** Filterable, sortable, exportable table (+ optional SQL console). */
  V.table = function (s) {
    const tbl = XQ.table(s.table);
    const root = el('div', 'xq-view');
    if (!tbl.n) { root.append(el('div', 'xq-empty', s.empty || 'This table is empty.')); return root; }
    if (s.intro) root.append(el('p', 'xq-intro', s.intro));
    const dt = XQ.dataTable(tbl, { columns: s.columns, sort: s.sort, desc: !!s.desc, pageSize: s.pageSize || 15, filename: (s.filename || s.table) + '.csv', labels: s.labels, maxChars: s.maxChars });
    if (s.filter !== false) {
      const fb = XQ.filterBox(s.sqlName || s.table, tbl, { examples: s.examples });
      fb.addEventListener('input', function () { dt.update(fb.value); });
      root.append(fb);
    }
    root.append(dt);
    if (s.sql !== false && s.filter !== false) {
      const det = el('details', 'xq-sqlbox');
      det.append(el('summary', null, 'Free SQL (DuckDB-Wasm): joins, GROUP BY, window functions'));
      const t = {}; t[s.sqlName || s.table] = tbl;
      let built = false;
      det.addEventListener('toggle', function () { if (det.open && !built) { built = true; det.append(XQ.sqlConsole(t, { sql: s.sqlExample })); } });
      root.append(det);
    }
    return root;
  };

  // ------------------------------------------------------------ DEG counts
  V.degCounts = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    return C.mount(function (w) {
      const lfcCut = tbl.has('log2FC_cutoff') ? tbl.num('log2FC_cutoff')[0] : null;
      const upC = tbl.first(['n_up_lfc', 'n_up']), dnC = tbl.first(['n_down_lfc', 'n_down']), totC = tbl.first(['n_sig_lfc', 'n_sig']);
      const pc = tbl.num('padj_cutoff');
      const rows = [];
      const noDir = s.lrt || !upC || !dnC;
      for (let i = 0; i < tbl.n; i++) {
        const lab = 'padj < ' + pc[i];
        if (noDir) { if (totC) rows.push({ lab: lab, kind: 'Significant', n: tbl.num(totC)[i], color: cfg.colors.up }); }
        else {
          rows.push({ lab: lab, kind: 'Up', n: tbl.num(upC)[i], color: cfg.colors.up });
          rows.push({ lab: lab, kind: 'Down', n: tbl.num(dnC)[i], color: cfg.colors.down });
        }
      }
      if (!rows.length) return C.empty('No DEG count table.');
      const P = g.XQV.Plot;
      const doms = Array.from(new Set(rows.map(function (r) { return r.kind; })));
      return P.plot({
        width: w, height: 300, marginLeft: 56, x: { label: null, axis: null, domain: doms }, fx: { label: null, domain: Array.from(new Set(rows.map(function (r) { return r.lab; }))) },
        y: { label: 'Number of genes', grid: true }, color: { domain: doms, range: doms.map(function (k) { return rows.find(function (r) { return r.kind === k; }).color; }), legend: true },
        marks: [P.barY(rows, { fx: 'lab', x: 'kind', y: 'n', fill: 'kind', tip: true, title: function (r) { return r.kind + ' (' + r.lab + (lfcCut != null ? ', |log2FC| > ' + lfcCut : '') + '): ' + F.int(r.n); } }), P.ruleY([0])]
      });
    }, { minHeight: 240 });
  };

  // --------------------------------------------- volcano / MA explorer (WebGL)
  V.deExplorer = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const root = el('div', 'xq-view');
    const sc = C.deScatter(tbl, { idCol: s.idCol, lfcCol: s.lfcCol, pCol: s.pCol, meanCol: s.meanCol, colors: s.colors || cfg.colors, tipCols: s.tipCols, height: s.height || 520, filename: s.filename || 'expressom_volcano' });
    if (!sc.api || !sc.api.set || sc.classList.contains('xq-empty')) { root.append(sc); return root; }
    const hasMean = !!tbl.resolve(s.meanCol || 'baseMean');
    const I = Inputs();
    const viewSel = hasMean ? I.radio(['Volcano', 'MA plot'], { value: 'Volcano', label: 'View' }) : null;
    const padj = select(padjChoices(s.padj), { value: s.padj, label: 'adj. p <' });
    const lfc = I.range([0, 5], { step: 0.05, value: s.lfc, label: '|log2FC| ≥' });
    const topN = I.range([0, 40], { step: 1, value: s.topN == null ? 10 : s.topN, label: 'Top labels' });
    const find = I.text({ placeholder: 'Find gene(s), e.g. TP53 MYC', label: 'Find' });
    const tiles = el('div', 'xq-tiles');
    const tableTitle = h4('');
    const tableHost = el('div');
    let dt = null;
    const cols = s.columns;
    function showTable(idx, title) {
      tableTitle.textContent = title;
      if (!dt) {
        dt = XQ.dataTable(tbl, { columns: cols, idx: idx, sort: s.sort || tbl.resolve(s.pCol || 'padj'), desc: false, pageSize: 10, filename: (s.filename || 'expressom_genes') + '_subset.csv' });
        tableHost.replaceChildren(dt);
      } else dt.update(idx);
    }
    let selecting = false;
    function refresh() {
      const c = sc.api.update({ view: viewSel && viewSel.value === 'MA plot' ? 'ma' : 'volcano', padj: +padj.value, lfc: +lfc.value, topN: +topN.value, gene: find.value });
      tiles.replaceChildren(XQ.stats([
        { label: s.lrt ? 'Significant' : 'Up', value: s.lrt ? F.int(c.sig) : F.int(c.up), cls: 'xq-up', sub: 'adj. p < ' + padj.value + (+lfc.value ? ' & |log2FC| ≥ ' + lfc.value : '') },
        s.lrt ? null : { label: 'Down', value: F.int(c.down), cls: 'xq-down' },
        { label: 'Not significant', value: F.int(c.ns) },
        { label: 'Features plotted', value: F.int(c.n) },
        find.value.trim() ? { label: 'Found', value: F.int(c.found), sub: 'ringed on the plot' } : null
      ].filter(Boolean)));
      if (!selecting) showTable(c.sigIdx, 'Significant features at these cut-offs (' + F.int(c.sig) + ') — lasso points on the plot to list them instead');
    }
    sc.addEventListener('input', function () {
      const v = sc.value;
      if (v && v.length) { selecting = true; showTable(v, 'Selected on the plot (' + F.int(v.length) + ')'); }
      else { selecting = false; refresh(); }
    });
    [viewSel, padj, lfc, topN, find].forEach(function (n) { if (n) onInput(n, refresh, 80); });
    root.append(controls(viewSel, padj, lfc, topN, find), tiles, fig(s.title || 'Volcano plot', sc, { download: false }), tableTitle, tableHost);
    refresh();
    return root;
  };

  // ----------------------------------------------------------- p diagnostics
  V.pDiag = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const nodes = [];
    (s.hists || []).forEach(function (h) {
      const c = tbl.resolve(h.col);
      if (c) nodes.push(fig(h.title, C.hist(tbl.num(c), { color: h.color }), { filename: h.col + '_hist' }));
    });
    if (s.qq) {
      const c = tbl.resolve(s.qq.col);
      if (c) nodes.push(fig(s.qq.title || 'p-value Q–Q plot', C.qq(tbl, c), { download: false }));
    }
    if (s.ranked) {
      const lfcC = tbl.resolve(s.ranked.lfcCol || 'log2FoldChange');
      const idC = tbl.resolve(s.ranked.idCol || 'gene');
      if (lfcC) {
        const l = tbl.num(lfcC), nm = labelsOf(tbl, idC);
        const ord = [];
        for (let i = 0; i < l.length; i++) if (l[i] === l[i] && isFinite(l[i])) ord.push(i);
        ord.sort(function (a, b) { return Math.abs(l[b]) - Math.abs(l[a]); });
        const top = ord.slice(0, s.ranked.n || 25).sort(function (a, b) { return l[a] - l[b]; });
        const data = top.map(function (i) { return { label: String(nm(i)), value: l[i], color: l[i] > 0 ? cfg.colors.up : cfg.colors.down }; });
        nodes.push(fig('Top ' + data.length + ' features by |log2 fold change|', C.bars({ data: data, xlabel: 'log2 fold change' }), { filename: 'top_effects' }));
      }
    }
    const root = el('div', 'xq-view');
    if (!nodes.length) { root.append(C.empty('No p-value columns in this table.')); return root; }
    // the ranked-effects bar chart is wide: it gets its own row below the grid
    const wide = s.ranked && nodes.length > 1 ? nodes.pop() : null;
    for (let k = 0; k < nodes.length; k += 2) root.append(grid(nodes[k], nodes[k + 1]));
    if (wide) root.append(wide);
    return root;
  };

  // -------------------------------------------------------- hierarchical order
  /** Average-linkage hierarchical clustering; returns the leaf order. */
  function hclustOrder(vecs) {
    const n = vecs.length;
    if (n < 3) return Array.from({ length: n }, function (_, i) { return i; });
    const D = [];
    for (let i = 0; i < n; i++) {
      D.push(new Float64Array(n));
    }
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      let d = 0; const a = vecs[i], b = vecs[j];
      for (let k = 0; k < a.length; k++) { const t = a[k] - b[k]; d += t * t; }
      D[i][j] = D[j][i] = Math.sqrt(d);
    }
    let cl = [];
    for (let i = 0; i < n; i++) cl.push({ id: i, ord: [i], size: 1 });
    // cluster-distance matrix between current clusters (indexed by position)
    let M = [];
    for (let i = 0; i < n; i++) { M.push(Array.from(D[i])); }
    while (cl.length > 1) {
      let bi = 0, bj = 1, bd = Infinity;
      for (let i = 0; i < cl.length; i++) for (let j = i + 1; j < cl.length; j++) if (M[i][j] < bd) { bd = M[i][j]; bi = i; bj = j; }
      const a = cl[bi], b = cl[bj];
      const merged = { ord: a.ord.concat(b.ord), size: a.size + b.size };
      const newRow = [];
      for (let k = 0; k < cl.length; k++) {
        if (k === bi || k === bj) continue;
        newRow.push((a.size * M[bi][k] + b.size * M[bj][k]) / (a.size + b.size));
      }
      const keep = [];
      for (let k = 0; k < cl.length; k++) if (k !== bi && k !== bj) keep.push(k);
      const nM = keep.map(function (r) { return keep.map(function (c) { return M[r][c]; }); });
      nM.forEach(function (row, r) { row.push(newRow[r]); });
      nM.push(newRow.concat([0]));
      M = nM;
      cl = keep.map(function (k) { return cl[k]; }).concat([merged]);
    }
    return cl[0].ord;
  }
  XQ.hclustOrder = hclustOrder;

  // ----------------------------------------------------------------- heatmap
  V.heatmap = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const I = Inputs();
    const samples = s.samples;
    const idC = tbl.resolve(s.idCol || 'gene');
    const ids = tbl.str(idC);
    const mats = samples.map(function (nm) { return tbl.num(nm); });
    const nMax = Math.min(tbl.n, s.nMax || 100);
    const choices = (s.choices || [10, 20, 30, 50, 100]).filter(function (n) { return n <= nMax; });
    if (choices.indexOf(Math.min(s.n || 20, nMax)) < 0) choices.push(Math.min(s.n || 20, nMax));
    choices.sort(function (a, b) { return a - b; });
    const nSel = select(choices, { value: Math.min(s.n || 20, nMax), label: 'Top genes' });
    const orderSel = I.radio(['Cluster samples', 'Group samples'], { value: s.groups ? 'Group samples' : 'Cluster samples', label: 'Columns' });
    const host = el('div');
    const groups = s.groups || {};
    const levels = Array.from(new Set(samples.map(function (x) { return groups[x]; }).filter(Boolean)));
    const gcol = C.groupColors(levels, { level: cfg.level, levelColor: cfg.pcaColors && cfg.pcaColors.level, base: cfg.base, baseColor: cfg.pcaColors && cfg.pcaColors.base });
    function draw() {
      const n = +nSel.value;
      const z = [];
      for (let r = 0; r < n; r++) {
        const row = samples.map(function (_, k) { const v = mats[k][r]; return v === v ? Math.log2(v + 1) : NaN; });
        const fin = row.filter(function (v) { return v === v; });
        const mu = fin.reduce(function (a, b) { return a + b; }, 0) / (fin.length || 1);
        const sd = Math.sqrt(fin.reduce(function (a, b) { return a + (b - mu) * (b - mu); }, 0) / Math.max(1, fin.length - 1));
        z.push(row.map(function (v) { return v === v && sd > 0 ? (v - mu) / sd : 0; }));
      }
      const rowOrd = n >= 3 ? hclustOrder(z) : z.map(function (_, i) { return i; });
      let colOrd;
      if (orderSel.value === 'Cluster samples' && samples.length >= 3) colOrd = hclustOrder(samples.map(function (_, k) { return z.map(function (row) { return row[k]; }); }));
      else if (s.groups) colOrd = samples.map(function (_, k) { return k; }).sort(function (a, b) { const ga = groups[samples[a]] || '', gb = groups[samples[b]] || ''; return ga < gb ? -1 : ga > gb ? 1 : a - b; });
      else colOrd = samples.map(function (_, k) { return k; });
      const colLabel = colOrd.map(function (k) { return samples[k]; });
      const vals = [];
      rowOrd.forEach(function (r) { colOrd.forEach(function (k) { vals.push(z[r][k]); }); });
      let zmax = 0; vals.forEach(function (v) { if (Math.abs(v) > zmax) zmax = Math.abs(v); });
      const hm = C.heatmap({
        rows: rowOrd.map(function (r) { return String(ids[r]); }), cols: colLabel, values: vals, label: 'row z-score',
        scheme: 'rdbu', pivot: 0, domain: [-(zmax || 1), zmax || 1],
        annot: s.groups ? colOrd.map(function (k) { return gcol[groups[samples[k]]] || '#999'; }) : null,
        annotLabels: s.groups ? colOrd.map(function (k) { return groups[samples[k]]; }) : null, annotName: s.groupName || 'group'
      });
      host.replaceChildren(fig('Top ' + n + ' differentially expressed genes: row z-score of log2(normalised count + 1)', hm, { filename: 'top_genes_heatmap' }));
      if (levels.length) {
        const lg = el('div', 'xq-legend');
        levels.forEach(function (lv) { const sp = el('span', 'xq-lg2'); const d = el('i'); d.style.background = gcol[lv]; sp.append(d, doc.createTextNode(lv)); lg.append(sp); });
        host.append(lg);
      }
    }
    [nSel, orderSel].forEach(function (n) { onInput(n, draw, 40); });
    const root = el('div', 'xq-view');
    root.append(controls(nSel, orderSel), host);
    draw();
    return root;
  };

  // --------------------------------------------------------------------- PCA
  V.pca = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const pcs = s.pcs;
    const root = el('div', 'xq-view');
    const xs = select(pcs, { value: pcs[0], label: 'X axis' });
    const ys = select(pcs, { value: pcs[1] || pcs[0], label: 'Y axis' });
    const host = el('div');
    const labC = tbl.resolve(s.labelCol || 'sample_label');
    const grpC = s.groupCol ? tbl.resolve(s.groupCol) : null;
    const labs = labC ? tbl.str(labC) : null;
    const grp = grpC ? tbl.str(grpC) : null;
    const levels = grp ? Array.from(new Set(grp.filter(function (x) { return x != null; }))) : [];
    const gcol = C.groupColors(levels, { level: cfg.level, levelColor: cfg.pcaColors && cfg.pcaColors.level, base: cfg.base, baseColor: cfg.pcaColors && cfg.pcaColors.base });
    const axisName = function (pc) { return pc + (s.variance && s.variance[pc] != null ? ' (' + s.variance[pc] + '% of variance)' : ''); };
    function draw() {
      const X = tbl.num(xs.value), Y = tbl.num(ys.value);
      const pts = [];
      for (let i = 0; i < tbl.n; i++) pts.push({ x: X[i], y: Y[i], label: labs ? String(labs[i]) : '#' + (i + 1), group: grp ? grp[i] : '', color: grp ? gcol[grp[i]] || '#888' : C.PALETTE[0] });
      const pl = C.points({ pts: pts, xlabel: axisName(xs.value), ylabel: axisName(ys.value), legend: levels.map(function (lv) { return { label: (grpC || 'group') + ': ' + lv, color: gcol[lv] }; }) });
      host.replaceChildren(fig(s.title || 'PCA', pl, { filename: 'pca_' + xs.value + '_' + ys.value }));
    }
    [xs, ys].forEach(function (n) { onInput(n, draw); });
    let scree = null;
    if (s.variance) {
      const names = Object.keys(s.variance).slice(0, 10);
      const data = names.map(function (k) { return { label: k, value: s.variance[k], color: '#0072B2' }; });
      scree = fig('Variance explained per component', C.bars({ data: data, xlabel: '% of variance' }), { filename: 'scree' });
    }
    root.append(controls(xs, ys), grid(host, scree), h4('Sample scores'), XQ.dataTable(tbl, { pageSize: 10, filename: (s.filename || 'pca_scores') + '.csv' }));
    draw();
    return root;
  };

  // ------------------------------------------------------ sample similarity
  V.similarity = function (s) {
    const cfg = XQ.cfg;
    const groups = s.groups || {};
    const levels = Array.from(new Set(Object.keys(groups).map(function (k) { return groups[k]; })));
    const gcol = C.groupColors(levels, { level: cfg.level, levelColor: cfg.pcaColors && cfg.pcaColors.level, base: cfg.base, baseColor: cfg.pcaColors && cfg.pcaColors.base });
    const node = [];
    if (s.cor) {
      let mn = 1; s.cor.values.forEach(function (v) { if (v < mn) mn = v; });
      const lo = Math.max(-1, Math.floor(Math.min(mn, 0.95) * 20) / 20);
      node.push(fig('Sample correlation (Pearson, log2 counts)', C.heatmap({ rows: s.cor.labels, cols: s.cor.labels, values: s.cor.values, scheme: 'ylgnbu', domain: [lo, 1], label: 'r' }), { filename: 'sample_correlation' }));
    }
    if (s.dist) node.push(fig('Sample-to-sample Euclidean distance', C.heatmap({ rows: s.dist.labels, cols: s.dist.labels, values: s.dist.values, scheme: 'viridis', label: 'distance' }), { filename: 'sample_distance' }));
    if (s.mds) {
      const pts = s.mds.map(function (p) { return { x: p.x, y: p.y, label: p.label, group: groups[p.label] || '', color: groups[p.label] ? gcol[groups[p.label]] : C.PALETTE[2] }; });
      node.push(fig('Multidimensional scaling of samples', C.points({ pts: pts, xlabel: 'MDS1', ylabel: 'MDS2', legend: levels.map(function (lv) { return { label: lv, color: gcol[lv] }; }) }), { filename: 'sample_mds' }));
    }
    if (s.dendro && s.dendro.segs && s.dendro.segs.length) node.push(fig('Hierarchical clustering of samples (complete linkage)', C.dendrogram({ segs: s.dendro.segs, leaves: s.dendro.leaves, ylabel: 'Euclidean distance' }), { filename: 'sample_dendrogram' }));
    const root = el('div', 'xq-view');
    if (!node.length) root.append(C.empty('At least three samples with expression values are needed.'));
    else root.append(grid.apply(null, node.slice(0, 2)), grid.apply(null, node.slice(2)));
    return root;
  };

  // ------------------------------------------------------ expression QC
  V.exprQC = function (s) {
    const cfg = XQ.cfg;
    const groups = s.groups || {};
    const levels = Array.from(new Set(Object.keys(groups).map(function (k) { return groups[k]; })));
    const gcol = C.groupColors(levels, { level: cfg.level, levelColor: cfg.pcaColors && cfg.pcaColors.level, base: cfg.base, baseColor: cfg.pcaColors && cfg.pcaColors.base });
    const colour = function (nm, k) { return groups[nm] ? gcol[groups[nm]] : C.PALETTE[k % C.PALETTE.length]; };
    const boxed = function (st) { return st.map(function (b, k) { return Object.assign({}, b, { color: colour(b.label, k) }); }); };
    const nodes = [];
    if (s.box) nodes.push(fig('Expression distribution by sample', C.boxes({ stats: boxed(s.box), ylabel: 'log2(normalised count + 1)' }), { filename: 'expression_boxplot' }));
    if (s.dens) nodes.push(fig('Expression density by sample', C.lines({ series: s.dens.map(function (d, k) { return { name: d.name, x: d.x, y: d.y, color: colour(d.name, k) }; }), xlabel: 'log2(normalised count + 1)', ylabel: 'Density', legend: true }), { filename: 'expression_density' }));
    if (s.rle) nodes.push(fig('Relative log expression (RLE)', C.boxes({ stats: boxed(s.rle), ylabel: 'log2 count + 1 minus gene median', zero: true }), { filename: 'rle' }));
    if (s.mv) {
      const t = XQ.table(s.mv);
      if (t.n) {
        const sc = XQ.scatter({ height: 380, selectable: false, pointSize: 3, tip: function (i) { return '<b>' + F.esc(t.str('gene')[i]) + '</b><br>log10(mean+1): ' + num(t.num('lmean')[i], 3) + '<br>log10(var+1): ' + num(t.num('lvar')[i], 3); } });
        sc.api.set({ x: t.num('lmean'), y: t.num('lvar'), xlabel: 'log10(mean normalised count + 1)', ylabel: 'log10(variance + 1)', palette: ['#8B6FC1'], sizes: [3], alphas: [0.55] });
        nodes.push(fig('Mean–variance relationship (most variable genes)', sc, { download: false }));
      }
    }
    const root = el('div', 'xq-view');
    if (!nodes.length) root.append(C.empty('No normalised expression matrix was found for this run.'));
    else root.append(grid.apply(null, nodes.slice(0, 2)), grid.apply(null, nodes.slice(2)));
    return root;
  };

  // --------------------------------------------------------------- enrichment
  const GENE_COLS = ['core_enrichment', 'geneID', 'gene_id', 'geneIDs', 'leadingEdge', 'leading_edge', 'Genes', 'genes'];
  const parseSet = function (x) { return x ? Array.from(new Set(String(x).split(/[\/;,| ]+/).filter(Boolean))) : []; };
  function jaccard(a, b) {
    if (!a.length || !b.length) return 0;
    const sb = new Set(b); let inter = 0;
    a.forEach(function (x) { if (sb.has(x)) inter++; });
    return inter / (a.length + b.length - inter);
  }
  /** Fruchterman-Reingold layout for <= ~40 nodes; weights = Jaccard. */
  function forceLayout(n, edges) {
    const pos = [];
    for (let i = 0; i < n; i++) pos.push({ x: Math.cos(2 * Math.PI * i / n), y: Math.sin(2 * Math.PI * i / n) });
    const k = Math.sqrt(4 / Math.max(1, n));
    let t = 0.25;
    for (let it = 0; it < 260; it++) {
      const disp = pos.map(function () { return { x: 0, y: 0 }; });
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
        let dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y;
        let d = Math.sqrt(dx * dx + dy * dy) || 0.01;
        const f = k * k / d; dx /= d; dy /= d;
        disp[i].x += dx * f; disp[i].y += dy * f; disp[j].x -= dx * f; disp[j].y -= dy * f;
      }
      edges.forEach(function (e) {
        let dx = pos[e.a].x - pos[e.b].x, dy = pos[e.a].y - pos[e.b].y;
        const d = Math.sqrt(dx * dx + dy * dy) || 0.01;
        const f = (d * d / k) * (0.5 + e.j) ; dx /= d; dy /= d;
        disp[e.a].x -= dx * f; disp[e.a].y -= dy * f; disp[e.b].x += dx * f; disp[e.b].y += dy * f;
      });
      for (let i = 0; i < n; i++) {
        const d = Math.sqrt(disp[i].x * disp[i].x + disp[i].y * disp[i].y) || 1e-6;
        const m = Math.min(d, t);
        pos[i].x += disp[i].x / d * m; pos[i].y += disp[i].y / d * m;
      }
      t *= 0.985;
    }
    return pos;
  }

  V.enrichment = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const I = Inputs();
    const root = el('div', 'xq-view');
    if (!tbl.n) { root.append(C.empty('This enrichment table is empty.')); return root; }
    const termC = tbl.resolve('term'), padjC = tbl.resolve('padj');
    if (!termC || !padjC) {
      root.append(el('p', 'xq-intro', 'No recognised term / adjusted p-value columns: showing the raw table only.'), XQ.dataTable(tbl, { pageSize: 10, filename: (s.filename || 'enrichment') + '.csv' }));
      return root;
    }
    const sizeC = tbl.resolve('size'), scoreC = tbl.resolve('score'), collC = tbl.resolve('collection');
    const geneC = firstOf(tbl, GENE_COLS);
    const term = tbl.str(termC), padjV = tbl.num(padjC), nlp = tbl.nlp(padjC);
    const size = sizeC ? tbl.num(sizeC) : null, score = scoreC ? tbl.num(scoreC) : null;
    const isGsea = s.kind === 'gsea' && score;
    const byColl = collC ? groupIdx(tbl, collC) : new Map([['', tbl.all()]]);
    const collections = s.collections && s.collections.length ? s.collections : Array.from(byColl.keys());
    const ALL = 'All collections';
    const collSel = collections.length > 1 ? select([ALL].concat(collections), { value: collections[0], label: 'Collection' }) : null;
    const topN = I.range([5, 40], { step: 1, value: s.topN || 15, label: 'Terms shown' });
    const cut = select(padjChoices(s.padj), { value: s.padj, label: 'adj. p <' });
    const stats = el('div', 'xq-tiles');
    const plots = el('div');
    const tableHost = el('div');
    const fb = XQ.filterBox(s.table, tbl, { examples: ["padj < 0.01", "size >= 10 AND size <= 300"] });
    let dt = null;
    const currentIdx = function () { return collSel && collSel.value !== ALL ? byColl.get(collSel.value) || new Uint32Array(0) : null; };
    let esSel = null, esHost = el('div');
    function topRows(idx) {
      const rows = tbl.order(padjC, false, idx || null);
      const out = [];
      for (let k = 0; k < rows.length && out.length < +topN.value; k++) { const i = rows[k]; if (padjV[i] === padjV[i]) out.push(i); }
      return out;
    }
    function drawES(rows) {
      esHost.replaceChildren();
      if (!isGsea || !s.rank || !geneC) return;
      const de = XQ.table(s.rank.table);
      const gC = de.resolve(s.rank.geneCol || 'gene'), stC = de.first(s.rank.statCols || ['contrast_stat', 'stat']);
      if (!gC || !stC) return;
      const names = rows.map(function (i) { return trunc(term[i], 70); });
      const sel = select(names.map(function (nm, k) { return k + '. ' + nm; }), { label: 'Gene set', value: names.length ? '0. ' + names[0] : undefined });
      const hostPlot = el('div');
      const draw = function () {
        const k = parseInt(String(sel.value).split('.')[0], 10) || 0;
        const i = rows[k];
        const set = new Set(parseSet(tbl.str(geneC)[i]));
        const gv = de.str(gC), sv = de.num(stC);
        const ranked = [];
        for (let r = 0; r < gv.length; r++) if (gv[r] && sv[r] === sv[r]) ranked.push({ g: gv[r], s: sv[r] });
        ranked.sort(function (a, b) { return b.s - a.s; });
        const seen = new Set(); const rk = [];
        ranked.forEach(function (r) { if (!seen.has(r.g)) { seen.add(r.g); rk.push(r); } });
        const hit = rk.map(function (r) { return set.has(r.g); });
        const nh = hit.reduce(function (a, b) { return a + (b ? 1 : 0); }, 0);
        if (nh < 2 || nh >= rk.length) { hostPlot.replaceChildren(C.empty('Fewer than two genes of this set are in the ranked list.')); return; }
        let wsum = 0; rk.forEach(function (r, j) { if (hit[j]) wsum += Math.abs(r.s); });
        const miss = 1 / (rk.length - nh);
        let es = 0, mx = -Infinity, mn = Infinity; const xs = [], ys = [], ticks = [];
        const stepK = Math.max(1, Math.floor(rk.length / 1200));
        rk.forEach(function (r, j) {
          es += hit[j] ? (wsum ? Math.abs(r.s) / wsum : 1 / nh) : -miss;
          if (es > mx) mx = es; if (es < mn) mn = es;
          if (hit[j]) ticks.push(j + 1);
          if (j % stepK === 0 || hit[j]) { xs.push(j + 1); ys.push(es); }
        });
        hostPlot.replaceChildren(fig('Running enrichment score: ' + trunc(term[i], 80), C.lines({ series: [{ name: 'Running ES', x: xs, y: ys, color: '#2b5c8f' }], ticks: ticks, tickFrom: mn - 0.07 * (mx - mn || 1), tickTo: mn - 0.01 * (mx - mn || 1), hline: 0, xlabel: 'Genes ranked by signed test statistic (left = up in ' + (cfg.level || 'test') + ')', ylabel: 'Running enrichment score', grid: false }), { filename: 'gsea_running_score' }));
      };
      onInput(sel, draw);
      esHost.append(h4('Running enrichment score of a gene set'), controls(sel), hostPlot);
      draw();
    }
    function draw() {
      const idx = currentIdx();
      const rows = topRows(idx);
      const nSig = (function () { let n = 0; const l = idx || tbl.all(); for (let k = 0; k < l.length; k++) { const p = padjV[l[k]]; if (p === p && p < +cut.value) n++; } return n; })();
      stats.replaceChildren(XQ.stats([
        { label: 'Terms in table', value: F.int(idx ? idx.length : tbl.n) },
        { label: 'adj. p < ' + cut.value, value: F.int(nSig), cls: 'xq-up' },
        { label: 'Genes per term (median)', value: (function () { if (!size) return '—'; const v = rows.map(function (i) { return size[i]; }).filter(function (x) { return x === x; }).sort(function (a, b) { return a - b; }); return v.length ? F.int(v[Math.floor(v.length / 2)]) : '—'; })() }
      ]));
      if (!rows.length) { plots.replaceChildren(C.empty('No term has an adjusted p-value in this collection.')); drawES(rows); return; }
      // ranked bars
      let bars;
      if (isGsea) {
        const r = rows.slice().sort(function (a, b) { return score[a] - score[b]; });
        bars = C.bars({ xlabel: 'Normalised enrichment score (NES)', data: r.map(function (i) { return { label: trunc(term[i], 60), value: score[i], color: score[i] >= 0 ? cfg.colors.up : cfg.colors.down, tip: term[i] + '\nNES ' + num(score[i], 3) + '\nadj. p ' + F.p(padjV[i]) + (size ? '\nsize ' + F.int(size[i]) : '') }; }) });
      } else {
        const r = rows.slice().reverse();
        bars = C.bars({ xlabel: '-log10(adjusted p-value)', data: r.map(function (i) { return { label: trunc(term[i], 60), value: nlp[i], color: '#0072B2', tip: term[i] + '\nadj. p ' + F.p(padjV[i]) + (size ? '\nsize ' + F.int(size[i]) : '') }; }) });
      }
      const bub = C.bubble({
        xlabel: isGsea ? 'Normalised enrichment score (NES)' : (size ? 'Gene count / overlap size' : 'Rank'), zero: isGsea,
        rows: rows.map(function (i, k) { return { label: trunc(term[i], 60), x: isGsea ? score[i] : (size && size[i] === size[i] ? size[i] : k + 1), nlp: nlp[i], size: size && size[i] === size[i] ? Math.sqrt(size[i]) : 3, tip: term[i] + '\nadj. p ' + F.p(padjV[i]) + (size ? '\nsize ' + F.int(size[i]) : '') + (score ? '\nscore ' + num(score[i], 3) : '') }; })
      });
      const parts = [grid(fig('Top terms by adjusted p-value', bars, { filename: 'enrichment_ranked' }), fig('Enrichment bubble plot', bub, { filename: 'enrichment_bubble' }))];
      // gene-overlap heat map + enrichment-map network (needs gene lists)
      if (geneC) {
        const gs = tbl.str(geneC);
        const o = rows.slice(0, Math.min(+topN.value, 20)).map(function (i) { return { i: i, set: parseSet(gs[i]) }; }).filter(function (x) { return x.set.length >= 2; });
        if (o.length >= 3) {
          const labs = o.map(function (x) { return trunc(term[x.i], 42); });
          const vals = []; o.forEach(function (a) { o.forEach(function (b) { vals.push(jaccard(a.set, b.set)); }); });
          const ov = C.heatmap({ rows: labs, cols: labs, values: vals, scheme: 'blues', domain: [0, 1], label: 'Jaccard' });
          const nn = o.slice(0, 14);
          const edges = [];
          for (let a = 0; a < nn.length; a++) for (let b = a + 1; b < nn.length; b++) { const j = jaccard(nn[a].set, nn[b].set); if (j >= 0.08) edges.push({ a: a, b: b, j: j }); }
          const pos = forceLayout(nn.length, edges);
          const net = C.network({
            nodes: nn.map(function (x, k) { return { label: term[x.i], x: pos[k].x, y: pos[k].y, size: Math.max(6, Math.min(22, Math.sqrt(x.set.length) * 2.4)), nlp: nlp[x.i], tip: term[x.i] + '\nadj. p ' + F.p(padjV[x.i]) + '\ngenes ' + x.set.length }; }),
            edges: edges
          });
          parts.push(grid(fig('Gene overlap between enriched terms (Jaccard)', ov, { filename: 'enrichment_overlap' }), fig('Enrichment map: terms linked by shared genes', net, { filename: 'enrichment_network' })));
        }
      }
      plots.replaceChildren.apply(plots, parts);
      drawES(rows);
    }
    function table() {
      const idx = currentIdx();
      const base = idx;
      const use = fb.value ? (base ? XQ._intersect(base, fb.value, tbl.n) : fb.value) : base;
      if (!dt) { dt = XQ.dataTable(tbl, { idx: use, sort: padjC, desc: false, pageSize: 10, filename: (s.filename || 'enrichment') + '.csv' }); tableHost.replaceChildren(dt); }
      else dt.update(use);
    }
    [collSel, topN, cut].forEach(function (n) { if (n) onInput(n, function () { draw(); if (n === collSel) table(); }, 60); });
    fb.addEventListener('input', table);
    root.append(controls(collSel, topN, cut), stats, plots, esHost, h4('Result table'), fb, tableHost);
    draw(); table();
    return root;
  };
  XQ._intersect = function (a, b, n) {
    const mask = new Uint8Array(n);
    for (let i = 0; i < b.length; i++) mask[b[i]] = 1;
    const out = []; for (let i = 0; i < a.length; i++) if (mask[a[i]]) out.push(a[i]);
    return Uint32Array.from(out);
  };

  // ---------------------------------------------------------------- decoupleR
  V.decouple = function (s) {
    const tbl = XQ.table(s.table);
    const cfg = XQ.cfg;
    const I = Inputs();
    const root = el('div', 'xq-view');
    const srcC = tbl.resolve('source'), scC = tbl.resolve('score'), pC = tbl.first(['p_value', 'padj', 'pval', 'p.value']);
    if (!tbl.n || !srcC || !scC) { root.append(C.empty('No decoupleR activity table.')); return root; }
    const src = tbl.str(srcC), sc = tbl.num(scC), pv = pC ? tbl.num(pC) : null;
    const topN = I.range([5, 50], { step: 1, value: s.topN || 15, label: 'Sources shown' });
    const cut = select(padjChoices(s.padj), { value: s.padj, label: 'p <' });
    const only = I.toggle({ label: 'Only p below cut-off', value: true });
    const host = el('div'), stats = el('div', 'xq-tiles');
    function draw() {
      const nSig = pv ? Array.prototype.filter.call(pv, function (p) { return p === p && p < +cut.value; }).length : null;
      const all = []; for (let i = 0; i < tbl.n; i++) if (sc[i] === sc[i]) all.push(i);
      let sel = all;
      if (only.value && pv) sel = all.filter(function (i) { return pv[i] < +cut.value; });
      const fallback = !sel.length;
      if (fallback) sel = all;
      sel = sel.slice().sort(function (a, b) { return Math.abs(sc[b]) - Math.abs(sc[a]); }).slice(0, +topN.value).sort(function (a, b) { return sc[a] - sc[b]; });
      stats.replaceChildren(XQ.stats([
        { label: 'Sources tested', value: F.int(tbl.n) },
        pv ? { label: 'p < ' + cut.value, value: F.int(nSig), cls: 'xq-up' } : null,
        { label: 'Most active', value: all.length ? String(src[all.reduce(function (a, b) { return sc[b] > sc[a] ? b : a; })]) : '—', sub: 'highest positive score' }
      ].filter(Boolean)));
      const bars = C.bars({ xlabel: 'Activity score (t-like)', data: sel.map(function (i) { return { label: trunc(src[i], 40), value: sc[i], color: sc[i] > 0 ? cfg.colors.up : cfg.colors.down, tip: src[i] + '\nscore ' + num(sc[i], 3) + (pv ? '\np ' + F.p(pv[i]) : '') }; }), empty: 'Nothing to show.' });
      host.replaceChildren(fallback && only.value && pv ? el('p', 'xq-intro', 'No source reaches p < ' + cut.value + ': showing the strongest scores anyway.') : null, fig((s.title || 'Activity') + ' (top ' + sel.length + ')', bars, { filename: 'decoupler_' + (s.resource || 'activity') }));
    }
    [topN, cut, only].forEach(function (n) { onInput(n, draw, 60); });
    const dt = XQ.dataTable(tbl, { sort: pC || scC, desc: !pC, pageSize: 15, filename: (s.filename || 'decouple') + '.csv' });
    const fb = XQ.filterBox(s.table, tbl, { examples: ['score > 2', 'p_value < 0.05 AND score < 0'] });
    fb.addEventListener('input', function () { dt.update(fb.value); });
    root.append(controls(topN, cut, only), stats, host, h4('All tested sources'), fb, dt);
    draw();
    return root;
  };

  // -------------------------------------------------------- TF regulon explorer
  V.regulon = function (s) {
    const reg = XQ.table(s.reg), le = s.le ? XQ.table(s.le) : null;
    const cfg = XQ.cfg;
    const root = el('div', 'xq-view');
    if (!reg.n) { root.append(C.empty('No regulon data.')); return root; }
    const byTf = groupIdx(reg, reg.resolve('tf'));
    const leTf = le && le.n ? groupIdx(le, le.resolve('tf')) : new Map();
    const tfs = (s.tfs && s.tfs.length ? s.tfs : Array.from(byTf.keys())).filter(function (t) { return byTf.has(t); });
    const tfSel = select(tfs, { value: tfs[0], label: 'Transcription factor / source' });
    const host = el('div');
    const CLS = ['Not significant', 'Other change', 'Activated (expected)', 'Repressed (expected)'];
    const PAL = ['#BDBDBD', '#ff7f00', '#1f78b4', '#e31a1c'];
    function draw() {
      const tf = tfSel.value;
      const idx = byTf.get(tf) || new Uint32Array(0);
      const lfc = reg.num('log2FoldChange'), nl = reg.nlp('padj'), regs = reg.str('Regulation'), tg = reg.str('target'), mor = reg.has('mor') ? reg.num('mor') : null;
      const m = idx.length;
      const X = new Float32Array(m), Y = new Float32Array(m), K = new Uint8Array(m);
      const cnt = [0, 0, 0, 0];
      for (let k = 0; k < m; k++) { const i = idx[k]; X[k] = lfc[i]; Y[k] = nl[i]; const c = Math.max(0, CLS.indexOf(regs[i])); K[k] = c; cnt[c]++; }
      const sc = XQ.scatter({ height: 420, pointSize: 6, tip: function (k) { const i = idx[k]; return '<b>' + F.esc(tg[i]) + '</b><br>log2FC: ' + num(lfc[i], 3) + '<br>adj. p: ' + F.nlp(nl[i]) + (mor ? '<br>mor: ' + num(mor[i], 2) : '') + '<br>' + F.esc(regs[i]); } });
      const thr = -Math.log10(Math.max(s.padj, 1e-300));
      sc.api.set({ x: X, y: Y, cls: K, palette: PAL, sizes: [5, 7, 8, 8], alphas: [0.6, 0.9, 0.95, 0.95], xlabel: 'log2 fold change', ylabel: '-log10(adjusted p)', hlines: [{ v: thr }], vlines: s.lfc > 0 ? [{ v: s.lfc }, { v: -s.lfc }] : [], legend: CLS.map(function (c, k) { return { label: c, color: PAL[k], n: cnt[k], cls: k }; }), labels: Array.from(idx).map(function (i, k) { return { i: k, text: tg[i] }; }).filter(function (l) { return K[l.i] >= 2; }).sort(function (a, b) { return Y[b.i] - Y[a.i]; }).slice(0, 8) });
      const parts = [fig('Targets of ' + tf + ' (volcano)', sc, { download: false })];
      const lidx = leTf.get(tf);
      if (le && lidx && lidx.length && s.nRanked) {
        const rank = le.num('rank');
        const pos = Array.from(lidx).map(function (i) { return rank[i]; }).filter(function (r) { return r >= 1 && r <= s.nRanked; }).sort(function (a, b) { return a - b; });
        if (pos.length) {
          const nh = pos.length, hw = 1 / nh, mw = -1 / Math.max(1, s.nRanked - nh);
          const hitSet = new Set(pos);
          let es = 0, mx = 0, mxAt = 1, mn = 0; const xs = [], ys = [];
          const stepK = Math.max(1, Math.floor(s.nRanked / 1200));
          for (let r = 1; r <= s.nRanked; r++) { es += hitSet.has(r) ? hw : mw; if (Math.abs(es) > Math.abs(mx)) { mx = es; mxAt = r; } if (es < mn) mn = es; if (r % stepK === 0 || hitSet.has(r)) { xs.push(r); ys.push(es); } }
          const span = Math.max(1e-6, Math.max.apply(null, ys) - mn);
          parts.push(fig('Leading edge of ' + tf + ' (running enrichment score)', C.lines({ series: [{ name: 'Running ES', x: xs, y: ys, color: '#1f78b4' }], ticks: pos, tickFrom: mn - 0.08 * span, tickTo: mn - 0.01 * span, hline: 0, xlabel: 'Rank by statistic (decreasing). Peak at rank ' + F.int(mxAt), ylabel: 'Running enrichment score', grid: false }), { filename: 'leading_edge_' + tf }));
        }
      }
      const rt = XQ.dataTable(reg, { idx: idx, sort: 'padj', pageSize: 10, filename: 'regulon_' + tf + '.csv', columns: reg.names.filter(function (c) { return c !== 'tf'; }) });
      const tabs = [grid.apply(null, parts), h4('Regulon targets of ' + tf + ' (' + F.int(m) + ' matched in the DE table; ' + cnt[2] + ' activated, ' + cnt[3] + ' repressed as expected)'), rt];
      if (le && lidx && lidx.length) tabs.push(h4('Leading-edge genes (rank order)'), XQ.dataTable(le, { idx: lidx, sort: 'rank', pageSize: 10, filename: 'leading_edge_' + tf + '.csv', columns: le.names.filter(function (c) { return c !== 'tf'; }) }));
      host.replaceChildren.apply(host, tabs);
    }
    onInput(tfSel, draw, 40);
    root.append(controls(tfSel), host);
    draw();
    return root;
  };

  // --------------------------------------------------------------- isoforms
  /** Generic abundance-vs-significance / LR scatter with two classes. */
  function twoClassScatter(o) {
    const sc = XQ.scatter({ height: o.height || 400, pointSize: 3, tip: function (i) { return '<b>' + F.esc(o.names(i)) + '</b><br>' + o.xname + ': ' + num(o.x[i], 3) + '<br>' + o.yname + ': ' + num(o.y[i], 3) + (o.q ? '<br>adj. p: ' + F.nlp(o.q[i]) : ''); } });
    const cls = new Uint8Array(o.x.length);
    let nsig = 0;
    if (o.q) for (let i = 0; i < cls.length; i++) if (o.q[i] >= o.thr) { cls[i] = 1; nsig++; }
    sc.api.set({ x: o.x, y: o.y, cls: cls, palette: ['#9aa0a6', o.color || '#7A6FAE'], sizes: [2.4, 4], alphas: [0.45, 0.9], xlabel: o.xlabel, ylabel: o.ylabel, xlog: !!o.xlog, hlines: o.hline == null ? [] : [{ v: o.hline }], vlines: [],
      legend: o.q ? [{ label: 'Not significant', color: '#9aa0a6', n: cls.length - nsig, cls: 0 }, { label: o.sigLabel || 'Below cut-off', color: o.color || '#7A6FAE', n: nsig, cls: 1 }] : [] });
    return sc;
  }

  function rankedBars(tbl, o) {
    const q = tbl.nlp(o.pCol), nm = labelsOf(tbl, o.idCol), lfc = o.lfcCol ? tbl.num(o.lfcCol) : null;
    const ord = []; for (let i = 0; i < q.length; i++) if (q[i] === q[i]) ord.push(i);
    ord.sort(function (a, b) { return q[b] - q[a]; });
    const top = ord.slice(0, o.n || 25).reverse();
    return C.bars({ xlabel: '-log10(adjusted p-value)', data: top.map(function (i) { return { label: trunc(String(nm(i)), 55), value: q[i], color: lfc ? (lfc[i] < 0 ? XQ.cfg.colors.down : XQ.cfg.colors.up) : (o.color || '#009E73'), tip: nm(i) + '\nadj. p ' + F.nlp(q[i]) + (lfc ? '\nlog2FC ' + num(lfc[i], 3) : '') }; }) });
  }

  V.dte = function (s) {
    const tbl = XQ.table(s.table);
    const root = el('div', 'xq-view');
    if (!tbl.n) { root.append(C.empty('No DTE transcript results.')); return root; }
    const idC = tbl.first(['gene_label', 'gene_symbol', 'gene', 'transcript_id', 'feature_id']);
    const lfcC = tbl.first(['log2FoldChange', 'log2FC', 'logFC', 'lfc']);
    const qC = tbl.first(['padj', 'adj_pvalue', 'FDR', 'qvalue', 'p.adjust']);
    const rawC = tbl.first(['pvalue', 'p_value']);
    const meanC = tbl.first(['baseMean', 'base_mean', 'mean', 'mean_expression']);
    const ex = V.deExplorer({ table: s.table, idCol: idC, lfcCol: lfcC, pCol: qC, meanCol: meanC, padj: s.padj, lfc: s.lfc, topN: 8, title: 'DTE volcano / MA (transcript level)', filename: 'dte', tipCols: s.tipCols, columns: s.columns });
    root.append(ex);
    const nodes = [];
    if (meanC && qC) {
      const mean = tbl.num(meanC), lm = new Float32Array(mean.length);
      for (let i = 0; i < mean.length; i++) lm[i] = mean[i] > 0 ? Math.log10(mean[i]) : NaN;
      nodes.push(fig('Transcript abundance vs adjusted significance', twoClassScatter({ x: lm, y: tbl.nlp(qC), q: tbl.nlp(qC), thr: -Math.log10(s.padj), names: labelsOf(tbl, idC), xname: 'log10(mean)', yname: '-log10(adj. p)', xlabel: 'log10(mean normalised abundance)', ylabel: '-log10(adjusted p)', hline: -Math.log10(s.padj), sigLabel: 'adj. p below cut-off' }), { download: false }));
    }
    if (qC) nodes.push(fig('Top 25 transcripts by adjusted p-value', rankedBars(tbl, { pCol: qC, idCol: idC, lfcCol: lfcC, n: 25 }), { filename: 'dte_top' }));
    if (rawC) nodes.push(fig('DTE raw p-value distribution', C.hist(tbl.num(rawC)), { filename: 'dte_raw_p' }));
    if (qC) nodes.push(fig('DTE adjusted p-value distribution', C.hist(tbl.num(qC), { color: '#56B4E9' }), { filename: 'dte_adj_p' }));
    if (nodes.length) root.append(h4('Diagnostics'), grid.apply(null, nodes.slice(0, 2)), grid.apply(null, nodes.slice(2)));
    return root;
  };

  V.dtu = function (s) {
    const tbl = XQ.table(s.table);
    const root = el('div', 'xq-view');
    if (!tbl.n) { root.append(C.empty('No DTU feature results.')); return root; }
    const idC = tbl.first(['gene_label', 'gene_symbol', 'gene_id', 'feature_id', 'transcript_id']);
    const qC = tbl.first(['adj_pvalue', 'padj', 'FDR', 'qvalue', 'p.adjust']);
    const rawC = tbl.first(['pvalue', 'p_value']);
    const lrC = tbl.first(['lr', 'LR', 'likelihood_ratio', 'statistic']);
    const dfC = tbl.first(['df', 'Df', 'degrees_of_freedom']);
    const geneC = tbl.first(['gene_symbol', 'gene', 'gene_id', 'gene_label']);
    const nodes = [];
    const pC = qC || rawC;
    if (pC) nodes.push(fig('Top 30 DTU features', rankedBars(tbl, { pCol: pC, idCol: idC, n: 30, color: '#009E73' }), { filename: 'dtu_top' }));
    if (lrC && qC) nodes.push(fig('Likelihood ratio vs adjusted significance', twoClassScatter({ x: tbl.num(lrC), y: tbl.nlp(qC), q: tbl.nlp(qC), thr: -Math.log10(s.padj), names: labelsOf(tbl, idC), xname: 'LR', yname: '-log10(adj. p)', xlabel: 'Likelihood-ratio statistic', ylabel: '-log10(adjusted p)', hline: -Math.log10(s.padj) }), { download: false }));
    if (lrC && dfC) nodes.push(fig('Likelihood ratio vs degrees of freedom', twoClassScatter({ x: tbl.num(dfC), y: tbl.num(lrC), q: null, names: labelsOf(tbl, idC), xname: 'df', yname: 'LR', xlabel: 'Degrees of freedom', ylabel: 'Likelihood-ratio statistic', color: '#4C78A8' }), { download: false }));
    if (geneC) {
      const gs = tbl.str(geneC), q = qC ? tbl.num(qC) : null, m = new Map();
      for (let i = 0; i < gs.length; i++) { const k = gs[i] || 'Unannotated'; const e = m.get(k) || { n: 0, q: Infinity }; e.n++; if (q && q[i] === q[i] && q[i] < e.q) e.q = q[i]; m.set(k, e); }
      const arr = Array.from(m.entries()).sort(function (a, b) { return b[1].n - a[1].n || a[1].q - b[1].q; }).slice(0, 30).reverse();
      nodes.push(fig('Genes with the most tested DTU features', C.bars({ xlabel: 'Number of DTU features', data: arr.map(function (e) { return { label: trunc(e[0], 40), value: e[1].n, color: '#7A6FAE', tip: e[0] + '\nfeatures ' + e[1].n + (isFinite(e[1].q) ? '\nbest adj. p ' + F.p(e[1].q) : '') }; }) }), { filename: 'dtu_gene_summary' }));
    }
    if (rawC) nodes.push(fig('DTU raw p-value distribution', C.hist(tbl.num(rawC)), { filename: 'dtu_raw_p' }));
    if (qC) nodes.push(fig('DTU adjusted p-value distribution', C.hist(tbl.num(qC), { color: '#56B4E9' }), { filename: 'dtu_adj_p' }));
    if (lrC) nodes.push(fig('Likelihood-ratio distribution', C.hist(tbl.num(lrC), { domain: 'auto', color: '#7A6FAE', xlabel: 'Likelihood-ratio statistic' }), { filename: 'dtu_lr_hist' }));
    for (let k = 0; k < nodes.length; k += 2) root.append(grid(nodes[k], nodes[k + 1]));
    return root;
  };

  V.switch = function (s) {
    const tbl = XQ.table(s.table);
    const I = Inputs();
    const root = el('div', 'xq-view');
    if (!tbl.n || !tbl.has('dIF') || !tbl.has('isoform_switch_q_value')) { root.append(C.empty('No switch table with dIF / isoform_switch_q_value.')); return root; }
    const dif = tbl.num('dIF'), q = tbl.nlp('isoform_switch_q_value');
    const idC = tbl.first(['gene_name', 'gene_id', 'isoform_id']);
    const nm = labelsOf(tbl, idC);
    const cutQ = select(padjChoices(s.padj), { value: s.padj, label: 'q <' });
    const cutD = I.range([0, 0.6], { step: 0.01, value: s.dif || 0.1, label: '|dIF| ≥' });
    const tiles = el('div', 'xq-tiles'), tableTitle = h4('');
    const sc = XQ.scatter({ height: 460, pointSize: 3, tip: function (i) { return '<b>' + F.esc(nm(i)) + '</b><br>dIF: ' + num(dif[i], 3) + '<br>q: ' + F.nlp(q[i]); } });
    const tableHost = el('div'); let dt = null;
    function refresh() {
      const thr = -Math.log10(+cutQ.value), d = +cutD.value;
      const cls = new Uint8Array(dif.length); const sig = [];
      for (let i = 0; i < dif.length; i++) if (q[i] >= thr && Math.abs(dif[i]) >= d) { cls[i] = 1; sig.push(i); }
      sig.sort(function (a, b) { return q[b] - q[a]; });
      sc.api.set({ x: dif, y: q, cls: cls, palette: ['#9aa0a6', '#c94a4a'], sizes: [2.6, 5], alphas: [0.45, 0.9], xlabel: 'dIF (change in isoform fraction)', ylabel: '-log10(q)', hlines: [{ v: thr }], vlines: d > 0 ? [{ v: d }, { v: -d }] : [], legend: [{ label: 'Not significant', color: '#9aa0a6', n: dif.length - sig.length, cls: 0 }, { label: 'Switch', color: '#c94a4a', n: sig.length, cls: 1 }], labels: sig.slice(0, 8).map(function (i) { return { i: i, text: String(nm(i)) }; }) });
      tiles.replaceChildren(XQ.stats([{ label: 'Isoforms', value: F.int(dif.length) }, { label: 'Significant switches', value: F.int(sig.length), cls: 'xq-up', sub: 'q < ' + cutQ.value + ' & |dIF| ≥ ' + d }]));
      tableTitle.textContent = 'Significant isoform switches (' + F.int(sig.length) + ')';
      const idx = Uint32Array.from(sig);
      if (!dt) { dt = XQ.dataTable(tbl, { idx: idx, sort: 'isoform_switch_q_value', pageSize: 10, filename: 'isoform_switches.csv' }); tableHost.replaceChildren(dt); } else dt.update(idx);
    }
    [cutQ, cutD].forEach(function (n) { onInput(n, refresh, 80); });
    root.append(controls(cutQ, cutD), tiles, fig('Isoform switch overview', sc, { download: false }), tableTitle, tableHost);
    refresh();
    return root;
  };

  /** Free SQL over several tables (DuckDB-Wasm). */
  V.sql = function (s) {
    const t = {};
    Object.keys(s.tables).forEach(function (name) { t[name] = XQ.table(s.tables[name]); });
    return stack(el('p', 'xq-intro', s.intro || ''), XQ.sqlConsole(t, { sql: s.sql }));
  };

  // ----------------------------------------------------------- hints toggle
  if (doc) {
    const KEY = 'expressom-hints';
    const get = function () { try { return g.localStorage.getItem(KEY); } catch (e) { return null; } };
    const put = function (v) { try { g.localStorage.setItem(KEY, v); } catch (e) { /* noop */ } };
    const init = function () {
      if (!doc.querySelector('.xq-help, script[id^="xqh-"]')) return;
      const b = el('button', 'xq-hint-toggle');
      b.type = 'button';
      const apply = function (hide) { doc.body.classList.toggle('xq-nohelp', hide); b.textContent = hide ? 'Show hints' : 'Hide hints'; b.setAttribute('aria-pressed', String(!hide)); };
      apply(get() === 'hide');
      b.addEventListener('click', function () {
        const hide = !doc.body.classList.contains('xq-nohelp');
        put(hide ? 'hide' : 'show'); apply(hide);
        if (!hide) doc.querySelectorAll('details.xq-help').forEach(function (d) { d.open = true; });
      });
      doc.body.append(b);
    };
    if (doc.readyState === 'loading') doc.addEventListener('DOMContentLoaded', init); else init();
  }
})(typeof window !== 'undefined' ? window : globalThis);
