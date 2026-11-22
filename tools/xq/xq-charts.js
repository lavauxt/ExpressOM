/* xq-charts.js -- chart recipes for the ExpressOM dashboard.
 *
 * Every recipe returns a DOM node that renders itself when it first becomes
 * visible (hidden tabs cost nothing) and re-renders when its width changes.
 * Large point clouds (volcano, MA, Q-Q) go through XQ.scatter (WebGL); small
 * summaries use Observable Plot (SVG). All input is plain arrays / XTables
 * that R embedded -- no chart ever recomputes statistics R already did.
 */
(function (g) {
  'use strict';
  const XQ = g.XQ;
  const doc = g.document;
  const el = XQ.el, F = XQ.fmt;
  const C = (XQ.charts = {});
  const PLOT = function () { return g.XQV.Plot; };

  // Okabe-Ito: colour-blind safe categorical palette
  C.PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#8B6FC1', '#999999', '#7A3E00', '#2A6F6F'];

  // ------------------------------------------------------------------ mounting
  /** Render `build(width)` when visible; re-render when the width changes. */
  function mount(build, o) {
    o = o || {};
    const host = el('div', 'xq-chart');
    host.style.minHeight = (o.minHeight || 120) + 'px';
    let lastW = 0, started = false;
    const render = function () {
      const w = Math.max(300, Math.floor(host.clientWidth || 0));
      lastW = w;
      try {
        const node = build(w);
        host.replaceChildren(node || el('div', 'xq-empty', 'Nothing to plot.'));
        host.style.minHeight = '';
      } catch (e) {
        console.error('XQ chart', e);
        host.replaceChildren(el('div', 'xq-err', 'Could not draw this chart: ' + (e && e.message ? e.message : e)));
      }
    };
    const start = function () { if (started) return; started = true; render(); };
    if (!('IntersectionObserver' in g)) { setTimeout(start, 0); return host; }
    const io = new IntersectionObserver(function (es) {
      if (es.some(function (x) { return x.isIntersecting; })) start();
    }, { rootMargin: '120px' });
    io.observe(host);
    let tmr = null;
    if ('ResizeObserver' in g) {
      const ro = new ResizeObserver(function () {
        if (!started) return;
        clearTimeout(tmr);
        tmr = setTimeout(function () {
          const w = Math.floor(host.clientWidth || 0);
          if (w && Math.abs(w - lastW) > 14) render();
        }, 220);
      });
      ro.observe(host);
      host._xqDispose = function () { io.disconnect(); ro.disconnect(); };
    } else host._xqDispose = function () { io.disconnect(); };
    return host;
  }
  C.mount = mount;

  const trunc = function (s, n) { s = String(s == null ? '' : s); return s.length > n ? s.slice(0, n - 1) + '…' : s; };
  const labelMargin = function (labels, max) {
    let m = 0;
    for (let i = 0; i < labels.length; i++) m = Math.max(m, String(labels[i]).length);
    return Math.min(max || 340, Math.round(m * 6.1 + 18));
  };
  const empty = function (msg) { return el('div', 'xq-empty', msg || 'No data for this view.'); };

  // ----------------------------------------------------------- DE scatter (GL)
  /**
   * Volcano / MA explorer on WebGL. `el.api.update({view, padj, lfc, gene, topN})`
   * re-classifies points and returns {up, down, ns, n, sig, shown}.
   */
  C.deScatter = function (tbl, o) {
    o = o || {};
    const idC = tbl.resolve(o.idCol) || tbl.first(['gene', 'gene_symbol', 'transcript_id', 'gene_id', 'feature_id']);
    const lfcC = tbl.resolve(o.lfcCol || 'log2FoldChange');
    const pC = tbl.resolve(o.pCol || 'padj');
    const mC = tbl.resolve(o.meanCol || 'baseMean');
    if (!tbl.n || !lfcC || !pC) { const e = empty('This table has no fold change / adjusted p-value columns.'); e.api = { update: function () { return { up: 0, down: 0, ns: 0, n: 0, sig: 0 }; } }; return e; }
    const lfc = tbl.num(lfcC);
    const nlp = tbl.nlp(pC);
    const names = idC ? tbl.str(idC) : null;
    const mean = mC ? tbl.num(mC) : null;
    let logMean = null;
    if (mean) { logMean = new Float32Array(mean.length); for (let i = 0; i < mean.length; i++) logMean[i] = mean[i] > 0 ? Math.log10(mean[i]) : NaN; }
    const col = Object.assign({ up: '#d62728', down: '#1f5fa8', ns: '#9aa0a6' }, o.colors || {});
    const extra = (o.tipCols || []).filter(function (c) { return tbl.has(c); });
    const lfcName = o.lfcLabel || 'log2FC', pName = o.pLabel || 'adj. p';
    const tip = function (i) {
      let h = '<b>' + F.esc(names ? names[i] : '#' + (i + 1)) + '</b>';
      h += '<br>' + F.esc(lfcName) + ': ' + F.num(lfc[i], 3) + '<br>' + F.esc(pName) + ': ' + F.nlp(nlp[i]);
      if (mean) h += '<br>mean: ' + F.num(mean[i], 3);
      extra.forEach(function (c) { const v = tbl.get(c, i); if (v != null) h += '<br>' + F.esc(c) + ': ' + F.esc(typeof v === 'number' ? F.num(v, 3) : v); });
      return h;
    };
    const sc = XQ.scatter({ height: o.height || 480, tip: tip, pointSize: o.pointSize || 3, filename: o.filename || 'expressom_volcano' });
    const palette = [col.ns, col.down, col.up];
    let last = {};
    sc.api.update = function (st) {
      st = Object.assign({ view: 'volcano', padj: 0.05, lfc: 1, gene: '', topN: 12 }, last, st || {});
      last = st;
      const n = tbl.n;
      const thr = -Math.log10(Math.max(st.padj, 1e-300));
      const noDir = !!o.noDirection;     // LRT / omnibus tests: no direction of change
      const lfcCut = noDir ? 0 : st.lfc;
      const cls = new Uint8Array(n);
      let up = 0, down = 0, sig = [];
      for (let i = 0; i < n; i++) {
        const l = lfc[i], p = nlp[i];
        if (p >= thr && (noDir || Math.abs(l) >= lfcCut)) {
          cls[i] = noDir || l > 0 ? 2 : 1;
          if (noDir || l > 0) up++; else down++;
          sig.push(i);
        }
      }
      sig.sort(function (a, b) { return nlp[b] - nlp[a]; });
      const counts = { up: noDir ? 0 : up, down: down, ns: n - up - down, n: n, sig: up + down, idx: Uint32Array.from(sig) };
      const labels = [];
      if (st.topN > 0 && names) {
        for (let k = 0; k < Math.min(st.topN, sig.length); k++) labels.push({ i: sig[k], text: names[sig[k]] });
      }
      const ring = [];
      const q = String(st.gene || '').trim();
      if (q && idC) {
        const seen = new Set();
        q.split(/[\s,;]+/).filter(Boolean).forEach(function (tok) {
          const hit = tbl.search(tok, [idC]);
          if (hit) for (let k = 0; k < hit.length && seen.size < 25; k++) {
            if (seen.has(hit[k])) continue;
            seen.add(hit[k]); ring.push(hit[k]); labels.push({ i: hit[k], text: names[hit[k]] });
          }
        });
      }
      counts.found = ring.length;
      counts.sigIdx = Uint32Array.from(sig);
      const legend = noDir
        ? [{ label: 'Not significant', color: col.ns, n: counts.ns, cls: 0 }, { label: 'Significant', color: col.up, n: up, cls: 2 }]
        : [
            { label: 'Down', color: col.down, n: down, cls: 1 },
            { label: 'Not significant', color: col.ns, n: counts.ns, cls: 0 },
            { label: 'Up', color: col.up, n: up, cls: 2 }
          ];
      const common = { cls: cls, palette: palette, sizes: [o.nsSize || 2, 3.6, 3.6], alphas: [0.4, 0.9, 0.9], legend: legend, labels: labels, ring: ring };
      if (st.view === 'ma' && logMean) {
        const hl = lfcCut > 0 ? [{ v: lfcCut }, { v: -lfcCut }] : [];
        hl.push({ v: 0, solid: true });
        sc.api.set(Object.assign({ x: logMean, y: lfc, xlog: true, xlabel: 'Mean of normalised counts (log scale)', ylabel: 'log2 fold change', hlines: hl, vlines: [], diag: false }, common));
      } else {
        sc.api.set(Object.assign({ x: lfc, y: nlp, xlog: false, xlabel: o.xlabel || 'log2 fold change', ylabel: o.ylabel || '-log10(adjusted p)',
          hlines: [{ v: thr }], vlines: lfcCut > 0 ? [{ v: lfcCut }, { v: -lfcCut }] : [], diag: false }, common));
      }
      return counts;
    };
    return sc;
  };

  /** p-value Q-Q plot on WebGL. */
  C.qq = function (tbl, pCol, o) {
    o = o || {};
    const pc = tbl.resolve(pCol || 'pvalue');
    if (!pc) return empty('No p-value column.');
    const nl = tbl.nlp(pc);
    const v = [];
    for (let i = 0; i < nl.length; i++) if (nl[i] === nl[i]) v.push(nl[i]);
    if (v.length < 3) return empty('At least three finite p-values are needed.');
    v.sort(function (a, b) { return b - a; }); // largest -log10(p) first = smallest p first
    const n = v.length;
    const x = new Float32Array(n), y = new Float32Array(n);
    for (let k = 0; k < n; k++) { x[k] = -Math.log10((k + 0.5) / n); y[k] = v[k]; }
    const sc = XQ.scatter({ height: o.height || 380, selectable: false, pointSize: 3, tip: function (i) { return 'Expected: ' + x[i].toFixed(2) + '<br>Observed: ' + y[i].toFixed(2); } });
    sc.api.set({ x: x, y: y, xlabel: 'Expected -log10(p)', ylabel: 'Observed -log10(p)', diag: true, palette: [o.color || '#0072B2'], sizes: [3], alphas: [0.6] });
    return sc;
  };

  // -------------------------------------------------------------- Plot recipes
  /** Horizontal bars: data = [{label, value, color, tip}]. */
  C.bars = function (o) {
    return mount(function (w) {
      const d = o.data || [];
      if (!d.length) return empty(o.empty);
      const P = PLOT();
      const ml = labelMargin(d.map(function (x) { return x.label; }), o.maxMargin || 300);
      return P.plot({
        width: w, height: Math.max(110, 24 * d.length + 56), marginLeft: ml, marginRight: 16, marginBottom: 38,
        x: { label: o.xlabel || null, grid: true }, y: { label: null, domain: d.map(function (x) { return x.label; }) },
        marks: [
          P.barX(d, { x: 'value', y: 'label', fill: { value: function (x) { return x.color || '#0072B2'; }, scale: null }, title: function (x) { return x.tip || (x.label + ': ' + F.num(x.value, 3)); }, tip: true }),
          P.ruleX([0])
        ].concat(o.ref != null ? [P.ruleX([o.ref], { stroke: '#666', strokeDasharray: '4 3' })] : [])
      });
    }, { minHeight: 140 });
  };

  /** Histogram. Default: p-values on [0,1] with the uniform-null reference line. */
  C.hist = function (values, o) {
    o = o || {};
    return mount(function (w) {
      const dom = o.domain || [0, 1];
      const isP = dom[0] === 0 && dom[1] === 1 && o.uniform !== false;
      const nb = o.bins || (isP ? 20 : 40);
      let lo = dom[0], hi = dom[1];
      if (o.domain === 'auto') {
        lo = Infinity; hi = -Infinity;
        for (let i = 0; i < values.length; i++) { const v = values[i]; if (v === v && isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } }
        if (!isFinite(lo) || lo === hi) { lo = 0; hi = (isFinite(hi) ? hi : 1) + 1; }
      }
      const cnt = new Array(nb).fill(0);
      let n = 0;
      for (let i = 0; i < values.length; i++) {
        const v = values[i];
        if (v === v && v >= lo && v <= hi) { cnt[Math.min(nb - 1, Math.floor((v - lo) / (hi - lo) * nb))]++; n++; }
      }
      if (!n) return empty('No finite values.');
      const P = PLOT();
      const step = (hi - lo) / nb;
      const bins = cnt.map(function (c, k) { return { x0: lo + k * step, x1: lo + (k + 1) * step, n: c }; });
      const marks = [P.rectY(bins, { x1: 'x0', x2: 'x1', y: 'n', fill: o.color || '#0072B2', inset: 0.5, tip: true, title: function (b) { return F.num(b.x0, 3) + ' – ' + F.num(b.x1, 3) + ': ' + b.n.toLocaleString('en-US'); } })];
      if (isP) marks.push(P.ruleY([n / nb], { stroke: '#d55e00', strokeDasharray: '5 4' }));
      marks.push(P.ruleY([0]));
      return P.plot({
        width: w, height: o.height || 300, marginLeft: 56,
        x: { label: o.xlabel || (isP ? 'p-value' : 'value'), domain: [lo, hi] }, y: { label: 'Number of features', grid: true },
        marks: marks
      });
    }, { minHeight: 200 });
  };

  /** Heatmap. o: rows[], cols[], values (row-major), scheme, domain, pivot, label, annot[colors per col] */
  C.heatmap = function (o) {
    return mount(function (w) {
      const rows = o.rows || [], cols = o.cols || [];
      if (!rows.length || !cols.length) return empty(o.empty);
      const P = PLOT();
      const cells = [];
      for (let r = 0; r < rows.length; r++) for (let c = 0; c < cols.length; c++) {
        const v = o.values[r * cols.length + c];
        cells.push({ row: rows[r], col: cols[c], v: v === v ? v : null });
      }
      const ml = labelMargin(rows, 220);
      const mb = Math.min(160, Math.round(labelMargin(cols, 200) * 0.75));
      const cs = Object.assign({ legend: true, label: o.label || 'value' }, o.scheme === 'rdbu' || !o.scheme
        ? { type: 'diverging', scheme: 'RdBu', reverse: true, pivot: o.pivot == null ? 0 : o.pivot }
        : { scheme: o.scheme });
      if (o.domain) cs.domain = o.domain;
      const cellH = Math.max(14, Math.min(26, Math.floor(520 / rows.length)));
      const marks = [P.cell(cells, { x: 'col', y: 'row', fill: 'v', inset: 0.4, tip: true, title: function (d) { return d.row + ' · ' + d.col + '\n' + (d.v == null ? 'NA' : F.num(d.v, 3)); } })];
      if (o.annot && o.annot.length === cols.length) {
        marks.push(P.cell(cols.map(function (c, k) { return { col: c, color: o.annot[k], g: o.annotLabels ? o.annotLabels[k] : '' }; }),
          { x: 'col', y: function () { return o.annotName || 'group'; }, fill: { value: 'color', scale: null }, inset: 0.4, tip: true, title: function (d) { return d.col + (d.g ? ' · ' + d.g : ''); } }));
      }
      const ydom = (o.annot && o.annot.length === cols.length ? [o.annotName || 'group'] : []).concat(rows);
      return P.plot({
        width: w, height: Math.max(180, ydom.length * cellH + mb + 50), marginLeft: ml, marginBottom: mb, marginTop: 24,
        x: { domain: cols, label: null, tickRotate: -40, axis: 'bottom' }, y: { domain: ydom, label: null },
        color: cs, marks: marks
      });
    }, { minHeight: 200 });
  };

  /** Box summary per label: stats = [{label, lo, q1, med, q3, hi, color}] */
  C.boxes = function (o) {
    return mount(function (w) {
      const d = o.stats || [];
      if (!d.length) return empty(o.empty);
      const P = PLOT();
      const mb = Math.min(140, Math.round(labelMargin(d.map(function (x) { return x.label; }), 200) * 0.7) + 8);
      const marks = [];
      if (o.zero) marks.push(P.ruleY([0], { stroke: '#777', strokeDasharray: '4 3' }));
      marks.push(
        P.ruleX(d, { x: 'label', y1: 'lo', y2: 'hi', stroke: 'currentColor', strokeOpacity: 0.6 }),
        P.barY(d, { x: 'label', y1: 'q1', y2: 'q3', fill: { value: function (x) { return x.color || '#0072B2'; }, scale: null }, fillOpacity: 0.8, tip: true, title: function (x) { return x.label + '\nmedian ' + F.num(x.med, 3) + '\nQ1–Q3 ' + F.num(x.q1, 3) + ' – ' + F.num(x.q3, 3) + '\nrange ' + F.num(x.lo, 3) + ' – ' + F.num(x.hi, 3); } }),
        P.tickY(d, { x: 'label', y: 'med', stroke: '#fff', strokeWidth: 2 })
      );
      return P.plot({ width: w, height: o.height || 340, marginBottom: mb, marginLeft: 56, x: { label: null, domain: d.map(function (x) { return x.label; }), tickRotate: -35 }, y: { label: o.ylabel || null, grid: true }, marks: marks });
    }, { minHeight: 200 });
  };

  /** Line series: series = [{name, x[], y[], color, points: bool}] */
  C.lines = function (o) {
    return mount(function (w) {
      const s = o.series || [];
      if (!s.length) return empty(o.empty);
      const P = PLOT();
      const rows = [];
      s.forEach(function (it, k) { for (let i = 0; i < it.x.length; i++) rows.push({ name: it.name, x: it.x[i], y: it.y[i], color: it.color || C.PALETTE[k % C.PALETTE.length] }); });
      const marks = [];
      if (o.hline != null) marks.push(P.ruleY([o.hline], { stroke: '#777', strokeDasharray: '4 3' }));
      marks.push(P.line(rows, { x: 'x', y: 'y', z: 'name', stroke: { value: 'color', scale: null }, strokeWidth: o.strokeWidth || 1.6, tip: o.tip !== false, title: function (d) { return d.name + '\n' + F.num(d.x, 3) + ', ' + F.num(d.y, 3); } }));
      if (o.ticks) marks.push(P.ruleX(o.ticks, { y1: o.tickFrom == null ? 0 : o.tickFrom, y2: o.tickTo == null ? 0.03 : o.tickTo, stroke: '#c94a4a', strokeOpacity: 0.7 }));
      if (o.legend && s.length > 1 && s.length <= 24) {
        // external legend (colour swatches)
        const lg = el('div', 'xq-legend');
        s.forEach(function (it, k) { const sp = el('span', 'xq-lg2'); const dot = el('i'); dot.style.background = it.color || C.PALETTE[k % C.PALETTE.length]; sp.append(dot, doc.createTextNode(it.name)); lg.append(sp); });
        const wrap = el('div');
        wrap.append(P.plot({ width: w, height: o.height || 320, marginLeft: 56, x: { label: o.xlabel, grid: o.grid !== false }, y: { label: o.ylabel, grid: true }, marks: marks }), lg);
        return wrap;
      }
      return P.plot({ width: w, height: o.height || 320, marginLeft: 56, x: { label: o.xlabel, grid: o.grid !== false, ticks: o.xticks }, y: { label: o.ylabel, grid: true }, marks: marks });
    }, { minHeight: 200 });
  };

  /** Dendrogram: segs = [{x1,y1,x2,y2}], leaves = [label...] (x positions 1..n) */
  C.dendrogram = function (o) {
    return mount(function (w) {
      const segs = o.segs || [];
      if (!segs.length) return empty(o.empty);
      const P = PLOT();
      const n = (o.leaves || []).length;
      const mb = Math.min(150, Math.round(labelMargin(o.leaves || [], 200) * 0.75) + 8);
      return P.plot({
        width: w, height: o.height || 340, marginBottom: mb, marginLeft: 56,
        x: { domain: [0.5, n + 0.5], ticks: (o.leaves || []).map(function (_, i) { return i + 1; }), tickFormat: function (i) { return o.leaves[i - 1]; }, tickRotate: -45, label: null },
        y: { label: o.ylabel || 'Height', grid: true, zero: true },
        marks: [P.link(segs, { x1: 'x1', y1: 'y1', x2: 'x2', y2: 'y2', stroke: '#496a81', strokeWidth: 1.8 })]
      });
    }, { minHeight: 200 });
  };

  /** Colour per group value; the main contrast uses the configured level/base colours. */
  C.groupColors = function (levels, o) {
    o = o || {};
    const out = {};
    let k = 0;
    levels.forEach(function (lv) {
      if (o.level != null && lv === o.level && o.levelColor) out[lv] = XQ.cssColor(o.levelColor);
      else if (o.base != null && lv === o.base && o.baseColor) out[lv] = XQ.cssColor(o.baseColor);
      else { while (o.level != null && C.PALETTE[k % C.PALETTE.length] === out[o.level]) k++; out[lv] = C.PALETTE[k++ % C.PALETTE.length]; }
    });
    return out;
  };

  /** Labelled point plot (PCA, MDS...). pts = [{x,y,label,group,color}] */
  C.points = function (o) {
    return mount(function (w) {
      const d = o.pts || [];
      if (!d.length) return empty(o.empty);
      const P = PLOT();
      const marks = [
        P.dot(d, { x: 'x', y: 'y', r: o.r || 7, fill: { value: 'color', scale: null }, stroke: '#fff', strokeWidth: 1.2, fillOpacity: 0.92, tip: true, title: function (p) { return p.label + (p.group ? '\n' + p.group : '') + '\n' + F.num(p.x, 3) + ', ' + F.num(p.y, 3); } })
      ];
      if (o.labels !== false) marks.push(P.text(d, { x: 'x', y: 'y', text: 'label', dy: -13, fontSize: 11, fill: 'currentColor', stroke: 'var(--xq-bg, #fff)', strokeWidth: 3, paintOrder: 'stroke' }));
      if (o.zero !== false) marks.push(P.ruleX([0], { stroke: '#999', strokeDasharray: '3 3' }), P.ruleY([0], { stroke: '#999', strokeDasharray: '3 3' }));
      const wrap = el('div');
      wrap.append(P.plot({ width: w, height: o.height || 430, marginLeft: 58, marginTop: 20, grid: true, x: { label: o.xlabel, nice: true }, y: { label: o.ylabel, nice: true }, marks: marks }));
      if (o.legend && o.legend.length) {
        const lg = el('div', 'xq-legend');
        o.legend.forEach(function (it) { const sp = el('span', 'xq-lg2'); const dot = el('i'); dot.style.background = it.color; sp.append(dot, doc.createTextNode(it.label)); lg.append(sp); });
        wrap.append(lg);
      }
      return wrap;
    }, { minHeight: 240 });
  };

  /** Enrichment bubble: rows = [{label, x, nlp, size, tip}] */
  C.bubble = function (o) {
    return mount(function (w) {
      const d = o.rows || [];
      if (!d.length) return empty(o.empty);
      const P = PLOT();
      const ml = labelMargin(d.map(function (x) { return x.label; }), 320);
      return P.plot({
        width: w, height: Math.max(160, 26 * d.length + 60), marginLeft: ml, marginRight: 16,
        x: { label: o.xlabel, grid: true }, y: { label: null, domain: d.map(function (x) { return x.label; }) },
        r: { range: [4, 15] }, color: { scheme: 'viridis', legend: true, label: '-log10(adj. p)' },
        marks: [P.dot(d, { x: 'x', y: 'label', r: 'size', fill: 'nlp', stroke: '#fff', strokeWidth: 0.8, tip: true, title: function (x) { return x.tip; } }), o.zero ? P.ruleX([0]) : null].filter(Boolean)
      });
    }, { minHeight: 180 });
  };

  /** Shared-gene network: nodes = [{label, x, y, size, nlp, tip}], edges = [{a, b, j}] */
  C.network = function (o) {
    return mount(function (w) {
      const nodes = o.nodes || [];
      if (nodes.length < 2) return empty(o.empty || 'Too few terms with shared genes to draw a network.');
      const P = PLOT();
      const link = (o.edges || []).map(function (e) { return { x1: nodes[e.a].x, y1: nodes[e.a].y, x2: nodes[e.b].x, y2: nodes[e.b].y, j: e.j }; });
      return P.plot({
        width: w, height: o.height || 520, margin: 30, marginTop: 24, axis: null, aspectRatio: 1,
        color: { scheme: 'viridis', legend: true, label: '-log10(adj. p)' }, r: { range: [6, 22] },
        marks: [
          P.link(link, { x1: 'x1', y1: 'y1', x2: 'x2', y2: 'y2', stroke: '#697d91', strokeOpacity: 0.35, strokeWidth: function (d) { return 1 + 6 * d.j; } }),
          P.dot(nodes, { x: 'x', y: 'y', r: 'size', fill: 'nlp', stroke: '#fff', strokeWidth: 1, tip: true, title: function (n) { return n.tip; } }),
          P.text(nodes, { x: 'x', y: 'y', text: function (n) { return trunc(n.label, 36); }, dy: -14, fontSize: 10.5, fill: 'currentColor', stroke: 'var(--xq-bg, #fff)', strokeWidth: 3, paintOrder: 'stroke' })
        ]
      });
    }, { minHeight: 300 });
  };

  C.trunc = trunc;
  C.empty = empty;
})(typeof window !== 'undefined' ? window : globalThis);
