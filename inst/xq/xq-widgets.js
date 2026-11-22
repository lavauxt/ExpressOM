/* xq-widgets.js -- ExpressOM Quarto dashboard runtime: UI building blocks.
 *
 *  XQ.lazy        mount a heavy node only when it first becomes visible
 *  XQ.note        "How to read this" callout shown at the top of every tab
 *  XQ.dataTable   virtualised table: sort, search, paging, CSV export
 *  XQ.filterBox   SQL WHERE box (JS engine, DuckDB-Wasm when loaded)
 *  XQ.sqlConsole  free SQL over the embedded tables (DuckDB-Wasm, lazy)
 *  XQ.frame       titled figure frame with SVG/PNG download
 *  XQ.stats       row of summary tiles
 */
(function (g) {
  'use strict';
  const XQ = g.XQ;
  const doc = g.document;
  const F = XQ.fmt;

  // ------------------------------------------------------------------ DOM kit
  function el(tag, cls, text) {
    const e = doc.createElement(tag);
    if (cls) e.className = cls;
    if (text != null) e.textContent = text;
    return e;
  }
  const SVGNS = 'http://www.w3.org/2000/svg';
  function sv(tag, attrs) {
    const e = doc.createElementNS(SVGNS, tag);
    if (attrs) for (const k in attrs) e.setAttribute(k, attrs[k]);
    return e;
  }
  XQ.el = el;
  XQ.sv = sv;

  // --------------------------------------------------------------------- lazy
  /**
   * Returns a placeholder that builds its content the first time it is on
   * screen (inside an active tab / page, within `rootMargin` of the viewport).
   * Hidden tabs therefore cost nothing at load. `build(host)` may return a
   * node, a promise of a node, or nothing (if it filled `host` itself).
   */
  XQ.lazy = function (build, opts) {
    opts = opts || {};
    const host = el('div', 'xq-lazy');
    if (opts.minHeight) host.style.minHeight = opts.minHeight + 'px';
    host.appendChild(el('div', 'xq-skel', opts.label || 'Loading…'));
    let done = false;
    const run = function () {
      if (done) return;
      done = true;
      Promise.resolve().then(function () { return build(host); }).then(function (node) {
        if (node && node !== host) host.replaceChildren(node); else if (node === undefined && host.firstChild && host.firstChild.className === 'xq-skel') host.firstChild.remove();
        if (opts.minHeight) host.style.minHeight = '';
        if (typeof opts.after === 'function') opts.after(host);
      }).catch(function (e) {
        console.error('XQ.lazy', e);
        const err = el('div', 'xq-err', 'Could not render this panel: ' + (e && e.message ? e.message : e));
        host.replaceChildren(err);
      });
    };
    if (!('IntersectionObserver' in g)) { setTimeout(run, 0); return host; }
    const io = new IntersectionObserver(function (es) {
      if (es.some(function (x) { return x.isIntersecting; })) { io.disconnect(); run(); }
    }, { rootMargin: opts.rootMargin || '120px' });
    io.observe(host);
    host._xqDispose = function () { io.disconnect(); };
    return host;
  };

  /** Dispose anything with destroy()/_xqDispose when an OJS cell re-runs. */
  XQ.own = function (invalidation, node) {
    if (invalidation && invalidation.then) {
      invalidation.then(function () {
        try { if (node && node.destroy) node.destroy(); if (node && node._xqDispose) node._xqDispose(); } catch (e) { /* noop */ }
      });
    }
    return node;
  };

  // --------------------------------------------------------------------- note
  /** Short "how to interpret" callout. `html` is trusted markup from R. */
  XQ.note = function (title, html, open) {
    const d = el('details', 'xq-help');
    if (open !== false) d.open = true;
    const s = el('summary');
    s.innerHTML = '<span class="xq-help-ic" aria-hidden="true">i</span>' + F.esc(title || 'How to read this');
    const b = el('div', 'xq-help-body');
    b.innerHTML = html;
    d.append(s, b);
    return d;
  };

  // ------------------------------------------------------------------ stats
  XQ.stats = function (items) {
    const row = el('div', 'xq-stats');
    items.forEach(function (it) {
      const t = el('div', 'xq-stat' + (it.cls ? ' ' + it.cls : ''));
      t.append(el('div', 'xq-stat-v', it.value == null ? '—' : String(it.value)), el('div', 'xq-stat-l', it.label));
      if (it.sub) t.append(el('div', 'xq-stat-s', it.sub));
      row.append(t);
    });
    return row;
  };

  // --------------------------------------------------------------- data table
  function intersect(a, b, n) {
    if (!a) return b;
    if (!b) return a;
    const mask = new Uint8Array(n);
    for (let i = 0; i < b.length; i++) mask[b[i]] = 1;
    const out = [];
    for (let i = 0; i < a.length; i++) if (mask[a[i]]) out.push(a[i]);
    return Uint32Array.from(out);
  }

  const URL_RE = /^(https?:\/\/|\.{0,2}\/|[A-Za-z0-9_\-./]+\.(html|pdf|png|csv|tsv|txt)$)/;

  /**
   * Virtualised table (only one page of rows is ever in the DOM).
   * opts: columns, idx (Uint32Array subset), sort, desc, pageSize, search,
   *       filename, labels {col: header}, onRow(i), maxChars, empty
   */
  XQ.dataTable = function (tbl, opts) {
    const o = opts || {};
    const cols = (o.columns || tbl.names).filter(function (c) { return tbl.has(c); });
    const st = { idx: o.idx || null, sort: o.sort || null, desc: !!o.desc, page: 0, size: o.pageSize || 15, q: '' };
    const root = el('div', 'xq-dt');
    if (!tbl.n || !cols.length) {
      root.append(el('div', 'xq-empty', o.empty || 'No rows to display.'));
      return root;
    }
    const bar = el('div', 'xq-dt-bar');
    let search = null;
    if (o.search !== false) {
      search = el('input', 'xq-dt-search');
      search.type = 'search';
      search.placeholder = 'Search text…';
      search.setAttribute('aria-label', 'Search the table');
      bar.append(search);
    }
    const info = el('span', 'xq-dt-info');
    const sizeSel = el('select', 'xq-dt-size');
    [10, 15, 25, 50, 100, 250].forEach(function (n) {
      const op = el('option', null, n + ' rows'); op.value = n; if (n === st.size) op.selected = true; sizeSel.append(op);
    });
    const csvBtn = el('button', 'xq-btn', 'CSV');
    csvBtn.type = 'button';
    csvBtn.title = 'Download the filtered rows (all pages) as CSV';
    bar.append(info, sizeSel, csvBtn);

    const wrap = el('div', 'xq-dt-wrap');
    const table = el('table', 'xq-dt-table');
    const thead = el('thead');
    const tr = el('tr');
    const ths = {};
    cols.forEach(function (c) {
      const th = el('th');
      th.textContent = (o.labels && o.labels[c]) || c;
      th.title = 'Click to sort';
      th.dataset.col = c;
      if (tbl.isNumeric(c)) th.className = 'xq-num';
      th.addEventListener('click', function () {
        if (st.sort === c) st.desc = !st.desc; else { st.sort = c; st.desc = tbl.isNumeric(c) && !tbl.isP(c); }
        st.page = 0; compute(); render();
      });
      ths[c] = th; tr.append(th);
    });
    thead.append(tr);
    const tbody = el('tbody');
    table.append(thead, tbody);
    wrap.append(table);
    const pager = el('div', 'xq-dt-pager');
    root.append(bar, wrap, pager);

    let view = null;
    function compute() {
      let base = st.idx;
      if (st.q) base = intersect(base, tbl.search(st.q, o.searchCols), tbl.n);
      view = st.sort ? tbl.order(st.sort, st.desc, base) : (base || tbl.all());
    }
    function fmtCell(td, c, i) {
      const t = tbl.type(c);
      if (t === 'nlp') { const v = tbl.nlp(c)[i]; td.textContent = F.nlp(v); td.className = 'xq-num'; return; }
      if (t === 'f32' || t === 'f64') { const v = tbl.num(c)[i]; td.textContent = F.num(v, 3); td.className = 'xq-num'; return; }
      if (t === 'i32') { const v = tbl.num(c)[i]; td.textContent = v !== v ? '' : String(v); td.className = 'xq-num'; return; }
      const v = tbl.str(c)[i];
      if (v == null) { td.className = 'xq-na'; td.textContent = 'NA'; return; }
      if (URL_RE.test(v) && v.length < 400) {
        const a = el('a', null, v.length > 60 ? '…' + v.slice(-58) : v);
        a.href = v; a.target = '_blank'; a.rel = 'noopener'; td.append(a); return;
      }
      const max = o.maxChars || 90;
      if (v.length > max) { td.textContent = v.slice(0, max - 1) + '…'; td.title = v.length > 1200 ? v.slice(0, 1200) + '…' : v; }
      else td.textContent = v;
    }
    function render() {
      const n = view.length;
      const pages = Math.max(1, Math.ceil(n / st.size));
      if (st.page >= pages) st.page = pages - 1;
      const from = st.page * st.size, to = Math.min(n, from + st.size);
      tbody.replaceChildren();
      const frag = doc.createDocumentFragment();
      for (let k = from; k < to; k++) {
        const i = view[k];
        const r = el('tr');
        if (o.onRow) { r.className = 'xq-click'; r.addEventListener('click', function () { o.onRow(i); }); }
        cols.forEach(function (c) { const td = el('td'); fmtCell(td, c, i); r.append(td); });
        frag.append(r);
      }
      if (!n) { const r = el('tr'); const td = el('td', 'xq-empty'); td.colSpan = cols.length; td.textContent = 'No matching rows.'; r.append(td); frag.append(r); }
      tbody.append(frag);
      Object.keys(ths).forEach(function (c) {
        ths[c].classList.toggle('xq-sorted', st.sort === c);
        ths[c].dataset.dir = st.sort === c ? (st.desc ? '▼' : '▲') : '';
      });
      info.textContent = n ? (from + 1).toLocaleString('en-US') + '–' + to.toLocaleString('en-US') + ' of ' + n.toLocaleString('en-US') + (n !== tbl.n ? ' (filtered from ' + tbl.n.toLocaleString('en-US') + ')' : '') : '0 rows';
      pager.replaceChildren();
      const mk = function (label, page, disabled) {
        const b = el('button', 'xq-btn', label); b.type = 'button'; b.disabled = disabled;
        b.addEventListener('click', function () { st.page = page; render(); });
        return b;
      };
      pager.append(mk('«', 0, st.page === 0), mk('‹', st.page - 1, st.page === 0),
        el('span', 'xq-dt-pg', 'Page ' + (st.page + 1) + ' / ' + pages),
        mk('›', st.page + 1, st.page >= pages - 1), mk('»', pages - 1, st.page >= pages - 1));
    }
    if (search) search.addEventListener('input', XQ.debounce(function () { st.q = search.value; st.page = 0; compute(); render(); }, 180));
    sizeSel.addEventListener('change', function () { st.size = +sizeSel.value; st.page = 0; render(); });
    csvBtn.addEventListener('click', function () { XQ.download(o.filename || 'expressom_table.csv', tbl.csv(view, cols)); });
    root.update = function (idx) { st.idx = idx || null; st.page = 0; compute(); render(); };
    root.count = function () { return view.length; };
    compute();
    render();
    return root;
  };

  // --------------------------------------------------------------- filter box
  /**
   * SQL WHERE box. `.value` is a Uint32Array of matching row indices (or null
   * when the box is empty), updated on Enter / button; fires an "input" event
   * so it works with Observable's `viewof`.
   */
  XQ.filterBox = function (name, tbl, opts) {
    const o = opts || {};
    const root = el('div', 'xq-fb');
    root.value = null;
    const row = el('div', 'xq-fb-row');
    const input = el('input', 'xq-fb-in');
    input.type = 'text';
    input.spellcheck = false;
    input.placeholder = o.placeholder || "SQL WHERE, e.g.  padj < 0.01 AND abs(log2FoldChange) > 1 AND gene LIKE 'HOX%'";
    input.setAttribute('aria-label', 'SQL WHERE filter');
    if (o.value) input.value = o.value;
    const run = el('button', 'xq-btn xq-btn-p', 'Filter'); run.type = 'button';
    const clr = el('button', 'xq-btn', 'Clear'); clr.type = 'button';
    const eng = el('select', 'xq-fb-eng');
    [['auto', 'Engine: auto'], ['js', 'Engine: JS (typed arrays)'], ['duckdb', 'Engine: DuckDB-Wasm']].forEach(function (p) {
      const op = el('option', null, p[1]); op.value = p[0]; eng.append(op);
    });
    row.append(input, run, clr, eng);
    const status = el('div', 'xq-fb-status');
    root.append(row);
    if (o.examples && o.examples.length) {
      const ex = el('div', 'xq-fb-ex');
      ex.append(el('span', null, 'Examples: '));
      o.examples.forEach(function (e) {
        const b = el('button', 'xq-chip', e); b.type = 'button';
        b.addEventListener('click', function () { input.value = e; go(); });
        ex.append(b);
      });
      root.append(ex);
    }
    root.append(status);
    const duckLabel = function () { return XQ.duck.state === 'ready' ? 'DuckDB-Wasm ready' : XQ.duck.state === 'loading' ? 'DuckDB-Wasm loading…' : XQ.duck.state === 'failed' ? 'DuckDB-Wasm unavailable (offline?) – JS engine in use' : 'DuckDB-Wasm loads on demand'; };
    let seq = 0, last = '';
    async function go() {
      const my = ++seq;
      const q = input.value.trim();
      status.className = 'xq-fb-status';
      try {
        const mode = eng.value;
        let r;
        if (!q) r = { idx: null, ms: 0, engine: 'js' };
        else if (mode === 'js') r = XQ.where(tbl, q);
        else if (mode === 'duckdb' || (mode === 'auto' && (XQ.duck.state === 'ready' || tbl.n > 400000))) {
          status.textContent = 'Running on DuckDB-Wasm…';
          const ok = await XQ.duck.ensure();
          r = ok ? await XQ.duck.where(name, tbl, q) : XQ.where(tbl, q);
        } else r = XQ.where(tbl, q);
        if (my !== seq) return;
        root.value = q ? r.idx : null;
        last = q ? r.engine : '';
        status.textContent = q
          ? r.idx.length.toLocaleString('en-US') + ' of ' + tbl.n.toLocaleString('en-US') + ' rows match · ' + r.ms.toFixed(r.ms < 10 ? 1 : 0) + ' ms · ' + (r.engine === 'duckdb' ? 'DuckDB-Wasm' : 'JS engine') + ' · ' + duckLabel()
          : 'Showing all ' + tbl.n.toLocaleString('en-US') + ' rows · ' + duckLabel();
        root.dispatchEvent(new Event('input', { bubbles: true }));
      } catch (e) {
        if (my !== seq) return;
        status.className = 'xq-fb-status xq-err';
        status.textContent = 'Filter error: ' + (e && e.message ? e.message : e);
      }
    }
    input.addEventListener('keydown', function (e) { if (e.key === 'Enter') go(); });
    run.addEventListener('click', go);
    clr.addEventListener('click', function () { input.value = ''; go(); });
    eng.addEventListener('change', function () { if (eng.value === 'duckdb') XQ.duck.ensure().then(go); else go(); });
    XQ.duck.on(function () { if (last || !input.value) status.textContent = status.textContent.replace(/ · DuckDB-Wasm[^·]*$| · JS engine.*$/, '') + (input.value ? '' : ' · ' + duckLabel()); });
    status.textContent = 'Showing all ' + tbl.n.toLocaleString('en-US') + ' rows · ' + duckLabel();
    if (o.value) setTimeout(go, 0);
    return root;
  };

  // -------------------------------------------------------------- SQL console
  XQ.sqlConsole = function (tables, opts) {
    const o = opts || {};
    const names = Object.keys(tables);
    const root = el('div', 'xq-sql');
    const ta = el('textarea', 'xq-sql-ta');
    ta.rows = 4;
    ta.spellcheck = false;
    ta.value = o.sql || ('SELECT * FROM ' + names[0] + ' LIMIT 50');
    const bar = el('div', 'xq-sql-bar');
    const run = el('button', 'xq-btn xq-btn-p', 'Run SQL'); run.type = 'button';
    const hint = el('span', 'xq-sql-hint', 'Tables: ' + names.join(', ') + ' · DuckDB-Wasm loads on first run (needs internet or a local duckdb_url)');
    bar.append(run, hint);
    const out = el('div', 'xq-sql-out');
    root.append(ta, bar, out);
    run.addEventListener('click', async function () {
      out.replaceChildren(el('div', 'xq-skel', 'Starting DuckDB-Wasm…'));
      try {
        const ok = await XQ.duck.ensure();
        if (!ok) { out.replaceChildren(el('div', 'xq-err', 'DuckDB-Wasm could not be loaded (' + (XQ.duck.error || 'offline') + '). The WHERE filters above keep working with the built-in JS engine.')); return; }
        const sql = ta.value;
        for (let k = 0; k < names.length; k++) if (new RegExp('\\b' + names[k] + '\\b').test(sql)) await XQ.duck.register(names[k], tables[names[k]]);
        const r = await XQ.duck.query(sql);
        const meta = el('div', 'xq-fb-status', r.table.n.toLocaleString('en-US') + ' rows · ' + r.ms.toFixed(0) + ' ms · DuckDB-Wasm');
        out.replaceChildren(meta, XQ.dataTable(r.table, { pageSize: 15, filename: 'sql_result.csv' }));
      } catch (e) {
        out.replaceChildren(el('div', 'xq-err', 'SQL error: ' + (e && e.message ? e.message : e)));
      }
    });
    return root;
  };

  // --------------------------------------------------------- figure frame/PNG
  function serializeSvg(svg) {
    const clone = svg.cloneNode(true);
    const cs = g.getComputedStyle(svg);
    clone.setAttribute('xmlns', SVGNS);
    clone.style.color = cs.color;
    clone.style.fontFamily = cs.fontFamily;
    if (!clone.getAttribute('width')) clone.setAttribute('width', svg.getBoundingClientRect().width);
    if (!clone.getAttribute('height')) clone.setAttribute('height', svg.getBoundingClientRect().height);
    return new XMLSerializer().serializeToString(clone);
  }
  XQ.svgToPng = async function (svg, filename, scale) {
    const w = Math.ceil(svg.getBoundingClientRect().width), h = Math.ceil(svg.getBoundingClientRect().height);
    const img = new Image();
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(serializeSvg(svg));
    await img.decode();
    const c = doc.createElement('canvas');
    const k = scale || 2;
    c.width = w * k; c.height = h * k;
    const ctx = c.getContext('2d');
    ctx.scale(k, k);
    const bg = g.getComputedStyle(doc.body).backgroundColor;
    ctx.fillStyle = bg && bg !== 'rgba(0, 0, 0, 0)' ? bg : '#ffffff';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);
    c.toBlob(function (b) {
      const a = doc.createElement('a'); a.href = URL.createObjectURL(b); a.download = filename || 'figure.png';
      doc.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
    });
  };

  /** Titled frame around a node. SVG charts get PNG/SVG buttons (the SVG is
   *  looked up at click time because charts render lazily). o.download=false
   *  hides them (WebGL scatterplots bring their own PNG button). */
  XQ.frame = function (node, o) {
    o = o || {};
    const fig = el('figure', 'xq-fig');
    if (o.title) fig.append(el('figcaption', 'xq-fig-t', o.title));
    const body = el('div', 'xq-fig-b');
    body.append(node);
    fig.append(body);
    if (o.download !== false) {
      const findSvg = function () {
        if (node.tagName && node.tagName.toLowerCase() === 'svg') return node;
        if (!node.querySelectorAll) return null;
        let best = null, bw = -1; // the plot itself, not a legend swatch
        node.querySelectorAll('svg:not(.xq-sc-ax)').forEach(function (v) { const w = v.getBoundingClientRect().width; if (w > bw) { bw = w; best = v; } });
        return best;
      };
      const tools = el('div', 'xq-tools xq-fig-tools');
      const png = el('button', 'xq-btn', 'PNG'); png.type = 'button';
      png.addEventListener('click', function () { const s = findSvg(); if (s) XQ.svgToPng(s, (o.filename || 'expressom_figure') + '.png'); });
      const sv2 = el('button', 'xq-btn', 'SVG'); sv2.type = 'button';
      sv2.addEventListener('click', function () { const s = findSvg(); if (s) XQ.download((o.filename || 'expressom_figure') + '.svg', serializeSvg(s), 'image/svg+xml'); });
      tools.append(png, sv2);
      fig.append(tools);
    }
    if (o.caption) fig.append(el('div', 'xq-fig-c', o.caption));
    return fig;
  };

  // ---------------------------------------------------- static figure lightbox
  // Clicking any <img class="xq-zoom"> opens it full size; Esc / click closes.
  if (doc) {
    doc.addEventListener('click', function (e) {
      const t = e.target;
      if (!t || !t.classList || !t.classList.contains('xq-zoom')) return;
      const ov = el('div', 'xq-lightbox');
      const im = el('img'); im.src = t.currentSrc || t.src; im.alt = t.alt || '';
      ov.append(im);
      const close = function () { ov.remove(); doc.removeEventListener('keydown', onKey); };
      const onKey = function (k) { if (k.key === 'Escape') close(); };
      ov.addEventListener('click', close);
      doc.addEventListener('keydown', onKey);
      doc.body.append(ov);
    });
  }
})(typeof window !== 'undefined' ? window : globalThis);
