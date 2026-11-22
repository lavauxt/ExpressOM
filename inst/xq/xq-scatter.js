/* xq-scatter.js -- WebGL scatter plot (regl-scatterplot) with SVG axes.
 *
 *   const el = XQ.scatter({height: 460, tip: i => 'html'});
 *   el.api.set({x, y, cls, palette, sizes, alphas, xlabel, ylabel, ...});
 *
 * 100k+ points pan/zoom smoothly because points are drawn by the GPU; the
 * SVG layer only holds axes, threshold lines and a handful of labels.
 * The GL context is created the first time the plot is visible and live
 * contexts are capped (browsers allow ~16 per page).
 *
 * `el.value` is the Uint32Array of lasso-selected points (original indices);
 * an "input" event is dispatched so `viewof` works in Observable JS.
 */
(function (g) {
  'use strict';
  const XQ = g.XQ;
  const doc = g.document;
  const el = XQ.el, sv = XQ.sv, F = XQ.fmt;
  const MAX_LIVE = 10;
  const live = [];

  function hexToRgba(h, a) {
    const c = doc.createElement('canvas').getContext('2d');
    c.fillStyle = '#000';
    c.fillStyle = h; // CSS colour names ("red2", "royalblue") work too
    const v = c.fillStyle; // normalised to #rrggbb or rgba(...)
    let r, gg, b;
    if (v[0] === '#') { r = parseInt(v.slice(1, 3), 16); gg = parseInt(v.slice(3, 5), 16); b = parseInt(v.slice(5, 7), 16); }
    else { const m = v.match(/[\d.]+/g); r = +m[0]; gg = +m[1]; b = +m[2]; }
    return [r / 255, gg / 255, b / 255, a == null ? 1 : a];
  }
  XQ.cssColor = function (h) { // -> '#rrggbb'
    const c = doc.createElement('canvas').getContext('2d');
    c.fillStyle = '#000'; c.fillStyle = h; const v = c.fillStyle;
    if (v[0] === '#') return v;
    const m = v.match(/[\d.]+/g);
    return '#' + [0, 1, 2].map(function (k) { return ('0' + Math.round(+m[k]).toString(16)).slice(-2); }).join('');
  };

  function niceExtent(arr, padFrac) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < arr.length; i++) { const v = arr[i]; if (v === v && isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } }
    if (!isFinite(lo)) { lo = 0; hi = 1; }
    if (lo === hi) { lo -= 0.5; hi += 0.5; }
    const p = (hi - lo) * (padFrac == null ? 0.04 : padFrac);
    return [lo - p, hi + p];
  }

  XQ.scatter = function (o) {
    o = o || {};
    const H = o.height || 460;
    const M = Object.assign({ l: 58, r: 16, t: 14, b: 46 }, o.margin || {});
    const root = el('div', 'xq-sc');
    root.style.height = H + 'px';
    root.value = new Uint32Array(0);
    const plot = el('div', 'xq-sc-plot');
    const canvas = el('canvas');
    plot.append(canvas);
    const svg = sv('svg', { class: 'xq-sc-ax' });
    const tip = el('div', 'xq-sc-tip');
    const legend = el('div', 'xq-sc-legend');
    const tools = el('div', 'xq-tools xq-sc-tools');
    const bReset = el('button', 'xq-btn', 'Reset view'); bReset.type = 'button';
    const bLasso = el('button', 'xq-btn', 'Lasso'); bLasso.type = 'button';
    bLasso.title = 'Toggle lasso selection (drag to select points; also Alt+drag)';
    const bClear = el('button', 'xq-btn', 'Clear'); bClear.type = 'button';
    const bPng = el('button', 'xq-btn', 'PNG'); bPng.type = 'button';
    const selInfo = el('span', 'xq-sc-sel');
    tools.append(selInfo, bLasso, bClear, bReset, bPng);
    if (o.selectable === false) { bLasso.style.display = 'none'; bClear.style.display = 'none'; }
    root.append(plot, svg, legend, tip, tools);
    const msg = el('div', 'xq-sc-msg');
    root.append(msg);

    // ---- state
    let spec = null;        // last full spec
    let sp = null;          // regl-scatterplot instance
    let visible = false;
    let W = 0, pw = 0, ph = 0;
    let perm = null;        // drawn position -> original index
    let xs = null, ys = null; // d3 scales owned by regl (data domain <-> pixels)
    let dom = null;         // {x:[a,b], y:[c,d]} data domains of the full view
    let drawn = null;       // arrays currently on the GPU
    let hiddenCls = {};
    let lassoOn = false;
    let mouse = null;
    let destroyed = false;
    const api = { el: root };

    function setSel(points) {
      const out = new Uint32Array(points.length);
      for (let k = 0; k < points.length; k++) out[k] = perm ? perm[points[k]] : points[k];
      root.value = out;
      selInfo.textContent = out.length ? out.length.toLocaleString('en-US') + ' selected' : '';
      root.dispatchEvent(new Event('input', { bubbles: true }));
    }

    // ---- axes
    function drawAxes() {
      if (!spec || !xs || !pw) return;
      const d = XQV().d3;
      const xd = xs.domain(), yd = ys.domain();
      const X = d.scaleLinear().domain(xd).range([M.l, M.l + pw]);
      const Y = d.scaleLinear().domain(yd).range([M.t + ph, M.t]);
      svg.replaceChildren();
      const defs = sv('defs');
      const cp = sv('clipPath', { id: 'xqclip' + (root._id = root._id || Math.random().toString(36).slice(2)) });
      cp.append(sv('rect', { x: M.l, y: M.t, width: pw, height: ph }));
      defs.append(cp);
      svg.append(defs);
      const gx = sv('g', { class: 'xq-ax-x', transform: 'translate(0,' + (M.t + ph) + ')' });
      const gy = sv('g', { class: 'xq-ax-y', transform: 'translate(' + M.l + ',0)' });
      svg.append(gx, gy);
      const xf = spec.xlog ? function (v) { return F.num(Math.pow(10, v), 2); } : null;
      const yf = spec.ylog ? function (v) { return F.num(Math.pow(10, v), 2); } : null;
      const ax = d.axisBottom(X).ticks(Math.max(3, Math.floor(pw / 90))).tickSize(-ph);
      const ay = d.axisLeft(Y).ticks(Math.max(3, Math.floor(ph / 55))).tickSize(-pw);
      if (xf) ax.tickFormat(xf);
      if (yf) ay.tickFormat(yf);
      d.select(gx).call(ax);
      d.select(gy).call(ay);
      svg.append(sv('rect', { x: M.l, y: M.t, width: pw, height: ph, class: 'xq-ax-box' }));
      const tx = sv('text', { class: 'xq-ax-t', x: M.l + pw / 2, y: H - 8, 'text-anchor': 'middle' }); tx.textContent = spec.xlabel || '';
      const ty = sv('text', { class: 'xq-ax-t', transform: 'translate(14,' + (M.t + ph / 2) + ') rotate(-90)', 'text-anchor': 'middle' }); ty.textContent = spec.ylabel || '';
      svg.append(tx, ty);
      const clip = sv('g', { 'clip-path': 'url(#' + cp.id + ')' });
      svg.append(clip);
      (spec.hlines || []).forEach(function (l) {
        const y = Y(l.v);
        clip.append(sv('line', { x1: M.l, x2: M.l + pw, y1: y, y2: y, class: 'xq-thr', 'stroke-dasharray': l.solid ? '' : '5 4' }));
      });
      (spec.vlines || []).forEach(function (l) {
        const x = X(l.v);
        clip.append(sv('line', { x1: x, x2: x, y1: M.t, y2: M.t + ph, class: 'xq-thr', 'stroke-dasharray': l.solid ? '' : '5 4' }));
      });
      if (spec.diag) {
        const a = Math.max(xd[0], yd[0]), b = Math.min(xd[1], yd[1]);
        if (a < b) clip.append(sv('line', { x1: X(a), y1: Y(a), x2: X(b), y2: Y(b), class: 'xq-thr', 'stroke-dasharray': '5 4' }));
      }
      // labelled points + highlighted ring
      (spec.ring || []).forEach(function (i) {
        const x = X(spec.x[i]), y = Y(spec.y[i]);
        if (x >= M.l && x <= M.l + pw && y >= M.t && y <= M.t + ph) clip.append(sv('circle', { cx: x, cy: y, r: 7, class: 'xq-ring' }));
      });
      const lab = sv('g', { 'clip-path': 'url(#' + cp.id + ')', class: 'xq-labels' });
      svg.append(lab);
      const placed = [];
      (spec.labels || []).forEach(function (m) {
        const x = X(spec.x[m.i]), y = Y(spec.y[m.i]);
        if (!(x >= M.l && x <= M.l + pw && y >= M.t && y <= M.t + ph)) return;
        const w = m.text.length * 6 + 4;
        for (let k = 0; k < placed.length; k++) { const p = placed[k]; if (Math.abs(p[0] - x) < (w + p[2]) / 2 && Math.abs(p[1] - y) < 12) return; }
        placed.push([x, y, w]);
        const t = sv('text', { x: x + 6, y: y - 5, class: 'xq-lab' });
        t.textContent = m.text;
        lab.append(t);
      });
    }

    function showTip(orig) {
      if (!spec || !spec.tip && !o.tip) { return; }
      const html = (spec.tip || o.tip)(orig);
      if (!html) { tip.style.display = 'none'; return; }
      tip.innerHTML = html;
      tip.style.display = 'block';
      placeTip();
    }
    function placeTip() {
      if (!mouse) return;
      const r = root.getBoundingClientRect();
      let x = mouse.clientX - r.left + 14, y = mouse.clientY - r.top + 14;
      const tw = tip.offsetWidth, th = tip.offsetHeight;
      if (x + tw > r.width - 4) x = mouse.clientX - r.left - tw - 14;
      if (y + th > r.height - 4) y = mouse.clientY - r.top - th - 14;
      tip.style.left = Math.max(4, x) + 'px';
      tip.style.top = Math.max(4, y) + 'px';
    }

    function drawLegend() {
      legend.replaceChildren();
      if (!spec || !spec.legend || !spec.legend.length) return;
      spec.legend.forEach(function (it, k) {
        const row = el('button', 'xq-lg' + (hiddenCls[it.cls != null ? it.cls : k] ? ' xq-off' : ''));
        row.type = 'button';
        const dot = el('span', 'xq-lg-dot'); dot.style.background = it.color;
        row.append(dot, el('span', null, it.label + (it.n != null ? ' (' + it.n.toLocaleString('en-US') + ')' : '')));
        row.title = 'Click to show / hide';
        row.addEventListener('click', function () {
          const c = it.cls != null ? it.cls : k;
          hiddenCls[c] = !hiddenCls[c];
          pushDraw(false);
          drawLegend();
        });
        legend.append(row);
      });
    }

    // ---- GL
    function init() {
      if (sp || destroyed || !spec) return;
      W = root.clientWidth;
      if (!W) return;
      pw = Math.max(120, W - M.l - M.r);
      ph = Math.max(100, H - M.t - M.b);
      plot.style.cssText = 'left:' + M.l + 'px;top:' + M.t + 'px;width:' + pw + 'px;height:' + ph + 'px';
      svg.setAttribute('width', W); svg.setAttribute('height', H);
      // free GL contexts of plots that are scrolled / tabbed away
      while (live.length >= MAX_LIVE) {
        const victim = live.find(function (a) { return !a._visible(); });
        if (!victim) break;
        victim._sleep();
      }
      const d = XQV().d3;
      xs = d.scaleLinear().domain(dom.x);
      ys = d.scaleLinear().domain(dom.y);
      try {
        sp = XQV().createScatterplot({
          canvas: canvas, width: pw, height: ph, xScale: xs, yScale: ys, aspectRatio: pw / ph,
          backgroundColor: [0, 0, 0, 0], pointSize: o.pointSize || 3, opacity: 0.85,
          pointColor: [[0.5, 0.5, 0.5, 1]], colorBy: 'valueA', sizeBy: 'valueA', opacityBy: 'valueA',
          pointScaleMode: 'constant', lassoColor: [0.1, 0.45, 0.85, 1], lassoMinDelay: 10,
          mouseMode: lassoOn ? 'lasso' : 'panZoom', lassoOnLongPress: false, actionKeyMap: { lasso: 'alt' }
        });
      } catch (e) {
        msg.textContent = 'WebGL is not available in this browser, so this plot cannot be drawn. (' + (e && e.message ? e.message : e) + ')';
        msg.style.display = 'block';
        sp = null;
        return;
      }
      live.push(api);
      sp.subscribe('view', function (ev) { if (ev.xScale) { xs = ev.xScale; ys = ev.yScale; } drawAxes(); });
      sp.subscribe('pointover', function (p) { showTip(perm ? perm[p] : p); });
      sp.subscribe('pointout', function () { tip.style.display = 'none'; });
      sp.subscribe('select', function (ev) { setSel(ev.points); });
      sp.subscribe('deselect', function () { setSel([]); });
      pushDraw(true);
    }

    function pushDraw(resetView) {
      if (!sp || !spec) return;
      const n = spec.x.length;
      const cls = spec.cls;
      const K = spec.palette ? spec.palette.length : 1;
      const sx0 = spec.x, sy0 = spec.y;
      // draw order: class 0 first (background), higher classes on top
      const cnt = new Uint32Array(K);
      const fin = function (i) { const a = sx0[i], b = sy0[i]; return a === a && b === b && isFinite(a) && isFinite(b); };
      for (let i = 0; i < n; i++) { const c = cls ? cls[i] : 0; if (!hiddenCls[c] && fin(i)) cnt[c]++; }
      const pos = new Uint32Array(K);
      let used = 0;
      for (let k = 0; k < K; k++) { pos[k] = used; used += cnt[k]; }
      const X = new Float32Array(used), Y = new Float32Array(used), V = new Float32Array(used);
      perm = new Uint32Array(used);
      const dx = dom.x, dy = dom.y;
      const sx = 2 / (dx[1] - dx[0]), sy = 2 / (dy[1] - dy[0]);
      for (let i = 0; i < n; i++) {
        const c = cls ? cls[i] : 0;
        if (hiddenCls[c] || !fin(i)) continue;
        const p = pos[c]++;
        perm[p] = i;
        X[p] = (sx0[i] - dx[0]) * sx - 1;
        Y[p] = (sy0[i] - dy[0]) * sy - 1;
        V[p] = c;
      }
      const colors = (spec.palette || ['#888888']).map(function (h) { return hexToRgba(h, 1); });
      const sizes = spec.sizes || colors.map(function () { return o.pointSize || 3; });
      const alphas = spec.alphas || colors.map(function () { return 0.8; });
      sp.set({ pointColor: colors, pointSize: sizes, opacity: alphas });
      drawn = { n: used };
      if (resetView) {
        // regl derives the visible domain from the scales it was given, so a new
        // geometry (e.g. volcano -> MA) needs fresh scales built on the new domain
        const d = XQV().d3;
        xs = d.scaleLinear().domain(dom.x);
        ys = d.scaleLinear().domain(dom.y);
        try { sp.set({ xScale: xs, yScale: ys }); sp.reset({ preventEvent: true }); } catch (e) { /* noop */ }
      }
      const maybe = sp.draw({ x: X, y: Y, valueA: V }, { zDataType: 'categorical', transition: false });
      Promise.resolve(maybe).then(function () { if (!destroyed) drawAxes(); });
      drawAxes();
    }

    api._visible = function () { return visible; };
    api._sleep = function () {
      if (!sp) return;
      try { sp.destroy(); } catch (e) { /* noop */ }
      sp = null;
      const k = live.indexOf(api); if (k >= 0) live.splice(k, 1);
    };

    /**
     * Update the plot. Pass a full spec once (x, y, ...); later calls may
     * pass only what changed (e.g. {cls, legend, labels, ring}).
     */
    api.set = function (s) {
      const prev = spec;
      spec = Object.assign({}, prev || {}, s);
      const geometry = !prev || s.x || s.y || s.xlog || s.ylog;
      if (geometry && spec.x) {
        dom = { x: s.xdomain || niceExtent(spec.x), y: s.ydomain || niceExtent(spec.y) };
        if (spec.diag) { const lo = Math.min(dom.x[0], dom.y[0]), hi = Math.max(dom.x[1], dom.y[1]); dom = { x: [lo, hi], y: [lo, hi] }; }
      }
      if (s.cls && !s.legend) { /* keep legend */ }
      drawLegend();
      if (!sp) { if (visible) init(); return api; }
      if (geometry) hiddenCls = hiddenCls; // keep hidden classes across views
      pushDraw(!!geometry);
      return api;
    };
    api.reset = function () { if (!sp) return; sp.reset(); drawAxes(); };
    api.destroy = function () {
      destroyed = true;
      if (io) io.disconnect();
      api._sleep();
    };
    api.selected = function () { return root.value; };
    api.clearSelection = function () { if (sp) sp.deselect(); setSel([]); };
    api.png = function (name) {
      if (!sp) return;
      const dpr = Math.max(1, g.devicePixelRatio || 1);
      const comp = doc.createElement('canvas');
      comp.width = W * 2; comp.height = H * 2;
      const ctx = comp.getContext('2d');
      ctx.scale(2, 2);
      const bg = g.getComputedStyle(doc.body).backgroundColor;
      ctx.fillStyle = bg && bg !== 'rgba(0, 0, 0, 0)' ? bg : '#ffffff';
      ctx.fillRect(0, 0, W, H);
      const done = function () {
        const clone = svg.cloneNode(true);
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        clone.style.color = g.getComputedStyle(root).color;
        clone.style.fontFamily = g.getComputedStyle(root).fontFamily;
        const img = new Image();
        img.onload = function () {
          ctx.drawImage(img, 0, 0, W, H);
          comp.toBlob(function (b) {
            const a = doc.createElement('a'); a.href = URL.createObjectURL(b); a.download = (name || o.filename || 'expressom_scatter') + '.png';
            doc.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
          });
        };
        img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(new XMLSerializer().serializeToString(clone));
      };
      try {
        const im = sp.export();
        const tmp = doc.createElement('canvas'); tmp.width = im.width; tmp.height = im.height;
        tmp.getContext('2d').putImageData(im, 0, 0);
        ctx.drawImage(tmp, M.l, M.t, pw, ph);
      } catch (e) { ctx.drawImage(canvas, M.l, M.t, pw, ph); }
      done();
    };

    bReset.addEventListener('click', function () { api.reset(); });
    bClear.addEventListener('click', function () { api.clearSelection(); });
    bPng.addEventListener('click', function () { api.png(); });
    bLasso.addEventListener('click', function () {
      lassoOn = !lassoOn;
      bLasso.classList.toggle('xq-on', lassoOn);
      if (sp) sp.set({ mouseMode: lassoOn ? 'lasso' : 'panZoom' });
    });
    root.addEventListener('mousemove', function (e) { mouse = e; if (tip.style.display === 'block') placeTip(); });
    root.addEventListener('mouseleave', function () { tip.style.display = 'none'; });

    let io = null;
    if ('IntersectionObserver' in g) {
      io = new IntersectionObserver(function (es) {
        visible = es.some(function (x) { return x.isIntersecting; });
        if (visible && !destroyed) { if (!sp) init(); else { const r = root.clientWidth; if (r && Math.abs(r - W) > 2) { api._sleep(); init(); } } }
      }, { rootMargin: '100px' });
      io.observe(root);
    } else { visible = true; setTimeout(init, 0); }
    root.api = api;
    root.destroy = api.destroy;
    root._xqDispose = api.destroy;
    return root;
  };

  function XQV() { return g.XQV; }
})(typeof window !== 'undefined' ? window : globalThis);
