/* Light / dark theme switch of the ExpressOM report.
   The placeholder assigned to DEFAULT below is replaced by "auto", "light" or "dark"
   when the report is built.
   Priority: the viewer's own choice (localStorage) > the `theme` report setting >
   the operating system preference (prefers-color-scheme). */
(function () {
  'use strict';
  var KEY = 'expressom-report-theme';
  var DEFAULT = '__XP_DEFAULT__';
  var root = document.documentElement;
  var mq = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;
  var printing = false;

  function stored() {
    try {
      var v = window.localStorage.getItem(KEY);
      return v === 'dark' || v === 'light' ? v : null;
    } catch (e) { return null; }
  }
  function remember(v) {
    try { window.localStorage.setItem(KEY, v); } catch (e) { /* private window, file:// ... */ }
  }
  function resolve() {
    var s = stored();
    if (s) return s;
    if (DEFAULT === 'dark' || DEFAULT === 'light') return DEFAULT;
    return mq && mq.matches ? 'dark' : 'light';
  }
  function cssVar(name) {
    return window.getComputedStyle(root).getPropertyValue(name).trim();
  }
  function getPath(obj, path) {
    return path.split('.').reduce(function (o, k) { return o == null ? undefined : o[k]; }, obj);
  }

  // ---- plotly -----------------------------------------------------------------
  // plotly.js paints its own SVG, so CSS cannot recolour it: restyle every plot
  // through Plotly.relayout / restyle. The light settings of a plot are recorded
  // the first time it is darkened and put back when switching to light again.
  function themePlot(gd, theme) {
    if (!window.Plotly || !gd || !gd.layout || !gd._fullLayout || !gd.data) return;
    if (gd.__xpTheme === undefined && theme === 'light') { gd.__xpTheme = 'light'; return; }
    if (gd.__xpTheme === theme) return;
    var dark = theme === 'dark';
    var keys = ['paper_bgcolor', 'plot_bgcolor', 'font.color'];
    var axes = Object.keys(gd._fullLayout).filter(function (k) { return /^[xy]axis\d*$/.test(k); });
    axes.forEach(function (a) {
      ['gridcolor', 'zerolinecolor', 'linecolor', 'tickcolor'].forEach(function (p) { keys.push(a + '.' + p); });
    });
    if (!gd.__xpOrig) {
      gd.__xpOrig = { layout: {}, shapes: [], annotations: [], traces: [] };
      keys.forEach(function (k) { gd.__xpOrig.layout[k] = getPath(gd.layout, k); });
      (gd.layout.shapes || []).forEach(function (s, i) { gd.__xpOrig.shapes[i] = getPath(s, 'line.color'); });
      (gd.layout.annotations || []).forEach(function (a, i) { gd.__xpOrig.annotations[i] = getPath(a, 'font.color'); });
    }
    var upd = {};
    if (dark) {
      var bg = cssVar('--xp-plot-bg'), fg = cssVar('--xp-fg'), grid = cssVar('--xp-plot-grid');
      var zero = cssVar('--xp-plot-zero'), ref = cssVar('--xp-ref-line'), muted = cssVar('--xp-muted');
      upd.paper_bgcolor = bg; upd.plot_bgcolor = bg; upd['font.color'] = fg;
      axes.forEach(function (a) {
        upd[a + '.gridcolor'] = grid; upd[a + '.zerolinecolor'] = zero;
        upd[a + '.linecolor'] = zero; upd[a + '.tickcolor'] = zero;
      });
      (gd.layout.shapes || []).forEach(function (s, i) { upd['shapes[' + i + '].line.color'] = ref; });
      (gd.layout.annotations || []).forEach(function (a, i) { upd['annotations[' + i + '].font.color'] = muted; });
    } else {
      keys.forEach(function (k) {
        var v = gd.__xpOrig.layout[k];
        upd[k] = v === undefined ? null : v;
      });
      gd.__xpOrig.shapes.forEach(function (c, i) { upd['shapes[' + i + '].line.color'] = c === undefined ? null : c; });
      gd.__xpOrig.annotations.forEach(function (c, i) { upd['annotations[' + i + '].font.color'] = c === undefined ? null : c; });
    }
    gd.__xpTheme = theme;
    try { window.Plotly.relayout(gd, upd); } catch (e) { /* plot not ready: the observer retries */ gd.__xpTheme = undefined; }
    // traces tagged meta = "xp-contrast" are drawn black on light, near-white on dark
    var idx = [];
    gd.data.forEach(function (t, i) {
      // R plotly recycles a scalar `meta` over the rows of the trace's data
      var m = Array.isArray(t.meta) ? t.meta[0] : t.meta;
      if (m === 'xp-contrast') idx.push(i);
    });
    if (idx.length) {
      try { window.Plotly.restyle(gd, { 'marker.color': dark ? cssVar('--xp-contrast') : 'black' }, idx); } catch (e) { /* ignore */ }
    }
  }

  function themeAllPlots(theme) {
    var plots = document.querySelectorAll('.js-plotly-plot');
    for (var i = 0; i < plots.length; i++) themePlot(plots[i], theme);
  }

  // plots draw themselves after the page loads (and tab by tab): theme them as
  // they appear. The class is added before plotly.js has finished its setup,
  // hence the short retry.
  function watchPlots() {
    if (!window.MutationObserver) return;
    var pending = typeof WeakSet !== 'undefined' ? new WeakSet() : null;
    function whenReady(gd, tries) {
      if (gd._fullLayout && gd.data && gd.layout) { themePlot(gd, resolve()); if (pending) pending['delete'](gd); return; }
      if (tries > 60) { if (pending) pending['delete'](gd); return; }
      window.setTimeout(function () { whenReady(gd, tries + 1); }, 50);
    }
    var mo = new MutationObserver(function (muts) {
      muts.forEach(function (m) {
        var t = m.target;
        if (t.classList && t.classList.contains('js-plotly-plot') && t.__xpTheme === undefined && !(pending && pending.has(t))) {
          if (pending) pending.add(t);
          whenReady(t, 0);
        }
      });
    });
    mo.observe(document.body, { subtree: true, attributes: true, attributeFilter: ['class'] });
  }

  // ---- page --------------------------------------------------------------------
  function updateButton(theme) {
    var b = document.getElementById('xp-theme-toggle');
    if (!b) return;
    var dark = theme === 'dark';
    b.setAttribute('aria-pressed', dark ? 'true' : 'false');
    b.setAttribute('title', dark ? 'Switch to light mode' : 'Switch to dark mode');
    var l = b.querySelector('.xp-theme-label');
    if (l) l.textContent = dark ? 'Light' : 'Dark';
  }

  function apply(theme) {
    root.setAttribute('data-theme', theme);
    updateButton(theme);
    themeAllPlots(theme);
  }

  apply(resolve());   // runs while the header is parsed: no flash of the wrong theme

  function init() {
    updateButton(resolve());
    var b = document.getElementById('xp-theme-toggle');
    if (b) {
      b.addEventListener('click', function () {
        var next = root.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
        remember(next);
        apply(next);
      });
    }
    watchPlots();
    themeAllPlots(resolve());
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
  window.addEventListener('load', function () { themeAllPlots(root.getAttribute('data-theme')); });

  // follow the operating system until the viewer makes a choice
  if (mq && DEFAULT === 'auto') {
    var onChange = function () { if (!stored()) apply(resolve()); };
    if (mq.addEventListener) mq.addEventListener('change', onChange); else if (mq.addListener) mq.addListener(onChange);
  }

  // print on a white page, then come back
  window.addEventListener('beforeprint', function () {
    if (root.getAttribute('data-theme') === 'dark') { printing = true; apply('light'); }
  });
  window.addEventListener('afterprint', function () {
    if (printing) { printing = false; apply(resolve()); }
  });
})();
