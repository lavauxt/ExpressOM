/* xq-core.js -- ExpressOM Quarto dashboard runtime: data layer.
 *
 *  - XQ.table(key)   lazily decodes a zlib-compressed, columnar payload that R
 *                    wrote in a <script type="application/x-xq-table"> block
 *                    (float32 / int32 / -log10(p) / dictionary / utf8 columns).
 *                    Nothing is inflated until a column is first read, so a
 *                    table that is never opened costs nothing but its bytes.
 *  - XQ.where()      SQL-style WHERE-clause engine over typed arrays (JS).
 *  - XQ.duck         DuckDB-Wasm, loaded lazily (first SQL use) with the JS
 *                    engine as an automatic fallback.
 *
 * Depends on window.XQV (xq-vendor.min.js) for unzlibSync only.
 */
(function (g) {
  'use strict';
  const XQ = (g.XQ = g.XQ || {});
  const doc = g.document;

  // ------------------------------------------------------------------ helpers
  const textDecoder = new TextDecoder('utf-8');
  const NA_I32 = -2147483648;
  const NA_STR = '\u0001';

  function b64ToBytes(s) {
    const bin = atob(s);
    const n = bin.length;
    const out = new Uint8Array(n);
    for (let i = 0; i < n; i++) out[i] = bin.charCodeAt(i);
    return out;
  }
  XQ._b64ToBytes = b64ToBytes;

  const jsonCache = {};
  /** Small JSON blobs written by R in <script type="application/json" id="xqj-KEY">. */
  XQ.json = function (key, fallback) {
    if (key in jsonCache) return jsonCache[key];
    let v = fallback === undefined ? null : fallback;
    if (doc) {
      const el = doc.getElementById('xqj-' + key);
      if (el) {
        try { v = JSON.parse(el.textContent); } catch (e) { console.error('XQ.json(' + key + ')', e); }
      }
    }
    jsonCache[key] = v;
    return v;
  };
  Object.defineProperty(XQ, 'cfg', { get: function () { return XQ.json('cfg', {}) || {}; } });

  // ---------------------------------------------------------------- formatting
  const F = (XQ.fmt = {});
  F.int = function (v) {
    return v == null || v !== v ? '' : Math.round(v).toLocaleString('en-US');
  };
  F.num = function (v, d) {
    if (v == null || v !== v) return '';
    if (v === Infinity) return 'Inf';
    if (v === -Infinity) return '-Inf';
    d = d == null ? 3 : d;
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(Math.max(1, d - 1));
    return String(+v.toPrecision(d + 1));
  };
  F.p = function (v) {
    if (v == null || v !== v) return '';
    if (v === 0) return '0';
    return v < 1e-3 ? v.toExponential(1) : String(+v.toPrecision(3));
  };
  F.esc = function (s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  };
  F.nlp = function (nl) { // -log10(p) -> p, as text
    return nl == null || nl !== nl ? '' : F.p(Math.pow(10, -nl));
  };

  // ------------------------------------------------------------------- XTable
  /**
   * Columnar table. Column types ("t"):
   *   f32 / i32 / f64   numbers (NaN = NA; i32 uses -2^31 as NA when spec.na)
   *   nlp               p-value stored as float32 -log10(p); num() returns p
   *   d                 dictionary-coded strings (uint16/uint32 codes + levels)
   *   s                 utf8 strings joined by "\n" (NA marker \u0001)
   *   b                 logical (uint8: 0/1/255)
   */
  function XTable(key) {
    this.key = key;
    this._spec = null; // column name -> spec
    this._names = null;
    this._cache = {};
    this._lc = {};
    this._n = 0;
    this._mem = null;
    this.empty = false;
  }

  XTable.fromElement = function (key, el) {
    const t = new XTable(key);
    t._el = el;
    let meta = {};
    try { meta = JSON.parse(el.getAttribute('data-meta') || '{}'); } catch (e) { /* noop */ }
    t._n = meta.n || 0;
    t._names = (meta.cols || []).map(function (c) { return c[0]; });
    t._types = {};
    (meta.cols || []).forEach(function (c) { t._types[c[0]] = c[1]; });
    return t;
  };

  XTable.fromBytes = function (key, bytes) { // bytes: Uint8Array of the zlib payload
    const t = new XTable(key);
    t._zbytes = bytes;
    t._init();
    return t;
  };

  /** In-memory table from {name: Array|TypedArray}. */
  XTable.mem = function (cols, key) {
    const t = new XTable(key || 'mem');
    const names = Object.keys(cols);
    t._mem = {};
    t._spec = {};
    t._names = names;
    t._types = {};
    let n = 0;
    names.forEach(function (nm) {
      const a = cols[nm];
      n = Math.max(n, a.length);
      let ty = 's';
      if (a instanceof Float32Array) ty = 'f32';
      else if (a instanceof Float64Array) ty = 'f64';
      else if (a instanceof Int32Array || a instanceof Uint32Array || a instanceof Uint16Array || a instanceof Uint8Array) ty = 'i32';
      else if (Array.isArray(a) && a.length && a.every(function (x) { return x == null || typeof x === 'number'; })) {
        ty = 'f64';
      }
      t._mem[nm] = ty === 'f64' && Array.isArray(a) ? Float64Array.from(a, function (x) { return x == null ? NaN : x; }) : a;
      t._spec[nm] = { name: nm, t: ty };
      t._types[nm] = ty;
    });
    t._n = n;
    t._ready = true;
    return t;
  };

  const P = XTable.prototype;

  P._init = function () {
    if (this._ready) return;
    if (!this._zbytes) {
      const txt = this._el.textContent.trim();
      this._zbytes = b64ToBytes(txt);
    }
    const unz = (g.XQV && g.XQV.unzlibSync) || null;
    if (!unz) throw new Error('XQ: unzlibSync not available (xq-vendor.min.js missing)');
    const raw = unz(this._zbytes);
    this._zbytes = null;
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const hl = dv.getUint32(0, true);
    const h = JSON.parse(textDecoder.decode(raw.subarray(4, 4 + hl)));
    this._raw = raw;
    this._base = 8 * Math.ceil((4 + hl) / 8);
    this._n = h.n;
    this._spec = {};
    this._names = [];
    this._types = {};
    for (let i = 0; i < h.cols.length; i++) {
      const c = h.cols[i];
      this._spec[c.name] = c;
      this._names.push(c.name);
      this._types[c.name] = c.t;
    }
    this._ready = true;
  };

  Object.defineProperty(P, 'n', { get: function () { return this._n; } });
  Object.defineProperty(P, 'names', { get: function () { return this._names || []; } });

  P.has = function (name) { return !!(this._types && name in this._types); };
  P.type = function (name) { return this._types ? this._types[name] : undefined; };
  P.isNumeric = function (name) {
    const t = this.type(name);
    return t === 'f32' || t === 'i32' || t === 'f64' || t === 'nlp';
  };
  P.isP = function (name) { return this.type(name) === 'nlp'; };

  /** Case-insensitive column lookup. */
  P.resolve = function (name) {
    if (this.has(name)) return name;
    const lc = String(name).toLowerCase();
    const names = this.names;
    for (let i = 0; i < names.length; i++) if (names[i].toLowerCase() === lc) return names[i];
    return null;
  };
  /** First existing column among candidates (case-insensitive) or null. */
  P.first = function (cands) {
    for (let i = 0; i < cands.length; i++) { const r = this.resolve(cands[i]); if (r) return r; }
    return null;
  };

  P._view = function (spec, Ctor, bytesPer) {
    const off = this._base + spec.o;
    const abs = this._raw.byteOffset + off;
    if (abs % bytesPer === 0) return new Ctor(this._raw.buffer, abs, this._n);
    // misaligned (should not happen with R's 8-byte padding): copy
    const copy = this._raw.slice(off, off + this._n * bytesPer);
    return new Ctor(copy.buffer, 0, this._n);
  };

  /** Raw stored column (typed array / string array). */
  P._col = function (name) {
    if (this._mem) return this._mem[name];
    if (name in this._cache) return this._cache[name];
    this._init();
    const s = this._spec[name];
    if (!s) return null;
    let v;
    switch (s.t) {
      case 'f32': case 'nlp': v = this._view(s, Float32Array, 4); break;
      case 'i32': v = this._view(s, Int32Array, 4); break;
      case 'f64': v = this._view(s, Float64Array, 8); break;
      case 'd': v = s.w === 4 ? this._view(s, Uint32Array, 4) : this._view(s, Uint16Array, 2); break;
      case 'b': v = this._raw.subarray(this._base + s.o, this._base + s.o + this._n); break;
      case 's': {
        const txt = textDecoder.decode(this._raw.subarray(this._base + s.o, this._base + s.o + s.l));
        v = txt.split('\n');
        if (v.length > this._n) v.length = this._n;
        for (let i = 0; i < v.length; i++) if (v[i] === NA_STR) v[i] = null;
        break;
      }
      default: v = null;
    }
    this._cache[name] = v;
    return v;
  };

  /** Numeric column as a typed array (NaN = NA); p-value columns return p. */
  P.num = function (name) {
    const key = '#num:' + name;
    if (key in this._cache) return this._cache[key];
    const t = this.type(name);
    let v = null;
    if (t === 'f32' || t === 'f64') v = this._col(name);
    else if (t === 'i32') {
      const a = this._col(name);
      const s = this._mem ? null : this._spec[name];
      if (s && s.na) {
        v = new Float64Array(a.length);
        for (let i = 0; i < a.length; i++) v[i] = a[i] === NA_I32 ? NaN : a[i];
      } else v = a;
    } else if (t === 'nlp') {
      const a = this._col(name);
      v = new Float64Array(a.length);
      for (let i = 0; i < a.length; i++) v[i] = Math.pow(10, -a[i]);
    } else if (t === 'b') {
      const a = this._col(name);
      v = new Float64Array(a.length);
      for (let i = 0; i < a.length; i++) v[i] = a[i] === 255 ? NaN : a[i];
    }
    this._cache[key] = v;
    return v;
  };

  /** -log10(p) as Float32Array. Exact for p-columns (no underflow). */
  P.nlp = function (name) {
    const key = '#nlp:' + name;
    if (key in this._cache) return this._cache[key];
    let v;
    if (this.type(name) === 'nlp') v = this._col(name);
    else {
      const a = this.num(name);
      v = new Float32Array(a.length);
      for (let i = 0; i < a.length; i++) {
        const x = a[i];
        v[i] = x !== x ? NaN : x <= 0 ? 320 : -Math.log10(x);
      }
    }
    this._cache[key] = v;
    return v;
  };

  /** String view of any column (array of string|null). */
  P.str = function (name) {
    const key = '#str:' + name;
    if (key in this._cache) return this._cache[key];
    const t = this.type(name);
    let v;
    if (t === 's') v = this._col(name);
    else if (t === 'd') {
      const codes = this._col(name);
      const lev = this._spec[name].k;
      const na = this._spec[name].w === 4 ? 4294967295 : 65535;
      v = new Array(codes.length);
      for (let i = 0; i < codes.length; i++) v[i] = codes[i] === na ? null : lev[codes[i]];
    } else if (t === 'b') {
      const a = this._col(name);
      v = Array.prototype.map.call(a, function (x) { return x === 255 ? null : x ? 'TRUE' : 'FALSE'; });
    } else if (t) {
      const a = this.num(name);
      v = new Array(a.length);
      for (let i = 0; i < a.length; i++) v[i] = a[i] !== a[i] ? null : String(a[i]);
    } else v = null;
    this._cache[key] = v;
    return v;
  };

  /** Dictionary codes + levels (for categorical columns). */
  P.categories = function (name) {
    const t = this.type(name);
    if (t === 'd') {
      const s = this._spec[name];
      return { codes: this._col(name), levels: s.k, na: s.w === 4 ? 4294967295 : 65535 };
    }
    const strs = this.str(name);
    if (!strs) return null;
    const map = new Map();
    const levels = [];
    const codes = new Uint32Array(strs.length);
    for (let i = 0; i < strs.length; i++) {
      const s = strs[i];
      if (s == null) { codes[i] = 4294967295; continue; }
      let c = map.get(s);
      if (c === undefined) { c = levels.length; map.set(s, c); levels.push(s); }
      codes[i] = c;
    }
    return { codes: codes, levels: levels, na: 4294967295 };
  };

  /** Distinct non-NA values of a categorical column (sorted). */
  P.unique = function (name) {
    const c = this.categories(name);
    if (!c) return [];
    const seen = new Uint8Array(c.levels.length);
    for (let i = 0; i < c.codes.length; i++) if (c.codes[i] !== c.na) seen[c.codes[i]] = 1;
    const out = [];
    for (let i = 0; i < seen.length; i++) if (seen[i]) out.push(c.levels[i]);
    return out.sort();
  };

  /** Plain value of cell (number | string | null). */
  P.get = function (name, i) {
    const t = this.type(name);
    if (t === 's' || t === 'd' || t === 'b') { const s = this.str(name); return s[i]; }
    const v = this.num(name)[i];
    return v !== v ? null : v;
  };

  /** Row object (display / tooltips only). */
  P.row = function (i, cols) {
    const o = {};
    const names = cols || this.names;
    for (let k = 0; k < names.length; k++) o[names[k]] = this.get(names[k], i);
    return o;
  };
  P.rows = function (idx, cols) {
    const names = cols || this.names;
    const out = [];
    const m = idx ? idx.length : this.n;
    for (let k = 0; k < m; k++) out.push(this.row(idx ? idx[k] : k, names));
    return out;
  };

  P.all = function () {
    const a = new Uint32Array(this.n);
    for (let i = 0; i < a.length; i++) a[i] = i;
    return a;
  };

  /** Row index by (string) key column, lazily built. */
  P.indexBy = function (name) {
    const key = '#idx:' + name;
    if (!(key in this._cache)) {
      const m = new Map();
      const s = this.str(name);
      for (let i = 0; i < s.length; i++) if (s[i] != null && !m.has(s[i])) m.set(s[i], i);
      this._cache[key] = m;
    }
    return this._cache[key];
  };

  /** Case-insensitive substring search across string columns -> Uint32Array. */
  P.search = function (query, cols) {
    const q = String(query || '').trim().toLowerCase();
    if (!q) return null;
    const self = this;
    const use = (cols || this.names.filter(function (c) { return !self.isNumeric(c); }));
    const arrays = use.map(function (c) {
      if (!(c in self._lc)) {
        const s = self.str(c);
        self._lc[c] = s ? s.map(function (x) { return x == null ? '' : x.toLowerCase(); }) : null;
      }
      return self._lc[c];
    }).filter(Boolean);
    const tokens = q.split(/\s+/).filter(Boolean);
    const out = [];
    for (let i = 0; i < this.n; i++) {
      let ok = true;
      for (let t = 0; t < tokens.length && ok; t++) {
        let hit = false;
        for (let a = 0; a < arrays.length; a++) { if (arrays[a][i].indexOf(tokens[t]) !== -1) { hit = true; break; } }
        ok = hit;
      }
      if (ok) out.push(i);
    }
    return Uint32Array.from(out);
  };

  /** Sort a set of row indices by a column (NA last). */
  P.order = function (name, desc, idx) {
    const base = idx ? Uint32Array.from(idx) : this.all();
    const t = this.type(name);
    if (!t) return base;
    const sgn = desc ? -1 : 1;
    if (t === 's' || t === 'd' || t === 'b') {
      const s = this.str(name);
      base.sort(function (a, b) {
        const x = s[a], y = s[b];
        if (x == null) return y == null ? a - b : 1;
        if (y == null) return -1;
        return x < y ? -sgn : x > y ? sgn : a - b;
      });
    } else {
      const v = t === 'nlp' ? this.nlp(name) : this.num(name);
      const nlp = t === 'nlp';
      base.sort(function (a, b) {
        const x = v[a], y = v[b];
        const xn = x !== x, yn = y !== y;
        if (xn) return yn ? a - b : 1;
        if (yn) return -1;
        // for p-columns the stored value is -log10(p): ascending p = descending stored
        return x === y ? a - b : (nlp ? (x < y ? sgn : -sgn) : (x < y ? -sgn : sgn));
      });
    }
    return base;
  };

  P.csv = function (idx, cols) {
    const names = (cols || this.names).filter(this.has, this);
    const esc = function (v) {
      if (v == null) return '';
      const s = typeof v === 'number' ? String(v) : String(v);
      return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    const lines = [names.map(esc).join(',')];
    const m = idx ? idx.length : this.n;
    for (let k = 0; k < m; k++) {
      const i = idx ? idx[k] : k;
      lines.push(names.map(function (c) { return esc(this.get(c, i)); }, this).join(','));
    }
    return lines.join('\n');
  };

  const tableCache = {};
  const emptyTable = function (key) { const t = XTable.mem({}, key); t.empty = true; return t; };
  /** Table embedded by R (lazy). Missing key -> empty table. */
  XQ.table = function (key) {
    if (tableCache[key]) return tableCache[key];
    const el = doc && doc.getElementById('xqt-' + key);
    const t = el ? XTable.fromElement(key, el) : emptyTable(key);
    tableCache[key] = t;
    return t;
  };
  XQ.XTable = XTable;
  XQ.mem = XTable.mem;
  XQ.hasTable = function (key) { return !!(doc && doc.getElementById('xqt-' + key)); };

  // ------------------------------------------------------------ WHERE engine
  // Grammar (a practical SQL subset; identical meaning in DuckDB):
  //   expr   := or ;  or := and (OR and)* ;  and := not (AND not)* ; not := NOT not | cmp
  //   cmp    := add [ (=|==|!=|<>|<|<=|>|>=) add | IS [NOT] NULL
  //                 | [NOT] (LIKE|ILIKE) add | [NOT] IN ( list ) | [NOT] BETWEEN add AND add ]
  //   add    := mul ((+|-) mul)* ; mul := unary ((*|/|%) unary)*
  //   unary  := - unary | primary
  //   primary:= number | 'string' | ident | "quoted ident" | fn(args) | ( expr )
  // Three-valued logic: comparisons with NULL/NaN are NULL; rows are kept
  // only when the whole predicate is TRUE (same as SQL).
  function tokenize(src) {
    const toks = [];
    let i = 0;
    const n = src.length;
    while (i < n) {
      const c = src[i];
      if (/\s/.test(c)) { i++; continue; }
      if (c === "'") {
        let j = i + 1, s = '';
        while (j < n) {
          if (src[j] === "'") { if (src[j + 1] === "'") { s += "'"; j += 2; continue; } break; }
          s += src[j++];
        }
        if (j >= n) throw new Error('unterminated string literal');
        toks.push({ t: 'str', v: s }); i = j + 1; continue;
      }
      if (c === '"') {
        const j = src.indexOf('"', i + 1);
        if (j < 0) throw new Error('unterminated quoted identifier');
        toks.push({ t: 'id', v: src.slice(i + 1, j), q: true }); i = j + 1; continue;
      }
      const m = /^(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)/.exec(src.slice(i));
      if (m) { toks.push({ t: 'num', v: parseFloat(m[0]) }); i += m[0].length; continue; }
      const id = /^[A-Za-z_][A-Za-z0-9_.]*/.exec(src.slice(i));
      if (id) { toks.push({ t: 'id', v: id[0] }); i += id[0].length; continue; }
      const op = /^(<=|>=|<>|!=|==|\|\||[=<>+\-*\/%(),])/.exec(src.slice(i));
      if (op) { toks.push({ t: 'op', v: op[0] }); i += op[0].length; continue; }
      throw new Error("unexpected character '" + c + "'");
    }
    return toks;
  }

  const KW = { and: 1, or: 1, not: 1, like: 1, ilike: 1, in: 1, is: 1, null: 1, between: 1, true: 1, false: 1 };

  function likeToRegex(pat, ci) {
    let re = '^';
    for (let i = 0; i < pat.length; i++) {
      const c = pat[i];
      if (c === '%') re += '[\\s\\S]*';
      else if (c === '_') re += '[\\s\\S]';
      else re += c.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    }
    return new RegExp(re + '$', ci ? 'i' : '');
  }

  const FUNCS = {
    abs: function (x) { return x == null ? null : Math.abs(x); },
    log10: function (x) { return x == null || x <= 0 ? null : Math.log10(x); },
    log2: function (x) { return x == null || x <= 0 ? null : Math.log2(x); },
    ln: function (x) { return x == null || x <= 0 ? null : Math.log(x); },
    log: function (x) { return x == null || x <= 0 ? null : Math.log10(x); },
    exp: function (x) { return x == null ? null : Math.exp(x); },
    sqrt: function (x) { return x == null || x < 0 ? null : Math.sqrt(x); },
    floor: function (x) { return x == null ? null : Math.floor(x); },
    ceil: function (x) { return x == null ? null : Math.ceil(x); },
    round: function (x, d) { if (x == null) return null; const f = Math.pow(10, d || 0); return Math.round(x * f) / f; },
    lower: function (s) { return s == null ? null : String(s).toLowerCase(); },
    upper: function (s) { return s == null ? null : String(s).toUpperCase(); },
    length: function (s) { return s == null ? null : String(s).length; },
    coalesce: function (a, b) { return a == null ? b : a; },
    contains: function (s, t) { return s == null || t == null ? null : String(s).indexOf(String(t)) !== -1; },
    starts_with: function (s, t) { return s == null || t == null ? null : String(s).indexOf(String(t)) === 0; },
    ends_with: function (s, t) { return s == null || t == null ? null : String(s).endsWith(String(t)); },
    regexp_matches: function (s, p) { if (s == null || p == null) return null; try { return new RegExp(p).test(String(s)); } catch (e) { return null; } }
  };

  function compileWhere(tbl, src) {
    const toks = tokenize(src);
    let pos = 0;
    const peek = function () { return toks[pos]; };
    const kw = function (w) { const t = toks[pos]; return t && t.t === 'id' && !t.q && t.v.toLowerCase() === w; };
    const eatKw = function (w) { if (kw(w)) { pos++; return true; } return false; };
    const eatOp = function (o) { const t = toks[pos]; if (t && t.t === 'op' && t.v === o) { pos++; return true; } return false; };
    const need = function (o) { if (!eatOp(o)) throw new Error("expected '" + o + "'"); };
    const colRef = function (name) {
      const c = tbl.resolve(name);
      if (!c) throw new Error('unknown column "' + name + '"');
      if (tbl.isNumeric(c)) { const a = tbl.num(c); return function (i) { const v = a[i]; return v !== v ? null : v; }; }
      if (tbl.type(c) === 'b') { const b = tbl._col(c); return function (i) { return b[i] === 255 ? null : b[i] === 1; }; }
      const s = tbl.str(c);
      return function (i) { return s[i]; };
    };

    function parseOr() {
      let l = parseAnd();
      while (eatKw('or')) { const a = l, b = parseAnd(); l = function (i) { const x = a(i); if (x === true) return true; const y = b(i); if (y === true) return true; return x == null || y == null ? null : false; }; }
      return l;
    }
    function parseAnd() {
      let l = parseNot();
      while (kw('and')) {
        pos++;
        const a = l, b = parseNot();
        l = function (i) { const x = a(i); if (x === false) return false; const y = b(i); if (y === false) return false; return x == null || y == null ? null : true; };
      }
      return l;
    }
    function parseNot() {
      if (eatKw('not')) { const e = parseNot(); return function (i) { const x = e(i); return x == null ? null : !x; }; }
      return parseCmp();
    }
    function cmpOp(op, a, b) {
      return function (i) {
        const x = a(i), y = b(i);
        if (x == null || y == null) return null;
        switch (op) {
          case '=': case '==': return x === y || (typeof x !== typeof y && String(x) === String(y));
          case '!=': case '<>': return !(x === y || (typeof x !== typeof y && String(x) === String(y)));
          case '<': return x < y; case '<=': return x <= y; case '>': return x > y; default: return x >= y;
        }
      };
    }
    function parseCmp() {
      const a = parseAdd();
      const t = peek();
      if (t && t.t === 'op' && /^(=|==|!=|<>|<|<=|>|>=)$/.test(t.v)) { pos++; return cmpOp(t.v, a, parseAdd()); }
      if (eatKw('is')) {
        const neg = eatKw('not');
        if (!eatKw('null')) throw new Error('expected NULL after IS');
        return function (i) { const x = a(i); return neg ? x != null : x == null; };
      }
      let neg = false;
      if (kw('not')) {
        const nx = toks[pos + 1];
        if (nx && nx.t === 'id' && /^(like|ilike|in|between)$/i.test(nx.v)) { neg = true; pos++; }
      }
      if (kw('like') || kw('ilike')) {
        const ci = kw('ilike'); pos++;
        const p = parseAdd();
        const cache = {};
        return function (i) {
          const x = a(i), pt = p(i);
          if (x == null || pt == null) return null;
          const k = String(pt);
          const re = cache[k] || (cache[k] = likeToRegex(k, ci));
          const r = re.test(String(x));
          return neg ? !r : r;
        };
      }
      if (eatKw('in')) {
        need('(');
        const items = [parseAdd()];
        while (eatOp(',')) items.push(parseAdd());
        need(')');
        return function (i) {
          const x = a(i);
          if (x == null) return null;
          let sawNull = false;
          for (let k = 0; k < items.length; k++) {
            const y = items[k](i);
            if (y == null) { sawNull = true; continue; }
            if (x === y || String(x) === String(y)) return !neg;
          }
          return sawNull ? null : neg;
        };
      }
      if (eatKw('between')) {
        const lo = parseAdd();
        if (!eatKw('and')) throw new Error('expected AND in BETWEEN');
        const hi = parseAdd();
        return function (i) {
          const x = a(i), l = lo(i), h = hi(i);
          if (x == null || l == null || h == null) return null;
          const r = x >= l && x <= h;
          return neg ? !r : r;
        };
      }
      return a;
    }
    function arith(op, a, b) {
      return function (i) {
        const x = a(i), y = b(i);
        if (x == null || y == null) return null;
        switch (op) {
          case '+': return x + y; case '-': return x - y; case '*': return x * y;
          case '/': return y === 0 ? null : x / y; default: return y === 0 ? null : x % y;
        }
      };
    }
    function parseAdd() {
      let l = parseMul();
      for (;;) {
        const t = peek();
        if (t && t.t === 'op' && (t.v === '+' || t.v === '-')) { pos++; l = arith(t.v, l, parseMul()); }
        else if (t && t.t === 'op' && t.v === '||') {
          pos++; const a = l, b = parseMul();
          l = function (i) { const x = a(i), y = b(i); return x == null || y == null ? null : String(x) + String(y); };
        } else return l;
      }
    }
    function parseMul() {
      let l = parseUnary();
      for (;;) {
        const t = peek();
        if (t && t.t === 'op' && (t.v === '*' || t.v === '/' || t.v === '%')) { pos++; l = arith(t.v, l, parseUnary()); }
        else return l;
      }
    }
    function parseUnary() {
      if (eatOp('-')) { const e = parseUnary(); return function (i) { const x = e(i); return x == null ? null : -x; }; }
      if (eatOp('+')) return parseUnary();
      return parsePrimary();
    }
    function parsePrimary() {
      const t = toks[pos++];
      if (!t) throw new Error('unexpected end of expression');
      if (t.t === 'num') { const v = t.v; return function () { return v; }; }
      if (t.t === 'str') { const v = t.v; return function () { return v; }; }
      if (t.t === 'op' && t.v === '(') { const e = parseOr(); need(')'); return e; }
      if (t.t === 'id') {
        const lc = t.v.toLowerCase();
        if (!t.q) {
          if (lc === 'null') return function () { return null; };
          if (lc === 'true') return function () { return true; };
          if (lc === 'false') return function () { return false; };
          if (toks[pos] && toks[pos].t === 'op' && toks[pos].v === '(') {
            const fn = FUNCS[lc];
            if (!fn) throw new Error('unsupported function ' + t.v + '() in the JS engine');
            pos++;
            const args = [];
            if (!eatOp(')')) { args.push(parseOr()); while (eatOp(',')) args.push(parseOr()); need(')'); }
            return function (i) { return fn.apply(null, args.map(function (a) { return a(i); })); };
          }
          if (KW[lc]) throw new Error("unexpected keyword '" + t.v + "'");
        }
        return colRef(t.v);
      }
      throw new Error("unexpected '" + (t.v) + "'");
    }

    const root = parseOr();
    if (pos < toks.length) throw new Error("unexpected '" + toks[pos].v + "'");
    return root;
  }

  /** Evaluate a WHERE clause with the JS engine -> {idx, ms, engine}. */
  XQ.where = function (tbl, clause, base) {
    const t0 = (g.performance || Date).now();
    const c = String(clause || '').trim();
    if (!c) return { idx: base || tbl.all(), ms: 0, engine: 'js' };
    const pred = compileWhere(tbl, c);
    const out = [];
    if (base) { for (let k = 0; k < base.length; k++) if (pred(base[k]) === true) out.push(base[k]); }
    else { const n = tbl.n; for (let i = 0; i < n; i++) if (pred(i) === true) out.push(i); }
    return { idx: Uint32Array.from(out), ms: (g.performance || Date).now() - t0, engine: 'js' };
  };

  // ------------------------------------------------------------- DuckDB-Wasm
  // Loaded lazily the first time a SQL box is used; never at page load. The
  // ESM + wasm come from the CDN (default jsDelivr, pinned) or from a URL the
  // user sets with expressom_report_settings(duckdb_url = ...). When the load
  // fails (offline) every box keeps working through the JS engine above.
  const duck = (XQ.duck = {
    state: 'idle', // idle | loading | ready | failed
    error: null,
    _listeners: [],
    _registered: {},
    on: function (fn) { this._listeners.push(fn); },
    _set: function (s, err) {
      this.state = s; this.error = err || null;
      this._listeners.forEach(function (f) { try { f(s); } catch (e) { /* noop */ } });
    },
    ensure: function () {
      if (this._p) return this._p;
      const cfg = XQ.cfg;
      if (cfg.duckdb === false) { this._set('failed', 'disabled'); return (this._p = Promise.resolve(false)); }
      this._set('loading');
      const url = cfg.duckdb_url || 'https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm';
      const self = this;
      this._p = (async function () {
        const mod = await import(url);
        const bundles = mod.getJsDelivrBundles();
        const bundle = await mod.selectBundle(bundles);
        const workerUrl = URL.createObjectURL(new Blob(['importScripts("' + bundle.mainWorker + '");'], { type: 'text/javascript' }));
        const worker = new Worker(workerUrl);
        const db = new mod.AsyncDuckDB(new mod.ConsoleLogger(mod.LogLevel.WARNING), worker);
        await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
        URL.revokeObjectURL(workerUrl);
        self._db = db;
        self._conn = await db.connect();
        self._set('ready');
        return true;
      })().catch(function (e) {
        console.warn('XQ: DuckDB-Wasm unavailable, using the JS filter engine.', e && e.message);
        self._set('failed', e && e.message);
        return false;
      });
      return this._p;
    },
    /** Register an XTable as a DuckDB table (with a __i row-index column). */
    register: async function (name, tbl) {
      if (this._registered[name]) return;
      const names = tbl.names;
      const esc = function (v) {
        if (v == null) return '';
        const s = String(v);
        return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
      };
      const cols = names.map(function (c) {
        const t = tbl.type(c);
        return t === 's' || t === 'd' ? tbl.str(c) : t === 'b' ? tbl.str(c) : tbl.num(c);
      });
      const lines = ['__i,' + names.map(esc).join(',')];
      for (let i = 0; i < tbl.n; i++) {
        const parts = new Array(cols.length + 1);
        parts[0] = i;
        for (let k = 0; k < cols.length; k++) {
          const v = cols[k][i];
          parts[k + 1] = v == null || v !== v ? '' : esc(typeof v === 'number' ? (Number.isInteger(v) ? v : +v.toPrecision(9)) : v);
        }
        lines.push(parts.join(','));
      }
      await this._db.registerFileText(name + '.csv', lines.join('\n'));
      const q = function (c) { return '"' + c.replace(/"/g, '""') + '"'; };
      const types = ['"__i": \'INTEGER\''].concat(names.map(function (c) {
        const t = tbl.type(c);
        const sql = t === 'i32' ? 'INTEGER' : t === 'f32' || t === 'f64' || t === 'nlp' ? 'DOUBLE' : t === 'b' ? 'VARCHAR' : 'VARCHAR';
        return q(c) + ": '" + sql + "'";
      }));
      await this._conn.query('DROP TABLE IF EXISTS ' + q(name));
      await this._conn.query('CREATE TABLE ' + q(name) + " AS SELECT * FROM read_csv('" + name + ".csv', header=true, columns={" + types.join(', ') + '})');
      this._registered[name] = true;
    },
    where: async function (name, tbl, clause) {
      const t0 = performance.now();
      await this.register(name, tbl);
      const q = '"' + name.replace(/"/g, '""') + '"';
      const sql = String(clause || '').trim() ? 'SELECT __i FROM ' + q + ' WHERE ' + clause + ' ORDER BY __i' : 'SELECT __i FROM ' + q + ' ORDER BY __i';
      const r = await this._conn.query(sql);
      const arr = r.getChildAt(0).toArray();
      return { idx: Uint32Array.from(arr), ms: performance.now() - t0, engine: 'duckdb' };
    },
    /** Free SQL -> {table, ms}. */
    query: async function (sql) {
      const t0 = performance.now();
      const r = await this._conn.query(sql);
      const names = r.schema.fields.map(function (f) { return f.name; });
      const cols = {};
      names.forEach(function (nm, k) {
        const arr = r.getChildAt(k).toArray();
        const first = arr.length ? arr[0] : null;
        cols[nm] = typeof first === 'bigint' ? Array.from(arr, Number)
          : (typeof first === 'number' || first == null) && (arr instanceof Float64Array || arr instanceof Float32Array || arr instanceof Int32Array)
            ? Float64Array.from(arr) : Array.from(arr, function (x) { return x == null ? null : typeof x === 'bigint' ? Number(x) : x; });
      });
      return { table: XTable.mem(cols, 'sql'), ms: performance.now() - t0, engine: 'duckdb' };
    }
  });

  /** Filter with the best available engine (DuckDB when ready, else JS). */
  XQ.filter = async function (name, tbl, clause, base) {
    const c = String(clause || '').trim();
    if (duck.state === 'ready' && !base) {
      try { return await duck.where(name, tbl, c); } catch (e) { /* fall through to JS to report a friendlier error */ }
    }
    return XQ.where(tbl, c, base);
  };

  // generic utilities used by widgets/charts
  XQ.download = function (filename, text, mime) {
    const blob = new Blob([text], { type: mime || 'text/csv;charset=utf-8' });
    const a = doc.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    doc.body.appendChild(a);
    a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  };
  XQ.debounce = function (fn, ms) {
    let h = null;
    return function () { const a = arguments, self = this; clearTimeout(h); h = setTimeout(function () { fn.apply(self, a); }, ms); };
  };
  XQ._compileWhere = compileWhere; // exposed for tests
})(typeof window !== 'undefined' ? window : globalThis);
