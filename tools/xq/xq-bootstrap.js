// xq-bootstrap.js -- the only runtime asset embedded in the report as plain
// text. It inflates the gzipped, base64-encoded asset payloads that R writes
// into <script type="application/x-xq-asset" id="xqa-<name>"> blocks (see
// R/xq-assets.R and tools/embed_xq_assets.py) and evaluates them in order,
// so the ~800 KB of dashboard JavaScript/CSS travels at ~1/3 of its size.
//
// Chicken-and-egg: the vendor bundle (fflate) is itself one of the compressed
// assets, so the bootstrap carries a tiny self-contained inflate (the classic
// public-domain "tinf" algorithm, ~2 KB) to decompress it; everything after
// the vendor bundle is inflated with fflate's unzlibSync via the same code
// path. Synchronous on purpose: the head is injected before </head>, so all
// assets are live before Quarto evaluates its Observable cells.
(function (g) {
  'use strict';
  var XQ = g.XQ = g.XQ || {};
  var done = false;

  // ---- tinf: minimal DEFLATE (RFC 1951) / gzip (RFC 1952) decoder ---------
  function tinf(src) {
    var LEN_BASE = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
    var LEN_EXTRA = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
    var DIST_BASE = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
    var DIST_EXTRA = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
    var CLCIDX = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
    var srcIdx = 0, dstIdx = 0, bitBuf = 0, bitCnt = 0;
    var out = new Uint8Array(src.length * 5 + 65536);
    function grow() { var n = new Uint8Array(out.length * 2); n.set(out); out = n; }
    function getbit() {
      if (!bitCnt) { bitBuf = src[srcIdx++]; bitCnt = 8; }
      var b = bitBuf & 1; bitBuf >>>= 1; bitCnt--; return b;
    }
    function readbits(n) {
      var v = 0;
      for (var i = 0; i < n; i++) v |= getbit() << i;
      return v;
    }
    // returns a tree: offs[l] = running count of codes shorter than l; the
    // decoder walks canonical codes bit by bit and needs exactly that
    function buildTree(lengths, num) {
      var counts = new Array(16), offs = new Array(16), symbols = new Array(num), i;
      for (i = 0; i < 16; i++) { counts[i] = 0; offs[i] = 0; }
      for (i = 0; i < num; i++) counts[lengths[i]]++;
      counts[0] = 0;
      var sum = 0;
      for (i = 1; i < 16; i++) { offs[i] = sum; sum += counts[i]; }
      for (i = 0; i < num; i++) if (lengths[i]) symbols[offs[lengths[i]]++] = i;
      // offs was consumed as write cursors; rebuild the running sums
      var run = new Array(16); sum = 0;
      for (i = 1; i < 16; i++) { run[i] = sum; sum += counts[i]; }
      run[0] = 0;
      return { counts: run, symbols: symbols };
    }
    function decodeSymbol(tree) {
      var code = 0, first = 0, index = 0, len;
      for (len = 1; len < 16; len++) {
        code |= getbit();
        var count = tree.counts[len + 1] - tree.counts[len];
        if (code - first < count) return tree.symbols[index + (code - first)];
        index += count; first += count; first <<= 1; code <<= 1;
      }
      throw new Error('tinf: bad Huffman code');
    }
    function inflateBlock(lit, dist) {
      while (true) {
        var sym = decodeSymbol(lit);
        if (sym === 256) return;
        if (sym < 256) {
          if (dstIdx >= out.length) grow();
          out[dstIdx++] = sym;
        } else {
          sym -= 257;
          var length = readbits(LEN_EXTRA[sym]) + LEN_BASE[sym];
          var dsym = decodeSymbol(dist);
          var d = readbits(DIST_EXTRA[dsym]) + DIST_BASE[dsym];
          if (dstIdx + length > out.length) grow();
          while (length--) { out[dstIdx] = out[dstIdx - d]; dstIdx++; }
        }
      }
    }
    function inflateStored() {
      bitCnt = 0;
      var len = src[srcIdx] | (src[srcIdx + 1] << 8);
      srcIdx += 4;
      if (dstIdx + len > out.length) grow();
      while (len--) out[dstIdx++] = src[srcIdx++];
    }
    function fixedTrees() {
      var litlen = new Array(288), i;
      for (i = 0; i < 144; i++) litlen[i] = 8;
      for (; i < 256; i++) litlen[i] = 9;
      for (; i < 280; i++) litlen[i] = 7;
      for (; i < 288; i++) litlen[i] = 8;
      var distlen = new Array(30);
      for (i = 0; i < 30; i++) distlen[i] = 5;
      return { lit: buildTree(litlen, 288), dist: buildTree(distlen, 30) };
    }
    function dynamicTrees() {
      var nlen = readbits(5) + 257, ndist = readbits(5) + 1, ncode = readbits(4) + 4;
      var clen = new Array(19), i;
      for (i = 0; i < ncode; i++) clen[CLCIDX[i]] = readbits(3);
      for (; i < 19; i++) clen[CLCIDX[i]] = 0;
      var ctree = buildTree(clen, 19);
      var lengths = new Array(nlen + ndist), k = 0;
      while (k < nlen + ndist) {
        var sym = decodeSymbol(ctree);
        if (sym < 16) { lengths[k++] = sym; }
        else {
          var n = 0, prev = 0;
          if (sym === 16) { prev = lengths[k - 1]; n = readbits(2) + 3; }
          else if (sym === 17) { n = readbits(3) + 3; }
          else { n = readbits(7) + 11; }
          while (n--) lengths[k++] = prev;
        }
      }
      return {
        lit: buildTree(lengths.slice(0, nlen), nlen),
        dist: buildTree(lengths.slice(nlen), ndist)
      };
    }
    // skip the gzip header (RFC 1952) or the 2-byte zlib header (RFC 1950)
    var p = 0;
    if (src[0] === 0x1f && src[1] === 0x8b) {
      p = 10; var flg = src[3];
      if (flg & 4) p += 2 + (src[p] | (src[p + 1] << 8));
      if (flg & 8) { while (src[p]) p++; p++; }
      if (flg & 16) { while (src[p]) p++; p++; }
      if (flg & 2) p += 2;
    } else if (src[0] === 0x78) {
      p = 2;  // zlib: 2-byte header, deflate data, 4-byte adler32 (ignored)
    }
    srcIdx = p;
    var t = null;
    // hard bound: corrupt data must fail fast, never loop
    var guard = 100000000;
    while (true) {
      if (--guard < 0) throw new Error('tinf: ran away (corrupt stream?)');
      var last = getbit(), type = readbits(2);
      if (type === 0) inflateStored();
      else {
        if (type === 1) t = fixedTrees();
        else if (type === 2) t = dynamicTrees();
        else throw new Error('tinf: bad block type');
        inflateBlock(t.lit, t.dist);
      }
      if (last) break;
    }
    return out.subarray(0, dstIdx);
  }

  function b64ToBytes(b64) {
    var bin = atob(b64);
    var out = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }

  function inflate(bytes) {
    var unz = g.XQV && g.XQV.unzlibSync;
    return unz ? unz(bytes) : tinf(bytes);
  }

  XQ.inflateAssets = function (doc, onLog) {
    if (done) return;  // Quarto may evaluate embedded cells more than once
    done = true;
    doc = doc || document;
    var nodes = doc.querySelectorAll('script[type="application/x-xq-asset"]');
    for (var i = 0; i < nodes.length; i++) {
      var el = nodes[i];
      var id = el.id.replace(/^xqa-/, '');
      var kind = el.getAttribute('data-kind') || 'js';
      try {
        var bytes = b64ToBytes(el.textContent.replace(/\s+/g, ''));
        var txt = new TextDecoder('utf-8').decode(inflate(bytes));
        if (kind === 'css') {
          var style = doc.createElement('style');
          style.id = el.id.replace(/^xqa-/, 'xqs-');
          style.textContent = txt;
          doc.head.appendChild(style);
        } else {
          (0, eval)(txt);
        }
      } catch (e) {
        (typeof console !== 'undefined' && console.error
          ? console.error : (onLog || function () {}))('ExpressOM asset ' + id + ' failed: ' + e);
      }
    }
  };

  if (typeof document !== 'undefined') {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', function () { XQ.inflateAssets(document); });
    } else {
      XQ.inflateAssets(document);
    }
  }
})(typeof window !== 'undefined' ? window : globalThis);
