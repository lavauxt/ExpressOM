#!/usr/bin/env python3
"""Verify R/xq-assets.R without R installed.

Tokenizes and parses the restricted R subset that tools/embed_xq_assets.py
emits (string literals, numbers, identifiers, list()/paste0()/paste(),
named arguments, commas, parens), rejecting exactly what made load_all()
fail: an empty argument -- ", )" -- which R's parser accepts but evaluation
rejects with "argument N is empty".

Then re-derives both constants semantically: inflates every b64 payload and
compares with the tools/xq source, and reconstructs .xq_bootstrap_js to the
bootstrap file byte-for-byte. Exit 0 only if everything holds.
"""
import base64
import pathlib
import re
import sys
import zlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "R" / "xq-assets.R"
XQ = ROOT / "tools" / "xq"

fails = []
def ok(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)

# ---- tokenizer -----------------------------------------------------------------
TOKEN_RE = re.compile(r"""
    (?P<comment>\#[^\n]*)
  | (?P<str>"(?:[^"\\]|\\.)*")
  | (?P<id>[A-Za-z_.][A-Za-z0-9_.]*)
  | (?P<num>\d+)
  | (?P<punct>\(|\)|,|=|<-)
  | (?P<nl>\n)
  | (?P<space>[ \t]+)
""", re.VERBOSE)

def tokenize(text):
    toks, i = [], 0
    while i < len(text):
        m = TOKEN_RE.match(text, i)
        if not m:
            fails.append(f"tokenizer stuck at offset {i}: {text[i:i+40]!r}")
            return None
        i = m.end()
        kind = m.lastgroup
        if kind in ("comment", "space", "nl"):
            continue
        toks.append((kind, m.group()))
    return toks

# ---- parser: expr := call | STR ; no empty args allowed ------------------------
class ParseError(Exception):
    pass

ESCAPES = {"n": "\n", "t": "\t", "r": "\r", '"': '"', "\\": "\\", "'": "'"}

def unescape_r(s):
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == "\\":
            if i + 1 >= len(s) or s[i + 1] not in ESCAPES:
                raise ParseError(f"unknown R escape {s[i:i+2]!r}")
            out.append(ESCAPES[s[i + 1]])
            i += 2
        else:
            out.append(c)
            i += 1
    return "".join(out)

def parse_call(toks, pos):
    """parse one expression starting at toks[pos]; return (value, next_pos).
    value: str for string literals; dict for calls {fn, args:[(name|None, val)]}"""
    kind, val = toks[pos]
    if kind == "str":
        return unescape_r(val[1:-1]), pos + 1
    if kind != "id":
        raise ParseError(f"expected identifier or string, got {kind} {val!r} at {pos}")
    fn, pos = val, pos + 1
    if toks[pos] != ("punct", "("):
        raise ParseError(f"expected ( after {fn}")
    pos += 1
    args = []
    while toks[pos] != ("punct", ")"):
        # optional named argument: id = expr
        name = None
        if toks[pos][0] == "id" and toks[pos + 1] == ("punct", "="):
            name = toks[pos][1]
            pos += 2
        elif toks[pos][0] == "str" and toks[pos + 1] == ("punct", "="):
            raise ParseError("string used as argument name")
        val2, pos = parse_call(toks, pos)
        args.append((name, val2))
        if toks[pos] == ("punct", ","):
            pos += 1
            if toks[pos] == ("punct", ")"):
                raise ParseError("empty trailing argument (the load_all() bug)")
    return {"fn": fn, "args": args}, pos + 1

def eval_r(ast):
    if isinstance(ast, str):
        return ast
    fn, args = ast["fn"], ast["args"]
    pos = [v for n, v in args if n is None]
    named = {n: v for n, v in args if n is not None}
    if fn == "paste0":
        if not pos:
            raise ParseError("paste0 with no positional argument")
        return "".join(eval_r(a) for a in pos)
    if fn == "paste":
        sep = eval_r(named["sep"]) if "sep" in named else " "
        if not pos:
            raise ParseError("paste with no positional argument")
        return sep.join(eval_r(a) for a in pos)
    if fn == "list":
        return {n: eval_r(v) for n, v in args}
    raise ParseError(f"unknown call {fn}() in generated file")

# ---- parse the whole file: comment header, then two assignments ---------------
text = SRC.read_text(encoding="utf-8")
toks = tokenize(text)
ok(toks is not None and not fails, "tokenizes cleanly")
if toks is None:
    sys.exit(1)

constants = {}
pos = 0
try:
    while pos < len(toks):
        if toks[pos][0] != "id":
            fails.append(f"unexpected top-level token {toks[pos]}")
            break
        if toks[pos + 1] != ("punct", "<-"):
            fails.append(f"expected <- after {toks[pos][1]}")
            break
        name = toks[pos][1]
        pos += 2
        ast, pos = parse_call(toks, pos)
        constants[name] = eval_r(ast)
except ParseError as e:
    fails.append(f"parse/eval error: {e}")
    print(f"FAIL parse/eval error: {e}")
ok(set(constants) == {".xq_asset_data", ".xq_bootstrap_js"},
   f"defines exactly the two constants (got {sorted(constants)})")

data = constants[".xq_asset_data"]
ok(isinstance(data, dict) and len(data) == 7, f".xq_asset_data has 7 records (got {len(data)})")

total = 0
for key, rec in data.items():
    ok(isinstance(rec, dict) and set(rec) == {"b64", "sha", "kind", "file"},
       f"{key}: record fields b64/sha/kind/file (got {sorted(rec)})")
    raw = zlib.decompress(base64.b64decode(rec["b64"]))
    src = (XQ / rec["file"]).read_bytes()
    ok(raw == src, f"{key}: payload inflates to tools/xq/{rec['file']} ({len(src)} B)")
    ok(rec["kind"] == ("css" if rec["file"].endswith(".css") else "js"), f"{key}: kind")
    total += len(raw)
ok(total == 795509, f"total inflated size matches the 7 assets ({total} B)")

boot = constants[".xq_bootstrap_js"]
boot_src = (XQ / "xq-bootstrap.js").read_text(encoding="utf-8")
# the generator drops the file's final newline when splitting lines; that is
# semantically irrelevant to JS, so compare with it normalised
ok(boot == boot_src or boot + "\n" == boot_src,
   f".xq_bootstrap_js == tools/xq/xq-bootstrap.js ({len(boot_src)} chars, one string)")
ok("\n" in boot, "bootstrap is a single string with real newlines (not a vector)")

print(("FAILED %d" % len(fails)) if fails else "ALL OK")
sys.exit(1 if fails else 0)
