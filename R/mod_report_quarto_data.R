# mod_report_quarto_data.R -- data layer of the Quarto dashboard engine.
#
# Every table the dashboard needs is written once, in R, as a zlib-compressed
# *columnar* binary container and embedded in a
# `<script type="application/x-xq-table" id="xqt-<key>">` block. The browser
# side (tools/xq/xq-core.js, embedded compressed in R/xq-assets.R) inflates a table
# first read, so tables that are never opened cost nothing but their bytes.
#
# Container layout (all little-endian), before zlib:
#
#   uint32 header_len | header JSON (utf8) | pad to 8 | column blocks (each 8-aligned)
#
# Header: {"n": rows, "cols": [{"name", "t", "o", "l", ...}]} where `o` is the
# block offset (relative to the end of the padded header) and `l` its length.
# Column types `t`:
#   f32  float32 (NaN = NA)                      i32  int32 (NA = -2^31, "na": true)
#   f64  float64 (huge integers)                 nlp  float32 of -log10(p): p-values
#                                                     keep their full dynamic range
#                                                     (no float32 underflow at 1e-38)
#   d    dictionary codes (uint16, or uint32 when "w": 4) + levels "k"
#   s    utf8 strings joined by "\n" (NA = "\001")   b  logical (uint8, 255 = NA)

# ---- binary helpers -------------------------------------------------------------

.xq_pad8 <- function(r) {
  pad <- (8L - length(r) %% 8L) %% 8L
  if (pad) c(r, raw(pad)) else r
}

.xq_bin <- function(x, size) {
  writeBin(x, raw(), size = size, endian = "little")
}

#' Does this numeric column hold p / q values (stored as -log10)?
#' @keywords internal
.xq_is_p_col <- function(name, v) {
  if (!is.numeric(v)) return(FALSE)
  if (grepl("cutoff|threshold|ratio|score|stat", name, ignore.case = TRUE)) return(FALSE)
  hit <- grepl(
    paste0("^(p|q)$|pval|p[._]?value|padj|p[._]?adj|adj[._ ]?p|fdr|qval|q[._]?val|",
           "^p(g|nde|pert)($|fdr|fwer)|^pg"),
    name, ignore.case = TRUE
  )
  if (!hit) return(FALSE)
  f <- v[is.finite(v)]
  length(f) > 0L && all(f >= 0 & f <= 1)
}

#' Encode one data.frame column
#' @return list(spec = list(...), raw = raw vector)
#' @keywords internal
.xq_encode_col <- function(name, v, n, allow_p = TRUE) {
  if (is.factor(v)) v <- as.character(v)
  if (is.list(v)) v <- vapply(v, function(z) paste(as.character(z), collapse = ";"), character(1))
  if (is.logical(v)) {
    code <- ifelse(is.na(v), 255L, as.integer(v))
    return(list(spec = list(t = "b"), raw = as.raw(code)))
  }
  if (is.numeric(v)) {
    v <- as.numeric(v)
    if (allow_p && .xq_is_p_col(name, v)) {
      nl <- -log10(pmax(v, 1e-320))
      nl[is.na(v)] <- NA_real_
      return(list(spec = list(t = "nlp"), raw = .xq_bin(nl, 4L)))
    }
    fin <- v[is.finite(v)]
    is_int <- length(fin) > 0L && all(fin == round(fin)) &&
      all(abs(fin) < 2147483647) && !any(is.infinite(v))
    if (is_int) {
      has_na <- anyNA(v)
      iv <- as.integer(v)   # NA_integer_ is exactly -2^31, the NA marker of the i32 type
      sp <- list(t = "i32")
      if (has_na) sp$na <- TRUE
      return(list(spec = sp, raw = .xq_bin(iv, 4L)))
    }
    # continuous values: float32, unless a huge integer would lose digits
    if (length(fin) > 0L && any(abs(fin) >= 2^24 & fin == round(fin))) {
      return(list(spec = list(t = "f64"), raw = .xq_bin(v, 8L)))
    }
    return(list(spec = list(t = "f32"), raw = .xq_bin(v, 4L)))
  }
  v <- as.character(v)
  v <- enc2utf8(v)
  lev <- unique(v[!is.na(v)])
  if (n > 8L && length(lev) <= 0.5 * n && length(lev) < 65535L) {
    code <- match(v, lev) - 1L
    code[is.na(code)] <- 65535L
    return(list(
      spec = list(t = "d", k = as.list(lev), w = 2L),
      raw = .xq_bin(code, 2L)
    ))
  }
  s <- gsub("[\r\n]+", " ", v)
  s[is.na(v)] <- "\001"
  list(spec = list(t = "s"), raw = charToRaw(paste(s, collapse = "\n")))
}

#' Encode a data.frame as the zlib-compressed columnar container (base64)
#'
#' @param df A data.frame (list columns are collapsed with `;`).
#' @param max_rows Optional hard row cap.
#' @param no_p Column names that must never be stored as -log10(p) (matrix
#'   columns whose names merely look like p-value names).
#' @return A list with `b64` (payload) and `meta` (n + `[name, type]` pairs).
#' @keywords internal
.xq_encode_table <- function(df, max_rows = Inf, no_p = character(0)) {
  if (is.null(df) || !is.data.frame(df) || ncol(df) == 0L) return(NULL)
  if (nrow(df) == 0L) return(NULL)
  if (nrow(df) > max_rows) df <- df[seq_len(as.integer(max_rows)), , drop = FALSE]
  n <- nrow(df)
  names_df <- colnames(df)
  if (is.null(names_df)) names_df <- paste0("V", seq_len(ncol(df)))
  specs <- vector("list", ncol(df))
  blocks <- vector("list", ncol(df))
  off <- 0
  for (j in seq_len(ncol(df))) {
    enc <- .xq_encode_col(names_df[[j]], df[[j]], n, allow_p = !(names_df[[j]] %in% no_p))
    r <- enc$raw
    sp <- enc$spec
    sp$name <- names_df[[j]]
    sp$o <- off
    sp$l <- length(r)
    r <- .xq_pad8(r)
    off <- off + length(r)
    specs[[j]] <- sp
    blocks[[j]] <- r
  }
  hdr <- charToRaw(enc2utf8(jsonlite::toJSON(
    list(n = n, cols = specs), auto_unbox = TRUE, digits = NA, null = "null", na = "null"
  )))
  pre <- .xq_pad8(c(.xq_bin(length(hdr), 4L), hdr))
  payload <- do.call(c, c(list(pre), blocks))
  comp <- memCompress(payload, type = "gzip")
  list(
    b64 = jsonlite::base64_enc(comp),
    meta = list(
      n = n,
      cols = lapply(specs, function(s) list(s$name, s$t))
    ),
    bytes = length(comp)
  )
}

#' HTML `<script>` block for one embedded table
#' @keywords internal
.xq_table_script <- function(key, enc) {
  meta <- jsonlite::toJSON(enc$meta, auto_unbox = TRUE, null = "null")
  meta <- gsub("'", "&#39;", as.character(meta), fixed = TRUE)
  paste0(
    "<script type=\"application/x-xq-table\" id=\"xqt-", key, "\" data-meta='",
    meta, "'>", enc$b64, "</script>"
  )
}

#' HTML `<script>` block for a small JSON payload (`XQ.json(key)`)
#' @keywords internal
.xq_json_script <- function(key, value) {
  json <- as.character(jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", na = "null", digits = NA
  ))
  # never let the payload close its own <script> element
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  json <- gsub("<!--", "<\\!--", json, fixed = TRUE)
  paste0("<script type=\"application/json\" id=\"xqj-", key, "\">", json, "</script>")
}

#' Safe element id fragment
#' @keywords internal
.xq_key <- function(...) {
  k <- paste(..., sep = "_")
  gsub("[^A-Za-z0-9_]+", "_", k)
}
