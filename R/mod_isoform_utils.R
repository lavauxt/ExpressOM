# mod_isoform_utils.R -- small shared utilities for isoform report/plot generation
#

#' Convert a PDF plot to a PNG alongside it
#' @keywords internal
convert_pdf_to_png <- function(pdf_file, dpi = 200) {
  if (!file.exists(pdf_file)) {
    warning("PDF not found: ", pdf_file)
    return(NULL)
  }

  sz <- file.info(pdf_file)$size

  if (is.na(sz) || sz < 2000) {
    warning("PDF is empty or too small (<2KB), skipping conversion: ", pdf_file)
    return(NULL)
  }

  png_file <- sub("\\.pdf$", ".png", pdf_file)

  if (requireNamespace("pdftools", quietly = TRUE)) {
    produced <- tryCatch(
      {
        pdftools::pdf_convert(
          pdf = pdf_file,
          filenames = png_file,
          dpi = dpi,
          verbose = FALSE
        )
      },
      error = function(e) {
        warning("pdftools failed to convert ", pdf_file, " to PNG: ", conditionMessage(e))
        character(0)
      }
    )

    existing <- produced[file.exists(produced)]

    if (length(existing) > 0) {
      if (length(existing) > 1) {
        message(
          "   -> ", basename(pdf_file), " has ", length(existing),
          " pages; using page 1 (", basename(existing[1]), ") for report embedding."
        )
      }

      return(normalizePath(existing[1], winslash = "/", mustWork = FALSE))
    }
  }

  if (requireNamespace("magick", quietly = TRUE)) {
    ok <- tryCatch(
      {
        img <- magick::image_read(
          pdf_file,
          density = sprintf("%dx%d", dpi, dpi)
        )

        magick::image_write(
          image = img,
          path = png_file,
          format = "png"
        )

        file.exists(png_file)
      },
      error = function(e) {
        warning("magick failed to convert ", pdf_file, " to PNG: ", conditionMessage(e))
        FALSE
      }
    )

    if (isTRUE(ok) && file.exists(png_file)) {
      return(normalizePath(png_file, winslash = "/", mustWork = FALSE))
    }
  }

  warning(
    "Could not convert ", pdf_file, " to PNG. ",
    "Install 'pdftools' (recommended) or 'magick' plus Ghostscript. ",
    "The PDF still exists and can be opened directly."
  )

  NULL
}

