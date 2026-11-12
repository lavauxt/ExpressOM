# mod_shiny.R -- the ExpressOM Shiny application.

# ---- parameter panel -----------------------------------------------------------------------------------

.app_spec <- function(name, label, type, group, choices = NULL, help = NULL) {
  list(name = name, label = label, type = type, group = group, choices = choices, help = help)
}

.app_param_spec <- function() {
  list(
    .app_spec("data_dir", "Quantification directory", "text", "Input", help = "Folder with one sub-folder per sample (quant.sf, abundance.tsv, ...)."),
    .app_spec("sample_table", "Sample table (CSV / TSV)", "text", "Input"),
    .app_spec("count_type", "Count type", "select", "Input", c("salmon", "kallisto", "rsem", "stringtie", "matrix")),
    .app_spec("matrix_file", "Count matrix file", "text", "Input", help = "Only for count type 'matrix'."),
    .app_spec("ensembl_package_name", "EnsDb package", "text", "Input", help = "Built automatically if it is not installed."),
    .app_spec("out_dir", "Output directory", "text", "Input"),
    .app_spec("subset_sample", "Sample filter", "text", "Input", help = "R expression on the sample table, e.g. cell_type == 'T_cells'."),
    .app_spec("remove_sample", "Samples to drop", "csv", "Input", help = "Comma-separated sample IDs."),

    .app_spec("level", "Level (numerator)", "text", "Design"),
    .app_spec("base", "Base (denominator)", "text", "Design"),
    .app_spec("model", "Design formula", "text", "Design"),
    .app_spec("test", "Test", "select", "Design", c("Wald", "LRT")),
    .app_spec("reduced", "Reduced formula (LRT)", "text", "Design"),
    .app_spec("batch_col", "Batch column", "text", "Design"),
    .app_spec("replicate_col", "Replicate column", "text", "Design"),
    .app_spec("group_col", "Group column (PCA colours)", "text", "Design"),
    .app_spec("shrink_method", "LFC shrinkage", "select", "Design", c("ashr", "apeglm", "normal")),
    .app_spec("run_dge", "Gene-level DGE", "check", "Design"),

    .app_spec("padj_cutoff", "Adjusted p-value cut-off", "number", "Cut-offs and figures"),
    .app_spec("deg_lfc_cutoff", "|log2FC| cut-off", "number", "Cut-offs and figures"),
    .app_spec("top_genes", "Top genes plotted", "number", "Cut-offs and figures"),
    .app_spec("pca_ntop", "Genes used for the PCA", "number", "Cut-offs and figures"),
    .app_spec("highlight_genes", "Genes to highlight", "csv", "Cut-offs and figures", help = "Comma-separated symbols marked on the volcano plot."),
    .app_spec("zscore_genes", "Genes for Z-score heatmaps", "csv", "Cut-offs and figures"),

    .app_spec("gsea_metric", "GSEA ranking metric", "select", "Functional analysis", c("stat", "signed_pval", "log2FoldChange")),
    .app_spec("gmt_file", "Gene-set collections / GMT files", "csv", "Functional analysis", help = "MSigDB codes (C2:CP:KEGG_LEGACY) or paths to GMT files."),
    .app_spec("ora_padj_cutoff", "ORA adjusted p-value cut-off", "number", "Functional analysis"),
    .app_spec("ora_lfc_cutoff", "ORA |log2FC| floor", "number", "Functional analysis"),
    .app_spec("ora_min_genes", "ORA minimum genes", "number", "Functional analysis"),
    .app_spec("run_spia", "SPIA pathway impact (slow)", "check", "Functional analysis"),
    .app_spec("run_decouple", "decoupleR TF / pathway activity", "check", "Functional analysis"),
    .app_spec("decouple_method", "decoupleR method", "select", "Functional analysis", c("ulm", "mlm", "wsum")),
    .app_spec("decouple_min_size", "decoupleR minimum targets", "number", "Functional analysis"),

    .app_spec("run_isoform", "Isoform analysis (DTE / DTU / switches)", "check", "Isoforms"),
    .app_spec("isoform_test_engine", "DTU engine", "select", "Isoforms", c("DEXSeq", "DRIMSeq", "satuRn")),
    .app_spec("run_dexseq", "Also run DEXSeq", "check", "Isoforms"),
    .app_spec("run_isoform_enrichment", "Enrichment on DTE / DTU", "check", "Isoforms"),
    .app_spec("run_predictors", "CPAT / Pfam / SignalP", "check", "Isoforms"),
    .app_spec("isoform_fasta", "Transcript FASTA", "text", "Isoforms"),
    .app_spec("isoform_gff", "Annotation GTF / GFF", "text", "Isoforms"),

    .app_spec("make_report", "Build the HTML report", "check", "Report")
  )
}

.app_to_input <- function(spec, value) {
  switch(
    spec$type,
    text = if (is.null(value) || !length(value) || is.na(value[[1L]])) "" else as.character(value)[[1L]],
    csv = if (is.null(value)) "" else paste(as.character(value), collapse = ", "),
    number = if (is.null(value) || !length(value)) NA_real_ else suppressWarnings(as.numeric(value)[[1L]]),
    select = if (!is.null(value) && as.character(value)[[1L]] %in% spec$choices) as.character(value)[[1L]] else spec$choices[[1L]],
    check = isTRUE(value)
  )
}

.app_from_input <- function(spec, x, default) {
  if (is.null(x)) return(default)
  switch(
    spec$type,
    text = {
      x <- trimws(as.character(x)[[1L]])
      if (nzchar(x)) x else NULL
    },
    csv = {
      parts <- trimws(unlist(strsplit(as.character(x)[[1L]], "[,;\n]")))
      parts <- parts[nzchar(parts)]
      if (length(parts)) parts else NULL
    },
    number = {
      x <- suppressWarnings(as.numeric(x)[[1L]])
      if (is.na(x)) NULL else x
    },
    select = as.character(x)[[1L]],
    check = isTRUE(x)
  )
}

.app_report_inputs <- function() {
  list(
    r_title = "Report title",
    r_max_points = "Max points in volcano / MA",
    r_top_terms = "Terms per enrichment plot",
    r_heatmap_genes = "Genes in the heatmap",
    r_static_figures = "Embed the static figures (PDF)"
  )
}

.app_report_settings_from_inputs <- function(values) {
  def <- expressom_report_settings()
  out <- list()
  num <- function(nm, key) {
    v <- suppressWarnings(as.numeric(values[[nm]]))
    if (length(v) == 1L && !is.na(v) && !identical(as.numeric(def[[key]]), v)) {
      out[[key]] <<- v
    }
  }
  title <- trimws(as.character(values$r_title %||% ""))
  if (nzchar(title) && !identical(title, def$title)) {
    out$title <- title
  }
  num("r_max_points", "max_points")
  num("r_top_terms", "top_terms")
  num("r_heatmap_genes", "heatmap_genes")
  if (!is.null(values$r_static_figures) && !identical(isTRUE(values$r_static_figures), isTRUE(def$static_figures))) {
    out$static_figures <- isTRUE(values$r_static_figures)
  }
  if (length(out)) out else NULL
}

.app_parse_extra <- function(text) {
  text <- trimws(text %||% "")
  if (!nzchar(text)) return(list())
  out <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = TRUE),
    error = function(e) {
      stop("The 'other parameters' box is not valid JSON: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.list(out) || is.null(names(out)) || any(!nzchar(names(out)))) {
    stop("The 'other parameters' box must be a JSON object, e.g. {\"zscore_genes\": [\"TP53\"]}.", call. = FALSE)
  }
  out
}

.app_build_params <- function(values, base, extra = "", spec = .app_param_spec()) {
  p <- base
  for (s in spec) {
    if (s$name %in% names(values)) {
      p[s$name] <- list(.app_from_input(s, values[[s$name]], base[[s$name]]))
    }
  }
  rs <- .app_report_settings_from_inputs(values)
  if (any(names(values) %in% names(.app_report_inputs()))) {
    p["report_settings"] <- list(rs)
  }
  ex <- .app_parse_extra(extra)
  unknown <- setdiff(names(ex), names(formals(expressom)))
  if (length(unknown)) {
    stop(
      "Unknown parameter(s) in the 'other parameters' box: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  for (nm in names(ex)) {
    p[nm] <- list(ex[[nm]])
  }
  .absolutize_paths(p)
}

.app_validate_params <- function(p) {
  bad <- character()
  need <- function(nm, label = nm) {
    v <- p[[nm]]
    if (is.null(v) || !length(v) || !nzchar(as.character(v)[[1L]])) {
      bad <<- c(bad, paste0("'", label, "' is empty."))
    }
  }
  need("out_dir", "Output directory")
  if (identical(p$count_type, "matrix")) {
    need("matrix_file", "Count matrix file")
  } else {
    need("data_dir", "Quantification directory")
  }
  need("sample_table", "Sample table")
  need("ensembl_package_name", "EnsDb package")
  full <- (!isTRUE(p$eda_only) && isTRUE(p$run_dge)) || isTRUE(p$run_isoform)
  if (full) {
    need("level", "Level")
    need("base", "Base")
    need("model", "Design formula")
    if (!is.null(p$level) && identical(p$level, p$base)) {
      bad <- c(bad, "Level and base must differ.")
    }
    if (!is.null(p$model) && is.null(tryCatch(stats::as.formula(p$model), error = function(e) NULL))) {
      bad <- c(bad, paste0("The design formula '", p$model, "' cannot be parsed."))
    }
    if (identical(p$test, "LRT") && is.null(p$reduced)) {
      bad <- c(bad, "The LRT needs a reduced formula.")
    }
  }
  a <- p$padj_cutoff
  if (is.null(a) || !is.numeric(a) || is.na(a) || a <= 0 || a > 1) {
    bad <- c(bad, "The adjusted p-value cut-off must be in (0, 1].")
  }
  if (!is.null(p$sample_table) && nzchar(p$sample_table) && !file.exists(p$sample_table)) {
    bad <- c(bad, paste0("Sample table not found: ", p$sample_table))
  }
  bad
}

.app_params_for_dir <- function(out_dir) {
  if (is.null(out_dir) || !nzchar(out_dir)) return(NULL)
  pf <- file.path(expressom_targets_dir(out_dir), "params.rds")
  saved <- .read_rds_safe(pf)
  if (is.list(saved)) return(saved)
  tryCatch(expressom_params_from_run(out_dir), error = function(e) NULL)
}

.app_initial_params <- function(out_dir = NULL, params = NULL) {
  base <- NULL
  if (!is.null(out_dir)) {
    base <- .app_params_for_dir(out_dir)
    if (is.null(base)) {
      base <- expressom_params(out_dir = out_dir)
    }
  }
  if (is.null(base)) {
    base <- expressom_default_params()
  }
  if (!is.null(params)) {
    base <- do.call(expressom_params, c(params, list(defaults = base)))
  }
  .absolutize_paths(base)
}

.app_param_widget <- function(spec, value, ns) {
  id <- ns(paste0("p_", spec$name))
  v <- .app_to_input(spec, value)
  w <- switch(
    spec$type,
    text = shiny::textInput(id, spec$label, value = v, width = "100%"),
    csv = shiny::textInput(id, spec$label, value = v, width = "100%"),
    number = shiny::numericInput(id, spec$label, value = v, step = NA, width = "100%"),
    select = shiny::selectInput(id, spec$label, choices = spec$choices, selected = v, width = "100%"),
    check = shiny::checkboxInput(id, spec$label, value = v)
  )
  if (!is.null(spec$help)) {
    w <- htmltools::tagAppendAttributes(w, title = spec$help)
  }
  w
}

.app_update_inputs <- function(session, p, spec = .app_param_spec()) {
  for (s in spec) {
    id <- paste0("p_", s$name)
    v <- .app_to_input(s, p[[s$name]])
    switch(
      s$type,
      text = shiny::updateTextInput(session, id, value = v),
      csv = shiny::updateTextInput(session, id, value = v),
      number = shiny::updateNumericInput(session, id, value = v),
      select = shiny::updateSelectInput(session, id, selected = v),
      check = shiny::updateCheckboxInput(session, id, value = v)
    )
  }
  rs <- utils::modifyList(expressom_report_settings(), as.list(p$report_settings))
  shiny::updateTextInput(session, "r_title", value = rs$title)
  shiny::updateNumericInput(session, "r_max_points", value = rs$max_points)
  shiny::updateNumericInput(session, "r_top_terms", value = rs$top_terms)
  shiny::updateNumericInput(session, "r_heatmap_genes", value = rs$heatmap_genes)
  shiny::updateCheckboxInput(session, "r_static_figures", value = isTRUE(rs$static_figures))
  shiny::updateTextAreaInput(session, "extra_json", value = "")
  invisible(NULL)
}

.app_params_ui <- function(id, init) {
  ns <- shiny::NS(id)
  tg <- htmltools::tags
  spec <- .app_param_spec()
  groups <- unique(vapply(spec, function(s) s$group, character(1)))
  rs <- utils::modifyList(expressom_report_settings(), as.list(init$report_settings))

  section <- function(title, ..., open = FALSE) {
    tg$details(class = "xp-section", open = if (open) NA else NULL, tg$summary(title), ...)
  }
  sections <- lapply(groups, function(g) {
    items <- Filter(function(s) identical(s$group, g), spec)
    extra <- if (identical(g, "Report")) {
      list(
        shiny::textInput(ns("r_title"), .app_report_inputs()$r_title, rs$title, width = "100%"),
        shiny::numericInput(ns("r_max_points"), .app_report_inputs()$r_max_points, rs$max_points, step = NA, width = "100%"),
        shiny::numericInput(ns("r_top_terms"), .app_report_inputs()$r_top_terms, rs$top_terms, step = NA, width = "100%"),
        shiny::numericInput(ns("r_heatmap_genes"), .app_report_inputs()$r_heatmap_genes, rs$heatmap_genes, step = NA, width = "100%"),
        shiny::checkboxInput(ns("r_static_figures"), .app_report_inputs()$r_static_figures, isTRUE(rs$static_figures))
      )
    }
    section(
      g,
      lapply(items, function(s) .app_param_widget(s, init[[s$name]], ns)),
      extra,
      open = g %in% c("Input", "Design")
    )
  })

  tg$div(
    class = "xp-side",
    tg$div(
      class = "xp-side-head",
      tg$strong("Parameters"),
      shiny::actionButton(
        ns("load_run"), "Load saved run", class = "btn-xs",
        title = "Fill the parameters from the run stored in the output directory"
      )
    ),
    tg$p(class = "xp-note", shiny::textOutput(ns("loaded"), inline = TRUE)),
    sections,
    section(
      "Other parameters",
      tg$p(class = "xp-note", "Any other argument of expressom() as a JSON object."),
      shiny::textAreaInput(
        ns("extra_json"), NULL, value = "", rows = 3, width = "100%",
        placeholder = '{"zscore_genes": ["TP53", "MYC"], "nBest": null}'
      )
    ),
    shiny::uiOutput(ns("problems"))
  )
}

.app_params_server <- function(id, init) {
  shiny::moduleServer(id, function(input, output, session) {
    spec <- .app_param_spec()
    base <- shiny::reactiveVal(init)
    loaded_msg <- shiny::reactiveVal("")

    params <- shiny::reactive({
      vals <- list()
      for (s in spec) {
        vals[s$name] <- list(input[[paste0("p_", s$name)]])
      }
      for (nm in names(.app_report_inputs())) {
        vals[nm] <- list(input[[nm]])
      }
      vals <- vals[!vapply(vals, is.null, logical(1))]
      tryCatch(
        .app_build_params(vals, base(), input$extra_json %||% ""),
        error = function(e) {
          structure(
            list(message = conditionMessage(e), call = NULL),
            class = c("xp_param_error", "error", "condition")
          )
        }
      )
    })

    reload <- function(p) {
      base(p)
      .app_update_inputs(session, p)
    }

    shiny::observeEvent(input$load_run, {
      od <- trimws(input$p_out_dir %||% "")
      found <- .app_params_for_dir(od)
      if (is.null(found)) {
        loaded_msg("")
        shiny::showNotification(
          paste0("No saved run found in '", od, "'."), type = "warning"
        )
      } else {
        found$out_dir <- .abs_path(od)
        reload(found)
        loaded_msg(paste0("Parameters loaded from ", found$out_dir))
      }
    })

    output$loaded <- shiny::renderText(loaded_msg())

    output$problems <- shiny::renderUI({
      p <- params()
      msgs <- if (inherits(p, "xp_param_error")) p$message else .app_validate_params(p)
      if (!length(msgs)) {
        return(NULL)
      }
      htmltools::tags$div(
        class = "xp-problems",
        htmltools::tags$strong("To fix before running:"),
        htmltools::tags$ul(lapply(msgs, htmltools::tags$li))
      )
    })

    list(params = params, reload = reload)
  })
}

# ---- pipeline tab -------------------------------------------------------------------------------------

.app_default_runner <- function(params, dir, stages, force) {
  if (length(force)) {
    ExpressOM::expressom_targets_invalidate(dir, force)
  }
  status <- ExpressOM::run_expressom_targets(
    params, dir = dir, stages = stages, callr = FALSE, reporter = "verbose"
  )
  as.character(status$label[status$state == "errored"])
}

.app_state_class <- function(state) {
  paste0("xp-pill xp-", gsub("[^a-z]+", "-", tolower(state)))
}

.app_stage_table <- function(status) {
  tg <- htmltools::tags
  if (is.null(status) || !nrow(status)) {
    return(tg$p(class = "xp-note", "Nothing to show yet."))
  }
  rows <- lapply(seq_len(nrow(status)), function(i) {
    secs <- status$seconds[[i]]
    tg$tr(
      tg$td(status$label[[i]]),
      tg$td(tg$span(class = .app_state_class(status$state[[i]]), status$state[[i]])),
      tg$td(if (is.na(secs)) "" else if (secs < 60) sprintf("%.1f s", secs) else sprintf("%.1f min", secs / 60)),
      tg$td(if (is.na(status$finished[[i]])) "" else status$finished[[i]])
    )
  })
  tg$table(
    class = "xp-stages",
    tg$thead(tg$tr(tg$th("Stage"), tg$th("State"), tg$th("Took"), tg$th("Last built"))),
    tg$tbody(rows)
  )
}

.app_strip_ansi <- function(x) gsub("\033\\[[0-9;]*[A-Za-z]", "", x)

.app_log_tail <- function(path, n = 400L) {
  if (is.null(path) || !file.exists(path)) {
    return("No run yet.")
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (!length(lines)) {
    return("(empty log)")
  }
  paste(.app_strip_ansi(utils::tail(lines, n)), collapse = "\n")
}

.app_pipeline_ui <- function(id) {
  ns <- shiny::NS(id)
  tg <- htmltools::tags
  tagList <- htmltools::tagList
  tagList(
    tg$p(
      class = "xp-note",
      "Each stage depends only on the parameters it uses: change a cut-off or a colour and only the stages that read it re-run; the report is rebuilt when something upstream changed."
    ),
    shiny::uiOutput(ns("buttons")),
    shiny::uiOutput(ns("stage_table")),
    tg$details(
      class = "xp-section",
      tg$summary("Choose stages"),
      shiny::checkboxGroupInput(
        ns("stages"), NULL,
        choiceNames = expressom_pipeline_stages()$label[-1L],
        choiceValues = expressom_pipeline_stages()$stage[-1L],
        inline = TRUE
      ),
      tg$p(
        class = "xp-note",
        "'Run selected' builds the chosen stages and what they need, nothing else. 'Force re-run' rebuilds them even if nothing changed, then refreshes everything downstream."
      )
    ),
    tg$h4("Log"),
    shiny::verbatimTextOutput(ns("log")),
    tg$script(htmltools::HTML(sprintf(
      "$(document).on('shiny:value', function(e) { if (e.name === '%s') { setTimeout(function() { var el = document.getElementById('%s'); if (el) el.scrollTop = el.scrollHeight; }, 50); } });",
      ns("log"), ns("log")
    )))
  )
}

.app_pipeline_server <- function(id, params, runner = NULL, poll_ms = 1000L) {
  runner <- runner %||% .app_default_runner
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    proc <- shiny::reactiveVal(NULL)
    version <- shiny::reactiveVal(0L)
    tick <- shiny::reactiveVal(0L)
    status_val <- shiny::reactiveVal(NULL)
    last_error <- shiny::reactiveVal(NULL)

    params_ok <- shiny::reactive({
      p <- params()
      if (inherits(p, "error")) NULL else p
    })
    pipe_dir <- shiny::reactive({
      p <- params_ok()
      if (is.null(p) || is.null(p$out_dir)) NULL else expressom_targets_dir(p$out_dir)
    })
    # The full pipeline run log lives under <out_dir>/Log/targets/ (not
    # under .expressom_targets/) so every diagnostic of a run is reachable
    # from the same Log/ folder as SessionInfo.txt, Warnings.txt,
    # DGE_params_*.json and WSL/predictor logs.
    log_file <- shiny::reactive({
      p <- params_ok()
      if (is.null(p) || is.null(p$out_dir)) {
        return(NULL)
      }
      file.path(p$out_dir, "Log", "targets", "last_run.log")
    })

    refresh_status <- function() {
      p <- shiny::isolate(params_ok())
      d <- shiny::isolate(pipe_dir())
      st <- if (is.null(p) || is.null(d)) {
        NULL
      } else {
        tryCatch(
          expressom_stage_status(p, dir = d, running = !is.null(shiny::isolate(proc()))),
          error = function(e) NULL
        )
      }
      if (!identical(st, shiny::isolate(status_val()))) {
        status_val(st)
      }
    }

    shiny::observe({
      params_ok()
      version()
      refresh_status()
    })

    finish <- function(pr) {
      res <- tryCatch(pr$get_result(), error = function(e) e)
      proc(NULL)
      if (inherits(res, "error")) {
        msg <- conditionMessage(res)
        last_error(msg)
        shiny::showNotification(
          paste0("The run stopped with an error (see the log). ", substr(msg, 1, 200)),
          type = "error", duration = NULL
        )
      } else if (is.character(res) && length(res)) {
        last_error(paste("Failed:", paste(res, collapse = ", ")))
        shiny::showNotification(
          paste0("Run finished, but these stages failed: ", paste(res, collapse = ", "), " (see the log)."),
          type = "warning", duration = NULL
        )
      } else {
        last_error(NULL)
        shiny::showNotification("Run finished.", type = "message")
      }
      version(shiny::isolate(version()) + 1L)
    }

    shiny::observe({
      pr <- proc()
      if (is.null(pr)) {
        return()
      }
      shiny::invalidateLater(poll_ms)
      tick(shiny::isolate(tick()) + 1L)
      refresh_status()
      if (!pr$is_alive()) {
        finish(pr)
      }
    })

    launch <- function(stages = NULL, force = character()) {
      p <- shiny::isolate(params_ok())
      d <- shiny::isolate(pipe_dir())
      if (is.null(p) || is.null(d)) {
        shiny::showNotification("Fix the parameters first.", type = "error")
        return(invisible(FALSE))
      }
      if (!is.null(shiny::isolate(proc()))) {
        shiny::showNotification("A run is already in progress.", type = "warning")
        return(invisible(FALSE))
      }
      problems <- .app_validate_params(p)
      if (length(problems)) {
        shiny::showNotification(paste(problems, collapse = " "), type = "error", duration = 10)
        return(invisible(FALSE))
      }
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      # targets store stays in .expressom_targets/, but the human-facing
      # run log goes to Log/targets/ next to the rest of the diagnostics.
      targets_log_dir <- file.path(p$out_dir, "Log", "targets")
      dir.create(targets_log_dir, recursive = TRUE, showWarnings = FALSE)
      log <- file.path(normalizePath(targets_log_dir, winslash = "/"), "last_run.log")
      pr <- callr::r_bg(
        runner,
        args = list(params = p, dir = d, stages = stages, force = force),
        stdout = log,
        stderr = "2>&1",
        supervise = TRUE
      )
      last_error(NULL)
      proc(pr)
      invisible(TRUE)
    }

    shiny::observeEvent(input$run_all, launch())
    shiny::observeEvent(input$run_sel, {
      if (!length(input$stages)) {
        shiny::showNotification("Tick at least one stage under 'Choose stages'.", type = "warning")
      } else {
        launch(stages = input$stages)
      }
    })
    shiny::observeEvent(input$force_sel, {
      if (!length(input$stages)) {
        shiny::showNotification("Tick at least one stage under 'Choose stages'.", type = "warning")
      } else {
        launch(force = input$stages)
      }
    })
    shiny::observeEvent(input$cancel, {
      pr <- proc()
      if (!is.null(pr)) {
        pr$kill()
        proc(NULL)
        last_error(NULL)
        shiny::showNotification("Run cancelled.", type = "warning")
        version(shiny::isolate(version()) + 1L)
      }
    })

    output$buttons <- shiny::renderUI({
      st <- status_val()
      tg <- htmltools::tags
      if (!is.null(proc())) {
        return(tg$div(
          class = "xp-run-row",
          shiny::actionButton(ns("cancel"), "Cancel run", class = "btn-danger"),
          tg$span(class = "xp-note", "Running in the background...")
        ))
      }
      if (is.null(st)) {
        return(tg$p(class = "xp-note", "Fix the parameters on the left to enable the pipeline."))
      }
      pending <- sum(st$state %in% c("stale", "never run", "errored"))
      main <- if (pending > 0L) {
        shiny::actionButton(
          ns("run_all"),
          sprintf("Run pipeline (%d stage%s out of date)", pending, if (pending == 1L) "" else "s"),
          class = "btn-danger"
        )
      } else {
        shiny::actionButton(ns("run_all"), "Up to date - run anyway", class = "btn-success")
      }
      tg$div(
        class = "xp-run-row",
        main,
        shiny::actionButton(ns("run_sel"), "Run selected", class = "btn-default"),
        shiny::actionButton(ns("force_sel"), "Force re-run selected", class = "btn-default"),
        if (!is.null(last_error())) tg$span(class = "xp-err", "Last run failed - see the log.")
      )
    })

    output$stage_table <- shiny::renderUI(.app_stage_table(status_val()))

    output$log <- shiny::renderText({
      tick()
      version()
      .app_log_tail(log_file())
    })

    list(
      version = version,
      running = shiny::reactive(!is.null(proc())),
      status = status_val
    )
  })
}

# ---- results tab --------------------------------------------------------------------------------------

.app_tree_ui <- function(nodes, prefix, ns = shiny::NS(NULL)) {
  leaves <- list()
  n <- 0L
  build <- function(nodes) {
    lapply(nodes, function(nd) {
      title <- .xp_clean_title(nd$title)
      if (.xp_is_leaf(nd)) {
        n <<- n + 1L
        lid <- paste0(prefix, "_", n)
        leaves[[lid]] <<- nd$content
        shiny::tabPanel(
          title,
          htmltools::tags$div(class = "xp-leaf", shiny::uiOutput(ns(lid)))
        )
      } else {
        shiny::tabPanel(
          title,
          do.call(shiny::tabsetPanel, c(list(type = "pills"), build(nd$children)))
        )
      }
    })
  }
  panels <- build(nodes)
  list(ui = do.call(shiny::tabsetPanel, c(list(type = "tabs"), panels)), leaves = leaves)
}

.app_safe_leaf <- function(fn) {
  tryCatch(
    fn(),
    error = function(e) {
      htmltools::tags$div(
        class = "xp-problems",
        "This tab could not be built: ", conditionMessage(e)
      )
    }
  )
}

.app_num_or_null <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || is.na(x)) NULL else x
}

.app_results_ui <- function(id) {
  ns <- shiny::NS(id)
  tg <- htmltools::tags
  tg$div(
    tg$div(
      class = "xp-viewer-bar",
      shiny::selectInput(ns("comparison"), "Comparison", choices = character(), width = "230px"),
      shiny::numericInput(ns("padj"), "padj cut-off (display)", value = NA, step = 0.001, width = "170px"),
      shiny::numericInput(ns("lfc"), "|log2FC| cut-off (display)", value = NA, step = 0.1, width = "170px"),
      shiny::textInput(ns("highlight"), "Highlight genes", value = "", width = "220px"),
      shiny::actionButton(ns("reload"), "Reload from disk", class = "btn-default")
    ),
    tg$p(
      class = "xp-note",
      "Display cut-offs only change what the figures colour and count. To recompute the differential expression with another cut-off, change it in the parameters and run the pipeline."
    ),
    shiny::uiOutput(ns("tabs"))
  )
}

.app_results_server <- function(id, out_dir, version = function() 0L) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    generation <- 0L

    comparisons <- shiny::reactive({
      version()
      input$reload
      od <- out_dir()
      if (is.null(od) || !dir.exists(od)) character() else list_expressom_comparisons(od)
    })

    shiny::observeEvent(comparisons(), {
      cmp <- comparisons()
      keep <- shiny::isolate(input$comparison)
      shiny::updateSelectInput(
        session, "comparison", choices = cmp,
        selected = if (!is.null(keep) && keep %in% cmp) keep else if (length(cmp)) cmp[[1L]] else NULL
      )
    }, ignoreNULL = FALSE)

    comparison <- shiny::reactive({
      cmp <- comparisons()
      sel <- input$comparison
      if (length(cmp) && !(isTRUE(nzchar(sel)) && sel %in% cmp)) {
        shiny::req(FALSE)
      }
      if (length(cmp)) sel else NULL
    })

    results <- shiny::reactive({
      version()
      input$reload
      od <- out_dir()
      shiny::validate(shiny::need(
        !is.null(od) && dir.exists(od),
        "No output directory yet: run the pipeline first."
      ))
      cmp <- comparison()
      shiny::withProgress(
        message = "Loading results...", value = 0.5,
        tryCatch(
          load_expressom_results(od, comparison = cmp),
          error = function(e) {
            shiny::validate(shiny::need(FALSE, paste0("No results to show yet: ", conditionMessage(e))))
          }
        )
      )
    })

    settings <- shiny::reactive({
      expressom_report_settings(
        padj_cutoff = .app_num_or_null(input$padj),
        lfc_cutoff = .app_num_or_null(input$lfc),
        highlight_genes = .app_from_input(
          .app_spec("h", "h", "csv", "x"), input$highlight, NULL
        )
      )
    })

    tree_ui <- shiny::reactive({
      res <- results()
      tree <- expressom_report_tabs(res, settings())
      generation <<- generation + 1L
      .app_tree_ui(tree, prefix = paste0("leaf", generation), ns = ns)
    })

    shiny::observeEvent(tree_ui(), {
      tu <- tree_ui()
      for (lid in names(tu$leaves)) {
        local({
          fn <- tu$leaves[[lid]]
          output[[lid]] <- shiny::renderUI(.app_safe_leaf(fn))
        })
      }
    })

    output$tabs <- shiny::renderUI(tree_ui()$ui)

    list(comparison = shiny::reactive(tryCatch(comparison(), error = function(e) NULL)), results = results)
  })
}

# ---- report tab -------------------------------------------------------------------------------------

# The report now lives at the ROOT of out_dir (ExpressOM_report_*.html), not
# under Report/. This finds the most recent one, optionally filtered by
# comparison.
.app_latest_report <- function(out_dir, comparison = NULL) {
  if (is.null(out_dir) || !dir.exists(out_dir)) return(NULL)
  files <- list.files(out_dir, pattern = "^ExpressOM_report.*\\.html$",
                      full.names = TRUE)
  if (!length(files)) return(NULL)
  if (!is.null(comparison)) {
    mine <- files[grepl(comparison, basename(files), fixed = TRUE)]
    if (length(mine)) files <- mine
  }
  files[which.max(file.mtime(files))]
}

.app_report_ui <- function(id) {
  ns <- shiny::NS(id)
  tg <- htmltools::tags
  tg$div(
    tg$div(
      class = "xp-run-row",
      shiny::actionButton(ns("build"), "Rebuild report now", class = "btn-primary"),
      shiny::downloadButton(ns("download"), "Download HTML", class = "btn-default"),
      shiny::textOutput(ns("note"), inline = TRUE)
    ),
    tg$p(
      class = "xp-note",
      "The report is a single self-contained HTML file written at the root of the output directory. Rebuilding only re-reads the results already on disk; it never re-runs an analysis."
    ),
    shiny::uiOutput(ns("frame"))
  )
}

.app_report_server <- function(id, params, out_dir, version = function() 0L, comparison = function() NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    proc <- shiny::reactiveVal(NULL)
    local_version <- shiny::reactiveVal(0L)
    note <- shiny::reactiveVal("")

    path <- shiny::reactive({
      version()
      local_version()
      .app_latest_report(out_dir(), comparison())
    })

    shiny::observeEvent(input$build, {
      od <- out_dir()
      if (is.null(od) || !dir.exists(od)) {
        shiny::showNotification("No results to put in a report yet.", type = "warning")
        return()
      }
      if (!is.null(proc())) {
        return()
      }
      p <- params()
      rs <- if (inherits(p, "error")) NULL else p$report_settings
      settings <- utils::modifyList(expressom_report_settings(), as.list(rs))
      note("Building the report...")
      proc(callr::r_bg(
        function(out_dir, comparison, settings) {
          ExpressOM::build_expressom_report(out_dir, comparison = comparison, settings = settings)
        },
        args = list(out_dir = od, comparison = comparison(), settings = settings),
        supervise = TRUE
      ))
    })

    shiny::observe({
      pr <- proc()
      if (is.null(pr)) {
        return()
      }
      shiny::invalidateLater(750L)
      if (!pr$is_alive()) {
        res <- tryCatch(pr$get_result(), error = function(e) e)
        proc(NULL)
        if (inherits(res, "error")) {
          note("The report could not be built.")
          shiny::showNotification(
            paste0("Report failed: ", substr(conditionMessage(res), 1, 300)),
            type = "error", duration = NULL
          )
        } else {
          note(paste("Report rebuilt at", format(Sys.time(), "%H:%M:%S")))
          local_version(shiny::isolate(local_version()) + 1L)
        }
      }
    })

    output$note <- shiny::renderText(note())

    output$frame <- shiny::renderUI({
      f <- path()
      if (is.null(f)) {
        return(htmltools::tags$p(
          class = "xp-note",
          "No report yet. Run the pipeline (with 'Build the HTML report' ticked) or press 'Rebuild report now'."
        ))
      }
      shiny::addResourcePath("xp_report", dirname(f))
      htmltools::tags$iframe(
        class = "xp-frame",
        src = paste0("xp_report/", utils::URLencode(basename(f)), "?v=", as.integer(file.mtime(f))),
        title = "ExpressOM report"
      )
    })

    output$download <- shiny::downloadHandler(
      filename = function() {
        f <- path()
        if (is.null(f)) "ExpressOM_report.html" else basename(f)
      },
      content = function(file) {
        f <- path()
        shiny::req(f)
        file.copy(f, file, overwrite = TRUE)
      }
    )

    list(path = path)
  })
}

# ---- the app ----------------------------------------------------------------------------------------

.app_css <- function() {
  paste0(
    .xp_css(),
    "\n",
    .xp_dark_css_text(),
    "
.xp-app { display:grid; grid-template-columns: minmax(300px, 360px) 1fr; gap: 22px; align-items: start; padding: 0 6px 30px; }
@media (max-width: 980px) { .xp-app { grid-template-columns: 1fr; } }
.xp-side { background:#f5f8fb; border:1px solid #d9e2ec; border-radius:10px; padding:12px 14px; position:sticky; top:8px; max-height: calc(100vh - 16px); overflow:auto; }
.xp-side-head { display:flex; justify-content:space-between; align-items:center; margin-bottom:4px; }
.xp-section { margin: 8px 0; border-top: 1px solid #d9e2ec; padding-top: 6px; }
.xp-section > summary { cursor:pointer; font-weight:650; color:#0b3d75; padding: 2px 0; }
.xp-section .form-group { margin-bottom: 8px; }
.xp-section label { font-weight: 500; font-size: 13px; margin-bottom: 2px; }
.xp-section .checkbox { margin: 4px 0; }
.xp-problems { background:#fdecea; border:1px solid #f5c2bd; color:#8a1c12; border-radius:8px; padding:8px 12px; margin-top:10px; font-size:13px; }
.xp-problems ul { margin:4px 0 0; padding-left:18px; }
.xp-run-row { display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin: 8px 0 14px; }
.xp-err { color:#b42318; font-weight:600; }
table.xp-stages { border-collapse:collapse; width:100%; margin-bottom:12px; font-size:14px; }
table.xp-stages th { text-align:left; color:#52606d; font-weight:600; border-bottom:2px solid #d9e2ec; padding:4px 10px; }
table.xp-stages td { padding:5px 10px; border-bottom:1px solid #e4ebf1; }
.xp-pill { display:inline-block; padding:1px 10px; border-radius:12px; font-size:12px; font-weight:650; }
.xp-current { background:#dcf3e6; color:#14663a; }
.xp-stale { background:#fdebd0; color:#9a5b00; }
.xp-never-run { background:#e4e7eb; color:#52606d; }
.xp-running { background:#d6e8fb; color:#0b5cad; }
.xp-errored { background:#fdecea; color:#b42318; }
#pipe-log { max-height: 380px; overflow:auto; font-size:12px; background:#101820; color:#d9e2ec; }
.xp-viewer-bar { display:flex; flex-wrap:wrap; gap:6px 14px; align-items:flex-end; }
.xp-viewer-bar .form-group { margin-bottom: 4px; }
.xp-leaf { min-height: 120px; }
iframe.xp-frame { width:100%; height: 82vh; border: 1px solid #d9e2ec; border-radius: 8px; }
a.xp-file-link { color: var(--xp-accent); text-decoration: underline; word-break: break-all; }
"
  )
}

#' Create the ExpressOM Shiny application
#'
#' Builds the app object without launching it (see [run_expressom_app()]).
#'
#' @inheritParams run_expressom_app
#' @param runner Function `(params, dir, stages, force)` executed in the
#'   background R process of each run. The default runs the `targets` pipeline
#'   ([run_expressom_targets()]); the tests use a stand-in.
#' @param poll_ms Polling interval (milliseconds) while a run is in progress.
#' @return A `shiny.appobj`.
#' @export
expressom_app <- function(out_dir = NULL, params = NULL, runner = NULL, poll_ms = 1000L) {
  for (pkg in c("shiny", "targets", "callr", "plotly", "DT", "htmltools", "htmlwidgets", "rmarkdown", "knitr")) {
    .rep_require(pkg, "to run the ExpressOM app")
  }
  init <- .app_initial_params(out_dir, params)
  tg <- htmltools::tags

  ui <- shiny::fluidPage(
    title = "ExpressOM",
    tg$head(tg$style(htmltools::HTML(.app_css()))),
    tg$div(
      class = "xp-header",
      tg$h1("ExpressOM"),
      tg$span(class = "xp-sub", "bulk RNA-seq: run, explore and report from one place")
    ),
    tg$div(
      class = "xp-app",
      .app_params_ui("params", init),
      tg$div(
        shiny::tabsetPanel(
          id = "main",
          type = "pills",
          shiny::tabPanel("Pipeline", .app_pipeline_ui("pipe")),
          shiny::tabPanel("Results", .app_results_ui("res")),
          shiny::tabPanel("Report", .app_report_ui("rep"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    pr <- .app_params_server("params", init)
    pl <- .app_pipeline_server("pipe", pr$params, runner = runner, poll_ms = poll_ms)
    out_dir_v <- shiny::reactiveVal(init$out_dir)
    shiny::observe({
      p <- pr$params()
      od <- if (inherits(p, "error")) NULL else p$out_dir
      if (!is.null(od) && !identical(od, shiny::isolate(out_dir_v()))) {
        out_dir_v(od)
      }
    })
    out_dir_r <- out_dir_v
    rs <- .app_results_server("res", out_dir_r, pl$version)
    .app_report_server("rep", pr$params, out_dir_r, pl$version, rs$comparison)
  }

  shiny::shinyApp(ui, server)
}

#' Launch the ExpressOM Shiny application
#'
#' An application to (re)run the analysis on the fly and explore it: set the
#' parameters, see which stages are out of date, run the `targets` pipeline in
#' the background (only the stages whose inputs changed are recomputed), follow
#' the log, browse the interactive results and rebuild or download the HTML
#' report.
#'
#' @param out_dir Output directory of an existing run: its parameters pre-fill
#'   the panel and its results are shown. `NULL` starts from the defaults.
#' @param params Optional named list of [expressom()] arguments (see
#'   [expressom_params()]) applied on top of `out_dir`'s saved run, or on top of
#'   the defaults when there is none.
#' @param launch.browser,port Passed to [shiny::runApp()].
#' @param ... Further arguments for [expressom_app()].
#' @return Called for its side effect (blocks until the app is closed).
#' @export
#' @examples
#' \dontrun{
#' run_expressom_app("results")
#' run_expressom_app(params = list(data_dir = "data", sample_table = "samples.csv",
#'                                 out_dir = "results", level = "treated", base = "control"))
#' }
run_expressom_app <- function(
  out_dir = NULL,
  params = NULL,
  launch.browser = interactive(),
  port = getOption("shiny.port"),
  ...
) {
  app <- expressom_app(out_dir = out_dir, params = params, ...)
  shiny::runApp(app, launch.browser = launch.browser, port = port)
}