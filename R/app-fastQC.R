###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastQC <- function(input = NA, output = NA, param = NA) {
  # support for ubam
  # isUBam <- input$readType() == "bam"
  # if (isTRUE(isUBam)) {
  #   if (isTRUE(param$perLibrary)) {
  #     fastqInput <- ezMethodBam2Fastq(
  #       input = input, param = param,
  #       OUTPUT_PER_RG = FALSE
  #     )
  #   } else {
  #     ## We only support one uBam when it's per cell mode
  #     stopifnot(input$getLength() == 1L)
  #     fastqInput <- ezMethodBam2Fastq(
  #       input = input, param = param,
  #       OUTPUT_PER_RG = TRUE
  #     )
  #   }
  #   input <- fastqInput$copy()
  # }

  ## trim the reads to maximum length. Useful for PacBio/ONT long reads with variable read length and very few ultralong reads
  if (ezIsSpecified(param$max_len1) && param$max_len1 > 0) {
    input <- ezMethodFastpTrim(input = input, param = param)
  }

  dataset <- input$meta
  samples <- rownames(dataset)

  if (any(grepl(',', dataset$Read1))) {
    ##local copy of data
    outputFiles <- c()
    for (j in 1:length(samples)) {
      files <- file.path(
        param$dataRoot,
        limma::strsplit2(dataset$Read1[j], ',')
      )
      outputFile <- paste0(samples[j], '_R1.fastq.gz')
      ezSystem(paste('touch', outputFile))
      sapply(files, function(x) system(paste("cat", x, " >>", outputFile)))
      outputFiles <- c(outputFiles, outputFile)
    }
    input$setColumn("Read1", file.path(getwd(), outputFiles))
    if (isTRUE(param$paired)) {
      outputFiles <- c()
      for (j in 1:length(samples)) {
        files <- file.path(
          param$dataRoot,
          limma::strsplit2(dataset$Read2[j], ',')
        )
        outputFile <- paste0(samples[j], '_R2.fastq.gz')
        ezSystem(paste('touch', outputFile))
        sapply(files, function(x) system(paste("cat", x, " >>", outputFile)))
        outputFiles <- c(outputFiles, outputFile)
      }
      input$setColumn("Read2", file.path(getwd(), outputFiles))
    }
    ##update input
    input <- EzDataset$new(meta = input$meta, dataRoot = '')
    dataset <- input$meta
  }

  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset

  if (is.null(dataset[['Read Count']])) {
    dataset[['Read Count']] <- countReadsInFastq(input$getFullPaths("Read1"))
  }
  readCount <- signif(dataset$"Read Count" / 1e6, digits = 3)
  names(readCount) <- rownames(dataset)

  ans4Report[["Read Counts"]] <- readCount

  if (sum(dataset$`Read Count`) > 1e9) {
    input <- ezMethodSubsampleFastq(input = input, param = param, n = 1e6)
    dataset <- input$meta
  }

  setwdNew(basename(output$getColumn("FastQC")))

  files <- c()
  for (sm in samples) {
    files[paste0(sm, "_R1")] <- input$getFullPaths("Read1")[sm]
    if (isTRUE(param$paired)) {
      files[paste0(sm, "_R2")] <- input$getFullPaths("Read2")[sm]
    }
  }
  nFiles <- length(files)

  ## guess the names of the report directories that will be creatd by fastqc
  reportDirs <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files))
  stopifnot(!any(duplicated(reportDirs)))
  filesUse <- files[!file.exists(reportDirs)]
  if (length(filesUse) > 0) {
    cmd <- paste(
      "fastqc",
      "--extract -o . -t",
      min(param$cores, 8),
      "-a",
      FASTQC_ADAPTERS,
      "--kmers 7",
      "--dir .",
      "-q",
      param$cmdOptions,
      paste(filesUse, collapse = " "),
      "> fastqc.out",
      "2> fastqc.err"
    )

    if (length(filesUse) > 384) {
      cat(cmd, file = 'fastqcCall.sh')
      result <- ezSystem('bash fastqcCall.sh')
    } else {
      result <- ezSystem(cmd)
    }
    gc()
  }

  if (ezIsSpecified(param$showNativeReports) && param$showNativeReports) {
    statusToPng <- c(
      PASS = "tick.png",
      WARN = "warning.png",
      FAIL = "error.png"
    )
    ## collect the overview table
    plots <- c(
      "Per base sequence quality" = "per_base_quality.png",
      "Per sequence quality scores" = "per_sequence_quality.png",
      "Per tile sequence quality" = "per_tile_quality.png",
      "Per base sequence content" = "per_base_sequence_content.png",
      "Per sequence GC content" = "per_sequence_gc_content.png",
      "Per base N content" = "per_base_n_content.png",
      "Sequence Length Distribution" = "sequence_length_distribution.png",
      "Sequence Duplication Levels" = "duplication_levels.png",
      "Adapter Content" = "adapter_content.png"
      #      "Kmer Content" = "kmer_profiles.png" ## kmers are sometimes mssing
    )

    ## make for each plot type an html report with all samples
    file.copy(
      system.file("templates/FastQC_overview.Rmd", package = "ezRun"),
      "FastQC_overview.Rmd"
    )
    plotPages <- sub(".png", ".html", plots)
    for (i in 1:length(plots)) {
      plotPage <- plotPages[i]
      pngs <- file.path(reportDirs, "Images", plots[i])
      rmarkdown::render(
        input = "FastQC_overview.Rmd",
        envir = new.env(),
        output_dir = ".",
        output_file = plotPage,
        quiet = TRUE
      )
    }

    ## Each sample can have different number of reports.
    ## Especially per tile sequence quality
    nrReports <- sapply(
      reportDirs,
      function(x) {
        smy <- ezRead.table(
          file.path(x, "summary.txt"),
          row.names = NULL,
          header = FALSE
        )
        nrow(smy)
      }
    )
    i <- which.max(nrReports)
    smy <- ezRead.table(
      file.path(reportDirs[i], "summary.txt"),
      row.names = NULL,
      header = FALSE
    )
    rowNames <- paste0(
      "<a href=",
      reportDirs,
      "/fastqc_report.html>",
      names(files),
      "</a>"
    )

    tbl <- ezMatrix("", rows = rowNames, cols = smy[[2]])
    for (i in 1:nFiles) {
      smy <- ezRead.table(
        file.path(reportDirs[i], "summary.txt"),
        row.names = NULL,
        header = FALSE
      )

      href <- paste0(
        reportDirs[i],
        "/fastqc_report.html#M",
        0:(ncol(tbl) - 1)
      )[colnames(tbl) %in% smy[[2]]]
      #img <- paste0(reportDirs[1], "/Icons/", statusToPng[smy[[1]]])
      tbl[i, colnames(tbl) %in% smy[[2]]] <- paste0(
        "<a href=",
        href,
        ">",
        smy[[1]],
        "</a>" #<img src=", img, "></a>"
      )
    }
    colnames(tbl) <- ifelse(
      colnames(tbl) %in% names(plotPages),
      paste0(
        "<a href=",
        plotPages[colnames(tbl)],
        ">",
        colnames(tbl),
        "</a>"
      ),
      colnames(tbl)
    )

    ans4Report[["Fastqc quality measures"]] <- tbl

    # gc()
    # qualMatrixList <- ezMclapply(files, getQualityMatrix, mc.cores = param$cores)
    #ans4Report[["Per Base Read Quality"]] <- qualMatrixList

    ## generate the main reports
    file.copy(
      system.file("templates/FastQC.Rmd", package = "ezRun"),
      "FastQC.Rmd"
    )
    rmarkdown::render(
      input = "FastQC.Rmd",
      envir = new.env(),
      output_dir = ".",
      output_file = basename(output$getColumn("FastQC Report")),
      quiet = TRUE
    )
    unlink(paste0(reportDirs, ".zip"), recursive = FALSE)
  } else {
    unlink(reportDirs, recursive = TRUE)
    unlink(paste0(reportDirs, ".html"), recursive = FALSE)
  }

  ## ==== AI feature flags + hardcoded LLM config ====
  ## The endpoint is sensitive: it is hardcoded here, kept out of all SLURM logs,
  ## stripped from MultiQC's stdout/log/HTML/JSON outputs, and never echoed by
  ## any message() call below. Only the model name is exposed publicly.
  AI_PROVIDER       <- "custom"
  AI_MODEL          <- "Qwen3.6-27B-FP8"
  AI_ENDPOINT       <- "http://fgcz-c-056:8081/v1/chat/completions"
  AI_CONTEXT_WINDOW <- 128000L

  gen_ai  <- isTRUE(as.logical(param$generate_ai_summary))
  per_sec <- isTRUE(as.logical(param$per_section_ai_summaries))

  ## Redaction helpers -- applied to any text that may be surfaced to the user
  ai_host_only   <- sub("^https?://([^/]+).*", "\\1", AI_ENDPOINT)   # fgcz-c-056:8081
  ai_scheme_host <- sub("^(https?://[^/]+).*", "\\1", AI_ENDPOINT)   # http://fgcz-c-056:8081
  redact_pairs <- list(
    c(AI_ENDPOINT,    "[REDACTED-LLM-ENDPOINT]"),
    c(ai_scheme_host, "[REDACTED-LLM-HOST]"),
    c(ai_host_only,   "[REDACTED-LLM-HOST]")
  )
  redact_string <- function(s) {
    if (is.null(s) || length(s) == 0L) return(s)
    for (p in redact_pairs) s <- gsub(p[1], p[2], s, fixed = TRUE)
    s
  }

  ## ==== Run multiqc ====
  ## We pass AI settings as proper CLI flags rather than via -c YAML, because:
  ##   * -c <(echo ...) process substitution → /dev/fd/N is a FIFO → MultiQC's
  ##     Path.is_file() check in load_config_file rejects it, so the YAML is
  ##     silently ignored and the seqera default from config_defaults.yaml wins.
  ##   * -c <real .yaml tempfile> loaded but didn't reliably propagate to
  ##     config.ai_provider in our setup either.
  ## CLI flags are applied unconditionally in update_config.py:228-235, so they
  ## are the only path that is guaranteed to land.
  ##
  ## TMPDIR is pointed at the current working dir so MultiQC's tempdir
  ## (tempfile.mkdtemp + shutil.rmtree at the end) lives somewhere we own and
  ## can rmdir. Default /tmp on SLURM compute nodes denied rmdir with EACCES,
  ## triggering MultiQC's quadratic-backoff retry loop (0+1+4+9+...+81 ≈ 285 s
  ## wasted at the end of each run).
  multiqc_dir <- file.path("..", basename(output$getColumn("MultiQC")))
  multiqc_tmpdir <- file.path(getwd(), ".multiqc_tmp")
  dir.create(multiqc_tmpdir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TMPDIR = multiqc_tmpdir)
  ## Use double-quotes (not shQuote, which emits single quotes) because the full
  ## command will be piped to sed for redaction, and ezSystem's pipe wrapper
  ## (bash -c '...') refuses any cmd that contains single quotes.
  multiqcCmd <- paste0('multiqc --outdir "', multiqc_dir, '" .')
  if (gen_ai) {
    multiqcCmd <- paste0(
      'OPENAI_API_KEY="dummy" TMPDIR="', multiqc_tmpdir,
      '" multiqc --outdir "', multiqc_dir, '" .',
      ' --ai-summary-full',
      ' --ai-provider ',              AI_PROVIDER,
      ' --ai-model "',                AI_MODEL, '"',
      ' --ai-custom-endpoint "',      AI_ENDPOINT, '"',
      ' --ai-custom-context-window ', AI_CONTEXT_WINDOW
    )
    redact_sed <- tempfile(pattern = "redact_", fileext = ".sed", tmpdir = multiqc_tmpdir)
    writeLines(vapply(redact_pairs, function(p) sprintf("s|%s|%s|g", p[1], p[2]), character(1)),
               redact_sed)
    multiqcCmdShell <- sprintf("%s 2>&1 | sed -f %s", multiqcCmd, redact_sed)
  } else {
    multiqcCmdShell <- multiqcCmd
  }

  t_mqc_start <- Sys.time()
  message(sprintf("[MultiQC] STARTED at %s (generate_ai_summary=%s, per_section_ai_summaries=%s, model=%s)",
                  format(t_mqc_start, "%Y-%m-%d %H:%M:%S"),
                  gen_ai, per_sec, AI_MODEL))
  ## echo = FALSE so ezSystem doesn't ezLog the full command (would leak the endpoint URL).
  ezSystem(multiqcCmdShell, echo = FALSE)
  t_mqc_end <- Sys.time()
  mqc_dur_s <- as.numeric(difftime(t_mqc_end, t_mqc_start, units = "secs"))
  message(sprintf("[MultiQC] FINISHED at %s -- duration: %.1f s",
                  format(t_mqc_end, "%Y-%m-%d %H:%M:%S"), mqc_dur_s))

  ## ==== Post-scrub: strip endpoint from any MultiQC-written files ====
  if (gen_ai) {
    scrub_files <- c(file.path(multiqc_dir, "multiqc_report.html"),
                     file.path(multiqc_dir, "multiqc_data", "multiqc.log"),
                     file.path(multiqc_dir, "multiqc_data", "multiqc_data.json"))
    for (f in scrub_files) {
      if (!file.exists(f)) next
      txt <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
      txt <- redact_string(txt)
      writeLines(txt, f, useBytes = TRUE)
    }
  }

  ## ==== Per-section AI summaries ====
  if (per_sec) {
    multiqc_html_path <- file.path(multiqc_dir, "multiqc_report.html")
    if (!file.exists(multiqc_html_path)) {
      message("[AI-Sections] multiqc_report.html missing; skipping per-section AI")
    } else {
      ## Markdown helpers (CommonMark-ish with MultiQC severity directives)
      style_directives <- function(txt) {
        if (is.null(txt) || !nzchar(txt)) return("")
        for (sev in c("red", "orange", "yellow", "green")) {
          txt <- gsub(sprintf(":span\\[([^]]*)\\]\\{\\.text-%s\\}", sev),
                      sprintf('<span class="text-%s">\\1</span>', sev), txt, perl = TRUE)
          txt <- gsub(sprintf(":sample\\[([^]]*)\\]\\{\\.text-%s\\}", sev),
                      sprintf('<span class="text-%s" style="font-weight:600;font-style:italic;">\\1</span>', sev),
                      txt, perl = TRUE)
        }
        txt <- gsub(":span\\[([^]]*)\\]\\{[^}]*\\}", "\\1", txt, perl = TRUE)
        txt <- gsub(":sample\\[([^]]*)\\]\\{[^}]*\\}",
                    '<span style="font-weight:600;font-style:italic;">\\1</span>', txt, perl = TRUE)
        txt
      }
      md_to_html <- function(txt) {
        txt <- style_directives(txt)
        lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
        out <- character(0)
        list_depth <- 0L
        for (l in lines) {
          indent_match <- regmatches(l, regexpr("^\\s*", l))
          indent <- nchar(indent_match)
          trimmed <- sub("^\\s+", "", l)
          if (grepl("^[-*] ", trimmed)) {
            target_depth <- (indent %/% 4L) + 1L
            while (list_depth < target_depth) { out <- c(out, "<ul>"); list_depth <- list_depth + 1L }
            while (list_depth > target_depth) { out <- c(out, "</ul>"); list_depth <- list_depth - 1L }
            item <- sub("^[-*] ", "", trimmed)
            out <- c(out, paste0("<li>", item, "</li>"))
          } else if (nzchar(trimmed)) {
            while (list_depth > 0L) { out <- c(out, "</ul>"); list_depth <- list_depth - 1L }
            out <- c(out, paste0("<p>", trimmed, "</p>"))
          }
        }
        while (list_depth > 0L) { out <- c(out, "</ul>"); list_depth <- list_depth - 1L }
        paste(out, collapse = "\n")
      }

      ## Data file -> MultiQC section anchor ID
      data_to_section <- function(fname) {
        base <- sub("\\.txt$", "", fname)
        base <- sub("_plot(_Counts|_Percentages)?$", "", base)
        base <- sub("_table$", "", base)
        base <- sub("[-_]heatmap$", "", base)
        base <- gsub("-", "_", base)
        if (base == "fastqc_status_check") base <- "fastqc_status_checks"
        if (base == "multiqc_general_stats") base <- "general_stats_table"
        base
      }

      ## Build candidate per-section data files
      data_dir <- file.path(multiqc_dir, "multiqc_data")
      section_files <- list.files(data_dir,
                                  pattern = "_plot.*\\.txt$|.*-heatmap\\.txt$|.*_table\\.txt$|^multiqc_general_stats\\.txt$",
                                  full.names = TRUE)
      section_files <- section_files[!grepl("/multiqc_(fastqc|sources|software_versions|citations)\\.txt$",
                                            section_files)]
      ## De-dup by mapped section id (a section may have multiple data files like
      ## _plot_Counts and _plot_Percentages; we only call the LLM once per section)
      sec_ids <- vapply(basename(section_files), data_to_section, character(1))
      keep <- !duplicated(sec_ids)
      section_files <- section_files[keep]
      sec_ids <- sec_ids[keep]
      N <- length(section_files)

      ## Read HTML once; we'll do all substitutions in-memory then write back
      html_text <- paste(readLines(multiqc_html_path, warn = FALSE, encoding = "UTF-8"),
                          collapse = "\n")

      sys_prompt <- paste(
        "You are an expert in bioinformatics, sequencing, and QC reports.",
        "You are given the data for a single section of a MultiQC FastQC report (a plot, table, or heatmap).",
        "Produce 1-2 concise bullet points summarising the key observations and any QC issues you see.",
        "If nothing concerning, say so in one bullet.",
        "Use markdown. Highlight severity with CommonMark directives like :span[39.2%]{.text-red}, .text-orange, .text-yellow, .text-green.",
        "Highlight sample names with :sample[name]{.text-red} etc.",
        "Use 4 spaces to indent nested lists. Do not add headers.",
        sep = "\n"
      )

      t_sec_total_start <- Sys.time()
      message(sprintf("[AI-Sections] STARTED at %s (N=%d sections, model=%s)",
                      format(t_sec_total_start, "%Y-%m-%d %H:%M:%S"), N, AI_MODEL))

      for (i in seq_along(section_files)) {
        sec_path <- section_files[i]
        sec_id   <- sec_ids[i]
        sec_data <- paste(readLines(sec_path, warn = FALSE), collapse = "\n")
        max_chars <- AI_CONTEXT_WINDOW * 3L
        if (nchar(sec_data) > max_chars) {
          sec_data <- paste0(substr(sec_data, 1L, max_chars), "\n[... truncated ...]")
        }
        user_msg <- sprintf("Section: %s\n\nData:\n%s", sec_id, sec_data)
        payload <- list(model = AI_MODEL,
                        messages = list(list(role = "system", content = sys_prompt),
                                        list(role = "user",   content = user_msg)))

        t_one_start <- Sys.time()
        resp <- tryCatch(
          httr2::request(AI_ENDPOINT) |>
            httr2::req_method("POST") |>
            httr2::req_headers(Authorization = "Bearer dummy") |>
            httr2::req_body_json(payload) |>
            httr2::req_timeout(600) |>
            httr2::req_error(is_error = function(r) FALSE) |>
            httr2::req_perform(),
          error = function(e) {
            message(sprintf("[AI-Sections] HTTP error on %s: %s",
                            sec_id, redact_string(conditionMessage(e))))
            NULL
          }
        )
        t_one_end <- Sys.time()
        one_dur <- as.numeric(difftime(t_one_end, t_one_start, units = "secs"))

        raw_body <- if (!is.null(resp)) httr2::resp_body_string(resp) else ""
        prompt_tokens     <- 0L
        completion_tokens <- 0L
        total_tokens      <- 0L
        response_text     <- ""
        if (nzchar(raw_body)) {
          parsed <- tryCatch(jsonlite::fromJSON(raw_body, simplifyVector = FALSE),
                              error = function(e) NULL)
          if (!is.null(parsed)) {
            response_text     <- tryCatch(parsed$choices[[1]]$message$content, error = function(e) "")
            prompt_tokens     <- as.integer(parsed$usage$prompt_tokens %||% 0L)
            completion_tokens <- as.integer(parsed$usage$completion_tokens %||% 0L)
            total_tokens      <- as.integer(parsed$usage$total_tokens %||%
                                              (prompt_tokens + completion_tokens))
          }
        }
        if (is.null(response_text) || !nzchar(response_text)) {
          response_text <- "_(no AI response)_"
        }
        response_text <- redact_string(response_text)

        tok_per_s <- if (one_dur > 0 && completion_tokens > 0) completion_tokens / one_dur else NA_real_
        footer_html <- sprintf(
          paste0('<div class="text-muted" style="font-size:0.85em;margin-top:8px;',
                 'border-top:1px solid #eee;padding-top:4px;">',
                 'Model: %s &middot; %.1f s &middot; %d tokens (%d in, %d out) &middot; %s tok/s',
                 '</div>'),
          AI_MODEL, one_dur, total_tokens, prompt_tokens, completion_tokens,
          if (is.na(tok_per_s)) "n/a" else sprintf("%.1f", tok_per_s)
        )
        summary_html <- paste0(md_to_html(response_text), "\n", footer_html)

        ## Inject: try MultiQC's pre-created empty div first, else inject into the section wrapper
        empty_div_re <- sprintf(
          '<div class="ai-summary-response" id="%s_ai_summary_response"[^>]*></div>',
          sec_id
        )
        empty_div_replacement <- sprintf(
          '<div class="ai-summary-response" id="%s_ai_summary_response" style="margin-bottom: -5px;">%s</div>',
          sec_id, summary_html
        )
        if (grepl(empty_div_re, html_text, perl = TRUE)) {
          html_text <- sub(empty_div_re, empty_div_replacement, html_text, perl = TRUE)
        } else {
          wrapper_re <- sprintf('(<div [^>]*id="mqc-section-wrapper-%s"[^>]*>)', sec_id)
          wrapper_replacement <- sprintf(
            '\\1\n<div class="ai-summary-response" style="margin:10px 0;padding:8px 12px;border-left:3px solid #4a90d9;background:#f5f9fd;">%s</div>',
            summary_html
          )
          if (grepl(wrapper_re, html_text, perl = TRUE)) {
            html_text <- sub(wrapper_re, wrapper_replacement, html_text, perl = TRUE)
          } else {
            message(sprintf("[AI-Sections] no injection point for %s; skipping", sec_id))
          }
        }

        message(sprintf("[AI-Sections] %d/%d %s: %.1f s, %d tokens, %s tok/s",
                        i, N, sec_id, one_dur, total_tokens,
                        if (is.na(tok_per_s)) "n/a" else sprintf("%.1f", tok_per_s)))
      }

      writeLines(html_text, multiqc_html_path, useBytes = TRUE)

      t_sec_total_end <- Sys.time()
      sec_total_dur_s <- as.numeric(difftime(t_sec_total_end, t_sec_total_start, units = "secs"))
      message(sprintf("[AI-Sections] FINISHED at %s -- total: %.1f s (mean %.1f s/section)",
                      format(t_sec_total_end, "%Y-%m-%d %H:%M:%S"),
                      sec_total_dur_s, sec_total_dur_s / max(1L, N)))
    }
  }

  unlink(c("fastqc.out", "fastqc.err"))
  unlink(multiqc_tmpdir, recursive = TRUE, force = TRUE)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }
EzAppFastqc <-
  setRefClass(
    "EzAppFastqc",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFastQC
        name <<- "EzAppFastqc"
        appDefaults <<- rbind(
          perLibrary = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Run FastQC per library or per cell for single cell experiment"
          ),
          showNativeReports = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Keep the original fastqc report"
          )
        )
      }
    )
  )

plotReadCountToLibConc <- function(dataset, colname) {
  if (colname %in% colnames(dataset) && nrow(dataset) > 1) {
    if (!all(dataset[[colname]] == 0) && !any(is.na(dataset[[colname]]))) {
      ## LibCon column can all be 0 or NA. Then don't plot.
      dataset <- dataset[order(dataset$"Read Count", decreasing = T), ]
      dataset$"Read Count" <- dataset$"Read Count" / 10^6
      corResult <- cor.test(
        dataset$"Read Count",
        dataset[[colname]],
        method = "spearman"
      )
      regressionResult <- lm(dataset[[colname]] ~ dataset$"Read Count")
      label <- sub(" \\[.*", "", colname)

      ## plotly
      require(plotly)
      # a function to calculate your abline
      xmin <- min(dataset$"Read Count") - 5
      xmax <- max(dataset$"Read Count") + 5
      intercept <- regressionResult$coefficients[1]
      slope <- regressionResult$coefficients[2]
      p_abline <- function(x, a, b) {
        y <- a * x + b
        return(y)
      }

      p <- plot_ly(
        x = dataset$"Read Count",
        y = dataset[[colname]],
        text = rownames(dataset)
      ) %>%
        add_markers() %>%
        add_text(textposition = "top right") %>%
        plotly::layout(showlegend = FALSE)
      a <- list(
        x = max(dataset$"Read Count"),
        y = max(dataset[[colname]]),
        text = paste0("r=", round(corResult$estimate, 2)),
        xref = "x",
        yref = "y",
        showarrow = FALSE
      )
      p <- p %>%
        plotly::layout(
          shapes = list(
            type = "line",
            line = list(dash = "dash"),
            x0 = xmin,
            x1 = xmax,
            y0 = p_abline(xmin, slope, intercept),
            y1 = p_abline(xmax, slope, intercept)
          ),
          annotations = a,
          title = label,
          yaxis = list(title = label),
          xaxis = list(title = "Counts [Mio]")
        )
      return(p)
    }
  }
}

plotQualityMatrixAsHeatmapGG2 <- function(
  qualMatrixList,
  isR2 = FALSE,
  xScale = 1,
  yScale = 1
) {
  colorsGray <- gray((30:256) / 256)
  minPercent <- 0
  maxPercent <- sqrt(40)
  minDiff <- -5
  maxDiff <- 5
  plotList <- list() # to store all the ggplot objects

  ## test if R2 exists
  index <- list("R1" = which(!isR2))

  if (any(isR2)) {
    stopifnot(sum(isR2) == sum(!isR2))
    index[["R2"]] <- which(isR2)
  }
  for (nm in names(index)) {
    plotList[[nm]] <- list()
    idx <- index[[nm]]
    ## Plot the color key for the average quality heatmap R1_1
    # colorKeyFile = paste0("averageReadsQuality-Key_", nm, ".png")
    by.label <- 1
    at <- seq(from = minPercent, to = maxPercent, by = by.label)
    p <- ezColorLegendGG2(
      colorRange = c(minPercent, maxPercent),
      colors = colorsGray,
      vertical = FALSE,
      by.label = by.label,
      at = at,
      labels = as.character(at^2)
    )
    plotList[[nm]][["Avg Qual Colors"]] <- p

    # result = ezMatrix(0, dim=dim(qualMatrixList[[idx[1]]]))
    result <- ezMatrix(0, dim = apply(sapply(qualMatrixList[idx], dim), 1, max))
    resultCount <- result
    for (i in idx) {
      qm <- qualMatrixList[[i]]
      result[1:nrow(qm), 1:ncol(qm)] <- result[1:nrow(qm), 1:ncol(qm)] + qm
      resultCount[1:nrow(qm), 1:ncol(qm)] <- resultCount[
        1:nrow(qm),
        1:ncol(qm)
      ] +
        1
    }
    result <- result / resultCount
    ## The hard way to deal with NaN in result
    result <- sweep(
      result,
      MARGIN = 2,
      STATS = colSums(result, na.rm = TRUE),
      FUN = "/"
    )
    avgQual <- signif(result * 100, digits = 3)
    p <- plotQualityHeatmapGG2(
      result = sqrt(avgQual),
      colorRange = c(minPercent, maxPercent),
      colors = colorsGray,
      main = paste("averageReadsQuality", nm, sep = "_"),
      xScale = xScale,
      yScale = yScale
    )
    plotList[[nm]][["Average"]] <- p

    ## plot the difference quality heatmap for R1_1
    at <- seq(from = minDiff, to = maxDiff, by = by.label)
    p <- ezColorLegendGG2(
      colorRange = c(minDiff, maxDiff),
      colors = getBlueRedScale(),
      vertical = FALSE,
      by.label = by.label,
      at = at,
      labels = as.character(at)
    )
    plotList[[nm]][["Diff Qual Colors"]] <- p

    for (sampleName in names(qualMatrixList[idx])) {
      qm <- qualMatrixList[[sampleName]]
      qm <- sweep(qm, MARGIN = 2, STATS = colSums(qm, na.rm = TRUE), FUN = "/")
      diffResult <- signif(qm * 100, digits = 3) -
        avgQual[1:nrow(qm), 1:ncol(qm)]
      p <- plotQualityHeatmapGG2(
        diffResult,
        colorRange = c(minDiff, maxDiff),
        colors = getBlueRedScale(),
        main = paste("diffReadsQuality", sampleName, sep = "_"),
        xScale = xScale,
        yScale = yScale
      )
      plotList[[nm]][[sampleName]] <- p
    }
  }
  return(plotList)
}

plotQualityHeatmapGG2 <- function(
  result,
  name = NULL,
  colorRange = c(0, sqrt(40)),
  colors = gray((1:256) / 256),
  main = NULL,
  xScale = 1,
  yScale = 1
) {
  require(reshape2)
  ## some ugly controls of labels
  labCol <- seq(0, ncol(result), by = 10)
  labCol[1] <- 1
  labRow <- seq(0, nrow(result) - 1, by = 5)

  result[result > colorRange[2]] <- colorRange[2]
  result[result < colorRange[1]] <- colorRange[1]
  toPlot <- melt(result)
  p <- ggplot(data = toPlot, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colours = colors, limits = colorRange) +
    theme_bw() +
    scale_y_continuous(
      breaks = seq(1, nrow(result), by = 5),
      labels = labRow,
      expand = c(0, 0)
    ) +
    scale_x_continuous(breaks = labCol, labels = labCol, expand = c(0, 0)) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5)
    ) +
    xlab("Read Position") +
    ylab("Read Quality") +
    ggtitle(main)
  return(p)
}

### Sample up to 300k reads from a fastq file or bam file.
### Calculate the quality matrix
getQualityMatrix <- function(fn) {
  ## This implementation is faster than the FastqStreamer.
  require(ShortRead)
  nReads <- 3e5

  if (grepl("\\.bam$", fn)) {
    cmd <- paste("samtools flagstat", fn)
    cmdOutput <- ezSystem(cmd, intern = TRUE)
    nrTotal <- eval(parse(text = sub(" in total.*", "", cmdOutput[1])))

    tempBamFn <- paste(Sys.getpid(), "temp.bam", sep = "-")
    on.exit(file.remove(tempBamFn), add = TRUE)
    cmd <- paste(
      "samtools view -s",
      nReads / nrTotal,
      "-b",
      fn,
      ">",
      tempBamFn
    )
    ezSystem(cmd)
    tempFastqFn <- paste(Sys.getpid(), "temp.fastq", sep = "-")
    on.exit(file.remove(tempFastqFn), add = TRUE)
    bam2fastq(
      bamFn = tempBamFn,
      OUTPUT_PER_RG = FALSE,
      fastqFns = tempFastqFn,
      paired = FALSE
    )
    qualMatrix <- as(quality(readFastq(tempFastqFn)), "matrix")
  } else {
    f <- FastqSampler(fn, nReads) ## we sample no more than 300k reads.
    qualMatrix <- as(quality(yield(f)), "matrix")
  }

  gc()
  maxQuality <- max(qualMatrix, na.rm = TRUE)
  qualCountMatrix <- ezMatrix(0, rows = 0:maxQuality, cols = 1:ncol(qualMatrix))
  for (basePos in 1:ncol(qualMatrix)) {
    qualCountByPos <- table(qualMatrix[, basePos])
    qualCountMatrix[names(qualCountByPos), basePos] <- qualCountByPos
  }
  return(qualCountMatrix)
}

plateStatistics <- function(
  dataset,
  colname = c(
    "Read Count",
    "LibConc_qPCR [Characteristic]",
    "LibConc_100_800bp [Characteristic]"
  )
) {
  colsExist <- colname %in% colnames(dataset)
  if (any(!colsExist)) {
    warning("No column ", colname[!colsExist], " in dataset!")
    colname <- colname[colsExist]
  }
  colsNumeric <- sapply(dataset[, colname, drop = FALSE], is, "numeric")
  if (any(!colsNumeric)) {
    warning("The column ", colname[!colsNumeric], " is non-numeric.")
    colname <- colname[colsNumeric]
  }
  colsNA <- is.na(colSums(dataset[, colname, drop = FALSE]))
  if (any(colsNA)) {
    ezLog("The column ", colname[colsNA], " has NA!")
    colname <- colname[!colsNA]
  }
  colsZero <- colSums(dataset[, colname, drop = FALSE]) == 0
  if (any(colsZero)) {
    ezLog("The column ", colname[colsZero], " is empty!")
    colname <- colname[!colsZero]
  }
  if (length(colname) == 0L) {
    warning("No suitable columns left in dataset!")
    return(NA)
  }

  if (
    !is.null(dataset$`PlatePosition [Characteristic]`) &&
      !any(is.na(dataset$`PlatePosition [Characteristic]`))
  ) {
    # `PlatePosition [Characteristic]` may not exist or is empty

    plateChar <- dataset$`PlatePosition [Characteristic]`
    ## Plate position should be in the format of *_C4;
    ## * is the plate number
    ## C is the column; 4 is the row
    require(stringr, quietly = TRUE)
    plateNumber <- sub("_[[:alpha:]]\\d+$", "", plateChar)
    datasetByPlate <- split(dataset, plateNumber)
    ans <- list()
    for (i in seq_len(length(datasetByPlate))) {
      platePos <- str_extract(
        datasetByPlate[[i]]$`PlatePosition [Characteristic]`,
        "_[[:alpha:]]\\d+$"
      )
      if (any(is.na(platePos))) {
        warning("The PlatePosition format is not supported!")
        return(NA)
      }
      plateRow <- str_extract(platePos, "[[:alpha:]]")
      plateCol <- as.numeric(str_extract(platePos, "\\d+"))
      ans[[names(datasetByPlate)[i]]] <- list()
      for (oneCol in colname) {
        ## always the entire plate should be shown which is either 8x12 or 16x24 ....
        counts <- datasetByPlate[[i]][[oneCol]]
        countMatrix <- ezMatrix(
          NA,
          rows = LETTERS[1:ifelse(max(plateRow) > "I", 16, 8)],
          cols = seq_len(ifelse(max(as.integer(plateCol)) > 12, 24, 12))
        )
        for (j in seq_len(length(counts))) {
          countMatrix[plateRow[j], plateCol[j]] <- counts[j]
        }
        ans[[names(datasetByPlate)[i]]][[oneCol]] <- countMatrix
      }
    }
    return(ans)
  } else {
    warning("PlatePosition [Characteristic] information is not available!")
    return(NA)
  }
}

heatmapPlate <- function(
  x,
  title = "unnamed",
  center = TRUE,
  log10 = TRUE,
  ...
) {
  require(plotly)
  ## do not plot if there are only NA or zeros
  if (all(x %in% c(NA, 0))) {
    return(NULL)
  }
  if (isTRUE(log10)) {
    ## shift zeros a bit
    isZero <- x == 0
    isZero[is.na(isZero)] <- FALSE
    x[isZero] <- min(0.25 * x[x > 0], na.rm = TRUE)

    x <- log10(x)
    if (isTRUE(center)) {
      medianX <- median(x, na.rm = TRUE)
      p <- plot_ly(
        z = x,
        x = colnames(x),
        y = rownames(x),
        type = "heatmap",
        zmin = medianX - log10(2),
        zmax = medianX + log10(2),
        hoverinfo = "text",
        text = matrix(
          paste0(
            "10^",
            format(x, digits = 3),
            "=",
            10^x
          ),
          ncol = ncol(x)
        ),
        ...
      )
    } else {
      p <- plot_ly(
        z = x,
        x = colnames(x),
        y = rownames(x),
        type = "heatmap",
        hoverinfo = "text",
        text = matrix(
          paste0(
            "10^",
            format(x, digits = 3),
            "=",
            10^x
          ),
          ncol = ncol(x)
        ),
        ...
      )
    }
  } else {
    if (isTRUE(center)) {
      medianX <- median(x, na.rm = TRUE)
      p <- plot_ly(
        z = x,
        x = colnames(x),
        y = rownames(x),
        type = "heatmap",
        zmin = 0,
        zmax = medianX * 2,
        ...
      )
    } else {
      p <- plot_ly(
        z = x,
        x = colnames(x),
        y = rownames(x),
        type = "heatmap",
        ...
      )
    }
  }

  p <- p %>%
    plotly::layout(
      xaxis = list(autotick = FALSE, dtick = 1),
      yaxis = list(autorange = "reversed"),
      title = title
    )
  p
}
