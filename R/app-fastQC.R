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

  ## ==== AI feature flags ====
  ai_mode <- all(nzchar(c(param$ai_provider, param$ai_model,
                          param$ai_custom_endpoint, param$ai_custom_context_window)))
  bake_sections <- isTRUE(param$bake_section_summaries) && ai_mode
  extra_call <- isTRUE(param$extra_ai_call) && nzchar(param$extra_ai_prompt) && ai_mode

  ## ==== Block 1: run multiqc (optionally with --ai-summary-full) ====
  multiqcCmd <- paste0("multiqc --outdir ../", basename(output$getColumn("MultiQC")), " .")
  if (ai_mode) {
    ai_yaml <- tempfile(pattern = "multiqc_ai_", fileext = ".yaml")
    writeLines(c(
      paste0("ai_provider: ",              param$ai_provider),
      paste0("ai_model: ",                 param$ai_model),
      paste0("ai_custom_endpoint: ",       param$ai_custom_endpoint),
      paste0("ai_custom_context_window: ", param$ai_custom_context_window)
    ), ai_yaml)
    Sys.setenv(OPENAI_API_KEY = "dummy")
    multiqcCmd <- paste0(multiqcCmd, " --ai-summary-full -c ", shQuote(ai_yaml))
  }

  t_mqc_start <- Sys.time()
  message(sprintf("[MultiQC] STARTED at %s (ai_mode=%s) -- cmd: %s",
                  format(t_mqc_start, "%Y-%m-%d %H:%M:%S"),
                  ai_mode, multiqcCmd))
  ezSystem(multiqcCmd)
  t_mqc_end <- Sys.time()
  mqc_dur_s <- as.numeric(difftime(t_mqc_end, t_mqc_start, units = "secs"))
  message(sprintf("[MultiQC] FINISHED at %s -- duration: %.1f s",
                  format(t_mqc_end, "%Y-%m-%d %H:%M:%S"), mqc_dur_s))

  ## ==== Block 2: extract MultiQC's global AI summary from the rendered HTML ====
  multiqc_dir <- file.path("..", basename(output$getColumn("MultiQC")))
  multiqc_html_path <- file.path(multiqc_dir, "multiqc_report.html")
  llms_full_path <- file.path(multiqc_dir, "multiqc_data", "llms-full.txt")
  mqc_global_summary_html <- ""
  mqc_detailed_html <- ""
  prompt_text <- ""
  if (ai_mode && file.exists(multiqc_html_path)) {
    html_doc <- tryCatch(xml2::read_html(multiqc_html_path), error = function(e) NULL)
    if (!is.null(html_doc)) {
      pull_div <- function(id) {
        node <- xml2::xml_find_first(html_doc, sprintf('//div[@id="%s"]', id))
        if (inherits(node, "xml_missing")) "" else as.character(node)
      }
      mqc_global_summary_html <- pull_div("global_ai_summary_response")
      mqc_detailed_html       <- pull_div("global_ai_summary_detailed_analysis_response")
    } else {
      mqc_global_summary_html <- "[MultiQC HTML schema changed -- see multiqc_report.html directly]"
    }
  }
  if (ai_mode && file.exists(llms_full_path)) {
    prompt_text <- paste(readLines(llms_full_path, warn = FALSE), collapse = "\n")
  }

  ## ==== Block 3: per-section AI calls (bake_section_summaries) ====
  section_summaries <- list()
  sec_total_dur_s <- NA_real_
  if (bake_sections) {
    section_files <- list.files(file.path(multiqc_dir, "multiqc_data"),
                                pattern = "_plot.*\\.txt$|.*-heatmap\\.txt$|.*_table\\.txt$",
                                full.names = TRUE)
    section_files <- section_files[!grepl("/multiqc_(general_stats|fastqc|sources|software_versions|citations)\\.txt$",
                                          section_files)]
    N <- length(section_files)
    t_sec_total_start <- Sys.time()
    message(sprintf("[AI-Sections] STARTED at %s (N=%d sections, model=%s)",
                    format(t_sec_total_start, "%Y-%m-%d %H:%M:%S"), N, param$ai_model))
    sys_prompt <- paste(
      "You are a bioinformatics expert reviewing a single section of a MultiQC FastQC report.",
      "In 2-3 concise bullet points, summarize the key observations from the data below and flag any concerning patterns.",
      "Use markdown formatting. Highlight severity with directives like :span[value]{.text-red}, .text-orange, .text-green.",
      "Highlight sample names with :sample[name]{.text-red} etc. Do not add headers.",
      sep = "\n"
    )
    for (i in seq_along(section_files)) {
      sec_path <- section_files[i]
      sec_name <- sub("\\.txt$", "", basename(sec_path))
      sec_data <- paste(readLines(sec_path, warn = FALSE), collapse = "\n")
      max_chars <- as.integer(param$ai_custom_context_window) * 3L
      if (!is.na(max_chars) && nchar(sec_data) > max_chars) {
        sec_data <- paste0(substr(sec_data, 1, max_chars), "\n[... truncated ...]")
      }
      user_msg <- sprintf("Section: %s\n\nData:\n%s", sec_name, sec_data)
      payload <- list(model = param$ai_model,
                      messages = list(list(role = "system", content = sys_prompt),
                                      list(role = "user",   content = user_msg)))
      t_one_start <- Sys.time()
      resp <- tryCatch(
        httr2::request(param$ai_custom_endpoint) |>
          httr2::req_method("POST") |>
          httr2::req_headers(Authorization = "Bearer dummy") |>
          httr2::req_body_json(payload) |>
          httr2::req_timeout(600) |>
          httr2::req_error(is_error = function(r) FALSE) |>
          httr2::req_perform(),
        error = function(e) { message("[AI-Sections] HTTP error on ", sec_name, ": ", conditionMessage(e)); NULL }
      )
      t_one_end <- Sys.time()
      one_dur <- as.numeric(difftime(t_one_end, t_one_start, units = "secs"))
      raw_body <- if (!is.null(resp)) httr2::resp_body_string(resp) else "[no response]"
      response_text <- tryCatch({
        parsed <- jsonlite::fromJSON(raw_body, simplifyVector = FALSE)
        parsed$choices[[1]]$message$content
      }, error = function(e) raw_body)
      if (is.null(response_text) || !nzchar(response_text)) response_text <- raw_body
      section_summaries[[sec_name]] <- list(response_text = response_text,
                                            raw_body = raw_body,
                                            dur_s = one_dur)
      message(sprintf("[AI-Sections] %d/%d %s: %.1f s", i, N, sec_name, one_dur))
    }
    t_sec_total_end <- Sys.time()
    sec_total_dur_s <- as.numeric(difftime(t_sec_total_end, t_sec_total_start, units = "secs"))
    message(sprintf("[AI-Sections] FINISHED at %s -- total: %.1f s (mean %.1f s/section)",
                    format(t_sec_total_end, "%Y-%m-%d %H:%M:%S"),
                    sec_total_dur_s, sec_total_dur_s / max(1L, N)))
  }

  ## ==== Block 4: FGCZ follow-up AI call (extra_ai_call + prompt) ====
  extra_response_text <- NULL
  extra_raw_body <- NULL
  extra_t_start <- NULL
  extra_t_end <- NULL
  extra_dur_s <- NA_real_
  if (extra_call && nzchar(prompt_text)) {
    composite_prompt <- paste0(prompt_text,
                               "\n\n----------------------\n",
                               "Additional instructions from the FGCZ user:\n",
                               param$extra_ai_prompt, "\n")
    payload <- list(model = param$ai_model,
                    messages = list(list(role = "user", content = composite_prompt)))
    extra_t_start <- Sys.time()
    message(sprintf("[AI] Extra LLM call STARTED at %s (model=%s, prompt_chars=%d)",
                    format(extra_t_start, "%Y-%m-%d %H:%M:%S"),
                    param$ai_model, nchar(composite_prompt)))
    resp <- tryCatch(
      httr2::request(param$ai_custom_endpoint) |>
        httr2::req_method("POST") |>
        httr2::req_headers(Authorization = "Bearer dummy") |>
        httr2::req_body_json(payload) |>
        httr2::req_timeout(600) |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform(),
      error = function(e) { message("[AI] HTTP error: ", conditionMessage(e)); NULL }
    )
    extra_t_end <- Sys.time()
    extra_dur_s <- as.numeric(difftime(extra_t_end, extra_t_start, units = "secs"))
    message(sprintf("[AI] Extra LLM call FINISHED at %s -- duration: %.1f s",
                    format(extra_t_end, "%Y-%m-%d %H:%M:%S"), extra_dur_s))
    extra_raw_body <- if (!is.null(resp)) httr2::resp_body_string(resp) else "[no response -- see SLURM log]"
    extra_response_text <- tryCatch({
      parsed <- jsonlite::fromJSON(extra_raw_body, simplifyVector = FALSE)
      parsed$choices[[1]]$message$content
    }, error = function(e) extra_raw_body)
    if (is.null(extra_response_text) || !nzchar(extra_response_text)) extra_response_text <- extra_raw_body
  }

  ## ==== Block 5: render the consolidated AI Interpretation report (whenever ai_mode) ====
  if (ai_mode) {
    saveRDS(list(
      flags = list(ai_mode = ai_mode, bake_sections = bake_sections, extra_call = extra_call),
      endpoint = param$ai_custom_endpoint, model = param$ai_model,
      t_mqc_start = t_mqc_start, t_mqc_end = t_mqc_end, mqc_dur_s = mqc_dur_s,
      mqc_global_summary_html = mqc_global_summary_html,
      mqc_detailed_html = mqc_detailed_html,
      prompt = prompt_text,
      section_summaries = section_summaries,
      sec_total_dur_s = sec_total_dur_s,
      extra_prompt = if (extra_call) param$extra_ai_prompt else "",
      extra_response_text = extra_response_text,
      extra_raw_body = extra_raw_body,
      extra_t_start = extra_t_start, extra_t_end = extra_t_end, extra_dur_s = extra_dur_s
    ), file.path(multiqc_dir, "ai_interpretation.rds"))

    file.copy(system.file("templates/FastQC_AI_Interpretation.Rmd", package = "ezRun"),
              "FastQC_AI_Interpretation.Rmd", overwrite = TRUE)
    rmarkdown::render(
      input       = "FastQC_AI_Interpretation.Rmd",
      envir       = new.env(),
      output_dir  = multiqc_dir,
      output_file = "ai_interpretation.html",
      quiet       = TRUE
    )
  }

  unlink(c("fastqc.out", "fastqc.err"))

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
