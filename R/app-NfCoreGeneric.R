EzAppNfCoreGeneric <- setRefClass(
  "EzAppNfCoreGeneric",
  contains = "EzApp",
  methods = list(
    initialize = function() {
      "Initializes the application"
      callSuper()
    },
    run = function(input, output, param) {
      pipeline <- param[["nfcorePipeline"]]
      version <- param[["pipelineVersion"]]
      dataRoot <- param[["dataRoot"]]
      
      # Output to current directory (scratch)
      # Results will be copied to gstore by the job script footer
      outdir <- "_result"
      
      # Ensure output directory exists
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      # Build sample sheet in current directory (scratch)
      sampleSheetPath <- "samplesheet.csv"
      writeNfCoreSampleSheet(input, sampleSheetPath, param)
      
      message("Sample sheet written to: ", sampleSheetPath)
      message("Sample sheet contents:")
      print(read.csv(sampleSheetPath))
      
      # Get resource limits from param (set by SUSHI, auto-detected from pipeline schema)
      maxCpus <- ifelse(is.null(param[["cores"]]), 16, as.integer(param[["cores"]]))
      maxMemory <- ifelse(is.null(param[["ram"]]), 128, as.integer(param[["ram"]]))
      
      customConfig <- "custom_resources.config"
      configContent <- paste0(
        "// Executor: total available resources for this job\n",
        "executor {\n",
        "  name = 'local'\n",
        "  cpus = ", maxCpus, "\n",
        "  memory = '", maxMemory, ".GB'\n",
        "}\n",
        "\n",
        "// Cap all process labels to fit within executor limits\n",
        "process {\n",
        "  resourceLimits = [\n",
        "    cpus: ", maxCpus, ",\n",
        "    memory: '", maxMemory, ".GB'\n",
        "  ]\n",
        "}\n"
      )
      writeLines(configContent, customConfig)
      message("Custom config written: ", customConfig)
      message(configContent)
      
      # Load Nextflow module via lmod (FGCZ environment)
      # Source the lmod R init script to get the module() function
      source("/usr/share/lmod/lmod/init/R")

      # Get module version from param (default to 25.10)
      nfModuleVersion <- param[["nextflowModuleVersion"]]
      if (is.null(nfModuleVersion) || nfModuleVersion == "") {
        nfModuleVersion <- "25.10"
      }

      # Load the nextflow module
      # Note: lmod's module() uses match.call() which captures variable names literally,
      # so we must use do.call() to pass evaluated values
      moduleName <- paste0("nextflow/", nfModuleVersion)
      message("Loading module: ", moduleName)
      do.call(module, list("load", moduleName))

      # Construct nextflow command
      # Output to local directory, will be copied to gstore by footer
      cmd <- paste(
        "nextflow run",
        paste0("nf-core/", pipeline),
        "-r", version,
        "--input", sampleSheetPath,
        "--outdir", outdir,
        "-profile singularity",
        "-c", customConfig
      )
      
      # Resolve reference genome from refBuild parameter
      refBuild <- param[["refBuild"]]
      if (!is.null(refBuild) && refBuild != "" && refBuild != "select") {
        genomesRoots <- strsplit(GENOMES_ROOT, ":")[[1]]
        refParts <- strsplit(refBuild, "/")[[1]]
        genomeBase <- paste(refParts[1:min(3, length(refParts))], collapse = "/")

        refRoot <- NULL
        for (root in genomesRoots) {
          if (dir.exists(file.path(root, genomeBase))) {
            refRoot <- root
            break
          }
        }

        if (!is.null(refRoot)) {
          fastaPath <- file.path(refRoot, genomeBase, "Sequence", "WholeGenomeFasta", "genome.fa")
          if (file.exists(fastaPath)) {
            cmd <- paste(cmd, "--fasta", fastaPath)
            message("Using FASTA: ", fastaPath)
          }

          gtfPath <- NULL
          if (length(refParts) > 3) {
            gtfPath <- file.path(refRoot, refBuild, "Genes", "genes.gtf")
          }
          if (is.null(gtfPath) || !file.exists(gtfPath)) {
            gtfPath <- file.path(refRoot, genomeBase, "Annotation", "Genes", "genes.gtf")
          }
          if (file.exists(gtfPath)) {
            cmd <- paste(cmd, "--gtf", gtfPath)
            message("Using GTF: ", gtfPath)
          }
        } else {
          message("WARNING: Reference genome not found for refBuild: ", refBuild)
        }
      }

      # Pass through pipeline-specific parameters to nextflow
      # These are SUSHI-internal params that should NOT be passed to nextflow
      sushiParams <- c(
        "cores", "ram", "scratch", "partition", "process_mode", "samples",
        "nfcorePipeline", "pipelineVersion", "inputType", "samplesheetMapping",
        "strandedness", "nextflowModuleVersion", "sushi_app",
        "refBuild", "dataRoot", "resultDir", "isLastJob",
        "name", "mail", "cmdOptions", "schemaDefaults"
      )
      # fasta and gtf are already resolved from refBuild above
      handledParams <- c("fasta", "gtf")
      skipParams <- c(sushiParams, handledParams)

      # Pipeline schema defaults ({param: default}) so we can skip params still at
      # their default: re-passing them is redundant and can break schema validation
      # (e.g. a string default like cf_ploidy="2" is re-cast to integer on the CLI).
      schemaDefaults <- tryCatch({
        sd <- param[["schemaDefaults"]]
        if (!is.null(sd) && length(sd) == 1 && !is.na(sd) && nzchar(sd)) {
          jsonlite::fromJSON(sd, simplifyVector = TRUE)
        } else list()
      }, error = function(e) list())

      for (key in names(param)) {
        if (key %in% skipParams) next
        value <- param[[key]]
        if (is.null(value) || is.na(value) || value == "") next
        # Skip params still at the pipeline schema default (case-insensitive).
        if (!is.null(schemaDefaults[[key]]) &&
            tolower(as.character(schemaDefaults[[key]])) == tolower(as.character(value))) next
        if (is.logical(value)) value <- tolower(as.character(value))
        cmd <- paste(cmd, paste0("--", key), shQuote(value))
      }

      message("Running command: ", cmd)

      ezSystem(cmd)
    },
    writeNfCoreSampleSheet = function(input, filePath, param) {
      # Handle input
      if (is.character(input) && file.exists(input)) {
        # Dataset mode: input is path to TSV
        ds <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
      } else {
        # Sample mode: input is list
        ds <- as.data.frame(input, stringsAsFactors=FALSE)
      }

      dataRoot <- param[["dataRoot"]]

      # Helper to resolve file paths
      resolve_path <- function(path) {
        if (is.na(path) || path == "") return(NA)
        if (startsWith(path, "/")) return(path)
        return(file.path(dataRoot, path))
      }

      # Resolve a mapped SUSHI source column name to an actual dataset column:
      # exact match first, then a "^prefix" fallback (e.g. "Read1" -> "Read1 [File]").
      matchSourceCol <- function(src) {
        if (src %in% names(ds)) return(src)
        pat <- paste0("^", gsub("([][(){}.^$*+?|\\\\])", "\\\\\\1", src))
        hit <- grep(pat, names(ds), value = TRUE)
        if (length(hit) >= 1) return(hit[1])
        return(NA_character_)
      }
      # A source column carries a file path (needs resolve_path) iff its name has [File].
      isFileCol <- function(colName) grepl("\\[File\\]", colName)

      # samplesheetMapping arrives as a JSON string {nf_core_col: SUSHI_col}.
      # Parse it; on nil/empty/malformed JSON fall back to the built-in mapping.
      mapping <- NULL
      mj <- param[["samplesheetMapping"]]
      if (!is.null(mj) && length(mj) == 1 && !is.na(mj) && nzchar(mj)) {
        mapping <- tryCatch(
          jsonlite::fromJSON(mj, simplifyVector = TRUE),
          error = function(e) {
            message("WARNING: could not parse samplesheetMapping as JSON (",
                    conditionMessage(e), "); using built-in default mapping.")
            NULL
          }
        )
      }

      strandedness <- param[["strandedness"]]

      if (is.null(mapping) || length(mapping) == 0) {
        # Built-in default mapping (nf-core/demo, rnaseq): sample, fastq_1, fastq_2
        message("samplesheetMapping not usable; using built-in sample/fastq_1/fastq_2 mapping.")
        read1Col <- grep("^Read1", names(ds), value = TRUE)[1]
        read2Col <- grep("^Read2", names(ds), value = TRUE)[1]

        samples <- data.frame(
          sample = ds$Name,
          fastq_1 = sapply(ds[[read1Col]], resolve_path),
          stringsAsFactors = FALSE
        )

        if (!is.null(read2Col) && !is.na(read2Col)) {
          read2Values <- ds[[read2Col]]
          if (!all(is.na(read2Values)) && !all(read2Values == "")) {
            samples$fastq_2 <- sapply(read2Values, resolve_path)
          }
        }

        # Inject strandedness from parameter if set (nf-core/rnaseq requires this column)
        if (!is.null(strandedness) && strandedness != "") {
          samples$strandedness <- rep(strandedness, nrow(samples))
          message("Added strandedness column: ", strandedness)
        }

        write.csv(samples, file=filePath, row.names=FALSE, quote=FALSE)
        return(invisible(filePath))
      }

      # Mapping-driven path: build each nf-core column (in mapping key order) from
      # its mapped SUSHI source column. File columns are path-resolved; plain
      # columns (patient/sex/status/...) are copied verbatim.
      samples <- data.frame(row.names = seq_len(nrow(ds)))
      for (nfCol in names(mapping)) {
        src <- mapping[[nfCol]]
        actual <- matchSourceCol(src)

        if (is.na(actual)) {
          # strandedness is a form param, not always a dataset column: fall back to it.
          if (nfCol == "strandedness" && !is.null(strandedness) && nzchar(strandedness)) {
            samples[[nfCol]] <- rep(strandedness, nrow(ds))
            message("Populated 'strandedness' from form param: ", strandedness)
            next
          }
          # Optional nf-core columns (e.g. sarek sex/status/lane) are often absent
          # from the dataset: skip with a loud message rather than failing the job.
          message("Skipping nf-core column '", nfCol, "' (source '", src,
                  "' not in dataset).")
          next
        }

        vals <- ds[[actual]]
        if (isFileCol(actual)) vals <- sapply(vals, resolve_path)
        if (all(is.na(vals) | vals == "")) {
          message("Skipping empty column '", nfCol, "' (source '", actual, "' all empty).")
          next
        }
        samples[[nfCol]] <- vals
      }

      # Inject strandedness from the form param only if the mapping neither produced
      # nor referenced it (avoid double-adding).
      if (!("strandedness" %in% names(samples)) &&
          !("strandedness" %in% names(mapping)) &&
          !is.null(strandedness) && nzchar(strandedness)) {
        samples$strandedness <- rep(strandedness, nrow(ds))
        message("Added strandedness column: ", strandedness)
      }

      write.csv(samples, file=filePath, row.names=FALSE, quote=FALSE)
      invisible(filePath)
    }
  )
)
