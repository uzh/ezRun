EzAppNfCoreGeneric <- setRefClass(
  "EzAppNfCoreGeneric",
  contains = "EzApp",
  methods = list(
    ## Shared by every nf-core app on this class (NfCoreDenovotranscriptApp, nf_core_app_factory.rb), so
    ## every entry is gated on param$nfcorePipeline; pipelines not listed here get no references.
    ## denovotranscript (checked against the dev branch, the app's default pipelineVersion): nf-core,
    ## Nextflow, the pipeline and fastp unconditional (its skip flag cannot be set: cmdOptions is not
    ## passed to nextflow). SortMeRNA with its default SILVA and Rfam rRNA sets only when
    ## remove_ribo_rna is true. Unless qc_only is true: Trinity when assemblers lists trinity or
    ## trinity_no_norm, rnaSPAdes when it lists rnaspades, and EvidentialGene tr2aacds, BUSCO, rnaQUAST,
    ## TransRate (Oyster River Protocol build; runs with the singularity profile) and Salmon always.
    ## FastQC (run QC), gawk, MultiQC and the containers are not cited.
    citation = function(param = list()) {
      if (!identical(param$nfcorePipeline, "denovotranscript")) return(character(0))
      assembly <- !isTRUE(as.logical(param$qc_only))
      assemblers <- if (ezIsSpecified(param$assemblers)) param$assemblers else "trinity,rnaspades"
      assemblers <- trimws(strsplit(assemblers, ",")[[1]])
      trinity <- assembly && any(c("trinity", "trinity_no_norm") %in% assemblers)
      spades <- assembly && "rnaspades" %in% assemblers
      ribo <- isTRUE(as.logical(param$remove_ribo_rna))
      c(
        "Ewels, P.A. et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38(3), 276-278 (2020). https://doi.org/10.1038/s41587-020-0439-x",
        "Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35(4), 316-319 (2017). https://doi.org/10.1038/nbt.3820",
        "nf-core/denovotranscript. nf-core. https://nf-co.re/denovotranscript",
        "Chen, S. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2(2), e107 (2023). https://doi.org/10.1002/imt2.107",
        if (ribo) "Kopylova, E., Noé, L. & Touzet, H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics 28(24), 3211-3217 (2012). https://doi.org/10.1093/bioinformatics/bts611",
        if (ribo) "Quast, C. et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research 41(D1), D590-D596 (2013). https://doi.org/10.1093/nar/gks1219",
        if (ribo) "Burge, S.W. et al. Rfam 11.0: 10 years of RNA families. Nucleic Acids Research 41(D1), D226-D232 (2013). https://doi.org/10.1093/nar/gks1005",
        if (trinity) "Grabherr, M.G. et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol 29(7), 644-652 (2011). https://doi.org/10.1038/nbt.1883",
        if (trinity) "Haas, B.J. et al. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature Protocols 8(8), 1494-1512 (2013). https://doi.org/10.1038/nprot.2013.084",
        if (spades) "Bushmanova, E., Antipov, D., Lapidus, A. & Prjibelski, A.D. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience 8(9), giz100 (2019). https://doi.org/10.1093/gigascience/giz100",
        if (spades) "Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. & Korobeynikov, A. Using SPAdes de novo assembler. Current Protocols in Bioinformatics 70(1), e102 (2020). https://doi.org/10.1002/cpbi.102",
        if (assembly) "Gilbert, D.G. Longest protein, longest transcript or most expression, for accurate gene reconstruction of transcriptomes? bioRxiv (2019) [preprint, not peer-reviewed]. https://doi.org/10.1101/829184",
        if (assembly) "Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A. & Zdobnov, E.M. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution 38(10), 4647-4654 (2021). https://doi.org/10.1093/molbev/msab199",
        if (assembly) "Bushmanova, E., Antipov, D., Lapidus, A., Suvorov, V. & Prjibelski, A.D. rnaQUAST: a quality assessment tool for de novo transcriptome assemblies. Bioinformatics 32(14), 2210-2212 (2016). https://doi.org/10.1093/bioinformatics/btw218",
        if (assembly) "Smith-Unna, R., Boursnell, C., Patro, R., Hibberd, J.M. & Kelly, S. TransRate: reference-free quality assessment of de novo transcriptome assemblies. Genome Research 26(8), 1134-1144 (2016). https://doi.org/10.1101/gr.196469.115",
        if (assembly) "MacManes, M.D. The Oyster River Protocol: a multi-assembler and kmer approach for de novo transcriptome assembly. PeerJ 6, e5428 (2018). https://doi.org/10.7717/peerj.5428",
        if (assembly) "Patro, R., Duggal, G., Love, M.I., Irizarry, R.A. & Kingsford, C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods 14(4), 417-419 (2017). https://doi.org/10.1038/nmeth.4197"
      )
    },
    initialize = function() {
      "Initializes the application"
      callSuper()
    },
    run = function(input, output, param) {
      succeeded <- FALSE
      on.exit(.self$sendExitMail(output, param, succeeded), add = TRUE)
      .self$runPipeline(input, output, param)
      succeeded <- TRUE
      invisible(NULL)
    },
    sendExitMail = function(output, param, succeeded) {
      if (!ezValidMail(param[["mail"]])) {
        return(invisible(NULL))
      }
      tryCatch(
        {
          outputDataset <- EzDataset$new(meta = output, dataRoot = param[["dataRoot"]])
          links <- .self$outputLinks(outputDataset, param)
          status <- if (succeeded) "done" else "failed"
          subject <- paste0("nf-core/", param[["nfcorePipeline"]], " ", param[["name"]], " ", status)
          text <- if (succeeded) links else c(links, " ", geterrmessage())
          ezMail(text = text, subject = subject, to = param[["mail"]])
        },
        error = function(e) {
          ezLog("could not send the notification e-mail: ", conditionMessage(e))
        }
      )
      invisible(NULL)
    },
    runPipeline = function(input, output, param) {
      pipeline <- param[["nfcorePipeline"]]
      version <- param[["pipelineVersion"]]
      dataRoot <- param[["dataRoot"]]
      
      # Output to current directory (scratch)
      # Results will be copied to gstore by the job script footer
      outdir <- "_result"
      if (!is.null(output[["Result [File]"]])) {
        outdir <- basename(output[["Result [File]"]])
      }

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
      pipelineConfig <- system.file(
        "extdata", "nfcore", paste0(pipeline, ".config"),
        package = "ezRun"
      )
      if (nzchar(pipelineConfig)) {
        cmd <- paste(cmd, "-c", pipelineConfig)
      }

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
        ds <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, comment.char="", quote="")
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
        mj <- trimws(mj)
        if (!startsWith(mj, "{")) {
          mj <- paste0("{", mj, "}")
        }
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
        if (startsWith(src, "=")) {
          samples[[nfCol]] <- rep(substring(src, 2), nrow(ds))
          next
        }
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
