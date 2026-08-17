###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ktIsYes <- function(x) {
  isTRUE(x) || identical(tolower(as.character(x)), "yes") ||
    identical(tolower(as.character(x)), "true")
}

ktTaxids <- function(x) {
  v <- unlist(strsplit(as.character(x), "[,;[:space:]]+"))
  v[nzchar(v)]
}

ktPaths <- function(x) {
  v <- unlist(strsplit(as.character(x), "[[:space:],]+"))
  v[nzchar(v)]
}

## KrakenTools cannot read a gzipped kraken per-read output: gzip.open() is used
## only for the -s/-s2 sequence files. Decompress to scratch first.
ktPlainAssignments <- function(gzPath, sampleName) {
  plain <- paste0(sampleName, ".assignments.txt")
  ezSystem(paste("pigz -dc", shQuote(gzPath), ">", shQuote(plain)))
  plain
}

ktCountSeqs <- function(file, format) {
  if (!file.exists(file)) return(0)
  cmd <- if (format == "fastq") {
    paste("pigz -dc", shQuote(file), "| wc -l")
  } else {
    paste("pigz -dc", shQuote(file), "| grep -c '^>' || true")
  }
  n <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)[1]))
  if (is.na(n)) return(0)
  if (format == "fastq") n / 4 else n
}

ezMethodKrakenTools <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  tool <- as.character(param$tool)
  sampleNames <- input$getNames()
  ezLog(paste0("===== KrakenToolsApp: tool=", tool,
               " | ", length(sampleNames), " sample(s) ====="))

  extraOpts <- if (!is.null(param$cmdOptions)) as.character(param$cmdOptions) else ""

  if (tool == "extract_kraken_reads") {
    sampleName <- sampleNames[1]
    paired <- ktIsYes(param$paired) || isTRUE(param$paired)
    taxids <- ktTaxids(param$extract_taxids)
    if (length(taxids) == 0L) {
      stop("extract_taxids is empty; nothing to extract.")
    }
    fmt <- if (identical(as.character(param$extract_outputFormat), "fasta")) "fasta" else "fastq"

    plain <- ktPlainAssignments(input$getFullPaths("KrakenAssignments"), sampleName)
    report <- input$getFullPaths("KrakenReport")
    read1 <- input$getFullPaths("Read1")

    if (paired) {
      out1 <- paste0(sampleName, "_extracted_R1.", fmt)
      out2 <- paste0(sampleName, "_extracted_R2.", fmt)
      read2 <- input$getFullPaths("Read2")
      seqArgs <- paste("-s", shQuote(read1), "-s2", shQuote(read2))
      outArgs <- paste("-o", shQuote(out1), "-o2", shQuote(out2))
      outFiles <- c(out1, out2)
    } else {
      out1 <- paste0(sampleName, "_extracted.", fmt)
      seqArgs <- paste("-s", shQuote(read1))
      outArgs <- paste("-o", shQuote(out1))
      outFiles <- out1
    }

    cmd <- paste(
      "extract_kraken_reads.py",
      "-k", shQuote(plain),
      "-r", shQuote(report),
      "--taxid", paste(taxids, collapse = " "),
      if (ktIsYes(param$extract_includeChildren)) "--include-children" else "",
      if (ktIsYes(param$extract_includeParents)) "--include-parents" else "",
      if (ktIsYes(param$extract_exclude)) "--exclude" else "",
      if (fmt == "fastq") "--fastq-output" else "",
      "--max", as.character(param$extract_max),
      seqArgs, outArgs, extraOpts,
      ">", shQuote(paste0(sampleName, ".extract.log")), "2>&1"
    )
    ezSystem(cmd)

    ## A taxid absent from this sample is a legitimate result, but the advertised
    ## output must exist or the g-req copy fails and kills the job.
    for (f in outFiles) {
      if (!file.exists(f)) {
        ezLog(paste("No reads matched; creating empty", f))
        file.create(f)
      }
    }
    ezSystem(paste("pigz -f --best -p", ezThreads(), paste(shQuote(outFiles), collapse = " ")))
    ezSystem(paste("rm -f", shQuote(plain)))

    written <- ktCountSeqs(paste0(outFiles[1], ".gz"), fmt)
    stats <- data.frame(
      sample = sampleName,
      taxids = paste(taxids, collapse = ","),
      include_children = ktIsYes(param$extract_includeChildren),
      include_parents = ktIsYes(param$extract_includeParents),
      exclude_mode = ktIsYes(param$extract_exclude),
      format = fmt,
      records_written = written,
      stringsAsFactors = FALSE
    )
    ezWrite.table(stats, file = paste0(sampleName, ".extract_stats.tsv"), row.names = FALSE)
    if (written == 0) {
      ezLog("WARNING: 0 records extracted. Check that the taxids appear in the report and that --include-children is on.")
    }

  } else if (tool == "kreport2krona") {
    sampleName <- sampleNames[1]
    report <- input$getFullPaths("KrakenReport")
    kronaTxt <- paste0(sampleName, ".krona.txt")
    cmd <- paste(
      "kreport2krona.py",
      "-r", shQuote(report),
      "-o", shQuote(kronaTxt),
      if (ktIsYes(param$k2krona_intermediateRanks)) "--intermediate-ranks" else "--no-intermediate-ranks",
      extraOpts
    )
    ezSystem(cmd)
    if (ktIsYes(param$k2krona_renderHtml)) {
      ezSystem(paste("ktImportText", shQuote(kronaTxt),
                     "-o", shQuote(paste0(sampleName, ".html")),
                     "-n", shQuote(sampleName),
                     ">", shQuote(paste0(sampleName, ".krona.log")), "2>&1"))
    }

  } else if (tool == "kreport2mpa") {
    sampleName <- sampleNames[1]
    report <- input$getFullPaths("KrakenReport")
    cmd <- paste(
      "kreport2mpa.py",
      "-r", shQuote(report),
      "-o", shQuote(paste0(sampleName, ".mpa.txt")),
      if (ktIsYes(param$k2mpa_displayHeader)) "--display-header" else "",
      if (identical(as.character(param$k2mpa_values), "percentages")) "--percentages" else "--read_count",
      if (ktIsYes(param$k2mpa_intermediateRanks)) "--intermediate-ranks" else "--no-intermediate-ranks",
      if (identical(as.character(param$k2mpa_spaces), "keep")) "--keep-spaces" else "--remove-spaces",
      extraOpts
    )
    ezSystem(cmd)

  } else if (tool == "combine_kreports") {
    reports <- input$getFullPaths("KrakenReport")
    cmd <- paste(
      "combine_kreports.py",
      "-r", paste(shQuote(reports), collapse = " "),
      "-o", shQuote("combined.kreports_combined.txt"),
      "--sample-names", paste(shQuote(sampleNames), collapse = " "),
      if (ktIsYes(param$combineK_displayHeaders)) "--display-headers" else "--no-headers",
      if (ktIsYes(param$combineK_onlyCombined)) "--only-combined" else "",
      extraOpts
    )
    ezSystem(cmd)

  } else if (tool == "combine_mpa") {
    profiles <- input$getFullPaths("MpaProfile")
    cmd <- paste(
      "combine_mpa.py",
      "-i", paste(shQuote(profiles), collapse = " "),
      "-o", shQuote("combined.mpa_combined.txt"),
      extraOpts
    )
    ezSystem(cmd)

  } else if (tool == "filter_bracken.out") {
    sampleName <- sampleNames[1]
    bracken <- input$getFullPaths("BrackenAbundance")
    include <- ktTaxids(param$filterB_include)
    exclude <- ktTaxids(param$filterB_exclude)
    cmd <- paste(
      "filter_bracken.out.py",
      "-i", shQuote(bracken),
      "-o", shQuote(paste0(sampleName, ".filtered.bracken")),
      if (length(include)) paste("--include", paste(include, collapse = " ")) else "",
      if (length(exclude)) paste("--exclude", paste(exclude, collapse = " ")) else "",
      extraOpts
    )
    ezSystem(cmd)

  } else if (tool == "alpha_diversity") {
    sampleName <- sampleNames[1]
    bracken <- input$getFullPaths("BrackenAbundance")
    metric <- as.character(param$alpha_metric)
    ## alpha_diversity.py prints to stdout and has no -o.
    cmd <- paste(
      "alpha_diversity.py",
      "-f", shQuote(bracken),
      "-a", metric,
      extraOpts,
      ">", shQuote(paste0(sampleName, ".alpha_", metric, ".txt")), "2>&1"
    )
    ezSystem(cmd)

  } else if (tool == "beta_diversity") {
    brackens <- input$getFullPaths("BrackenAbundance")
    ## beta_diversity.py prints the matrix to stdout and has no -o.
    cmd <- paste(
      "beta_diversity.py",
      "-i", paste(shQuote(brackens), collapse = " "),
      "--type", as.character(param$beta_type),
      "--level", as.character(param$beta_level),
      "--cols", as.character(param$beta_cols),
      extraOpts,
      ">", shQuote("combined.beta_diversity.txt"), "2>&1"
    )
    ezSystem(cmd)

  } else if (tool == "make_kreport") {
    sampleName <- sampleNames[1]
    taxonomy <- as.character(param$mkreport_taxonomy)
    if (!file.exists(taxonomy)) {
      stop(paste("mkreport_taxonomy not found:", taxonomy))
    }
    plain <- ktPlainAssignments(input$getFullPaths("KrakenAssignments"), sampleName)
    cmd <- paste(
      "make_kreport.py",
      "-i", shQuote(plain),
      "-t", shQuote(taxonomy),
      "-o", shQuote(paste0(sampleName, ".report.txt")),
      if (ktIsYes(param$mkreport_useReadLen)) "--use-read-len" else "",
      extraOpts
    )
    ezSystem(cmd)
    ezSystem(paste("rm -f", shQuote(plain)))

  } else if (tool == "make_ktaxonomy") {
    for (p in c("ktax_nodes", "ktax_names", "ktax_seqid2taxid")) {
      path <- as.character(param[[p]])
      if (!file.exists(path)) stop(paste0(p, " not found: ", path))
    }
    cmd <- paste(
      "make_ktaxonomy.py",
      "--nodes", shQuote(as.character(param$ktax_nodes)),
      "--names", shQuote(as.character(param$ktax_names)),
      "--seqid2taxid", shQuote(as.character(param$ktax_seqid2taxid)),
      "-o", shQuote("ktaxonomy.taxonomy.txt"),
      extraOpts
    )
    ezSystem(cmd)

  } else if (tool == "fix_unmapped") {
    inFile <- as.character(param$fixun_input)
    if (!file.exists(inFile)) stop(paste("fixun_input not found:", inFile))
    refs <- ktPaths(param$fixun_accession2taxid)
    missing <- refs[!file.exists(refs)]
    if (length(missing)) {
      stop(paste("fixun_accession2taxid not found:", paste(missing, collapse = ", ")))
    }
    cmd <- paste(
      "fix_unmapped.py",
      "-i", shQuote(inFile),
      "--accession2taxid", paste(shQuote(refs), collapse = " "),
      "-o", shQuote("fix_unmapped.mapped.txt"),
      "-r", shQuote("fix_unmapped.still_unmapped.txt"),
      extraOpts
    )
    ezSystem(cmd)

  } else {
    stop(paste("Unsupported tool:", tool))
  }

  ezLog(paste0("===== KrakenToolsApp: ", tool, " finished ====="))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodKrakenTools()
##' @templateVar htmlArg )
##' @description Use this reference class to run KrakenTools scripts on the
##' output of a Kraken or Bracken run.
EzAppKrakenTools <-
  setRefClass(
    "EzAppKrakenTools",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodKrakenTools
        name <<- "EzAppKrakenTools"
        appDefaults <<- rbind(
          tool = ezFrame(
            Type = "character",
            DefaultValue = "extract_kraken_reads",
            Description = "which KrakenTools script to run"
          ),
          extract_taxids = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "comma/space separated NCBI taxids for extract_kraken_reads"
          ),
          extract_includeChildren = ezFrame(
            Type = "character",
            DefaultValue = "yes",
            Description = "pass --include-children (extract the whole clade)"
          ),
          extract_outputFormat = ezFrame(
            Type = "character",
            DefaultValue = "fastq",
            Description = "fastq (--fastq-output) or fasta"
          ),
          extract_max = ezFrame(
            Type = "character",
            DefaultValue = "100000000",
            Description = "--max ceiling on records written (KrakenTools default)"
          ),
          alpha_metric = ezFrame(
            Type = "character",
            DefaultValue = "Sh",
            Description = "alpha_diversity -a metric: BP, Sh, F, Si, ISi"
          ),
          beta_type = ezFrame(
            Type = "character",
            DefaultValue = "bracken",
            Description = "beta_diversity --type input format"
          )
        )
      }
    )
  )
