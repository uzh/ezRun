###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodSplitPipe <- function(input = NA, output = NA, param = NA) {
  ## Parse Biosciences split-pipe pipeline (Evercode WT).
  ## Runs in DATASET mode: every row of the input is one sublibrary of the
  ## same experiment. Each sublibrary is processed with 'split-pipe --mode all'
  ## and, when more than one sublibrary is present, they are merged with
  ## 'split-pipe --mode comb' into a single per-sample result.

  sublibNames <- input$getNames()
  nSublib <- length(sublibNames)

  ## 1. Obtain (or build) the Parse genome reference
  refDir <- getParseReference(param)

  ## 2. Sample -> well mapping, shared across all sublibraries
  sampleArgs <- getParseSampleArgs(param)

  ## 3. Run split-pipe --mode all for each sublibrary
  read1List <- input$getFullPathsList("Read1")
  read2List <- input$getFullPathsList("Read2")
  sublibDirs <- character(0)
  for (nm in sublibNames) {
    fq1 <- read1List[[nm]]
    fq2 <- read2List[[nm]]
    if (length(fq1) != length(fq2)) {
      stop(sprintf(
        "Unequal number of Read1 (%d) and Read2 (%d) files for sublibrary %s",
        length(fq1),
        length(fq2),
        nm
      ))
    }
    outDir <- file.path("sublibs", nm)
    cmd <- paste(
      "split-pipe --mode all",
      paste0("--chemistry ", param$chemistry),
      paste0("--kit ", param$kit),
      paste0("--nthreads ", param$cores),
      "--fq1",
      paste(fq1, collapse = " "),
      "--fq2",
      paste(fq2, collapse = " "),
      paste0("--genome_dir ", refDir),
      paste0("--output_dir ", outDir),
      sampleArgs,
      if (ezIsSpecified(param$saveAnndata) && param$saveAnndata) {
        "--save_anndata"
      } else {
        ""
      },
      if (ezIsSpecified(param$cmdOptions)) param$cmdOptions else ""
    )
    ezSystem(cmd)
    sublibDirs <- c(sublibDirs, outDir)
  }

  ## 4. Combine sublibraries (or promote the single sublibrary to the result)
  resultDir <- param$name
  if (nSublib > 1) {
    ## Keep sublibs/ afterwards: each sublibrary's STAR Log.final.out (genome
    ## mapping rate for the report) lives under sublibs/<Name>/process/.
    cmd <- paste(
      "split-pipe --mode comb",
      paste0("--nthreads ", param$cores),
      "--sublibraries",
      paste(sublibDirs, collapse = " "),
      paste0("--output_dir ", resultDir)
    )
    ezSystem(cmd)
  } else {
    ezSystem(paste("mv", sublibDirs[1], resultDir))
    unlink("sublibs", recursive = TRUE)
  }

  ## 5. Convert each biological-sample DGE matrix to 10x format so the output is
  ## readable by the downstream ScSeurat / CellBender apps (which use Read10X /
  ## read10xCounts). The native Parse DGE_* directories are left in place.
  writeTenxMatrices(resultDir)

  ## 6. Report (stable FGCZ-Quarto 00index.html inside the result directory)
  makeSplitPipeReport(resultDir, param)

  return("Success")
}

##' Convert every Parse split-pipe DGE matrix under a result directory to 10x
##' format. For each biological-sample dir `<resultDir>/<sample>/` it writes
##' `filtered_feature_bc_matrix/` (from DGE_filtered) and `raw_feature_bc_matrix/`
##' (from DGE_unfiltered) using DropletUtils::write10xCounts.
writeTenxMatrices <- function(resultDir) {
  sampleDirs <- list.dirs(resultDir, recursive = FALSE)
  sampleDirs <- sampleDirs[dir.exists(file.path(sampleDirs, "DGE_filtered"))]
  for (sampleDir in sampleDirs) {
    convertParseToTenx(
      file.path(sampleDir, "DGE_filtered"),
      file.path(sampleDir, "filtered_feature_bc_matrix")
    )
    if (dir.exists(file.path(sampleDir, "DGE_unfiltered"))) {
      convertParseToTenx(
        file.path(sampleDir, "DGE_unfiltered"),
        file.path(sampleDir, "raw_feature_bc_matrix")
      )
    }
  }
}

##' Convert a single Parse DGE directory (count_matrix.mtx [cells x genes],
##' all_genes.csv, cell_metadata.csv) into a 10x v3 matrix directory
##' (matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz). The matrix is transposed
##' to the 10x genes x cells orientation; gene_id becomes the feature id (so
##' Read10X(gene.column = 1) yields Ensembl-ID rownames) and gene_name the symbol.
convertParseToTenx <- function(dgeDir, outDir) {
  ezLoadPackage("DropletUtils")
  m <- Matrix::readMM(file.path(dgeDir, "count_matrix.mtx")) ## cells x genes
  genes <- ezRead.table(
    file.path(dgeDir, "all_genes.csv"),
    sep = ",",
    row.names = NULL
  )
  cells <- ezRead.table(
    file.path(dgeDir, "cell_metadata.csv"),
    sep = ",",
    row.names = NULL
  )
  m <- Matrix::t(m) ## -> genes x cells (10x orientation)
  rownames(m) <- genes$gene_id
  colnames(m) <- cells$bc_wells
  unlink(outDir, recursive = TRUE)
  DropletUtils::write10xCounts(
    path = outDir,
    x = m,
    barcodes = cells$bc_wells,
    gene.id = genes$gene_id,
    gene.symbol = genes$gene_name,
    gene.type = "Gene Expression",
    version = "3",
    overwrite = TRUE
  )
}

##' Build the '--sample'/'--samp_list'/'--samp_sltab' arguments for split-pipe.
##' A sample loading table takes precedence; otherwise sampleWells is parsed as
##' a ';'-separated list of "<name> <wells>" specs (e.g. "all-well A1-A12" or
##' "sampleA A1-A6; sampleB A7-A12"). A ';'-separator is used rather than
##' newlines because SUSHI serialises each parameter onto a single R line.
getParseSampleArgs <- function(param) {
  if (ezIsSpecified(param$sampleLoadingTable)) {
    tbl <- param$sampleLoadingTable
    if (!file.exists(tbl)) {
      stop("sampleLoadingTable not found: ", tbl)
    }
    ext <- tolower(tools::file_ext(tbl))
    flag <- if (ext %in% c("xls", "xlsx")) "--samp_sltab" else "--samp_list"
    return(paste(flag, tbl))
  }
  if (!ezIsSpecified(param$sampleWells)) {
    ## Let split-pipe add the default all-well sample.
    return("--yes_allwell")
  }
  specs <- trimws(strsplit(param$sampleWells, ";")[[1]])
  specs <- specs[nzchar(specs)]
  if (length(specs) == 0) {
    return("--yes_allwell")
  }
  paste(paste("--sample", specs), collapse = " ")
}

##' Obtain the Parse (split-pipe mkref) genome directory, building it on the fly
##' when it does not yet exist. Modeled on getCellRangerGEXReference: the index
##' directory is derived from the reference GTF path and the selected transcript
##' types, so an existing index with the same transcript-type set is reused.
getParseReference <- function(param) {
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)

  if (ezIsSpecified(param$transcriptTypes)) {
    parseBase <- paste(sort(param$transcriptTypes), collapse = "-")
  } else {
    parseBase <- ""
  }
  refDir <- sub(
    "\\.gtf$",
    paste0("_Parse_SC_", parseBase, "_Index"),
    param$ezRef["refFeatureFile"]
  )

  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste(
      "reference building still in progress after",
      INDEX_BUILD_TIMEOUT,
      "min"
    ))
  }
  ## there is no lock file
  if (file.exists(refDir)) {
    ## we assume the index is built and complete
    return(refDir)
  }

  ## we have to build the reference
  setwd(dirname(refDir))
  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)

  job <- ezJobStart("Parse split-pipe mkref build")

  ## GTF filtered by the selected transcript types
  gtfFile <- tempfile(pattern = "genes", tmpdir = getwd(), fileext = ".gtf")
  on.exit(
    {
      if (file.exists(gtfFile)) file.remove(gtfFile)
    },
    add = TRUE
  )
  if (ezIsSpecified(param$transcriptTypes)) {
    export.gff2(gtfByTxTypes(param, param$transcriptTypes), con = gtfFile)
  } else {
    file.copy(from = param$ezRef@refFeatureFile, to = gtfFile)
  }

  ## Assembly name for --genome_name (split-pipe sanitises dots internally)
  refFields <- strsplit(param$refBuild, "/")[[1]]
  genomeName <- refFields[3]

  cmd <- paste(
    "split-pipe --mode mkref",
    paste0("--genome_name ", genomeName),
    paste0("--nthreads ", param$cores),
    paste0("--fasta ", param$ezRef@refFastaFile),
    paste0("--genes ", gtfFile),
    paste0("--output_dir ", refDir)
  )
  ezSystem(cmd)

  return(refDir)
}

##' Render a stable FGCZ-Quarto 00index.html report inside the result directory.
makeSplitPipeReport <- function(resultDir, param) {
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  setwd(resultDir)
  makeQuartoReport(
    param = param,
    qmdFile = "SplitPipe.qmd",
    reportTitle = paste("Parse split-pipe -", param$name),
    buttons = TRUE,
    number = TRUE,
    colour = FALSE
  )
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodSplitPipe(input=NA, output=NA, param=NA)
##' @description Use this reference class to run the Parse Biosciences split-pipe
##'   pipeline.
EzAppSplitPipe <-
  setRefClass(
    "EzAppSplitPipe",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSplitPipe
        name <<- "EzAppSplitPipe"
        appDefaults <<- rbind(
          kit = ezFrame(
            Type = "character",
            DefaultValue = "WT",
            Description = "Parse Evercode WT kit."
          ),
          chemistry = ezFrame(
            Type = "character",
            DefaultValue = "v3",
            Description = "Parse chemistry version."
          ),
          sampleWells = ezFrame(
            Type = "character",
            DefaultValue = "all-well A1-A12",
            Description = "';'-separated '<name> <wells>' specs for --sample."
          ),
          sampleLoadingTable = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Optional Parse sample loading table (overrides sampleWells)."
          ),
          transcriptTypes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "Transcript types to keep when building the reference."
          ),
          saveAnndata = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Also write an AnnData (.h5ad) output."
          )
        )
      }
    )
  )
