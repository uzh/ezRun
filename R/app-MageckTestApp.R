ezMethodMageckTest = function(input = NA, output = NA, param = NA) {
  require(Herper)
  require(stringr)
  require(limma)

  # Loading the variables
  dataset <- input$meta
  dir.create(param$comparison, showWarnings = FALSE)
  outputPrefix <- file.path(param$comparison, output$getNames())

  mergedCountFileName <- paste0(output$getNames(), ".merged.count.tsv")
  mergedCountFileLoc <- file.path(param$comparison, mergedCountFileName)
  sampleNames <- rownames(dataset)

  # Combining the count files into a single count file
  countFilePaths <- file.path(param$dataRoot, dataset$`Count [File]`)
  counts <- lapply(countFilePaths, data.table::fread)
  mergeCounts <- counts %>% reduce(inner_join, by = "sgRNA")
  mergeCounts[['Gene']] = NULL
  mergeCounts <- mergeCounts %>%
    rename("Gene" = "Gene.x") %>%
    select(!starts_with("Gene."))

  # Rename the sample columns to the actual sample names
  sampleColumns <- !(colnames(mergeCounts) %in% c("sgRNA", "Gene"))
  colnames(mergeCounts)[sampleColumns] <- sampleNames
  mergeCounts$Gene <- gsub(' ', '_', mergeCounts$Gene) #bug in Mageck regarding Spaces in GeneNames

  ezWrite.table(mergeCounts, file = mergedCountFileLoc, row.names = FALSE)

  # We give the design of the experiment as indices corresponding to the
  # columns (skipping the first 2 positions) which are sample vs ref groups
  fullColumnName <- paste(param$grouping, "[Factor]")
  rId <- paste(
    which(dataset[[fullColumnName]] == param$refGroup) - 1,
    collapse = ","
  )
  sId <- paste(
    which(dataset[[fullColumnName]] == param$sampleGroup) - 1,
    collapse = ","
  )

  # ---- Optional CRISPRcleanR copy-number-bias correction (before hit calling) ----
  # Corrects gene-independent CN artifacts by flattening biased genomic segments,
  # then feeds the corrected counts to mageck2. Off by default; only for human
  # screens with a coordinate-annotated library. Falls back to raw counts on any
  # failure so the app never breaks because of the optional step.
  ccrApplied <- FALSE
  ccrInfo <- NULL
  if (isTRUE(as.logical(param$useCRISPRcleanR))) {
    ccrRes <- tryCatch(
      ccrCorrectMergedCounts(
        mergeCounts,
        param,
        sampleNames,
        dataset[[fullColumnName]],
        mergedCountFileLoc
      ),
      error = function(e) {
        ezLog(paste(
          "CRISPRcleanR correction failed, using raw counts:",
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(ccrRes)) {
      mergedCountFileLoc <- ccrRes$correctedFile # mageck2 now reads corrected counts
      ccrInfo <- ccrRes$info
      ccrApplied <- TRUE
      ezLog(paste(
        "CRISPRcleanR applied:",
        ccrInfo$nSegments,
        "segments,",
        ccrInfo$nGuides,
        "guides (baseline:",
        ccrInfo$baseline,
        ")"
      ))
    }
  }

  opt <- c(
    "test",
    "-k",
    mergedCountFileLoc,
    "-t",
    sId,
    "-c",
    rId,
    "-n",
    outputPrefix,
    "--pdf-report",
    as.vector(str_split(param$cmdOptions, "\ +", simplify = TRUE))
  )

  # Load the conda environment
  local_CondaEnv("gi_mageck2", pathToMiniConda = "/usr/local/ngseq/miniforge3")

  ctrlFile <- list.files(
    param$libName,
    pattern = 'MAGeCK_Ctrl.csv$',
    full.names = TRUE
  )
  # CRISPRcleanR removes control (non-targeting) guides (no genomic locus), so the
  # corrected count table has none -- skip --control-sgrna and fall back to median
  # normalization when correction was applied.
  if (length(ctrlFile) == 1L && param$useControls && !ccrApplied) {
    opt <- c(opt, "--control-sgrna", ctrlFile)
  }
  if ((!param$useControls || ccrApplied) && param$normalizationMethod == 'control') {
    warning(
      "No usable control sgRNAs, switching normalization method to 'median'"
    )
    param$normalizationMethod <- 'median'
  }
  # Additional options based on parameters
  if (param$normalizationMethod != 'median') {
    opt <- c(opt, "--norm-method", param$normalizationMethod)
  }

  if (param$geneLFCMethod != 'median') {
    opt <- c(opt, "--gene-lfc-method", param$geneLFCMethod)
  }
  # Execute the command
  system2("mageck2", args = opt)

  geneSummaryFile <- paste0(outputPrefix, ".gene_summary.txt")
  sgrnaSummaryFile <- paste0(outputPrefix, ".sgrna_summary.txt")
  if (!file.exists(geneSummaryFile) || !file.exists(sgrnaSummaryFile)) {
    stop(
      "mageck test did not produce the expected summary files:\n",
      geneSummaryFile,
      "\n",
      sgrnaSummaryFile
    )
  }

  # add official gene symbol to gene_summary file for human/mouse samples
  if (param$species %in% c('hsa', 'mmu')) {
    dat <- ezRead.table(geneSummaryFile, row.names = NULL)
    dat[['GeneSymbol_Addgene']] = dat[['id']]
    for (j in 1:nrow(dat)) {
      gene <- c()
      if (param$species == 'hsa') {
        gene <- alias2Symbol(dat$id[j], species = "Hs")
      } else if (param$species == 'mmu') {
        gene <- alias2Symbol(dat$id[j], species = "Mm")
      }
      if (length(gene) == 1L) {
        dat[['id']][j] <- gene
      }
    }
    dat <- dat[!duplicated(dat$id), ]
    ezWrite.table(dat, geneSummaryFile, row.names = FALSE)
  }

  # We convert the raw outputs to xlsx files
  lapply(c(".sgrna_summary", ".gene_summary"), function(fileComp) {
    dat <- ezRead.table(
      paste0(outputPrefix, fileComp, ".txt"),
      row.names = NULL
    )
    writexl::write_xlsx(dat, paste0(outputPrefix, fileComp, ".xlsx"))
  })

  # Read the (symbol-annotated) results for the report
  geneRes <- ezRead.table(geneSummaryFile, row.names = NULL)
  sgrnaRes <- ezRead.table(sgrnaSummaryFile, row.names = NULL)

  # Optional MAGeCK MLE for the cross-condition nine-square (SquareView).
  # RRA gives a single pairwise contrast; the nine-square needs comparable
  # per-condition effect sizes (beta scores) relative to a shared baseline, which
  # is exactly what MLE estimates. Only runs when a day0/baseline condition is set.
  mleGeneRes <- NULL
  mleAxes <- NULL
  if (ezIsSpecified(param$day0Label)) {
    cond <- dataset[[fullColumnName]] # per-sample condition, in count-column order
    day0 <- param$day0Label
    safe <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    if (!(day0 %in% cond)) {
      ezLog(paste0(
        "day0Label '", day0, "' not found among '", param$grouping,
        "' values; skipping MLE nine-square."
      ))
    } else {
      nonDay0 <- setdiff(unique(cond), day0)
      dm <- data.frame(Samples = sampleNames, baseline = 1L, check.names = FALSE)
      for (cc in nonDay0) dm[[safe(cc)]] <- as.integer(cond == cc)
      ## MAGeCK2 requires the FIRST design-matrix row to be a baseline sample
      ## (all-zero condition columns); reorder day0 rows to the top. Rows are
      ## name-matched to the count columns via the Samples column, so reordering
      ## is safe.
      dm <- dm[order(cond != day0), , drop = FALSE]
      designFile <- file.path(param$comparison, "mle_design.txt")
      ezWrite.table(dm, designFile, row.names = FALSE)
      mlePrefix <- file.path(param$comparison, paste0(output$getNames(), "_mle"))
      mleOpt <- c("mle", "-k", mergedCountFileLoc, "-d", designFile, "-n", mlePrefix)
      if (length(ctrlFile) == 1L && param$useControls && !ccrApplied) {
        mleOpt <- c(mleOpt, "--control-sgrna", ctrlFile)
      }
      system2("mageck2", args = mleOpt)
      mleFile <- paste0(mlePrefix, ".gene_summary.txt")
      if (file.exists(mleFile)) {
        mleGeneRes <- ezRead.table(mleFile, row.names = NULL)
        # Nine-square axes: the two contrast conditions, each vs the day0 baseline.
        mleAxes <- list(ctrl = safe(param$refGroup), treat = safe(param$sampleGroup))
      } else {
        ezLog("MAGeCK MLE did not produce a gene_summary; skipping nine-square.")
      }
    }
  }

  # The report is rendered inside the comparison folder (00index.html + xlsx links).
  setwd(param$comparison)
  makeQuartoReport(
    geneRes = geneRes,
    sgrnaRes = sgrnaRes,
    mleGeneRes = mleGeneRes,
    mleAxes = mleAxes,
    ccrInfo = ccrInfo,
    param = param,
    output = output,
    qmdFile = "MageckTest.qmd",
    reportTitle = param$comparison,
    number = TRUE
  )
  return("Success")
}

##' @title Build a CRISPRcleanR sgRNA-coordinate annotation for a library
##' @description Returns a data.frame (rownames = the library's sgRNA IDs, columns
##'   CODE/GENES/CHRM/STARTpos/ENDpos/STRAND/seq) usable as CRISPRcleanR's
##'   \code{libraryAnnotation}. FGCZ libraries carry no genomic coordinates, so
##'   these are obtained on the fly: first by transferring them from CRISPRcleanR's
##'   built-in annotations (Brunello/KY/AVANA/GeCKO) matched by 20nt guide sequence,
##'   otherwise by aligning the guides to the genome. Cached per library dir.
getCRISPRcleanRlibAnnotation <- function(param) {
  require(CRISPRcleanR)
  libDir <- param$libName
  cacheFile <- file.path(
    libDir,
    paste0(basename(libDir), "_CRISPRcleanR_annotation.rds")
  )
  if (file.exists(cacheFile)) {
    return(readRDS(cacheFile))
  }
  dictFile <- list.files(libDir, pattern = "_MAGeCK\\.csv$", full.names = TRUE)
  if (length(dictFile) != 1L) {
    return(NULL)
  }
  fg <- ezRead.table(dictFile, header = FALSE, row.names = NULL, sep = ",")
  colnames(fg)[1:3] <- c("ID", "seq20", "gene")
  fg <- fg[
    !grepl("control", fg$ID, ignore.case = TRUE) &
      !grepl("non.?targeting|^control$", fg$gene, ignore.case = TRUE),
  ]
  fg$seq20 <- toupper(trimws(fg$seq20))

  ## Built-in path: map the library to a bundled annotation and transfer coords.
  builtins <- c(
    brunello = "Brunello_Library",
    ky = "KY_Library_v1.1",
    avana = "AVANA_Library",
    gecko = "GeCKO_Library_v2"
  )
  hit <- names(builtins)[vapply(
    names(builtins),
    function(k) grepl(k, libDir, ignore.case = TRUE),
    logical(1)
  )]
  annot <- NULL
  source <- NA_character_
  if (length(hit) >= 1L) {
    dname <- builtins[[hit[1]]]
    e <- new.env()
    utils::data(list = dname, package = "CRISPRcleanR", envir = e)
    b <- e[[dname]]
    ## auto-detect the spacer offset within the built-in context sequence
    offs <- vapply(
      1:8,
      function(o) mean(toupper(substr(b$seq, o, o + 19)) %in% fg$seq20),
      numeric(1)
    )
    off <- which.max(offs)
    if (max(offs) >= 0.5) {
      idx <- match(fg$seq20, toupper(substr(b$seq, off, off + 19)))
      ok <- !is.na(idx)
      annot <- data.frame(
        CODE = fg$ID[ok],
        GENES = fg$gene[ok],
        CHRM = b$CHRM[idx[ok]],
        STARTpos = b$STARTpos[idx[ok]],
        ENDpos = b$ENDpos[idx[ok]],
        STRAND = b$STRAND[idx[ok]],
        seq = fg$seq20[ok],
        stringsAsFactors = FALSE
      )
      source <- paste0("built-in:", dname)
    }
  }

  ## Fallback: align guides to the genome (only for unknown libraries).
  if (is.null(annot)) {
    annot <- ccrAlignGuidesToGenome(fg, param)
    if (!is.null(annot)) {
      source <- "genome-alignment"
    }
  }

  if (is.null(annot) || nrow(annot) < 100) {
    return(NULL)
  }
  rownames(annot) <- annot$CODE
  attr(annot, "source") <- source
  tmp <- paste0(cacheFile, ".tmp.", Sys.getpid())
  ok <- tryCatch(
    {
      saveRDS(annot, tmp)
      file.rename(tmp, cacheFile)
    },
    error = function(e) FALSE
  )
  return(annot)
}

##' @title Fallback: align sgRNA guides to the genome to get coordinates
##' @description Used only when the library has no CRISPRcleanR built-in annotation.
##'   Aligns the 20nt guides with bowtie2 to the species genome and keeps unique
##'   hits. Returns a coordinate data.frame or NULL. Human-focused; mouse needs a
##'   refBuild that points at a mouse genome index.
ccrAlignGuidesToGenome <- function(fg, param) {
  idxBase <- if (!is.null(param$refBuild) && nzchar(param$refBuild)) {
    file.path(param$dataRoot, dirname(param$refBuild), "Sequence", "BOWTIE2Index", "genome")
  } else if (identical(param$species, "hsa")) {
    "/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Bowtie2MainChromosomeIndex/coreChromosomes"
  } else {
    NA_character_
  }
  if (is.na(idxBase) || !length(Sys.glob(paste0(idxBase, "*.bt2")))) {
    ezLog("CRISPRcleanR fallback: no bowtie2 genome index found; cannot derive coordinates.")
    return(NULL)
  }
  wd <- tempfile("ccrAlign")
  dir.create(wd)
  faFile <- file.path(wd, "guides.fa")
  samFile <- file.path(wd, "guides.sam")
  writeLines(paste0(">", fg$ID, "\n", fg$seq20), faFile)
  cmd <- paste(
    "bowtie2 -x", shQuote(idxBase), "-f -U", shQuote(faFile),
    "-L 18 -N 0 --no-unal -k 2 -p 4 -S", shQuote(samFile)
  )
  res <- tryCatch(ezSystem(cmd), error = function(e) 1L)
  if (!file.exists(samFile)) {
    return(NULL)
  }
  sam <- tryCatch(
    data.table::fread(
      cmd = paste("grep -v '^@'", shQuote(samFile)),
      header = FALSE,
      sep = "\t",
      fill = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(sam) || !nrow(sam)) {
    return(NULL)
  }
  colnames(sam)[1:6] <- c("qname", "flag", "rname", "pos", "mapq", "cigar")
  sam <- sam[!bitwAnd(sam$flag, 4L) & sam$mapq >= 20, ] # mapped + ~unique
  sam <- sam[!duplicated(sam$qname), ]
  m <- match(fg$ID, sam$qname)
  ok <- !is.na(m)
  if (sum(ok) < 100) {
    return(NULL)
  }
  data.frame(
    CODE = fg$ID[ok],
    GENES = fg$gene[ok],
    CHRM = sub("^chr", "", sam$rname[m[ok]]),
    STARTpos = sam$pos[m[ok]],
    ENDpos = sam$pos[m[ok]] + 20L,
    STRAND = ifelse(bitwAnd(sam$flag[m[ok]], 16L) > 0, "-", "+"),
    seq = fg$seq20[ok],
    stringsAsFactors = FALSE
  )
}

##' @title Run CRISPRcleanR copy-number correction on a merged count table
##' @description Normalises, genome-sorts, segments (CBS) and corrects sgRNA fold
##'   changes, then writes corrected counts (columns restored to the original
##'   sample order so the mageck2 -t/-c indices stay valid). Returns a list with
##'   the corrected-count file path and before/after info for the report, or NULL
##'   if correction is not applicable.
ccrCorrectMergedCounts <- function(
  mergeCounts,
  param,
  sampleNames,
  cond,
  mergedCountFileLoc
) {
  require(CRISPRcleanR)
  if (!identical(param$species, "hsa")) {
    ezLog(
      "CRISPRcleanR: built-in annotations are human; skipping (species != hsa)."
    )
    return(NULL)
  }
  annot <- getCRISPRcleanRlibAnnotation(param)
  if (is.null(annot)) {
    ezLog("CRISPRcleanR: no coordinate annotation available; skipping correction.")
    return(NULL)
  }
  baseline <- if (ezIsSpecified(param$day0Label)) {
    param$day0Label
  } else {
    param$refGroup
  }
  isBase <- cond == baseline
  if (!any(isBase) || all(isBase)) {
    ezLog(paste0(
      "CRISPRcleanR: baseline condition '",
      baseline,
      "' unusable; skipping."
    ))
    return(NULL)
  }
  baseSamples <- sampleNames[isBase]
  otherSamples <- sampleNames[!isBase]
  ncontrols <- length(baseSamples)

  Dframe <- data.frame(
    sgRNA = mergeCounts$sgRNA,
    gene = mergeCounts$Gene,
    mergeCounts[, c(baseSamples, otherSamples), drop = FALSE],
    check.names = FALSE
  )
  normANDfcs <- ccr.NormfoldChanges(
    Dframe = Dframe,
    min_reads = 30,
    EXPname = param$comparison,
    libraryAnnotation = annot,
    ncontrols = ncontrols,
    saveToFig = FALSE,
    display = FALSE
  )
  gw <- ccr.logFCs2chromPos(normANDfcs$logFCs, annot)
  cc <- ccr.GWclean(gw, display = FALSE, label = param$comparison)
  corr <- as.data.frame(ccr.correctCounts(
    param$comparison,
    normANDfcs$norm_counts,
    cc,
    annot,
    minTargetedGenes = 3,
    ncontrols = ncontrols
  ))
  colnames(corr)[colnames(corr) == "gene"] <- "Gene"
  ## Restore the original sample-column order so mageck2 -t/-c indices stay valid.
  corrOut <- data.frame(
    sgRNA = corr$sgRNA,
    Gene = corr$Gene,
    corr[, sampleNames, drop = FALSE],
    check.names = FALSE
  )
  correctedFile <- sub("\\.tsv$", ".CCR.tsv", mergedCountFileLoc)
  ezWrite.table(corrOut, correctedFile, row.names = FALSE)

  cl <- cc$corrected_logFCs
  info <- list(
    rawGeneFC = tapply(cl$avgFC, cl$genes, mean),
    corrGeneFC = tapply(cl$correctedFC, cl$genes, mean),
    nSegments = nrow(cc$segments),
    nGuides = nrow(cl),
    nGuidesCorrected = sum(abs(cl$correctedFC - cl$avgFC) > 1e-6),
    annotationSource = attr(annot, "source"),
    baseline = baseline,
    ncontrols = ncontrols
  )
  list(correctedFile = correctedFile, info = info)
}

##' @template app-template
##' @templateVar method ezMethodMageckTest(input=NA, output=NA, param=NA)
##' @description Use this reference class to run Mageck Test
##' @author Falko Noé
EzAppMageckTest <-
  setRefClass(
    "EzAppMageckTest",
    contains = "EzApp",
    methods = list(
      ## mageck test unconditional. limma::alias2Symbol gated on species hsa/mmu.
      ## Downstream biology (essential-gene QC via bundled CEGv2/NEGv1, GO/KEGG
      ## over-representation, MSigDB GSEA, KEGG pathview) runs in MageckTest.qmd,
      ## not here. MAGeCKFlute intentionally NOT used: its useful outputs (RankView,
      ## essential-gene depletion, enrichment) are reproduced natively for styling
      ## control + CVD-safe palettes; FluteRRA is a validation-only reference.
      citation = function() {
        c(
          "Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4",
          "Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141",
          "Korotkevich, G. et al. Fast gene set enrichment analysis. bioRxiv (2021). https://doi.org/10.1101/060012",
          "Hart, T. et al. Evaluation and Design of Genome-Wide CRISPR/SpCas9 Knockout Screens. G3 7(8), 2719-2727 (2017). https://doi.org/10.1534/g3.117.041277",
          "Ritchie, M.E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47 (2015). https://doi.org/10.1093/nar/gkv007"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMageckTest
        name <<- "EzAppMageckTest"
        appDefaults <<- rbind(
          outputDir = ezFrame(
            Type = "character",
            DefaultValue = ".",
            Description = "Output directory"
          ),
          normalizationMethod = ezFrame(
            Type = "character",
            DefaultValue = "median",
            Description = "Normalization method"
          ),
          geneLFCMethod = ezFrame(
            Type = "character",
            DefaultValue = "median",
            Description = "Gene LFC method"
          ),
          useControls = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Use control sgRNAs"
          ),
          fdrThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "FDR cutoff for calling gene hits and labelling plots"
          ),
          nTopGenes = ezFrame(
            Type = "numeric",
            DefaultValue = 15,
            Description = "number of top genes per direction to label in plots"
          ),
          runEnrichment = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "run GO/KEGG over-representation on the hit sets"
          ),
          runGSEA = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "run MSigDB Hallmark/C2 GSEA on the ranked gene list"
          ),
          runPathview = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "draw KEGG pathview maps for top pathways (needs KEGG network access)"
          ),
          positiveControlGenes = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "optional comma-separated known positive-control genes to highlight"
          ),
          day0Label = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "baseline/day0 condition (plasmid/T0) enabling the MAGeCK MLE cross-condition nine-square (SquareView); empty = RRA test only"
          ),
          useCRISPRcleanR = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "optional CRISPRcleanR copy-number-bias correction of counts before RRA/MLE (human; best for fitness/dropout screens)"
          ),
          refBuild = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "genome for the CRISPRcleanR coordinate alignment fallback (only used for libraries without a built-in annotation)"
          )
        )
      }
    )
  )
