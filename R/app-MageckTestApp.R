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
  local_CondaEnv("gi_mageck", pathToMiniConda = "/usr/local/ngseq/miniforge3")

  ctrlFile <- list.files(
    param$libName,
    pattern = 'MAGeCK_Ctrl.csv$',
    full.names = TRUE
  )
  if (length(ctrlFile) == 1L && param$useControls) {
    opt <- c(opt, "--control-sgrna", ctrlFile)
  }
  if (!param$useControls && param$normalizationMethod == 'control') {
    warning(
      "No control file provided, switching normalization method to 'median'"
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
  system2("mageck", args = opt)

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
      designFile <- file.path(param$comparison, "mle_design.txt")
      ezWrite.table(dm, designFile, row.names = FALSE)
      mlePrefix <- file.path(param$comparison, paste0(output$getNames(), "_mle"))
      mleOpt <- c("mle", "-k", mergedCountFileLoc, "-d", designFile, "-n", mlePrefix)
      if (length(ctrlFile) == 1L && param$useControls) {
        mleOpt <- c(mleOpt, "--control-sgrna", ctrlFile)
      }
      system2("mageck", args = mleOpt)
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
    param = param,
    output = output,
    qmdFile = "MageckTest.qmd",
    reportTitle = param$comparison,
    number = TRUE
  )
  return("Success")
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
          )
        )
      }
    )
  )
