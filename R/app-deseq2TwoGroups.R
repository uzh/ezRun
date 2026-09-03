###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDeseq2 = function(input = NA, output = NA, param = NA) {
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  stopifnot(param$sampleGroup != param$refGroup)

  input = cleanupTwoGroupsInput(input, param)
  param$groupingName <- param$grouping
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1) {
    param$grouping2Name <- param$grouping2
    param$grouping2 = input$getColumn(param$grouping2)
  }

  rawData = loadCountDataset(input, param)
  if (isError(rawData)) {
    writeErrorReport("00index.html", param = param, error = rawData$error)
    return("Error")
  }

  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)) {
    writeErrorReport("00index.html", param = param, error = deResult$error)
    return("Error")
  }
  dds = metadata(deResult)$nativeResult$dds
  dataset <- data.frame(colData(deResult), check.names = FALSE)
  dataset <- dataset[rownames(dataset) %in% rownames(dds@colData), ]
  seqAnno <- data.frame(
    rowData(deResult),
    row.names = rownames(deResult),
    check.names = FALSE
  )

  ## -- exploreDE-compatible .h5ad (opt-in; must never fail the SUSHI job) ----
  if (ezIsSpecified(param$writeAnnData) && param$writeAnnData) {
    if (requireNamespace("exploreDE", quietly = TRUE)) {
      tryCatch(
        {
          factorCols <- paste(param$groupingName, "[Factor]")
          if (ezIsSpecified(param$grouping2Name)) {
            factorCols <- c(factorCols, paste(param$grouping2Name, "[Factor]"))
          }
          ## exploreDE's own organism inference (eds_infer_organism()) is
          ## package-internal / not exported -- reproduce its one-line logic
          ## here rather than reach into exploreDE:::.
          organism <- NA_character_
          if (ezIsSpecified(param$refBuild)) {
            organism <- gsub("_", " ", strsplit(param$refBuild, "/", fixed = TRUE)[[1]][1], fixed = TRUE)
          }
          enrichResult <- metadata(deResult)$enrichResult
          enrichResultGSEA <- metadata(deResult)$enrichResultGSEA
          pathwaysArg <- if (!is.null(enrichResult)) {
            setNames(list(list(ora = enrichResult, gsea = enrichResultGSEA)), param$comparison)
          } else {
            list()
          }
          exploreDE::build_explore_h5ad(
            x             = list(raw = assays(deResult)$counts, normalised = assays(deResult)$xNorm),
            out           = paste0("result--", param$comparison, ".h5ad"),
            omics         = "rnaseq",
            feature_level = param$featureLevel,
            organism      = organism,
            ref_build     = param$refBuild,
            coldata       = colData(deResult),
            rowdata       = rowData(deResult),
            factors       = factorCols,
            de            = setNames(list(rowData(deResult)), param$comparison),
            pathways      = pathwaysArg
          )
        },
        error = function(e) {
          ezLog("EzAppDeseq2: failed to write exploreDE .h5ad for ", param$comparison, ": ", conditionMessage(e))
        }
      )
    } else {
      ezLog("EzAppDeseq2: exploreDE package not available -- skipping .h5ad output")
    }
  }

  makeRmdReport(
    output = output,
    param = param,
    deResult = deResult,
    rmdFile = "twoGroups.Rmd",
    reportTitle = param$comparison
  )
  rmStatus <- file.remove(list.files(pattern = "enrichr-.*rds"))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeseq2(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppDeseq2 <-
  setRefClass(
    "EzAppDeseq2",
    contains = "EzApp",
    methods = list(
      ## DESeq2 unconditional; ashr gated on useLfcShrink; RUVSeq gated on runRUV;
      ## clusterProfiler/GO.db/Enrichr gated on runGO (Enrichr additionally on
      ## doEnrichr's organism/featureLevel check) -- listed regardless of gating.
      citation = function() {
        c(
          "Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8",
          "Stephens, M. False discovery rates: a new deal. Biostatistics 18(2), 275-294 (2017). https://doi.org/10.1093/biostatistics/kxw041",
          "Risso, D., Ngai, J., Speed, T.P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32(9), 896-902 (2014). https://doi.org/10.1038/nbt.2931",
          "Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141",
          "Bioconductor. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.23.1. https://doi.org/10.18129/B9.bioc.GO.db",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDeseq2
        name <<- "EzAppDeseq2"
        appDefaults <<- rbind(
          testMethod = ezFrame(
            Type = "character",
            DefaultValue = "deseq2",
            Description = "which test method in DESeq to use: deseq2"
          ),
          normMethod = ezFrame(
            Type = "character",
            DefaultValue = "DESeq2_MedianRatio",
            Description = "Deseq2's default norm method; this is actually not read"
          ),
          useRefGroupAsBaseline = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "should the log-ratios be centered at the reference samples"
          ),
          onlyCompGroupsHeatmap = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Only show the samples from comparison groups in heatmap"
          ),
          useLfcShrink = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "apply lfcShrink"
          ),
          rankMetric = ezFrame(
              Type = "character",
              DefaultValue = 'log2Ratio',
              Description = "how to rank genes for GSEA"
          ),
          writeAnnData = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Also write an exploreDE-compatible .h5ad file alongside the xlsx/Rmd report"
          )
        )
      }
    )
  )
