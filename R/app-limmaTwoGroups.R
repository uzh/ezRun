###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## for now only gene set testing from limma

EzAppLimma <-
  setRefClass(
    "EzAppLimma",
    contains = "EzApp",
    methods = list(
      ## limma, Smyth 2004 (lmFit/eBayes) and edgeR DGEList/TMM/cpm unconditional; voom only for
      ## modelMethod voom (default limma-trend); duplicateCorrelation (Smyth 2005) only when
      ## grouping2 is set; clusterProfiler/GO.db unless runGO is false (default TRUE): the
      ## GO step is clusterProfiler ezEnricher/ezGSEA (twoGroupCountComparison), as in DESeq2/edgeR;
      ## Enrichr only when twoGroups.Rmd precomputes it (doPrecomputeEnrichr, default TRUE).
      citation = function(param = list()) {
        runGO <- !isFALSE(as.logical(param$runGO))
        enrichr <- runGO && !isFALSE(as.logical(param$doPrecomputeEnrichr))
        c(
          "Ritchie, M.E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47 (2015). https://doi.org/10.1093/nar/gkv007",
          "Smyth, G.K. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3, 1-25 (2004). https://doi.org/10.2202/1544-6115.1027",
          if (identical(param$modelMethod, "voom")) "Law, C.W. et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29 (2014). https://doi.org/10.1186/gb-2014-15-2-r29",
          if (ezIsSpecified(param$grouping2)) "Smyth, G.K. et al. Use of within-array replicate spots for assessing differential expression in microarray experiments. Bioinformatics 21, 2067-2075 (2005). https://doi.org/10.1093/bioinformatics/bti270",
          "Robinson, M.D. et al. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140 (2010). https://doi.org/10.1093/bioinformatics/btp616",
          "Chen, Y. et al. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. Nucleic Acids Research 53, gkaf018 (2025). https://doi.org/10.1093/nar/gkaf018",
          "Robinson, M.D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25 (2010). https://doi.org/10.1186/gb-2010-11-3-r25",
          if (runGO) "Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141",
          if (runGO) "Carlson, M. GO.db: A set of annotation maps describing the entire Gene Ontology. R package. https://doi.org/10.18129/B9.bioc.GO.db",
          if (enrichr) "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          if (enrichr) "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodLimma
        name <<- "EzAppLimma"
        appDefaults <<- rbind(
          testMethod = ezFrame(
            Type = "character",
            DefaultValue = "limma",
            Description = "which test method in limma to use: limma"
          ),
          modelMethod = ezFrame(
            Type = "character",
            DefaultValue = "limma-trend",
            Description = "which mean-variance relationship model method in limma to use: limma-trend or voom"
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
          priorCount = ezFrame(
            Type = "numeric",
            DefaultValue = 10,
            Description = "prior count to be added to shrink the log-fold-changes"
          )
        )
      }
    )
  )

ezMethodLimma = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))

  stopifnot(param$sampleGroup != param$refGroup)

  input = cleanupTwoGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1) {
    param$grouping2 = input$getColumn(param$grouping2)
  }

  rawData = loadCountDataset(input, param)
  if (isError(rawData)) {
    writeErrorReport(htmlFile, param = param, error = rawData$error)
    return("Error")
  }

  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)) {
    writeErrorReport(htmlFile, param = param, error = deResult$error)
    return("Error")
  }

  makeRmdReport(
    output = output,
    param = param,
    deResult = deResult,
    rmdFile = "twoGroups.Rmd"
  )
  rmStatus <- file.remove(list.files(pattern = "enrich-.*rds"))
  return("Success")
}
