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
      ## limma/edgeR defaults quoted here were checked against limma 3.68.4 and
      ## edgeR 4.10.1 formals() (R 4.6.0 system lib).
      ## limma + edgeR (DGEList, calcNormFactors TMM, cpm) unconditional; eBayes(trend=TRUE) for
      ## modelMethod limma-trend; voom gated on modelMethod = voom; duplicateCorrelation gated on
      ## grouping2; goseq + GO.db gated on runGO + GO annotation; Enrichr gated on runGO + doEnrichr.
      citation = function() {
        c(
          "Ritchie, M.E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47 (2015). https://doi.org/10.1093/nar/gkv007",
          "Smyth, G.K. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3, 1-25 (2004). https://doi.org/10.2202/1544-6115.1027",
          "Law, C.W. et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29 (2014). https://doi.org/10.1186/gb-2014-15-2-r29",
          "Smyth, G.K. et al. Use of within-array replicate spots for assessing differential expression in microarray experiments. Bioinformatics 21, 2067-2075 (2005). https://doi.org/10.1093/bioinformatics/bti270",
          "Robinson, M.D. et al. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140 (2010). https://doi.org/10.1093/bioinformatics/btp616",
          "Chen, Y. et al. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. Nucleic Acids Research 53, gkaf018 (2025). https://doi.org/10.1093/nar/gkaf018",
          "Robinson, M.D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25 (2010). https://doi.org/10.1186/gb-2010-11-3-r25",
          "Young, M.D. et al. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biology 11, R14 (2010). https://doi.org/10.1186/gb-2010-11-2-r14",
          "Bioconductor. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.23.1. https://doi.org/10.18129/B9.bioc.GO.db",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## ngsio.R:117-127; twoGroupCountComparison (twoGroups.R:87-95, 143-147)
          "A gene was called present in a sample when its count exceeded sigThresh (ezRun default 10); all genes were fitted, but only genes present in at least half of the samples of the sample group or of the reference group were counted as tested, and the reported FDR is the Benjamini-Hochberg adjustment of the limma p-values over these genes, computed by ezRun rather than taken from limma's adj.P.Val.",
          ## runLimma (twoGroups.R:421-431); calcNormFactors default method TMM
          "Only the samples of the sample and reference groups were fitted, with TMM normalization factors (edgeR calcNormFactors default) computed on those samples.",
          ## runLimma (twoGroups.R:433-455)
          "The linear model used the design ~ group with the reference group as baseline, and the reported log2 fold change and moderated t-test p-value are for the sample-group coefficient.",
          if (ezIsSpecified(param$grouping2)) "grouping2 was treated as a blocking factor, with a consensus within-block correlation from duplicateCorrelation, not as a fixed covariate.",
          ## runLimma limma-trend branch (twoGroups.R:436-450); eBayes robust default FALSE
          if (identical(param$modelMethod, "limma-trend")) "Counts were converted to log2 CPM with cpm (TMM-normalized library sizes, prior.count = priorCount), fitted with lmFit and moderated with eBayes(trend = TRUE), without robust estimation.",
          ## runLimma voom branch (twoGroups.R:451-460); voom normalize.method "none", eBayes trend/robust FALSE
          if (identical(param$modelMethod, "voom")) "Counts were transformed with voom using the TMM-normalized library sizes (no further between-array normalization), fitted with lmFit using the voom precision weights and moderated with eBayes without trend or robust estimation; priorCount was not used."
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
