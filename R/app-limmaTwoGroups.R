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
      methods_facts = function() {
        c(
          ## ngsio.R:117-127; twoGroupCountComparison (twoGroups.R:87-95, 143-147)
          "A gene was called present in a sample when its count exceeded sigThresh (ezRun default 10); all genes were fitted, but only genes present in at least half of the samples of the sample group or of the reference group were counted as tested, and the reported FDR is the Benjamini-Hochberg adjustment of the limma p-values over these genes, computed by ezRun rather than taken from limma's adj.P.Val.",
          ## runLimma (twoGroups.R:421-431); calcNormFactors default method TMM
          "Only the samples of the sample and reference groups were fitted, with TMM normalization factors (edgeR calcNormFactors default) computed on those samples.",
          ## runLimma (twoGroups.R:433-455)
          "The linear model used the design ~ group with the reference group as baseline, and the reported log2 fold change and moderated t-test p-value are for the sample-group coefficient; when grouping2 is set it was treated as a blocking factor, with a consensus within-block correlation from duplicateCorrelation, not as a fixed covariate.",
          ## runLimma limma-trend branch (twoGroups.R:436-450); eBayes robust default FALSE
          "When modelMethod is limma-trend, counts were converted to log2 CPM with cpm (TMM-normalized library sizes, prior.count = priorCount), fitted with lmFit and moderated with eBayes(trend = TRUE), without robust estimation.",
          ## runLimma voom branch (twoGroups.R:451-460); voom normalize.method "none", eBayes trend/robust FALSE
          "When modelMethod is voom, counts were transformed with voom using the TMM-normalized library sizes (no further between-array normalization), fitted with lmFit using the voom precision weights and moderated with eBayes without trend or robust estimation; priorCount was not used."
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
