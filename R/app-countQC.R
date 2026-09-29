###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountQC = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  if (param$useFactorsAsSampleName) {
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(
      ezDesignFromDataset(dataset),
      1,
      paste,
      collapse = "_"
    ))
  }
  input$meta <- dataset
  rawData <- loadCountDataset(input, param)

  if (isError(rawData)) {
    writeErrorReport(htmlFile, param = param, error = rawData$error)
    return("Error")
  }

  metadata(rawData)$output <- output
  makeQuartoReport(
    rawData = rawData,
    qmdFile = "CountQC.qmd",
    reportTitle = "CountQC",
    colour = isTRUE(param$colour),
    number = isTRUE(param$number)
  )

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppCountQC <-
  setRefClass(
    "EzAppCountQC",
    contains = "EzApp",
    methods = list(
      ## ezRun unconditional; goseq + GO.db unless runGO is false (default TRUE; doGo also needs
      ## GO annotation, data-dependent); RUVSeq only when runRUV is true (loadCountDataset; no
      ## Ruby app declares it). pheatmap/WGCNA-dendrogram/cluster-validation stats deliberately
      ## excluded: no citable paper, or only a plotting helper is used rather than the method.
      citation = function(param = list()) {
        runGO <- !isFALSE(as.logical(param$runGO))
        c(
          "Rehrauer, H. et al. ezRun: An R meta-package for the analysis of Next Generation Sequencing Data. https://github.com/uzh/ezRun",
          if (runGO) "Young, M.D., Wakefield, M.J., Smyth, G.K. & Oshlack, A. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biology 11, R14 (2010). https://doi.org/10.1186/gb-2010-11-2-r14",
          if (runGO) "Carlson, M. GO.db: A set of annotation maps describing the entire Gene Ontology. R package. https://doi.org/10.18129/B9.bioc.GO.db",
          if (isTRUE(as.logical(param$runRUV))) "Risso, D., Ngai, J., Speed, T.P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32(9), 896-902 (2014). https://doi.org/10.1038/nbt.2931"
        )
      },
      methods_facts = function(param = list()) {
        known <- length(param) > 0
        ## the run's value, or the ezRun default (EZ_PARAM_DEFAULTS / appDefaults) when unknown
        val <- function(name, default) {
          v <- param[[name]]
          if (known && length(v) == 1 && !is.na(v)) v else paste("ezRun default", default)
        }
        runGO <- known && isTRUE(as.logical(param$runGO))
        ftest <- isTRUE(as.logical(param$selectByFtest))
        noSig <- known && isFALSE(as.logical(param$useSigThresh))
        norm <- if (known && ezIsSpecified(param$normMethod)) param$normMethod else NULL
        ## CountQC.qmd enrichr-markers: eval unique(dataset$featureLevel) != 'smRNA'
        fl <- methodsInput(param, "featureLevel")
        enrichr <- is.null(fl) || !all(fl == "smRNA")
        topN <- val("topGeneSize", 100)
        sel <- if (!known) {
          "highest standard deviation of the log2 signal across samples (or, with selectByFtest true, the smallest one-way ANOVA F-test p-value across conditions)"
        } else if (ftest) {
          "smallest one-way ANOVA F-test p-value across conditions (selectByFtest true, one lm per gene)"
        } else {
          "highest standard deviation of the log2 signal across samples (selectByFtest false)"
        }
        c(
          ## loadCountDataset (ngsio.R) presentFlag = counts > sigThresh (0 without useSigThresh);
          ## CountQC.qmd prepare-signal isValid, prepare-data-correlation (IsControl), filter-valid-features
          paste0("A gene was called present in a sample when its count exceeded ",
                 if (noSig) "0 (useSigThresh false)" else paste0("sigThresh (", val("sigThresh", 10), ")"),
                 "; the correlation, clustering and MDS plots used only the genes present in at least half of the samples of at least half of the conditions (combinations of the first two factors) and not flagged as control sequences, while the data files keep every gene."),
          ## CountQC.qmd prepare-signal: ezNorm(counts, method = normMethod) + minSignal; ezLogmeanScalingFactor (util.R)
          paste0("The normalized signal is the counts normalized with ",
                 if (is.null(norm)) "the normMethod method" else paste0("the ", norm, " method (normMethod)"),
                 " plus minSignal (", val("minSignal", 5), ")",
                 if (identical(norm, "logMean")) "; logMean scales each sample so that the geometric mean of its counts over the genes present in all samples equals the geometric mean of these per-sample values",
                 "."),
          ## CountQC.qmd filter-valid-features: log2signal <- log2(signal + backgroundExpression)
          paste0("backgroundExpression (", val("backgroundExpression", 10), ") is not a threshold and removed no genes: it is an offset added to the normalized signal before the log2 transform used for the correlation, clustering and MDS plots; gene presence is the sigThresh rule."),
          ## CountQC.qmd filter-valid-features (topGenes), plot-sample-correlation, sample-clustering, mds-plots
          paste0("topGeneSize (", topN, ") only sets the 'top genes' panels: the genes with the ", sel,
                 " feed the top-gene sample correlation plot, sample dendrogram and MDS plots; the gene-clustering heatmap",
                 if (runGO) " and the GO analysis do not use them." else " does not use them."),
          ## CountQC.qmd plot-sample-correlation, sample-clustering: hclust(as.dist(1 - cor(...)), "ward.D2")
          "Sample correlations are Pearson correlations of the log2 signal; with more than 3 samples, samples were clustered hierarchically (hclust, ward.D2) on 1 minus the Pearson correlation, and the sample tree was not cut into clusters.",
          ## CountQC.qmd clustering-high-variance; clusterPheatmap (go-analysis.R): pheatmap ward.D2, scale none, cutree(tree_row)
          paste0("The clustering heatmap clusters genes: with more than 3 samples, the up to ", val("maxGenesForClustering", 2000),
                 " genes (maxGenesForClustering) with the highest standard deviation of the log2 signal, minus those with a standard deviation at or below ",
                 val("highVarThreshold", 0.5), " (highVarThreshold), were row-centred without scaling and clustered with pheatmap (ward.D2), and the gene tree was cut into ",
                 val("nSampleClusters", 6), " gene clusters; nSampleClusters sets this number of gene clusters, not a number of sample clusters, and the step needs more than ",
                 val("minGenesForClustering", 30), " such genes (minGenesForClustering)."),
          ## goClusterResults / ezGoseq (go-analysis.R): goseq method "Hypergeometric", normalizedAvgSignal NULL, p.adjust fdr;
          ## universeProbeIds = rownames(seqAnno) (all genes of the table); .getGoTermsAsTd (go-reports.R) maxNumberOfTerms 40
          if (runGO) paste0("When the reference had GO annotation (doGo), each gene cluster of the heatmap was tested for GO over-representation (BP, MF, CC) with goseq using the hypergeometric method and no gene-length bias correction, against all genes of the count table annotated in that ontology and over GO terms with at least ",
                            val("minCountFisher", 3), " annotated genes (minCountFisher); p-values were Benjamini-Hochberg adjusted, and the tables list at most 40 terms per cluster with p below ",
                            val("pValThreshFisher", 1e-4), " (pValThreshFisher) and at least ", val("minCountFisher", 3), " cluster genes."),
          ## ezMdsPlotly (ngsPlots.R): plotMDS(logSignal, plot = FALSE); CountQC.qmd mds-plots
          "With more than 3 samples, the MDS plots are limma plotMDS coordinates (leading log2 fold change dimensions, plotMDS default settings) of the log2 signal of all present genes and of the top genes.",
          ## CountQC.qmd enrichr-markers: top 500 TPM genes per sample minus HRT housekeeping genes, JavaScript links only
          if (enrichr) "The Enrichr table only offers links that send each sample's top 500 genes by TPM, minus the housekeeping genes of the HRT atlas, and an all-sample list of the first 250 of each, to the Enrichr website; no enrichment is computed in the report.",
          ## ezMethodCountQC / CountQC.qmd: no DE model in the template
          "CountQC runs no differential expression test between conditions.",
          ## ezMethodCountQC: makeQuartoReport(qmdFile = "CountQC.qmd")
          "The HTML report is rendered by ezRun from its CountQC Quarto template."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCountQC
        name <<- "EzAppCountQC"
        appDefaults <<- rbind(
          runGO = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "whether to run the GO analysis"
          ),
          nSampleClusters = ezFrame(
            Type = "numeric",
            DefaultValue = 6,
            Description = "Number of SampleClusters, default value 6"
          ),
          selectByFtest = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "select topGenes by Test instead of SD"
          ),
          topGeneSize = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "number of genes to consider in gene clustering, mds etc"
          ),
          colour = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "coloured tab hierarchy in the report"
          ),
          number = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "numbered tab hierarchy in the report"
          )
        )
      }
    )
  )
