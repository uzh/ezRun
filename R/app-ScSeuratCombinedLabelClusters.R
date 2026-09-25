###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCombinedLabelClusters <-
  setRefClass(
    "EzAppScSeuratCombinedLabelClusters",
    contains = "EzApp",
    methods = list(
      ## Seurat defaults quoted here were checked against Seurat 5.5.1 formals() (R 4.6 system lib).
      ## Seurat v5 unconditional; Enrichr, AUCell/CellMarker 2.0, SingleR/celldex gated on human/mouse
      ## (AUCell on tissue, SingleR on SingleR != none); decoupleR/DoRothEA/PROGENy gated on
      ## computePathwayTFActivity; Azimuth gated on Azimuth inherited from the upstream run.
      citation = function() {
        c(
          "Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377",
          "Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 14, 1083-1086 (2017). https://doi.org/10.1038/nmeth.4463",
          "Hu, C. et al. CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data. Nucleic Acids Research 51, D870-D876 (2023). https://doi.org/10.1093/nar/gkac947",
          "Aran, D. et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology 20(2), 163-172 (2019). https://doi.org/10.1038/s41590-018-0276-y",
          "Aran, D. et al. celldex: Reference Index for Cell Types. R package. https://doi.org/10.18129/B9.bioc.celldex",
          "Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016",
          "Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. & Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research 29, 1363-1375 (2019). https://doi.org/10.1101/gr.240663.118",
          "Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6",
          "Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587 (2021). https://doi.org/10.1016/j.cell.2021.04.048"
        )
      },
      methods_facts = function(param = list()) {
        ## refBuild is not on this app's form: the job reads it from the upstream run
        ## (param.rds), so species stays conditional prose unless the parameters carry it.
        speciesKnown <- ezIsSpecified(param$refBuild)
        annotate <- !speciesKnown || methodsSpeciesIs(param, c("Human", "Mouse"))
        onlyHM <- if (speciesKnown) "." else " (human and mouse data only)."
        c(
          ## ezMethodScSeuratCombinedLabelClusters loads the saved object and only relabels (app-ScSeuratCombinedLabelClusters.R:98, 123-126)
          "No normalisation, integration, embedding or clustering was recomputed; the integrated object from the upstream ScSeuratCombine run was reused with its clusters and embeddings.",
          ## ClusterAnnotationFile parsing (app-ScSeuratCombinedLabelClusters.R:116-126)
          "Cluster labels were read from the second (cluster) and third (label) columns of the ClusterAnnotationFile, ignoring the first column; clusters given the same label were merged into one group, and clusters not listed in the file were left without a label.",
          ## ezUpdateMissingParam(param, oldParams) + refBuild override (app-ScSeuratCombinedLabelClusters.R:99-101; util.R:1097-1109)
          "Parameters absent from both this job's form and ezRun's global defaults (for example normalizationMethod, integrationMethod, npcs and resolution) were taken from the upstream ScSeuratCombine run's saved parameters, and refBuild was always taken from that run.",
          ## getSeuratMarkers via getSeuratMarkersAndAnnotate (app-ScSeuratCombinedLabelClusters.R:129; seuratUtils.R:795-825); the assay depends on the upstream run, so it stays conditional
          sprintf("Markers were recomputed for the new labels with Seurat FindAllMarkers on the object's default assay (SCT when the upstream run used SCTransform, RNA for LogNormalize), using the test in DE.method and only positive markers; p-values were Bonferroni-adjusted over all genes (Seurat p_val_adj), and markers with p_val_adj below %s (pvalue_allMarkers, ezRun default, not on the form) were kept.",
                  if (ezIsSpecified(param$pvalue_allMarkers)) param$pvalue_allMarkers else 0.05),
          ## getSeuratMarkers passes no latent.vars (seuratUtils.R:799-809)
          if ("LR" %in% param$DE.method) "Markers were tested with the LR test without latent variables; DE.regress was not applied.",
          ## getSeuratMarkersAndAnnotate (seuratUtils.R:847-877); Enrichr skipped without a database (app-ScSeurat.R)
          if (annotate && ezIsSpecified(param$enrichrDatabase)) paste0("The labelled groups were annotated with Enrichr on each group's significant markers (terms with adjusted p below 0.001 and more than 3 overlapping genes, top 5 per database)", onlyHM),
          ## cellsLabelsWithAUC skipped without a tissue (scTools.R:521-560)
          if (annotate && ezIsSpecified(param$tissue)) paste0("The labelled groups were annotated with AUCell using CellMarker 2.0 (2023-09-27 release) gene sets of at least 3 genes for the selected tissue", onlyHM),
          ## cellsLabelsWithSingleR skipped when SingleR is unset or none (scTools.R:630-660)
          if (annotate && ezIsSpecified(param$SingleR) && !("none" %in% param$SingleR)) paste0("Cells and labelled groups were annotated with SingleR on the fine labels (label.fine) of the celldex reference named in SingleR, per cell and per group", onlyHM),
          ## computeTFActivityAnalysis / computePathwayActivityAnalysis (app-ScSeurat.R)
          if (annotate && isTRUE(as.logical(param$computePathwayTFActivity))) paste0(sprintf(
            "Transcription-factor and pathway activities were inferred with decoupleR run_wmean (100 permutations, minsize 5, seed 42) using DoRothEA regulons of confidence A-C and PROGENy models (%s)",
            if (methodsSpeciesIs(param, "Human")) "top 500 genes per pathway"
            else if (methodsSpeciesIs(param, "Mouse")) "top 100 genes per pathway"
            else "top 500 genes per pathway for human, top 100 for mouse"), onlyHM)
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSeuratCombinedLabelClusters
        name <<- "EzAppScSeuratCombinedLabelClusters"
        appDefaults <<- rbind(
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          computePathwayTFActivity = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Whether we should compute pathway and TF activities."
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcox",
            Description = "Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"
          ),
          DE.regress = ezFrame(
            Type = "charVector",
            DefaultValue = "Batch",
            Description = "Variables to regress out if the test LR is chosen"
          ),
          min.pct = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
          ),
          logfc.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.25,
            Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
          )
        )
      }
    )
  )

ezMethodScSeuratCombinedLabelClusters = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(AUCell)
  library(decoupleR)
  library(qs2)
  library(BiocParallel)

  BPPARAM <- MulticoreParam(workers = param$cores)
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock (e.g. AUCell labeling).
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  register(BPPARAM)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  reportCwd <- getwd()

  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  filePath <- input$getFullPaths("SeuratObject")
  filePath_course <- file.path(
    "/srv/GT/analysis/course_sushi/public/projects",
    input$getColumn("SeuratObject")
  )

  if (!file.exists(filePath[1])) {
    filePath <- filePath_course
  }
  names(filePath) <- input$getNames()
  stopifnot(
    "App only supports single integrated dataset!" = length(input$getNames()) ==
      1
  )

  # Load previous dataset
  scData <- ezLoadRobj(filePath, nthreads = param$cores)
  oldParams <- ezLoadRobj(file.path(input$getFullPaths("Report"), "param.rds"))
  param <- ezUpdateMissingParam(param, oldParams)
  param$refBuild <- oldParams$refBuild

  # load cluster annotation file
  clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
  if (ezIsSpecified(param$ClusterAnnotationFile)) {
    clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
    stopifnot(
      "The cluster annotation file does not exist or is not an .xlsx file!" = file.exists(
        clusterAnnoFn
      ) &&
        str_ends(clusterAnnoFn, ".xlsx$")
    )
  } else {
    stop("Must supply cluster annotation file path.", call. = FALSE)
  }
  clusterAnno <- readxl::read_xlsx(clusterAnnoFn) %>%
    as_tibble() %>%
    dplyr::select(1:3) %>% # remove all other columns
    dplyr::rename(c("_" = 1, "Cluster" = 2, "ClusterLabel" = 3)) # we don't use the first column
  labelMap <- as.character(clusterAnno$ClusterLabel)
  names(labelMap) <- as.character(clusterAnno$Cluster)

  # Do the renaming
  scData$cellTypeIntegrated <- unname(labelMap[as.character(Idents(scData))])
  Idents(scData) <- scData$cellTypeIntegrated
  scData$ident <- Idents(scData)

  # perform all of the analysis
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)

  # save the markers
  writexl::write_xlsx(anno$markers, path = "posMarkers.xlsx")
  qs2::qs_save(scData, "scData.qs2", nthreads = param$cores)

  # Save some results in external files
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(
    param = param,
    output = output,
    scData = scData,
    enrichRout = anno$enrichRout,
    TFActivity = anno$TFActivity,
    pathwayActivity = anno$pathwayActivity,
    aziResults = anno$aziResults,
    cells.AUC = anno$cells.AUC,
    singler.results = anno$singler.results,
    rmdFile = "ScSeuratCombine.Rmd",
    reportTitle = reportTitle
  )
  return("Success")
}
