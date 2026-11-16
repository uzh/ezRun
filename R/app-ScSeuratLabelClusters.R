###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratLabelClusters <-
  setRefClass(
    "EzAppScSeuratLabelClusters",
    contains = "EzApp",
    methods = list(
      ## Seurat v5: every job (FindAllMarkers on the new labels; nothing is re-normalised or
      ## re-clustered, so no SCTransform). In getSeuratMarkersAndAnnotate, human/mouse only (getSpecies(refBuild)):
      ## Enrichr when enrichrDatabase is set; AUCell (SCENIC paper) + CellMarker 2.0 when tissue is set;
      ## SingleR + celldex when SingleR is not none; decoupleR + DoRothEA + PROGENy when
      ## computePathwayTFActivity (appDefault TRUE when absent). Azimuth (any species) when not none.
      ## Note: the input loading in ezMethodScSeuratLabelClusters (undefined `nparam` for a .qs2
      ## SC Seurat; allCellsMeta.rds and param.rds where makeRmdReport now writes .qs2) stops every
      ## current run before the analysis; the list is for the code path after loading.
      citation = function(param = list()) {
        annotate <- getSpecies(param$refBuild) %in% c("Human", "Mouse")
        enrichr <- annotate && ezIsSpecified(param$enrichrDatabase)
        aucell <- annotate && ezIsSpecified(param$tissue)
        singleR <- annotate && ezIsSpecified(param$SingleR) && !identical(param$SingleR, "none")
        tf <- annotate && (is.null(param$computePathwayTFActivity) || isTRUE(as.logical(param$computePathwayTFActivity)))
        azimuth <- ezIsSpecified(param$Azimuth) && !identical(param$Azimuth, "none")
        c(
          "Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y",
          if (enrichr) "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          if (enrichr) "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377",
          if (aucell) "Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 14, 1083-1086 (2017). https://doi.org/10.1038/nmeth.4463",
          if (aucell) "Hu, C. et al. CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data. Nucleic Acids Research 51, D870-D876 (2023). https://doi.org/10.1093/nar/gkac947",
          if (singleR) "Aran, D. et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology 20(2), 163-172 (2019). https://doi.org/10.1038/s41590-018-0276-y",
          if (singleR) "Aran, D. et al. celldex: Reference Index for Cell Types. R package. https://doi.org/10.18129/B9.bioc.celldex",
          if (tf) "Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016",
          if (tf) "Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. & Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research 29, 1363-1375 (2019). https://doi.org/10.1101/gr.240663.118",
          if (tf) "Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6",
          if (azimuth) "Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587 (2021). https://doi.org/10.1016/j.cell.2021.04.048"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSeuratLabelClusters
        name <<- "EzAppScSeuratLabelClusters"
        appDefaults <<- rbind(
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcoxon",
            Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
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
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          ),
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          computePathwayTFActivity = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Whether we should compute pathway and TF activities."
          )
        )
      }
    )
  )

ezMethodScSeuratLabelClusters <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(HDF5Array)
  library(stringr)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(scDblFinder)
  library(BiocParallel)
  library(scuttle)
  library(DropletUtils)
  library(decoupleR)
  library(tidyverse)

  if (param$cores > 1) {
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    ## scDblFinder fails with many cells and MulticoreParam
    BPPARAM <- SerialParam()
  }
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock (e.g. AUCell labeling).
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  register(BPPARAM)
  require(future)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = param$ram * 1024^3)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Cluster Report")))
  on.exit(setwd(cwd), add = TRUE)

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
    dplyr::rename(c("Sample" = 1, "Cluster" = 2, "ClusterLabel" = 3)) %>%
    dplyr::filter(Sample == input$getNames())
  labelMap <- as.character(clusterAnno$ClusterLabel)
  names(labelMap) <- as.character(clusterAnno$Cluster)

  # load cell data
  allCellsMeta <- ezLoadRobj(file.path(
    input$getFullPaths("SC Cluster Report"),
    "allCellsMeta.rds"
  ))
  scData <- ezLoadRobj(input$getFullPaths("SC Seurat"), nthreads = nparam$cores)
  param_scSeurat <- ezLoadRobj(file.path(
    dirname(input$getFullPaths("SC Seurat")),
    'param.rds'
  ))
  param$npcs <- param_scSeurat$npcs
  param$nfeatures <- param_scSeurat$nfeatures
  # change labels and store in a variable
  scData$cellType <- unname(labelMap[as.character(Idents(scData))])
  Idents(scData) <- scData$cellType

  # get markers and annotations
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)

  # save markers
  markers <- anno$markers
  writexl::write_xlsx(markers, path = "posMarkers.xlsx")

  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(
    Sample = input$getNames(),
    Cluster = levels(Idents(scData)),
    ClusterLabel = ""
  )
  if (!is.null(anno$singler.results)) {
    clusterInfos$SinglerCellType <- anno$singler.results$singler.results.cluster[
      clusterInfos$Cluster,
      "pruned.labels"
    ]
  }
  nTopMarkers <- 10
  topMarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(
    split(topMarkers$gene, topMarkers$cluster),
    paste,
    collapse = ", "
  )
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path = clusterInfoFile)
  qs2::qs_save(scData, "scData.qs2", nthreads = param$cores)

  makeRmdReport(
    param = param,
    output = output,
    scData = scData,
    allCellsMeta = allCellsMeta,
    enrichRout = anno$enrichRout,
    cells.AUC = anno$cells.AUC,
    singler.results = anno$singler.results,
    aziResults = anno$aziResults,
    pathwayActivity = anno$pathwayActivity,
    TFActivity = anno$TFActivity,
    rmdFile = "ScSeurat.Rmd",
    reportTitle = paste0(param$name, ": ", input$getNames())
  )
  #remove no longer used objects
  rm(scData)
  gc()
  return("Success")
}
