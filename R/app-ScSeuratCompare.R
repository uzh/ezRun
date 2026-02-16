###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCompare <-
  setRefClass(
    "EzAppScSeuratCompare",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSeuratCompare
        name <<- "EzAppScSeuratCompare"
        appDefaults <<- rbind(
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
          sccomp.variability = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "Whether to test for differential variability in sccomp"
          )
        )
      }
    )
  )

ezMethodScSeuratCompare = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  # writableRPackageDir <- file.path('/srv/GT/databases/writable_R_package', strsplit(version[['version.string']], ' ')[[1]][3])
  # .libPaths(c(writableRPackageDir, .libPaths()[length(.libPaths())]))

  library(Seurat)
  library(HDF5Array)
  library(SingleCellExperiment)
  library(qs2)
  library(tidyverse)
  library(sccomp)
  library(ComplexHeatmap)
  library(clusterProfiler)

  cache_stan_model <- system.file("stan", package = "sccomp", mustWork = TRUE)
  ## matches path in installation script
  cmdstanr::set_cmdstan_path(
    "/misc/ngseq12/src/CmdStan/cmdstan-2.36.0/cmdstan-2.36.0"
  )

  # dir.create((dirname(sccomp:::sccomp_stan_models_cache_dir)), showWarnings=FALSE)
  # ## remove the link if it is there
  # if (file.exists(sccomp:::sccomp_stan_models_cache_dir)){
  #   file.remove(sccomp:::sccomp_stan_models_cache_dir)
  # }
  # ## recreate the link to the current installed library
  # ezSystem(paste("ln -s", system.file("stan", package="sccomp", mustWork = TRUE), sccomp:::sccomp_stan_models_cache_dir))

  # ## link to the compiled sccomp models
  # ## This workaround is needed, since we do not have a writeable library folder
  # ## We have all compiled models in the library folder cached but because of a bug/feature
  # ## sccomp does not let you specify that cache folder but uses instead a hardcoded folder sccomp:::sccomp_stan_models_cache_dir
  # ## that can't be modified, so we:
  # ## create the expected directory
  # dir.create((dirname(sccomp:::sccomp_stan_models_cache_dir)), showWarnings=FALSE)
  # ## remove the link if it is there
  # if (file.exists(sccomp:::sccomp_stan_models_cache_dir)){
  #   file.remove(sccomp:::sccomp_stan_models_cache_dir)
  # }
  # ## recreate the link to the current installed library
  # ezSystem(paste("ln -s", system.file("stan", package="sccomp", mustWork = TRUE), sccomp:::sccomp_stan_models_cache_dir))
  # library(sccomp)

  # library(cmdstanr)
  # cmdstanr::set_cmdstan_path(file.path(writableRPackageDir, 'cmdstanr/cmdstan-2.36.0'))

  ## Setup sccomp
  # You need to install first sccomp and cmdstanr R packages in '/srv/GT/databases/writable_R_package'
  # Then you need to explicitely install cmdstan
  # install_cmdstan(dir = "/srv/GT/databases/writable_R_package/cmdstanr")

  # Create scratch directory
  # scratch_dir <- "/scratch/sccomp_output"
  # dir.create(scratch_dir, recursive = TRUE, mode = "0777", showWarnings = FALSE)
  #
  # # Load sccomp and set up cmdstan
  # library(sccomp)

  # # Load model
  # mod <- sccomp:::load_model(
  #   name = "glm_multi_beta_binomial",
  #   threads = 4,
  #   cache_dir = file.path(writableRPackageDir, "sccomp")
  # )
  #

  # Determine pseudobulk mode for DEG analysis
  # Note: This is independent of sccomp (which runs based on sample counts)
  pseudoBulkMode <- ezIsSpecified(param$replicateGrouping) &&
    param$pseudoBulkMode == "true"

  ###
  set.seed(38)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)

  scData <- ezLoadRobj(
    input$getFullPaths("SeuratObject"),
    nthreads = param$cores
  )

  # Auto-detect species/refBuild if not set
  if (!ezIsSpecified(param$refBuild) || param$refBuild == "") {
    # Try to get from scData misc slot
    if ("refBuild" %in% names(scData@misc) && !is.null(scData@misc$refBuild)) {
      param$refBuild <- scData@misc$refBuild
      ezLog(paste("Using refBuild from scData:", param$refBuild))
    } else {
      # Infer from gene names
      genes <- rownames(scData)
      # Sample 100 genes to check capitalization
      sample_genes <- head(genes[grepl("^[A-Z]", genes)], 100)

      # Human genes are typically ALL CAPS (ACTB, GAPDH, etc.)
      # Mouse genes are Title Case (Actb, Gapdh, etc.)
      uppercase_count <- sum(grepl("^[A-Z][A-Z]", sample_genes))
      titlecase_count <- sum(grepl("^[A-Z][a-z]", sample_genes))

      if (uppercase_count > titlecase_count) {
        param$refBuild <- "Homo_sapiens/Ensembl/GRCh38/Annotation/Release_110-2023-10-30"
        ezLog("Inferred Human species from gene names")
      } else {
        param$refBuild <- "Mus_musculus/Ensembl/GRCm39/Annotation/Release_109-2023-06-29"
        ezLog("Inferred Mouse species from gene names")
      }
    }
  } else {
    ezLog(paste("Using refBuild from parameters:", param$refBuild))
  }

  # Auto-detect best cell identity column
  # Prioritize: celltype > celltypeintegrated > manualAnnot > ident > seurat_clusters
  available_cols <- colnames(scData@meta.data)
  priority_cols <- c(
    "celltype",
    "celltypeintegrated",
    "cellTypeIntegrated",
    "manualAnnot",
    "ident"
  )

  for (col in priority_cols) {
    if (col %in% available_cols) {
      # Check if column has meaningful values (not all NA/empty)
      values <- scData@meta.data[[col]]
      if (!all(is.na(values)) && length(unique(values)) > 1) {
        ezLog(paste("Using", col, "as CellIdentity (auto-detected)"))
        param$CellIdentity <- col
        break
      }
    }
  }

  # If no priority column found and CellIdentity not set, default to seurat_clusters
  if (
    !ezIsSpecified(param$CellIdentity) ||
      !(param$CellIdentity %in% available_cols)
  ) {
    if ("seurat_clusters" %in% available_cols) {
      ezLog("Using seurat_clusters as CellIdentity (fallback)")
      param$CellIdentity <- "seurat_clusters"
    } else {
      ezLog("Using ident as CellIdentity (default fallback)")
      param$CellIdentity <- "ident"
    }
  } else {
    ezLog(paste(
      "Using",
      param$CellIdentity,
      "as CellIdentity (from parameters)"
    ))
  }

  DefaultAssay(scData) = "SCT"
  #subset the object to only contain the conditions we are interested in
  Idents(scData) <- scData@meta.data[[param$grouping]]
  stopifnot(c(param$sampleGroup, param$refGroup) %in% Idents(scData))
  scData <- subset(scData, idents = c(param$sampleGroup, param$refGroup))

  # Run sccomp if we have biological replicates (>=3 samples per condition)
  run_sccomp <- FALSE
  if (
    ezIsSpecified(param$replicateGrouping) &&
      param$replicateGrouping != ""
  ) {
    # Check if replicateGrouping column exists in metadata
    if (param$replicateGrouping %in% colnames(scData@meta.data)) {
      # Count samples per condition (simple base R approach)
      sample_condition_pairs <- unique(
        scData@meta.data[, c(param$replicateGrouping, param$grouping)]
      )
      sample_counts <- table(sample_condition_pairs[[param$grouping]])

      ezLog("Sample counts per condition:")
      print(sample_counts)

      # Need at least 3 samples per condition for sccomp
      if (all(sample_counts >= 3)) {
        run_sccomp <- TRUE
        ezLog(paste(
          "Running sccomp compositional analysis",
          "(>=3 samples per condition)"
        ))
      } else {
        ezLog(paste(
          "Skipping sccomp: Need >=3 samples per condition,",
          "found:",
          paste(sample_counts, collapse = ", ")
        ))
      }
    } else {
      ezLog(paste(
        "Skipping sccomp: replicateGrouping column",
        param$replicateGrouping,
        "not found in metadata"
      ))
    }
  } else {
    ezLog("Skipping sccomp: replicateGrouping parameter not specified")
  }

  if (run_sccomp) {
    # Run sccomp analysis
    sccomp_res <- scData |>
      sccomp_estimate(
        formula_composition = as.formula(paste("~", param$grouping)),
        sample = param$replicateGrouping,
        cell_group = param$CellIdentity,
        cores = as.integer(param$cores),
        output_directory = ".",
        cache_stan_model = cache_stan_model,
        verbose = TRUE
      )

    sccomp_res <- sccomp_res |>
      sccomp_remove_outliers(
        cores = as.integer(param$cores),
        cache_stan_model = cache_stan_model
      ) |>
      sccomp_test()

    # Save sccomp results
    saveRDS(sccomp_res, "sccomp_results.rds")
    ezLog("sccomp analysis completed and saved to sccomp_results.rds")
  }

  pvalue_allMarkers <- 0.05

  #Before calculating the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
  Idents(scData) <- scData@meta.data[[param$CellIdentity]]
  clusters_freq <- table(
    grouping = scData@meta.data[[param$grouping]],
    cellIdent = Idents(scData)
  ) %>%
    data.frame()
  small_clusters <- clusters_freq[clusters_freq$Freq < 10, "cellIdent"] %>%
    as.character() %>%
    unique()
  big_clusters <- setdiff(Idents(scData), small_clusters)

  if (length(slot(scData[['SCT']], "SCTModel.list")) > 2) {
    toKeep <- which(
      sapply(SCTResults(scData[['SCT']], slot = "cell.attributes"), nrow) != 0
    )
    slot(scData[['SCT']], "SCTModel.list") = slot(
      scData[['SCT']],
      "SCTModel.list"
    )[toKeep]
  }
  if (pseudoBulkMode) {
    scData_agg <- AggregateExpression(
      scData,
      assays = "RNA",
      return.seurat = TRUE,
      group.by = c(param$grouping, param$replicateGrouping, param$CellIdentity)
    )
    # Special case: If the param$CellIdentity is "ident", the resulting object
    # will not have "ident" in the metadata anymore, will just be "orig.ident"
    if (param$CellIdentity == "ident") {
      scData_agg$ident <- scData_agg$orig.ident
    }
    #Fix for strange bug in Seurat: it replaces '_' by '-' in the metadata columns
    scData_agg[[]][param$grouping] <- gsub(
      "-",
      "_",
      scData_agg[[param$grouping]][, 1]
    )
    Idents(scData_agg) <- scData_agg@meta.data[[param$CellIdentity]]
    consMarkers <- conservedMarkers(
      scData_agg,
      grouping.var = param$grouping,
      pseudoBulkMode = pseudoBulkMode
    )
    diffGenes <- diffExpressedGenes(
      scData_agg,
      param,
      grouping.var = param$grouping
    )
  } else {
    scData <- PrepSCTFindMarkers(scData)
    consMarkers <- conservedMarkers(
      scData,
      grouping.var = param$grouping,
      pseudoBulkMode = pseudoBulkMode
    )
    diffGenes <- diffExpressedGenes(
      scData,
      param,
      grouping.var = param$grouping
    )
  }

  # Save the files for the report
  writexl::write_xlsx(consMarkers, path = "consMarkers.xlsx")
  writexl::write_xlsx(diffGenes, path = "diffGenes.xlsx")
  qs2::qs_save(scData, "scData.qs2", nthreads = as.integer(param$cores))
  makeRmdReport(
    param = param,
    output = output,
    scData = scData,
    rmdFile = "ScSeuratCompare.Rmd",
    reportTitle = paste0(param$name)
  )
  return("Success")
}
