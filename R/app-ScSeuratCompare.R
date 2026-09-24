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
      ## Defaults quoted here were checked on the R version this app loads (Dev/R/4.5.0:
      ## Seurat 5.4.0, clusterProfiler 4.16.0, sccomp 2.1.30) and are identical on R 4.6 (Seurat 5.5.1).
      methods_facts = function() {
        c(
          ## set.seed(38) (app-ScSeuratCompare.R:81); sccomp mcmc_seed default sample_seed() = sample(1e5, 1)
          "The R random seed was set to 38 at the start, and sccomp's sampling seed was drawn from it.",
          ## subset to the two groups (app-ScSeuratCompare.R:164-168)
          "All analyses used only the cells whose grouping value is sampleGroup or refGroup.",
          ## CellIdentity auto-detection loop (app-ScSeuratCompare.R:121-162)
          "The cell identity used for per-group tests was chosen automatically as the first of the metadata columns celltype, celltypeintegrated, cellTypeIntegrated, manualAnnot or ident that has more than one value, overriding the CellIdentity parameter; if none qualified, CellIdentity was used when it names a metadata column, otherwise seurat_clusters.",
          ## refBuild inference from gene-name case (app-ScSeuratCompare.R:92-119); the app declares no refBuild parameter
          "Species for GO and KEGG enrichment was inferred from gene-symbol capitalisation (human when most symbols are all upper case, otherwise mouse), as the app takes no reference parameter.",
          ## sccomp block (app-ScSeuratCompare.R:170-235); sccomp 2.1.30 defaults
          "When replicateGrouping names a metadata column with at least 3 replicates in every condition, cell-type composition was tested with sccomp (formula ~ grouping, pathfinder inference, variability formula ~1 so no differential-variability test), after removing outliers with sccomp_remove_outliers; sccomp_test called effects above a logit fold change of 0.1 at 5% false positives (sccomp defaults).",
          ## non-pseudobulk branch: PrepSCTFindMarkers + diffExpressedGenes (app-ScSeuratCompare.R:289-301; seuratUtils.R:706-779)
          "Without pseudobulk mode, PrepSCTFindMarkers was run on the subset, and within each cell identity genes were tested between sampleGroup and refGroup cells with Seurat FindMarkers on the SCT assay (RNA as fallback if SCT failed), using DE.method with DE.regress as latent variables when DE.method is LR, the Seurat defaults logfc.threshold 0.1 and min.pct 0.01, and both directions; p-values were Bonferroni-adjusted over all genes (Seurat p_val_adj).",
          ## min.cells.group = 3 default; small_clusters computed but unused (app-ScSeuratCompare.R:239-249)
          "A cell identity was skipped when FindMarkers failed for it, in particular when either condition had fewer than 3 cells (fewer than 3 replicate pseudobulks in pseudobulk mode); no other minimum size was applied.",
          ## pseudobulk branch (app-ScSeuratCompare.R:260-288; seuratUtils.R:647-651, 715-719); Seurat DESeq2DETest
          "With pseudoBulkMode true and replicateGrouping set, RNA counts were summed per condition, replicate and cell identity (Seurat AggregateExpression), and differential expression and conserved markers used the DESeq2 Wald test through Seurat FindMarkers (local dispersion fit), overriding DE.method, with Seurat's Bonferroni adjustment over all genes.",
          ## conservedMarkers (seuratUtils.R:635-704); FindConservedMarkers meta.method default metap::minimump
          "Conserved markers of each cell identity across the two conditions were found with Seurat FindConservedMarkers (positive markers only; Wilcoxon test regardless of DE.method, or DESeq2 in pseudobulk mode), combining per-condition p-values with the minimum-p method (metap minimump, Seurat default) and ordering by the combined p-value; a condition with fewer than 3 cells of that identity was skipped, so some conserved markers rest on one condition only.",
          ## ScSeuratCompare.Rmd diff-genes table and enrichment chunks; enrichGO/enrichKEGG defaults
          "The reported differential expression table kept genes with p_val_adj below 0.05 (at most 1000, by absolute log2 fold change); for GO and KEGG, genes with p_val_adj below 0.05 and log2 fold change above 0.25 (up) or below -0.25 (down) in each cell identity were tested when at least 5, with clusterProfiler enrichGO (Biological Process, gene symbols, Benjamini-Hochberg, p cutoff 0.05, q cutoff 0.2, gene sets of 10-500 genes, all annotated genes of the organism database as background) and enrichKEGG after symbol-to-Entrez mapping (same cutoffs, KEGG data downloaded at run time).",
          ## ScSeuratCompare.Rmd run_pseudobulk_pca and MSE chunks
          "With at least 3 samples, the report's sample-level PCA used RNA counts summed per sample, log-normalised, on the top 2000 variable genes (vst) with up to 10 components, and per cell identity only samples with at least 5 cells, when at least 3 such samples existed; MSE distances used per-sample mean log-normalised RNA expression over up to 2000 variable features."
        )
      },
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
  # Load sccomp (+ deps) from the shared writable package lib, where the Stan model is
  # precompiled and world-readable. Avoids the system-lib copy whose model-cache dir is
  # frozen to the build user's home (unwritable by trxcopy on SLURM).
  writableRPackageDir <- file.path('/srv/GT/databases/writable_R_package', strsplit(version[['version.string']], ' ')[[1]][3])
  .libPaths(c(writableRPackageDir, .libPaths()))

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

  # sccomp freezes its model-cache dir to the build user's ~/.sccomp_models (a known sccomp
  # limitation - the path is a top-level constant baked at install time). Re-point it to the
  # SHARED, world-readable precompiled models inside the writable package lib so every user
  # (incl. trxcopy on SLURM) gets a cache hit and never compiles into the read-only pkg dir.
  utils::assignInNamespace(
    "sccomp_stan_models_cache_dir",
    system.file("stan", package = "sccomp", mustWork = TRUE),
    ns = "sccomp"
  )

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
