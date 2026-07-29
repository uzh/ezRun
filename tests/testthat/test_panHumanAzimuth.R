context("Pan-Human Azimuth species gate")

## Regression guard. Pan-Human Azimuth (CloudAzimuth) is human-only, but for a
## long time NOTHING enforced that: ScSeuratApp.rb carried "HUMAN DATASETS ONLY"
## only in a description string, and app-ScSeurat.R branched on the toggle alone.
## With the toggle default-on, an unguarded mouse job ships its matrix to a model
## keyed on human symbols (Gpr4 vs GPR4), where almost nothing matches and the
## annotation comes back confident and meaningless. These tests are cheap and
## need no network, no Seurat object and no AzimuthAPI.

test_that("the gate is off unless the toggle is explicitly on", {
  human <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"

  expect_false(ezRun:::panHumanAzimuthPlan(list(refBuild = human))$run)
  expect_false(
    ezRun:::panHumanAzimuthPlan(
      list(AzimuthPanHuman = FALSE, refBuild = human)
    )$run
  )
  expect_false(
    ezRun:::panHumanAzimuthPlan(
      list(AzimuthPanHuman = "false", refBuild = human)
    )$run
  )
  expect_false(
    ezRun:::panHumanAzimuthPlan(list(AzimuthPanHuman = "", refBuild = human))$run
  )
})

test_that("a human dataset with the toggle on runs, as logical or as string", {
  human <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"

  ## SUSHI passes params through an eval'd parameterset, so the toggle arrives
  ## as the STRING "true" via the Ruby path and as TRUE from a direct R call.
  ## Both must be accepted or the feature is silently off in production.
  expect_true(
    ezRun:::panHumanAzimuthPlan(
      list(AzimuthPanHuman = TRUE, refBuild = human)
    )$run
  )
  expect_true(
    ezRun:::panHumanAzimuthPlan(
      list(AzimuthPanHuman = "true", refBuild = human)
    )$run
  )
  expect_true(
    ezRun:::panHumanAzimuthPlan(
      list(AzimuthPanHuman = "TRUE", refBuild = human)
    )$run
  )
})

test_that("non-human datasets are skipped even when the toggle is on", {
  mouse <- "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"
  other <- "Danio_rerio/Ensembl/GRCz11/Annotation/Release_110"

  mousePlan <- ezRun:::panHumanAzimuthPlan(
    list(AzimuthPanHuman = "true", refBuild = mouse)
  )
  expect_false(mousePlan$run)
  expect_match(mousePlan$reason, "human-only")
  expect_match(mousePlan$reason, "Mouse")

  otherPlan <- ezRun:::panHumanAzimuthPlan(
    list(AzimuthPanHuman = TRUE, refBuild = other)
  )
  expect_false(otherPlan$run)
  expect_match(otherPlan$reason, "human-only")
})

test_that("a missing refBuild is a skip, not an unchecked run", {
  ## Fail closed: without a refBuild the species cannot be established, so the
  ## step must not run on the assumption that it is human.
  plan <- ezRun:::panHumanAzimuthPlan(list(AzimuthPanHuman = "true"))
  expect_false(plan$run)
  expect_match(plan$reason, "refBuild")

  plan2 <- ezRun:::panHumanAzimuthPlan(
    list(AzimuthPanHuman = "true", refBuild = "")
  )
  expect_false(plan2$run)
})

test_that("every plan reports a non-empty reason", {
  ## The reason is what gets logged, and a silent skip is the failure mode this
  ## whole change exists to remove.
  params <- list(
    list(AzimuthPanHuman = TRUE, refBuild = "Homo_sapiens/x"),
    list(AzimuthPanHuman = TRUE, refBuild = "Mus_musculus/x"),
    list(AzimuthPanHuman = FALSE, refBuild = "Homo_sapiens/x"),
    list(AzimuthPanHuman = TRUE)
  )
  for (p in params) {
    plan <- ezRun:::panHumanAzimuthPlan(p)
    expect_true(is.logical(plan$run) && length(plan$run) == 1L)
    expect_true(nzchar(plan$reason))
  }
})

context("Pan-Human Azimuth sentinel fill")

## panhumanpy's refiner returns a SENTINEL, not a label, whenever its 8 hierarchy
## heads contradict each other (annotate_tools.py:1503 returns the string "False").
## CellRanger fills that sentinel before writing its CSV; AzimuthAPI::CloudAzimuth
## does not, so the raw sentinel lands in the delivered Seurat object. Measured on
## p42258/o42614 (2026-07-29): 3,504 of 17,777 cells (19.7%) read "False" in
## azimuth_broad/medium/fine, making it the THIRD most common value in all three
## columns -- so exploreSC lists "False" as a browsable cell type and every DimPlot
## grouped by a tier column renders a legend entry called "False".
##
## The sentinel arrives in four forms (annotate.py:545 masks None/False/"False"/
## "false"); after the R round-trip we observed both the string and NA (130 NAs in
## azimuth_medium). All four must be caught.

test_that("every sentinel form is replaced, in all three tier columns", {
  meta <- data.frame(
    full_hierarchical_labels = c(
      "Immune cell|Lymphoid cell|T/NK cell|T cell|CD4 T cell|Naive CD4 T cell",
      "Immune cell|Lymphoid cell|B cell",
      "Immune cell",
      "Immune cell|Myeloid cell"
    ),
    azimuth_broad  = c("Immune cell", "False", NA, "false"),
    azimuth_medium = c("T cell", "False", "False", NA),
    azimuth_fine   = c("Naive CD4 T cell", NA, "False", "False"),
    stringsAsFactors = FALSE
  )

  out <- ezRun:::fillAzimuthSentinels(meta)

  sentinel <- function(x) is.na(x) | x %in% c("False", "false", "None")
  for (cc in c("azimuth_broad", "azimuth_medium", "azimuth_fine")) {
    expect_false(any(sentinel(out[[cc]])), info = cc)
  }
})

test_that("real labels are never touched", {
  meta <- data.frame(
    full_hierarchical_labels = "Immune cell|Lymphoid cell|T/NK cell|T cell",
    azimuth_broad = "Immune cell", azimuth_medium = "T cell",
    azimuth_fine = "GZMB CD8 T cell", stringsAsFactors = FALSE
  )
  out <- ezRun:::fillAzimuthSentinels(meta)
  expect_identical(out$azimuth_broad, "Immune cell")
  expect_identical(out$azimuth_medium, "T cell")
  expect_identical(out$azimuth_fine, "GZMB CD8 T cell")
})

test_that("a sentinel is filled from its PARENT tier, not from a placeholder", {
  ## Forward-fill is annotate.py:551's fill_unannotated=FALSE branch
  ## (medium_series[mask] <- broad_series[mask]). Chosen over the
  ## "Not Confidently Annotated" string so every tier keeps a real cell type,
  ## just a coarser one, and no cell drops out of a plot legend.
  meta <- data.frame(
    full_hierarchical_labels = "Immune cell|Lymphoid cell|T/NK cell",
    azimuth_broad = "Immune cell", azimuth_medium = "False",
    azimuth_fine = "False", stringsAsFactors = FALSE
  )
  out <- ezRun:::fillAzimuthSentinels(meta)
  expect_identical(out$azimuth_medium, "Immune cell")
  ## fine inherits the ALREADY-FILLED medium, matching the Python's order
  expect_identical(out$azimuth_fine, "Immune cell")
})

test_that("a sentinel broad falls back to the first hierarchy component", {
  ## Verified on 15,163 real cells: the first component of
  ## full_hierarchical_labels equals panhumanpy's level_1_labels AND
  ## CellRanger's broad_cell_type on 15,163/15,163 cells. So this fallback
  ## reproduces exactly what CellRanger's broad tier shows.
  meta <- data.frame(
    full_hierarchical_labels = c("Immune cell|Lymphoid cell", "Unassigned"),
    azimuth_broad = c("False", "False"),
    azimuth_medium = c("Lymphoid cell", "Unassigned"),
    azimuth_fine = c("Lymphoid cell", "Unassigned"), stringsAsFactors = FALSE
  )
  out <- ezRun:::fillAzimuthSentinels(meta)
  expect_identical(out$azimuth_broad, c("Immune cell", "Unassigned"))
})

test_that("Unassigned is preserved, not treated as a sentinel", {
  ## "Unassigned" is a TRAINED reject class (one of 13 level-0 output neurons,
  ## inference_encoders_level0.txt:13), not a low-confidence marker -- those
  ## cells carry a median confidence of 0.90. Rewriting it would be wrong.
  meta <- data.frame(
    full_hierarchical_labels = "Unassigned", azimuth_broad = "Unassigned",
    azimuth_medium = "Unassigned", azimuth_fine = "Unassigned",
    stringsAsFactors = FALSE
  )
  out <- ezRun:::fillAzimuthSentinels(meta)
  expect_identical(unname(unlist(out[1, c("azimuth_broad", "azimuth_medium", "azimuth_fine")])),
                   rep("Unassigned", 3))
})

test_that("a missing tier column is tolerated", {
  meta <- data.frame(full_hierarchical_labels = "Immune cell|Lymphoid cell",
                     azimuth_medium = "False", stringsAsFactors = FALSE)
  out <- ezRun:::fillAzimuthSentinels(meta)
  expect_identical(out$azimuth_medium, "Immune cell")
  expect_false("azimuth_broad" %in% colnames(out))
})

context("Pan-Human Azimuth: reuse CellRanger's local annotation")

## CellRanger >= 10.1.0 runs the SAME Pan-Human model locally and by default,
## writing outs/cell_types/Azimuth/cell_types.csv. Verified on p42258/o42614
## (15,163 cells, identical matrix): CellRanger's vendored ONNX fork and upstream
## panhumanpy 1.0.0 agree BIT-FOR-BIT on confidence (Pearson r = 1.0000, max
## |diff| = 0.0000) and on full_hierarchical_labels (15,163/15,163). So when that
## file exists there is nothing to gain by calling the remote API -- and a concrete
## cost: CloudAzimuth POSTs the expression matrix to azimuthapi.satijalab.org.

cellTypesFixture <- function(root) {
  d <- file.path(root, "outs", "cell_types", "Azimuth")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    data.frame(
      barcode = c("AAACCATTCACAGGCT-1", "AAACCCGCAATGTCGG-1", "AAACCAAAGCCCTACC-1"),
      broad_cell_type = c("Immune cell", "Immune cell", "Unassigned"),
      coarse_cell_type = c("T cell", "Coarse Not Confidently Annotated", "Unassigned"),
      fine_cell_type = c("Naive CD4 T cell", "Fine Not Confidently Annotated", "Unassigned"),
      full_hierarchical_labels = c(
        "Immune cell|Lymphoid cell|T/NK cell|T cell|CD4 T cell|Naive CD4 T cell",
        "Immune cell|Lymphoid cell|T/NK cell|T cell", "Unassigned"),
      final_level_softmax_prob = c(0.7774416, 0.6047713, 0.9525460),
      coarse_cell_type_unfiltered = c("T cell", "Coarse Not Confidently Annotated", "Unassigned"),
      umi_count = c(6656L, 1125L, 636L), stringsAsFactors = FALSE
    ),
    file.path(d, "cell_types.csv"), row.names = FALSE, quote = FALSE
  )
  file.path(root, "outs", "filtered_feature_bc_matrix")
}

test_that("the annotation is found from the CountMatrix path and mapped to ezRun columns", {
  root <- file.path(tempdir(), paste0("cr_", as.integer(Sys.time())))
  cmDir <- cellTypesFixture(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  ann <- ezRun:::readCellRangerPanHuman(cmDir)
  expect_false(is.null(ann))
  for (cc in c("azimuth_broad", "azimuth_medium", "azimuth_fine", "azimuth_label",
               "final_level_labels", "final_level_confidence",
               "full_hierarchical_labels")) {
    expect_true(cc %in% colnames(ann), info = cc)
  }
  expect_identical(rownames(ann)[1], "AAACCATTCACAGGCT-1")
  expect_identical(ann$azimuth_medium[1], "T cell")
  expect_equal(ann$final_level_confidence[1], 0.7774416, tolerance = 1e-6)
})

test_that("final_level_labels is the deepest hierarchy component", {
  ## Verified on 15,163 real cells: the last "|"-separated component of
  ## full_hierarchical_labels equals panhumanpy's final_level_labels on
  ## 15,163/15,163. ScSeurat.Rmd:109 gates the whole report section on this
  ## column, so it must be present and correct.
  root <- file.path(tempdir(), paste0("cr2_", as.integer(Sys.time())))
  cmDir <- cellTypesFixture(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  ann <- ezRun:::readCellRangerPanHuman(cmDir)
  expect_identical(ann$final_level_labels,
                   c("Naive CD4 T cell", "T cell", "Unassigned"))
  expect_identical(ann$azimuth_label, ann$final_level_labels)
})

test_that("an h5 CountMatrix path resolves to the same outs directory", {
  root <- file.path(tempdir(), paste0("cr3_", as.integer(Sys.time())))
  cellTypesFixture(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  h5 <- file.path(root, "outs", "filtered_feature_bc_matrix.h5")
  file.create(h5)
  expect_false(is.null(ezRun:::readCellRangerPanHuman(h5)))
})

test_that("absent or unreadable annotation returns NULL, never a partial object", {
  ## A silent half-populated return would put the report in exactly the state
  ## app-ScSeurat.R:807 exists to prevent.
  expect_null(ezRun:::readCellRangerPanHuman(file.path(tempdir(), "no_such_dir")))
  expect_null(ezRun:::readCellRangerPanHuman(NULL))
  expect_null(ezRun:::readCellRangerPanHuman(character(0)))

  root <- file.path(tempdir(), paste0("cr4_", as.integer(Sys.time())))
  d <- file.path(root, "outs", "cell_types", "Azimuth")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ## present but missing the columns we need
  write.csv(data.frame(barcode = "AAA-1", something_else = 1),
            file.path(d, "cell_types.csv"), row.names = FALSE, quote = FALSE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_null(ezRun:::readCellRangerPanHuman(file.path(root, "outs", "filtered_feature_bc_matrix")))
})

test_that("the CellRanger route carries no sentinel", {
  ## CellRanger already applied fill_unannotated=TRUE, so its columns are clean.
  ## Passing them through fillAzimuthSentinels must be a no-op.
  root <- file.path(tempdir(), paste0("cr5_", as.integer(Sys.time())))
  cmDir <- cellTypesFixture(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  ann <- ezRun:::readCellRangerPanHuman(cmDir)
  expect_identical(ezRun:::fillAzimuthSentinels(ann), ann)
})
