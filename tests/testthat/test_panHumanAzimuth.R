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
