## Guards on the LLM-written Methods text (R/methods-guards.R). Synthetic text only.

guardConfig <- paste(
  "param[['resolution']] = '0.6'", "param[['npcs']] = '20'", "param[['nfeatures']] = '3000'",
  "emptyDrops niters 100000", "window 2000", "threshold 0.5", sep = "\n")
guardAll <- paste(guardConfig, "Seurat_5.1.0", "R_4.6.0", "12,345 cells passed", "85.3% mapped", sep = "\n")
checkNum <- function(text, all = guardAll, n = 24)
  methods_check_numbers(text, guardConfig, all, sample_count = n)

test_that("planted result numbers, versions and parameters are flagged", {
  expect_identical(checkNum("In total 12,345 cells passed QC."), "12,345")
  expect_identical(checkNum("Of the reads, 85.3% of reads mapped."), "85.3%")
  expect_identical(checkNum("Clustering used Seurat 4.9.9."), "4.9.9")
  expect_identical(checkNum("Clusters were found at resolution 0.8."), "0.8")
})

test_that("configured values in other spellings, identifiers and small integers pass", {
  expect_length(checkNum("SCTransform selected 3,000 variable features."), 0)
  expect_length(checkNum("emptyDrops ran with 1e5 iterations."), 0)
  expect_length(checkNum("Windows of 2 kb were used."), 0)
  expect_length(checkNum("The threshold was [0.5], i.e. c(0.5)."), 0)
  expect_length(checkNum("Analyses ran in R version 4.6.0."), 0)
  expect_length(checkNum("R version 4.6.0 was used.", all = paste(guardConfig, "R version 4.6.0")), 0)
  expect_length(checkNum("The neighbour graph used PCs 1 to 20 (1:20, 1-20, 1–20)."), 0)
  expect_length(checkNum("Reads were aligned to GRCh38 from 10x Genomics libraries with log2FC and CD45."), 0)
  expect_length(checkNum("All 24 samples were processed in 3 batches."), 0)
})

test_that("number normalisation covers units, percentages and scientific forms", {
  cfg <- "a 0.05 b 1e-5 c 100000 d 3000000"
  expect_length(methods_check_numbers("5% and 10^-5 and 1E-5, 100k and 3 Mb.", cfg, cfg, 1), 0)
  expect_identical(methods_check_numbers("17% of 2e6", cfg, cfg, 1), c("17%", "2e6"))
  ## a year is a plain number: allowed only when the citation text (config) has it
  expect_identical(methods_check_numbers("Smith et al. (2019)", cfg, cfg, 1), "2019")
  expect_length(methods_check_numbers("Smith et al. (2019)", paste(cfg, "Smith 2019"), cfg, 1), 0)
})
