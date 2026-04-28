context("multiOmics utils")

test_that("detectModalities flags only RNA for plain CellRanger output", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "filtered_feature_bc_matrix.h5"))
  m <- detectModalities(d)
  expect_true(m$hasRNA)
  expect_false(m$hasADT)
  expect_false(m$hasVDJ_T)
  expect_false(m$hasVDJ_B)
  expect_false(m$hasATAC)
})

test_that("detectModalities flags VDJ-T and VDJ-B siblings (CellRanger Multi)", {
  root <- tempfile(); dir.create(file.path(root, "count"), recursive = TRUE)
  file.create(file.path(root, "count", "sample_filtered_feature_bc_matrix.h5"))
  dir.create(file.path(root, "vdj_t")); file.create(file.path(root, "vdj_t", "filtered_contig_annotations.csv"))
  dir.create(file.path(root, "vdj_b")); file.create(file.path(root, "vdj_b", "filtered_contig_annotations.csv"))
  m <- detectModalities(file.path(root, "count"))
  expect_true(m$hasVDJ_T)
  expect_true(m$hasVDJ_B)
})

test_that("detectModalities flags ATAC siblings (CellRanger ARC)", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "filtered_feature_bc_matrix.h5"))
  file.create(file.path(d, "atac_fragments.tsv.gz"))
  file.create(file.path(d, "atac_peaks.bed"))
  m <- detectModalities(d)
  expect_true(m$hasATAC)
})
