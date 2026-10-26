context("SpaceRanger: probes lost to the reference annotation")

## The app filters the Visium probe set to gene_ids in the reference geneInfo.tab.
## DEPRECATED_ and included=FALSE rows are harmless to drop; an included probe on
## a live gene missing from the annotation is a real loss and must be reported.

probes <- data.frame(
  gene_id = c("ENSG1", "ENSG1", "ENSG2", "DEPRECATED_ENSG3", "ENSG4", "ENSG5", "ENSG5"),
  probe_id = c("ENSG1|AAA|h1", "ENSG1|AAA|h2", "ENSG2|BBB|h3", "DEPRECATED_ENSG3|CCC|h4",
               "ENSG4|DDD|h5", "ENSG5|EEE|h6", "ENSG5|EEE|h7"),
  included = c("TRUE", "TRUE", "TRUE", "TRUE", "FALSE", "TRUE", "TRUE")
)
refIds <- c("ENSG1", "ENSG9")

test_that("only included, non-deprecated probes of missing genes count as dropped", {
  res <- findProbesMissingFromReference(probes, refIds)
  expect_equal(res$keep, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(res$dropped$probe_id, c("ENSG2|BBB|h3", "ENSG5|EEE|h6", "ENSG5|EEE|h7"))
  expect_equal(res$dropped$gene_name, c("BBB", "EEE", "EEE"))
  expect_equal(length(unique(res$dropped$gene_id)), 2)
})

test_that("a gene_name column wins over parsing probe_id, logical included works", {
  p <- probes
  p$gene_name <- paste0("sym", seq_len(nrow(p)))
  p$included <- as.logical(p$included)
  res <- findProbesMissingFromReference(p, refIds)
  expect_equal(res$dropped$gene_name, c("sym3", "sym6", "sym7"))
})

test_that("a matched annotation drops nothing real", {
  res <- findProbesMissingFromReference(probes, c("ENSG1", "ENSG2", "ENSG5"))
  expect_equal(nrow(res$dropped), 0)
  expect_equal(sum(!res$keep), 2) # DEPRECATED_ and included=FALSE still filtered
})

test_that("the warning names the annotation the probe set was built on", {
  d <- findProbesMissingFromReference(probes, refIds)$dropped
  expect_match(probesetAnnotationMismatchMessage(
    "Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv", d), "Release_32")
  expect_match(probesetAnnotationMismatchMessage(
    "Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv", d), "Release_42")
  expect_match(probesetAnnotationMismatchMessage(
    "Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv", d), "M31")
  expect_match(probesetAnnotationMismatchMessage(
    "Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv", d), "mm10-2020-A")
  msg <- probesetAnnotationMismatchMessage("custom.csv", d)
  expect_match(msg, "^3 probes for 2 genes")
  expect_false(grepl("this probe set matches", msg))
  expect_match(msg, "transcriptTypes")
})
