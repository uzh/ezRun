context("mLLMCelltype annotation on the FGCZ-internal vLLM")

VLLM_URL <- "http://fgcz-c-056:8000/v1/chat/completions"

vllmReachable <- function(url) {
  isTRUE(tryCatch(
    {
      jsonlite::fromJSON(sub("/chat/completions/*$", "/models", url))
      TRUE
    },
    error = function(e) FALSE
  ))
}

## Cluster ids are deliberately NOT in sorted order and include "10", so a
## regression that lets mLLMCelltype re-sort the clusters (which it does for
## bare numeric ids passed as a data.frame) shifts the labels and fails here.
topMarkers <- rbind(
  data.frame(cluster = "2", gene = c("LYZ", "CD14", "S100A8", "S100A9", "VCAN")),
  data.frame(cluster = "10", gene = c("CD3D", "CD3E", "IL7R", "TRAC", "CD2")),
  data.frame(cluster = "0", gene = c("MS4A1", "CD79A", "CD79B", "IGHM", "TCL1A")),
  stringsAsFactors = FALSE
)
topMarkers$cluster <- factor(topMarkers$cluster, levels = c("2", "10", "0"))

test_that("annotateClustersWithMLLMCelltype labels each cluster from its own markers", {
  skip_if_not_installed("mLLMCelltype")
  skip_if_not(vllmReachable(VLLM_URL), "FGCZ vLLM endpoint not reachable")

  res <- ezRun:::annotateClustersWithMLLMCelltype(
    topMarkers = topMarkers,
    tissueName = "human PBMC",
    vllmUrl = VLLM_URL
  )

  expect_setequal(names(res$clusterLabels), c("2", "10", "0"))
  expect_true(nzchar(res$model))
  ## labels are keyed by cluster id, and each carries the markers the model cited
  expect_setequal(names(res$clusterSupport), c("2", "10", "0"))
  expect_true(all(nzchar(res$clusterSupport)))
  ## Each label must match ITS OWN cluster's markers, not a neighbour's.
  expect_match(res$clusterLabels[["2"]], "[Mm]ono")
  expect_match(res$clusterLabels[["10"]], "T")
  expect_match(res$clusterLabels[["0"]], "B")
})

## Regression: at ~15 clusters the served reasoning model spends its whole token
## budget on a hidden reasoning trace and returns content = NULL for roughly one
## request in three, which mLLMCelltype reports as "Unexpected response format".
## registerFgczVllmProvider() turns thinking off; this test fails if that
## regresses. Three consecutive calls, because the failure is intermittent.
test_that("a realistic cluster count returns a label for every cluster", {
  skip_if_not_installed("mLLMCelltype")
  skip_if_not(vllmReachable(VLLM_URL), "FGCZ vLLM endpoint not reachable")

  panel <- c(
    "IL7R,CCR7,LEF1,TCF7,SELL,MAL,CD3D,CD3E,TRAC,LDHB",
    "CD14,LYZ,S100A8,S100A9,VCAN,FCN1,MNDA,CST3,CTSS,TYROBP",
    "CCL5,NKG7,GZMK,CST7,IL32,CD3D,TRAC,KLRG1,GZMA,CD8A",
    "MS4A1,CD79A,CD79B,IGHM,IGHD,TCL1A,BANK1,CD74,HLA-DRA,LINC00926",
    "GNLY,NKG7,GZMB,PRF1,KLRD1,FGFBP2,SPON2,CTSW,HOPX,CST7",
    "FCGR3A,MS4A7,CDKN1C,LST1,AIF1,IFITM3,SERPINA1,COTL1,SAT1,RHOC",
    "CD8A,CD8B,CCL5,GZMH,NKG7,GZMB,KLRD1,CST7,TRGC2,LAG3",
    "IGKC,IGHG1,JCHAIN,MZB1,DERL3,XBP1,TNFRSF17,CD38,PRDM1,SEC11C",
    "FCER1A,CST3,CLEC10A,CD1C,HLA-DQA1,ENHO,PLD4,GSN,IL1R2,CD74",
    "PPBP,PF4,NRGN,GNG11,SDPR,TUBB1,SPARC,CLU,HIST1H2AC,GP9",
    "LILRA4,CLEC4C,IL3RA,SERPINF1,ITM2C,IRF7,TCF4,GZMB,JCHAIN,PLD4",
    "TRDC,TRGC1,TRGC2,KLRC1,NKG7,GNLY,CD7,CTSW,KLRD1,IL2RB",
    "MKI67,TOP2A,STMN1,TUBB,HMGB2,TYMS,PCNA,CENPF,UBE2C,BIRC5",
    "HBB,HBA1,HBA2,ALAS2,SLC4A1,AHSP,CA1,HBD,SNCA,BNIP3L",
    "CD34,SPINK2,PRSS57,SOX4,EGFL7,CYTL1,NPR3,MYB,KIAA0125,STMN1"
  )
  big <- do.call(rbind, lapply(seq_along(panel), function(i) {
    data.frame(
      cluster = as.character(i - 1),
      gene = strsplit(panel[i], ",")[[1]],
      stringsAsFactors = FALSE
    )
  }))
  big$cluster <- factor(big$cluster, levels = as.character(seq_along(panel) - 1))

  for (i in 1:3) {
    res <- ezRun:::annotateClustersWithMLLMCelltype(
      topMarkers = big, tissueName = "human PBMC", vllmUrl = VLLM_URL
    )
    expect_length(res$clusterLabels, length(panel))
    expect_true(all(nzchar(res$clusterLabels)))
  }
})

test_that("an unreachable endpoint errors instead of returning empty labels", {
  skip_if_not_installed("mLLMCelltype")

  expect_error(
    suppressWarnings(ezRun:::annotateClustersWithMLLMCelltype(
      topMarkers = topMarkers,
      tissueName = "human PBMC",
      vllmUrl = "http://fgcz-c-056:9/v1/chat/completions"
    ))
  )
})
