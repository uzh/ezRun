#!/usr/bin/env Rscript
## Verification loop for the XeniumSeurat SUSHI app.
##
## Runs the REAL ezMethodXeniumSeurat end to end on a synthetic Xenium bundle,
## then asserts the properties that past bugs violated. Every assertion here
## corresponds to a defect that actually shipped - it is a regression net, not a
## style checklist.
##
## It also runs NEGATIVE controls (each guard must block a bad config) and a
## POSITIVE control (a valid config must still run). A guard that never fires
## looks exactly like a guard that found nothing, so both directions are tested.
##
## Usage:
##   Rscript inst/verification/verify_XeniumSeurat.R [pkgRoot] [workDir]
##
##   pkgRoot  ezRun checkout whose R/ and inst/templates/ are under test
##            (default: two levels up from this script)
##   workDir  scratch dir on /srv/GT (NOT /tmp - node-local, invisible to a
##            compute node). Default /srv/GT/analysis/scratch/$USER/xenium_verify
##
## Exit code 0 = all checks passed, 1 = at least one failed.

args <- commandArgs(trailingOnly = TRUE)
PKG <- if (length(args) >= 1) normalizePath(args[1]) else
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
WORK <- if (length(args) >= 2) args[2] else
  file.path("/srv/GT/analysis/scratch", Sys.info()[["user"]], "xenium_verify")

`%||%` <- function(a, b) if (is.null(a)) b else a
suppressMessages({
  library(ezRun); library(Seurat); library(Matrix); library(data.table)
  library(qs2); library(writexl); library(ggplot2); library(dplyr)
})

## ---- load the code UNDER TEST, over whatever ezRun is installed ------------
for (f in c("R/xeniumQcMetrics.R", "R/xeniumCooccurrence.R", "R/app-XeniumSeurat.R")) {
  source(file.path(PKG, f))
}
cat("verifying:", PKG, "\n")
cat("R:", R.version.string, "| Seurat", as.character(packageVersion("Seurat")),
    "| spacexr", as.character(packageVersion("spacexr")), "\n")
cat("workdir:", WORK, "\n\n")

RESULTS <- new.env(); RESULTS$pass <- 0L; RESULTS$fail <- 0L; RESULTS$failed <- character()
ok <- function(label, cond) {
  good <- isTRUE(cond)
  if (good) RESULTS$pass <- RESULTS$pass + 1L else {
    RESULTS$fail <- RESULTS$fail + 1L; RESULTS$failed <- c(RESULTS$failed, label)
  }
  cat(sprintf("  [%s] %s\n", if (good) "PASS" else "FAIL", label))
}

## ---- build the synthetic input --------------------------------------------
IN <- file.path(WORK, "input"); RUN <- file.path(WORK, "run")
unlink(WORK, recursive = TRUE); dir.create(RUN, recursive = TRUE)
gen <- file.path(PKG, "inst/verification/make_synthetic_xenium.R")
stopifnot(file.exists(gen))
system2("Rscript", c(shQuote(gen), shQuote(IN)), stdout = TRUE)

## ---- point makeRmdReport at the template UNDER TEST ------------------------
localTmpl <- file.path(RUN, "templates"); dir.create(localTmpl)
file.copy(list.files(system.file("templates", package = "ezRun"), full.names = TRUE), localTmpl)
file.copy(file.path(PKG, "inst/templates/XeniumSeurat.Rmd"), localTmpl, overwrite = TRUE)
.orig <- base::system.file    # capture BEFORE shimming, else infinite recursion
assignInNamespace("system.file", function(..., package = "base", lib.loc = NULL, mustWork = FALSE) {
  a <- list(...)
  if (identical(package, "ezRun") && length(a) && identical(a[[1]], "templates")) {
    return(if (length(a) > 1) file.path(localTmpl, a[[2]]) else localTmpl)
  }
  .orig(..., package = package, lib.loc = lib.loc, mustWork = mustWork)
}, ns = "base")

mkInput <- function(nm = "VerifyRun") list(getNames = function() nm,
                                           getColumn = function(x) IN)
baseParam <- list(
  dataRoot = "", cores = 4,
  minCounts = 10, minFeatures = 5, qcNmads = 3,
  clusterResolution = 0.5, lambda = 0.8, nicheResolution = 0.5,
  banksyDims = 30, banksyKgeom = 30,
  rctdReference = "None", rctdFile = "", rctdClassFile = "",
  rctdUMImin = 10, rctdUMIminSigma = 300,
  doSPLIT = FALSE, splitMode = "neighborhood", splitNeighborThreshold = 0.05,
  coocRadius = 30, coocNperm = 100, coocFdr = TRUE,
  resultDir = RUN, name = "VerifyRun"
)

## ---- 1. the happy path must run and produce correct outputs ----------------
cat("\n=== 1. end-to-end run ===\n")
setwd(RUN)
t0 <- Sys.time()
res <- ezMethodXeniumSeurat(input = mkInput(), output = NA, param = baseParam)
cat("  returned:", res, "in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n\n")
setwd(file.path(RUN, "VerifyRun"))

cat("=== 2. output assertions ===\n")
ok("app returned Success", identical(res, "Success"))

## log.txt: ezWrite has no `append` arg and `con` only matches by name, so the
## old idiom printed corrupted text to stdout and never created the file.
ok("log.txt exists", file.exists("log.txt"))
if (file.exists("log.txt")) {
  lg <- readLines("log.txt")
  ok("log.txt accumulates (not truncated per call)", length(lg) > 5)
  ok("log.txt free of corrupted 'log.txtTRUE' text", !any(grepl("log.txtTRUE", lg)))
}

sc <- qs2::qs_read("scData.qs2", nthreads = 1)

## THE critical one: FindClusters(cluster.name=) overwrites seurat_clusters
## unconditionally on Seurat >= 5.5.1, so BANKSY niches replaced the clusters.
if ("banksy_cluster" %in% colnames(sc@meta.data)) {
  ok("seurat_clusters NOT overwritten by banksy niches",
     !identical(as.character(sc$seurat_clusters), as.character(sc$banksy_cluster)))
  cat(sprintf("       %d clusters vs %d niches\n",
              nlevels(factor(sc$seurat_clusters)), nlevels(factor(sc$banksy_cluster))))
}

## Xenium Explorer matches the RAW cell_id; RenameCells adds a sample prefix.
cl <- read.csv("clusters_for_explorer.csv")
ok("explorer cell_id has no sample prefix", !any(grepl("^VerifyRun_", cl$cell_id)))
if (file.exists("niches_for_explorer.csv")) {
  ok("clusters CSV differs from niches CSV",
     !identical(cl$group, read.csv("niches_for_explorer.csv")$group))
}

## Xenium is a targeted panel: HVG selection silently dropped 3001/5001 genes on Prime 5K.
ok("all panel genes are variable features",
   length(VariableFeatures(sc)) == nrow(sc))

## BANKSY assay is a dense scale.data block nothing reads; labels must survive.
ok("BANKSY assay dropped from the saved object", !("BANKSY" %in% SeuratObject::Assays(sc)))
ok("banksy_cluster labels retained", "banksy_cluster" %in% colnames(sc@meta.data))
ok("counts layer intact after DietSeurat",
   "counts" %in% SeuratObject::Layers(sc[["Xenium"]]))
ok("no reduction left pointing at a dropped assay",
   !("pca.banksy" %in% names(sc@reductions)))

## exploreSC only lists character/factor meta.data in "colour by".
ok("no logical meta.data columns (invisible in exploreSC)",
   !any(vapply(sc@meta.data, is.logical, logical(1))))

## Genomic Control is a real Xenium assay the QC helper used to ignore.
ok("GenomicControl counted in QC",
   "nCount_GenomicControl" %in% colnames(sc@meta.data) ||
     "GenomicControl" %in% names(sc@assays))

ok("scData.unfiltered.qs2 cleaned up", !file.exists("scData.unfiltered.qs2"))
ok("no junk use.qs.qs2 (makeRmdReport arg typo)", !file.exists("use.qs.qs2"))
ok("cell-type / cluster markers exported", file.exists("posMarkers.xlsx"))

ok("report 00index.html produced", file.exists("00index.html"))
if (file.exists("00index.html")) {
  h <- paste(readLines("00index.html", warn = FALSE), collapse = "\n")
  ok("report links to exploreSC", grepl("exploreSC", h))
  ok("report has no dead Vitessce link", !grepl("exploreVitessceXenium", h))
  ok("segmentation-health section rendered", grepl("Segmentation health", h))
  ok("no unrendered placeholder text", !grepl("figure pending", h, ignore.case = TRUE))
}

## ---- 3. negative controls: every guard must BLOCK --------------------------
cat("\n=== 3. guards must block bad configs ===\n")
guard <- function(label, p, expect) {
  d <- file.path(WORK, paste0("guard_", abs(as.integer(runif(1, 1, 1e6)))))
  dir.create(d, recursive = TRUE); old <- getwd(); setwd(d)
  msg <- tryCatch({ ezMethodXeniumSeurat(mkInput("G"), NA, p); "NO ERROR" },
                  error = function(e) conditionMessage(e))
  setwd(old); unlink(d, recursive = TRUE)
  ok(paste("blocks:", label), grepl(expect, msg, ignore.case = TRUE))
}
guard("splitMode='shift' without rctdClassFile",
      modifyList(baseParam, list(doSPLIT = TRUE, splitMode = "shift",
                                 rctdReference = "allen/allen_cortex_rctd.rds")),
      "requires rctdClassFile")
guard("doSPLIT without any RCTD reference",
      modifyList(baseParam, list(doSPLIT = TRUE)), "requires an RCTD reference")
guard("rctdFile that does not exist",
      modifyList(baseParam, list(rctdFile = "/nonexistent/ref.rds")), "rctdFile not found")
guard("rctdClassFile that does not exist",
      modifyList(baseParam, list(rctdClassFile = "/nonexistent/c.tsv")), "rctdClassFile not found")
guard("QC thresholds that remove every cell",
      modifyList(baseParam, list(minCounts = 1e9)), "QC removed all")

## ---- 4. positive control: a valid config must NOT be blocked ---------------
## Without this, "all guards blocked" is unfalsifiable - a function that always
## threw would score a perfect result above.
cat("\n=== 4. positive control ===\n")
d <- file.path(WORK, "poscontrol"); dir.create(d, recursive = TRUE)
old <- getwd(); setwd(d)
pc <- tryCatch({ ezMethodXeniumSeurat(mkInput("PosCtl"), NA,
                                      modifyList(baseParam, list(name = "PosCtl"))); "ran" },
               error = function(e) paste("ERROR:", conditionMessage(e)))
setwd(old)
ok("a valid config still runs to completion", identical(pc, "ran"))

## ---- summary ---------------------------------------------------------------
cat(sprintf("\n=== %d passed, %d failed ===\n", RESULTS$pass, RESULTS$fail))
if (RESULTS$fail > 0) {
  cat("failed checks:\n"); cat(paste0("  - ", RESULTS$failed, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("report for inspection:", file.path(RUN, "VerifyRun", "00index.html"), "\n")
quit(status = 0)
