#!/usr/bin/env Rscript
## Generate a fully SYNTHETIC Xenium output bundle.
##
## Why synthetic: verifying the XeniumSeurat app needs a real Xenium directory
## (LoadXenium reads cell_feature_matrix + cells.csv.gz), but real Xenium runs
## are patient tissue. Anything above "public" in the UZH/ETH classification
## must not leave FGCZ, which rules out pasting real report figures into a
## ticket, a chat, or a cloud agent's context. This generator has no provenance
## question at all, so the resulting report can be shared freely.
##
## It plants K cell types into spatial domains with distinct marker programs, so
## clustering, BANKSY niches and co-occurrence all have real structure to find -
## which is what makes the assertions in verify_XeniumSeurat.R meaningful rather
## than merely "it did not crash".
##
## Usage:
##   Rscript make_synthetic_xenium.R <outDir> [nCells] [nGenes] [seed]

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else stop("usage: make_synthetic_xenium.R <outDir> [nCells] [nGenes] [seed]")
nCells <- if (length(args) >= 2) as.integer(args[2]) else 6000L
nGenes <- if (length(args) >= 3) as.integer(args[3]) else 300L
seed <- if (length(args) >= 4) as.integer(args[4]) else 7L

suppressMessages({library(Matrix); library(data.table)})
set.seed(seed)
dir.create(file.path(OUT, "cell_feature_matrix"), recursive = TRUE, showWarnings = FALSE)

K <- 6
typeNames <- c("Tumour", "Fibroblast", "Tcell", "Macrophage", "Endothelial", "Epithelial")
nMark <- max(5L, floor(nGenes / (K * 2)))

## --- spatial layout: one blob per type, immune types deliberately interleaved
centres <- data.frame(x = c(1200, 2600, 1500, 1800, 3000, 600),
                      y = c(1200, 1400, 2400, 1300, 2600, 2200))
prob <- c(.34, .18, .13, .13, .10, .12)
ctIdx <- sample.int(K, nCells, replace = TRUE, prob = prob)
spread <- c(520, 430, 300, 340, 380, 400)[ctIdx]
x <- rnorm(nCells, centres$x[ctIdx], spread)
y <- rnorm(nCells, centres$y[ctIdx], spread)

## --- expression: per-type marker block over a shared background
genes <- sprintf("GENE%03d", seq_len(nGenes))
base <- rgamma(nGenes, shape = 1.1, rate = 3.5)
prog <- matrix(base, nGenes, K, dimnames = list(genes, typeNames))
for (k in seq_len(K)) {
  mk <- ((k - 1) * nMark + 1):min(k * nMark, nGenes)
  prog[mk, k] <- prog[mk, k] + rgamma(length(mk), shape = 7, rate = 1.4)
}
depth <- rnbinom(nCells, mu = 190, size = 4) + 12   # Xenium-like counts/cell
lambda <- prog[, ctIdx, drop = FALSE]
lambda <- sweep(lambda, 2, Matrix::colSums(lambda), "/")
lambda <- sweep(lambda, 2, depth, "*")
gex <- Matrix(matrix(rpois(length(lambda), as.vector(lambda)), nGenes, nCells,
                     dimnames = list(genes, NULL)), sparse = TRUE)

## --- negative controls INCLUDING Genomic Control (21 probes, as on Prime 5K).
## The QC helper must count this assay; omitting it under-reports the
## negative-control rate on exactly the panel where controls matter most.
mkCtrl <- function(n, pref, mu) {
  Matrix(matrix(rpois(n * nCells, mu), n, nCells,
                dimnames = list(sprintf("%s%02d", pref, seq_len(n)), NULL)), sparse = TRUE)
}
ctrlProbe <- mkCtrl(20, "NegControlProbe_", 0.012)
ctrlCode  <- mkCtrl(30, "NegControlCodeword_", 0.008)
genomic   <- mkCtrl(21, "GenomicControl_", 0.010)
unassign  <- mkCtrl(15, "UnassignedCodeword_", 0.006)

full <- rbind(gex, ctrlProbe, ctrlCode, genomic, unassign)
ftype <- c(rep("Gene Expression", nrow(gex)),
           rep("Negative Control Probe", nrow(ctrlProbe)),
           rep("Negative Control Codeword", nrow(ctrlCode)),
           rep("Genomic Control", nrow(genomic)),
           rep("Unassigned Codeword", nrow(unassign)))
## Xenium-shaped cell ids ("aaaddlda-1") so the Explorer-export prefix check is real
bc <- make.unique(sprintf("%s-1", replicate(nCells, paste(sample(letters, 8, TRUE), collapse = ""))),
                  sep = "x")
colnames(full) <- bc

writeMM(full, file.path(OUT, "cell_feature_matrix", "matrix.mtx"))
system2("gzip", c("-f", shQuote(file.path(OUT, "cell_feature_matrix", "matrix.mtx"))))
write.table(data.frame(id = rownames(full), name = rownames(full), type = ftype),
            gzfile(file.path(OUT, "cell_feature_matrix", "features.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bc, gzfile(file.path(OUT, "cell_feature_matrix", "barcodes.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## --- cells.csv.gz with morphology, so the QC and segmentation-health panels populate
cellArea <- rlnorm(nCells, log(c(140, 110, 45, 95, 80, 130)[ctIdx]), 0.42)
nucFrac <- rbeta(nCells, 6, 8)                 # median ~0.43, stain-segmented-like
nucArea <- cellArea * nucFrac
nucCount <- sample(c(0, 1, 2, 3), nCells, TRUE, prob = c(.05, .86, .07, .02))
nucArea[nucCount == 0] <- NA                   # anucleate -> NA, as real runs do
cells <- data.table(
  cell_id = bc, x_centroid = x, y_centroid = y,
  transcript_counts = Matrix::colSums(gex),
  control_probe_counts = Matrix::colSums(ctrlProbe),
  genomic_control_counts = Matrix::colSums(genomic),
  control_codeword_counts = Matrix::colSums(ctrlCode),
  unassigned_codeword_counts = Matrix::colSums(unassign),
  deprecated_codeword_counts = 0L,
  total_counts = Matrix::colSums(full),
  cell_area = cellArea, nucleus_area = nucArea, nucleus_count = nucCount,
  segmentation_method = sample(c("segmented by boundary stain",
                                 "segmented by interior stain",
                                 "segmented by nucleus expansion"),
                               nCells, TRUE, prob = c(.55, .30, .15))
)
data.table::fwrite(cells, file.path(OUT, "cells.csv.gz"))

cat(sprintf("synthetic Xenium bundle -> %s\n  %d cells, %d features (%d panel genes), %d planted cell types\n",
            OUT, nCells, nrow(full), nrow(gex), K))
cat(sprintf("  median counts/cell %g | median nucleus:cell ratio %.2f\n",
            stats::median(cells$transcript_counts), stats::median(nucFrac)))
