context("splitPipe: Parse DGE -> 10x conversion")

## convertParseToTenx() turns split-pipe's cells x genes count_matrix.mtx into a
## 10x genes x cells directory. It used to load the matrix with Matrix::readMM,
## which cannot read a matrix with more than 2^31 - 1 non-zeros: the o43379
## all-sample matrix (1.36M cells, 4.2e9 non-zeros) failed with
## "scan() expected 'an integer', got '4198638981'" after an 11 h combine.

makeParseDge <- function(dir, header, entries, cells, genes) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "%%MatrixMarket matrix coordinate integer general                                ",
      sprintf("%%Rows=cells (%s), Cols=genes (%s)", header[1], header[2]),
      paste(header, collapse = " "),
      entries
    ),
    file.path(dir, "count_matrix.mtx")
  )
  write.csv(
    data.frame(gene_id = genes, gene_name = paste0("sym", seq_along(genes)), genome = "g"),
    file.path(dir, "all_genes.csv"),
    row.names = FALSE, quote = FALSE
  )
  write.csv(
    data.frame(bc_wells = cells, sample = "A1"),
    file.path(dir, "cell_metadata.csv"),
    row.names = FALSE, quote = FALSE
  )
}

test_that("the 10x matrix is the transpose, with barcodes and features", {
  d <- tempfile("dge")
  cells <- c("c1", "c2", "c3")
  genes <- c("G1", "G2", "G3", "G4")
  ## cells x genes: (1,2)=5, (2,4)=1, (3,1)=2
  makeParseDge(d, c(3, 4, 3), c("1 2 5", "2 4 1", "3 1 2"), cells, genes)
  out <- tempfile("tenx")
  convertParseToTenx(d, out)

  m <- Matrix::readMM(gzfile(file.path(out, "matrix.mtx.gz")))
  expect_equal(dim(m), c(4L, 3L))
  expect_equal(m[2, 1], 5)
  expect_equal(m[4, 2], 1)
  expect_equal(m[1, 3], 2)
  expect_equal(sum(m), 8)
  expect_equal(readLines(gzfile(file.path(out, "barcodes.tsv.gz"))), cells)
  feat <- read.delim(gzfile(file.path(out, "features.tsv.gz")), header = FALSE)
  expect_equal(feat$V1, genes)
  expect_equal(feat$V2, paste0("sym", 1:4))
  expect_true(all(feat$V3 == "Gene Expression"))
})

test_that("a matrix past 2^31 non-zeros converts (header is not parsed as integer)", {
  d <- tempfile("dge")
  ## the real o43379 all-sample header, with one entry: only the header is at stake
  makeParseDge(d, c(1363240, 65965, 4198638981), "1 1 1", c("c1"), c("G1"))
  out <- tempfile("tenx")
  expect_error(convertParseToTenx(d, out), NA)
  hdr <- readLines(gzfile(file.path(out, "matrix.mtx.gz")), n = 2)
  expect_match(hdr[1], "^%%MatrixMarket matrix coordinate integer general")
  expect_equal(trimws(hdr[2]), "65965 1363240 4198638981")
})
