context("NfCoreAtacSeq sample sheet and QC-mode subsampling")

writeFakeFastq <- function(file, ids) {
  con <- gzfile(file, "w")
  writeLines(
    as.vector(rbind(paste0("@", ids), "ACGT", "+", "IIII")),
    con
  )
  close(con)
  file
}

makeAtacDataset <- function(condition, readCount = NULL, lanes = 1) {
  root <- tempfile("atacds")
  dir.create(root)
  names <- paste0("S", seq_along(condition))
  r1 <- r2 <- character(length(names))
  for (i in seq_along(names)) {
    f1 <- f2 <- character(lanes)
    for (l in seq_len(lanes)) {
      ids <- paste0(names[i], "_L", l, "_read", 1:10)
      f1[l] <- basename(writeFakeFastq(
        file.path(root, paste0(names[i], "_L", l, "_R1.fastq.gz")),
        ids
      ))
      f2[l] <- basename(writeFakeFastq(
        file.path(root, paste0(names[i], "_L", l, "_R2.fastq.gz")),
        ids
      ))
    }
    r1[i] <- paste(f1, collapse = ",")
    r2[i] <- paste(f2, collapse = ",")
  }
  meta <- data.frame(
    "Condition [Factor]" = condition,
    "Read1 [File]" = r1,
    "Read2 [File]" = r2,
    check.names = FALSE,
    row.names = names
  )
  if (!is.null(readCount)) {
    meta[["Read Count"]] <- readCount
  }
  EzDataset$new(meta = meta, dataRoot = root)
}

atacParam <- list(grouping = "Condition", paired = TRUE, cores = 1)

test_that("missing conditions fall back to the sample name", {
  sheet <- getAtacSampleSheet(makeAtacDataset(c(NA, NA, NA)), atacParam)
  expect_equal(sheet$sample, c("S1", "S2", "S3"))
  expect_equal(sheet$replicate, c(1, 1, 1))
  expect_equal(sheet$sid, c("S1", "S2", "S3"))

  sheet <- getAtacSampleSheet(makeAtacDataset(c("", "NA", "")), atacParam)
  expect_equal(sheet$sample, c("S1", "S2", "S3"))
})

test_that("partially filled conditions keep groups and replicates", {
  sheet <- getAtacSampleSheet(makeAtacDataset(c("WT", "", "WT")), atacParam)
  expect_equal(sheet$sample, c("WT", "S2", "WT"))
  expect_equal(sheet$replicate, c(1, 1, 2))
})

test_that("missing or unset grouping column is tolerated", {
  ds <- makeAtacDataset(c("A", "B"))
  sheet <- getAtacSampleSheet(ds, modifyList(atacParam, list(grouping = "")))
  expect_equal(sheet$sample, c("S1", "S2"))
  sheet <- getAtacSampleSheet(
    ds,
    modifyList(atacParam, list(grouping = "Genotype"))
  )
  expect_equal(sheet$sample, c("S1", "S2"))
})

test_that("group names are made nf-core safe", {
  sheet <- getAtacSampleSheet(makeAtacDataset(c("a+b", "a+b")), atacParam)
  expect_true(all(grepl("^\\S+$", sheet$sample)))
  expect_equal(sheet$replicate, c(1, 2))
})

test_that("multi-lane samples give one row per fastq file", {
  sheet <- getAtacSampleSheet(makeAtacDataset(c("A", "B"), lanes = 2), atacParam)
  expect_equal(nrow(sheet), 4)
  expect_equal(sheet$sid, c("S1", "S1", "S2", "S2"))
  expect_equal(sheet$replicate, c(1, 1, 1, 1))
})

readIds <- function(file) {
  x <- readLines(gzfile(file))
  x[seq(1, length(x), by = 4)]
}

test_that("QC mode keeps the first n read pairs", {
  withr::local_dir(withr::local_tempdir())
  ds <- makeAtacDataset(c("A", "B"), readCount = c(10, 10))
  param <- modifyList(atacParam, list(qcMode = TRUE, qcReadsPerSample = 4))
  sheet <- subsampleAtacFastqs(getAtacSampleSheet(ds, param), ds, param)
  expect_equal(nrow(sheet), 2)
  expect_equal(readIds(sheet$fastq_1[1]), paste0("@S1_L1_read", 1:4))
  expect_identical(readIds(sheet$fastq_1[1]), readIds(sheet$fastq_2[1]))
  expect_equal(sheet$sample, c("A", "B"))
})

test_that("QC mode uses small samples completely and concatenates lanes", {
  withr::local_dir(withr::local_tempdir())
  ds <- makeAtacDataset(c("A", "B"), readCount = c(20, 5), lanes = 2)
  param <- modifyList(atacParam, list(qcMode = TRUE, qcReadsPerSample = 15))
  fullSheet <- getAtacSampleSheet(ds, param)
  sheet <- subsampleAtacFastqs(fullSheet, ds, param)
  ## S1 is subsampled across both lanes into a single file
  s1 <- sheet[sheet$sid == "S1", ]
  expect_equal(nrow(s1), 1)
  expect_equal(
    readIds(s1$fastq_1),
    c(paste0("@S1_L1_read", 1:10), paste0("@S1_L2_read", 1:5))
  )
  ## S2 has fewer reads than requested and keeps its original files
  s2 <- sheet[sheet$sid == "S2", ]
  expect_equal(s2$fastq_1, fullSheet$fastq_1[fullSheet$sid == "S2"])
  expect_equal(sheet$sid, c("S1", "S2", "S2"))
})
