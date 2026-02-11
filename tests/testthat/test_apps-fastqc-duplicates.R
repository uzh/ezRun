context("Tests for ezMethodFastQC duplicate detection and resolution")

test_that("Duplicate detection logic works correctly", {
  # Test 1: No duplicates should not trigger any special handling
  files_no_dup <- c(
    "sample1_R1" = "path/to/sample1.fastq.gz",
    "sample2_R1" = "path/to/sample2.fastq.gz"
  )
  reportDirs_no_dup <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files_no_dup))
  expect_false(any(duplicated(reportDirs_no_dup)))
  
  # Test 2: Duplicates that can be resolved with sample names
  files_dup_resolvable <- c(
    "sample1_R1" = "path1/reads.fastq.gz",
    "sample2_R1" = "path2/reads.fastq.gz"
  )
  reportDirs_dup <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files_dup_resolvable))
  expect_true(any(duplicated(reportDirs_dup)))
  
  # Check if renaming to sample names resolves it
  reportDirs_renamed <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", names(files_dup_resolvable))
  expect_false(any(duplicated(reportDirs_renamed)))
  
  # Test 3: Duplicates that cannot be resolved (both basename and sample names duplicate)
  files_dup_unresolvable <- c(
    "sample1_R1" = "path1/reads.fastq.gz",
    "sample1_R1" = "path2/reads.fastq.gz"
  )
  reportDirs_dup_unres <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files_dup_unresolvable))
  expect_true(any(duplicated(reportDirs_dup_unres)))
  
  reportDirs_renamed_unres <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", names(files_dup_unresolvable))
  expect_true(any(duplicated(reportDirs_renamed_unres)))
})

test_that("File extension handling works correctly", {
  # Test various file extensions
  files <- c(
    "sample1" = "path/to/file.fastq.gz",
    "sample2" = "path/to/file.fastq",
    "sample3" = "path/to/file.fq.gz",
    "sample4" = "path/to/file.fq",
    "sample5" = "path/to/file.bam"
  )
  
  reportDirs <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files))
  
  expect_equal(reportDirs["sample1"], "file_fastqc")
  expect_equal(reportDirs["sample2"], "file_fastqc")
  expect_equal(reportDirs["sample3"], "file_fastqc")
  expect_equal(reportDirs["sample4"], "file_fastqc")
  expect_equal(reportDirs["sample5"], "file_fastqc")
})
