fakeFtpPath <- "ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/001/SRR0000001/SRR0000001.fastq.gz"

test_that(".enaDownloadFastq reuses an existing file with matching md5", {
  withr::local_dir(withr::local_tempdir())
  writeLines("complete", basename(fakeFtpPath))
  md5 <- unname(tools::md5sum(basename(fakeFtpPath)))
  local_mocked_bindings(ezSystem = function(...) stop("should not download"))
  expect_equal(.enaDownloadFastq(fakeFtpPath, md5), basename(fakeFtpPath))
})

test_that(".enaDownloadFastq replaces a partial file and retries failed downloads", {
  withr::local_dir(withr::local_tempdir())
  writeLines("partial", basename(fakeFtpPath))
  nCalls <- 0
  local_mocked_bindings(
    ezSystem = function(cmd, ...) {
      nCalls <<- nCalls + 1
      if (nCalls == 1) {
        stop("No such directory")
      }
      writeLines("complete", basename(fakeFtpPath))
    }
  )
  writeLines("complete", "ref")
  md5 <- unname(tools::md5sum("ref"))
  expect_equal(
    .enaDownloadFastq(fakeFtpPath, md5, retryWait = 0),
    basename(fakeFtpPath)
  )
  expect_equal(nCalls, 2)
})

test_that(".enaDownloadFastq gives up after the given number of tries", {
  withr::local_dir(withr::local_tempdir())
  local_mocked_bindings(ezSystem = function(...) stop("No such directory"))
  expect_error(
    .enaDownloadFastq(fakeFtpPath, NA, tries = 2, retryWait = 0),
    "Could not download"
  )
})

test_that(".enaSushiRegisterCmd separates the ssh host from the remote command", {
  cmd <- .enaSushiRegisterCmd(
    "/srv/GT/analysis/datasets/ENA_App/dataset_42.tsv", "42", "ENA_App_PRJ1", "777"
  )
  expect_match(cmd, "^ssh -o BatchMode=yes trxcopy@fgcz-h-082 \"bash -lc '", fixed = FALSE)
  expect_match(cmd, "--dataset /srv/GT/analysis/datasets/ENA_App/dataset_42.tsv", fixed = TRUE)
  expect_match(cmd, "--project 42 --dataset_name ENA_App_PRJ1 --dataset_id 777'\"$")
})

test_that(".enaNormaliseFastqFiles keeps the mates of a run with an unpaired file", {
  ftp <- paste0(
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz;",
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_1.fastq.gz;",
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_2.fastq.gz"
  )
  res <- .enaNormaliseFastqFiles(ftp, "aaa;bbb;ccc")
  expect_equal(
    res$ftp,
    paste0(
      "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_1.fastq.gz;",
      "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1_2.fastq.gz"
    )
  )
  expect_equal(res$md5, "bbb;ccc")
})

test_that(".enaNormaliseFastqFiles leaves 1- and 2-file runs alone and rejects unknown triples", {
  expect_equal(.enaNormaliseFastqFiles("a/x.fastq.gz", "m1")$ftp, "a/x.fastq.gz")
  two <- .enaNormaliseFastqFiles("a/x_1.fastq.gz;a/x_2.fastq.gz", "m1;m2")
  expect_equal(two$md5, "m1;m2")
  expect_error(
    .enaNormaliseFastqFiles("a/x.fq.gz;a/y.fq.gz;a/z.fq.gz", "1;2;3"),
    "Cannot identify"
  )
})
