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
