## methods_facts(): fixed behaviour of an app's code that its parameter form does
## not show (seed, clustering algorithm, multiple-testing correction, ...). The
## Methods LLM only reads the job script and logs, so without this it has to write
## [not recorded] for these or guess them. methods_description() hands the facts to
## llm_write_methods as one more evidence file.

## Stand-in for the AI/llm_methods_caller binary: records its arguments and writes
## a minimal methods.md, so the plumbing is tested without calling the LLM.
withStubWriter <- function(code) {
  bin <- tempfile("bin"); dir.create(bin)
  argsFile <- file.path(bin, "args.txt")
  writeLines(c(
    "#!/bin/sh",
    paste0("printf '%s\\n' \"$@\" > ", shQuote(argsFile)),
    "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
    "printf 'Methods text.\\n' > \"$out\""
  ), file.path(bin, "llm_write_methods"))
  Sys.chmod(file.path(bin, "llm_write_methods"), "755")
  oldPath <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, oldPath, sep = ":"))
  on.exit(Sys.setenv(PATH = oldPath))
  code(argsFile)
}

test_that("the base app has no facts and every drafted app has some", {
  expect_identical(EzApp$new()$methods_facts(), character(0))
  for (cls in list(EzAppScSeurat, EzAppDeseq2)) {
    facts <- cls$new()$methods_facts()
    expect_type(facts, "character")
    expect_gt(length(facts), 0)
    expect_false(any(is.na(facts) | !nzchar(facts)))
  }
})

test_that("methods_description passes app_facts.txt only when there are facts", {
  withStubWriter(function(argsFile) {
    script <- tempfile(fileext = ".sh"); writeLines("echo job", script)

    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$methods_description(script, character(0), 1, out)
    args <- readLines(argsFile)
    factsFile <- file.path(out, "app_facts.txt")
    expect_true(file.exists(factsFile))
    scriptsAt <- which(args == "--scripts")
    expect_true(factsFile %in% args[-seq_len(scriptsAt)])
    expect_true(all(EzAppScSeurat$new()$methods_facts() %in% readLines(factsFile)))
    expect_match(readLines(factsFile)[1], "EzAppScSeurat")

    out2 <- tempfile("out"); dir.create(out2)
    EzApp$new()$methods_description(script, character(0), 1, out2)
    expect_false(file.exists(file.path(out2, "app_facts.txt")))
    expect_false(any(grepl("app_facts.txt", readLines(argsFile), fixed = TRUE)))
  })
})
