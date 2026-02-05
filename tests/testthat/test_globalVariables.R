context("Checks EZ_GLOBAL_VARIABLES for valid file paths.")

globalVarsFile = system.file(
  "extdata/EZ_GLOBAL_VARIABLES.txt",
  package = "ezRun",
  mustWork = TRUE
)
globalVars = read.table(globalVarsFile, sep = "=", row.names = 1)
varNames = rownames(globalVars)
globalVars = sapply(globalVars, trimWhiteSpace)

containsWhiteSpace = sapply(globalVars, grepl, pattern = " ")
isPythonPath = sapply(globalVars, grepl, pattern = "PYTHONPATH")
isJavaCall = sapply(globalVars, grepl, pattern = "java")
isInternetAddress = sapply(globalVars, grepl, pattern = "http")
hasSlash = sapply(globalVars, grepl, pattern = "/")

test_that("Test file paths", {
  globalPaths = globalVars[
    -which(containsWhiteSpace | isInternetAddress | !hasSlash),
  ]
  names(globalPaths) = varNames[
    !(containsWhiteSpace | isInternetAddress | !hasSlash)
  ]
  exception = globalPaths["REFSEQ_mRNA_REF"]
  expect_true(file.exists(sub("/transcriptome", "", exception)))
  globalPaths = globalPaths[-which(globalPaths == exception)]
  fileExists = sapply(globalPaths, file.exists)
  if (any(!fileExists)) {
    print("Some paths do not exist:")
    print(globalPaths[which(!fileExists)])
  }
  #expect_true(sum(fileExists) == length(fileExists))
})
