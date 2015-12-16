
# delete entire "/scratch" directory
ezSystem("rm -fr /scratch/*")

# save() save.image()

# open report files by providing character(s) to grepl the path for (case is ignored)
editReportFile = function(patterns=""){
  reportFiles = paste0("/scratch/", list.files("/scratch", "00index", recursive = TRUE))
  select = ezMatrix(FALSE, rows = patterns, cols = 1:length(reportFiles))
  for (i in 1:length(patterns)){
    select[i, ] = grepl(patterns[i], reportFiles, ignore.case = TRUE)
  }
  use = apply(select, 2, any)
  file.edit(reportFiles[use])
}

# open all report files
editReportFile()

# specific
patterns = c("deseq", "edger")
editReportFile(patterns)

# delete left-overs from devtools::run_examples()
rm0 = "tests/testthat/run_examples/*"
rm1 = "DESCRIPTION_head"
rm2 = "inst/extdata/genes.bed"
rm3 = "inst/extdata/genesWithPrespliced.gtf"
cmd = paste("rm -fr", rm0, rm1, rm2, rm3)
ezSystem(cmd)





