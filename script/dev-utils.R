
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







