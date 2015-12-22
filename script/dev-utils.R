
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
rm2 = "inst/extdata/genesWithPrespliced.gtf"
cmd = paste("rm -fr", rm0, rm1, rm2)
ezSystem(cmd)

# create bed dummy for teqc
randStart = round(runif(500), digits=3)*1000
randEnd = randStart + round(runif(500), digits=2)*100+50
dat = data.frame(rep("I",500), randStart*(1:500), randEnd*(1:500))
ezWrite.table(dat, file="./inst/extdata/genes.bed", head="track dummy", row.names = FALSE, col.names = FALSE)


