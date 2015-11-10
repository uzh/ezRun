

require(ezRun)
require(ShortRead)
ds = my.read.table("inst/extdata/ventricles_100k/dataset.tsv")

ds$"Read1 [File]" = sub("p1001/ILLUMINA/ventricles_100k", "inst/extdata/ventricles_10k", ds$"Read1 [File]")
ds$"Read2 [File]" = sub("p1001/ILLUMINA/ventricles_100k", "inst/extdata/ventricles_10k", ds$"Read2 [File]")
for (fq in c(ds$"Read1 [File]", ds$"Read2 [File]")){ 
  reads = readFastq(fq)
  #idx = sample(1:length(reads), 100000, replace=FALSE)
  writeFastq(reads[1:10000], file=sub(".gz$", "", fq), compress=FALSE)
}
my.write.table(ds, file="inst/extdata/ventricles_10k/dataset.txt", head="Name")
