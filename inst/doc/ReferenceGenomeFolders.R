## ---- warning=FALSE,message=FALSE,split=TRUE-----------------------------
library(ezRun)

## ----split=TRUE,eval=TRUE------------------------------------------------
refBuild = "Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18"
genomesRoot = "./refExample"
system("rm -rf refExample")
param = ezParam(list(refBuild=refBuild, genomesRoot=genomesRoot))
gtf = system.file("extdata/genes.gtf", package="ezRun", mustWork = TRUE)
fp = system.file("extdata/genome.fa", package="ezRun", mustWork = TRUE)
buildRefDir(param$ezRef, fp, gtf)
buildIgvGenome(param$ezRef)
seqAnno = writeAnnotationFromGtf(param=param)# featureFile=param$ezRef["refFeatureFile"], featAnnoFile=myRef["refAnnotationFile"])

