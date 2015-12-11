

# Make additionally to the GO stuff some links to the enrichr web page.
# - on separate html pages in a sub-directory enrichr
# - only for human and mouse
# - send the gene_names to enrichr
# - http://amp.pharm.mssm.edu/Enrichr/#help

genes = c("Nsun3","Polrmt","Nlrx1","Sfxn5","Zc3h12c","Slc25a39","Arsg","Defb29","Ndufb6","Zfand1")


if (grepl("musculus", param$refBuild) | grepl("sapiens", param$refBuild)){
  setwdNew("./enrichr")
  on.exit(setwd(".."))
  
  genesList = paste(genes, collapse="\n")
  jsCall = paste0("enrich({list: '", genesList, "', popup: true});")
  jsFile = system.file("extdata/enrichrFoo.js", package="ezRun", mustWork=TRUE)
  jsFunction = paste(scan(javaFile, what = "character", sep = "\n", quiet = TRUE), collapse = "\n")
  wholeJS = paste(jsFunction, "\n", jsCall, collapse=" ")
  
  enrichrDoc = openBsdocReport()
  addJavascript(enrichrDoc, text=wholeJS)
  closeBsdocReport(enrichrDoc, "enrichr.html")
}

# maybe with cmd line
cmd = paste0('javac ', enrichJavaFile, ' "enrich({list: ', genesList, ', popup: true});"')
ezSystem(cmd)


ezSystem(paste0('java -jar js.jar ', enrichJavaFile))





