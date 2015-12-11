

# Make additionally to the GO stuff some links to the enrichr web page.
# - on separate html pages in a sub-directory enrichr
# - only for human and mouse
# - send the gene_names to enrichr
# - http://amp.pharm.mssm.edu/Enrichr/#help

genes = c("Nsun3","Polrmt","Nlrx1","Sfxn5","Zc3h12c","Slc25a39","Arsg","Defb29","Ndufb6","Zfand1")


if (grepl("musculus", param$refBuild) | grepl("sapiens", param$refBuild)){
  genesList = paste(genes, collapse="\\n")
  jsCall = paste0('enrich({list: "', genesList, '", popup: true});')
  jsFile = system.file("extdata/enrichr.js", package="ezRun", mustWork=TRUE)
  jsFunction = paste(scan(jsFile, what = "character", sep = "\n", quiet = TRUE))
  wholeJS = c(jsFunction, jsCall)
  
  doc = openBsdocReport()
#   setwdNew("./enrichr")
#   enrichrDoc = openBsdocReport()
#   enrichrDoc = addJavascript(enrichrDoc, jsFile)
#   addParagraph(enrichrDoc, as.html(pot("Go back", hyperlink="../testEnrichrLink.html")))
#   closeBsdocReport(enrichrDoc, "enrichr.html")
#   setwd("..")
  addJavascript(doc, jsFile)
  addParagraph(doc, as.html(pot("Enrichr", hyperlink="enrichr/enrichr.html")))
  addParagraph(doc, pot(paste0("<a href='javascript:void(0)' onClick='", jsCall, "'>blabla</a>")))
  closeBsdocReport(doc, "testEnrichrLink.html")
}


