

# Make additionally to the GO stuff some links to the enrichr web page.
# - on separate html pages in a sub-directory enrichr
# - only for human and mouse
# - send the gene_names to enrichr
# - http://amp.pharm.mssm.edu/Enrichr/#help

genes = c("Nsun3","Polrmt","Nlrx1","Sfxn5","Zc3h12c","Slc25a39","Arsg","Defb29","Ndufb6","Zfand1")


if (grepl("musculus", param$refBuild) | grepl("sapiens", param$refBuild)){
  setwdNew("./enrichr")
  on.exit(setwd(".."))
  
  genesList = paste(genes, collapse="\\n")
  jsCall = paste0("enrich({list: '", genesList, "', popup: true});")
  jsFile = system.file("extdata/enrichrFoo.js", package="ezRun", mustWork=TRUE)
  jsFunction = paste(scan(jsFile, what = "character", sep = "\n", quiet = TRUE))
  wholeJS = c(jsFunction, jsCall)
  
  enrichrDoc = openBsdocReport()
  .enrichrLink = function(document, js){
    ezLegend(title="Enrichr link")
    add.plot.interactivity(text, labels="click me", x=1, y=0.8, click.actions=addJavascript(document, text=js))
  }
  # addJavascript(enrichrDoc, text=wholeJS)
  addPlot(enrichrDoc, .enrichrLink, fontname="serif", par.properties=parLeft(), document=enrichrDoc, js=wholeJS)
  closeBsdocReport(enrichrDoc, "enrichr.html")
}





# maybe with cmd line
cmd = paste0('javac ', enrichJavaFile, ' "enrich({list: ', genesList, ', popup: true});"')
ezSystem(cmd)


ezSystem(paste0('java -jar js.jar ', enrichJavaFile))





