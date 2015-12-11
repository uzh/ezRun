

# Make additionally to the GO stuff some links to the enrichr web page.
# - on separate html pages in a sub-directory enrichr
# - only for human and mouse
# - send the gene_names to enrichr
# - http://amp.pharm.mssm.edu/Enrichr/#help

genes = c("Nsun3","Polrmt","Nlrx1","Sfxn5","Zc3h12c","Slc25a39","Arsg","Defb29","Ndufb6","Zfand1")


if (grepl("musculus", param$refBuild) | grepl("sapiens", param$refBuild)){
  setwdNew("./enrichr")
  enrichrDoc = openBsdocReport()
  javaFile = system.file("extdata/enrichrFoo.js", package="ezRun", mustWork=TRUE)
  addJavascript(enrichrDoc, javaFile)
#   javaObject = rJava::.jnew("java/awt/Frame", "Hello")
#   rJava::.jcall(javaObject, "V", "addJavascript", javaFile)
  genesList = paste(genes, collapse="\n")
  .enrichrLink = function(){
    ezLegend(title="Enrichr link")
    add.plot.interactivity(text, labels="click me", x=1, y=0.8, click.actions=paste0("enrich('{list: ", genesList, ", popup: true}');"))
  }
  addPlot(enrichrDoc, .enrichrLink, fontname="serif", par.properties=parLeft())
  enrichrLink = ezImageFileLink(plotCmd, file="enrichrLink.png", width=70, height=40, addPdfLink=FALSE)
  addParagraph(enrichrDoc, enrichrLink)
  closeBsdocReport(enrichrDoc, "enrichr.html")
  setwd("..")
}

# maybe with cmd line
cmd = paste0('javac ', enrichJavaFile, ' "enrich({list: ', genesList, ', popup: true});"')
ezSystem(cmd)


ezSystem(paste0('java -jar js.jar ', enrichJavaFile))





