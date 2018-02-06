###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodHomerDiffPeaks = function(input=NA, output=NA, param=NA, 
                                  htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  
  stopifnot(param$sampleGroup != param$refGroup)
  useSamples = input$getNames()[input$getColumn(param$grouping) %in% 
                                  c(param$sampleGroup, param$refGroup)]
  input <- input$subset(useSamples)
  
  bamFiles <- input$getFullPaths("BAM")
  localBamFiles <- sapply(bamFiles, getBamLocally)
  
  mcmapply(makeTagDirectory, inBam=localBamFiles, 
           outputDir=names(localBamFiles), mc.cores=param$cores)
  
  if(all(localBamFiles != bamFiles))
    file.remove(c(localBamFiles, paste0(localBamFiles, ".bai")))
  
  
  
}


makeTagDirectory <- function(inBam, outputDir, isAntisense=FALSE,
                             strandedPaired=FALSE){
  setEnvironments("HOMER")
  cmd <- paste("makeTagDirectory", outputDir, paste(inBam, collapse=" "),
               "-format sam")
  if(isTRUE(isAntisense))
    cmd <- paste(cmd, "-flip")
  if(isTRUE(strandedPaired))
    cmd <- paste(cmd, "-sspe")
  ezSystem(cmd)
}
