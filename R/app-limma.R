###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## for now only gene set testing from limma

twoGrouplsRoast = function(param, testResult, seqAnno){
  job = ezJobStart("twoGroupsRoast")
  library(GOstats, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  library(annotate, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  if (param$featureLevel != "gene"){
    stop("only feature level: gene is supported")
  }
  ontologies = c("BP", "MF", "CC")
  #goResults = list()
  #for (onto in ontologies){
  goResults = ezMclapply(ontologies, function(onto){
    gene2goList = goStringsToList(seqAnno[[paste("GO", onto)]], listNames=rownames(seqAnno))[presentGenes]
    if (param$includeGoParentAnnotation){
    }
  })
}
