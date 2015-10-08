###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


getReferenceFeaturesBed = function(param){
  bedFile = sub(".gtf$", ".bed", param$ezRef["refFeatureFile"])
  if (!file.exists(bedFile)){
    ezSystem(paste(GTF2BED, param$ezRef["refFeatureFile"], ">", bedFile))
    ezSystem(paste("chmod", "g+w", bedFile))
  }
  return(bedFile)
}

getRefChromSizesFile = function(param){
  if (!file.exists(param$ezRef@refChromSizesFile)){
    faFiles = list.files(path = param$ezRef@refChromDir, pattern = ".fa$", full.names = TRUE)
    chromSizes = character()
    for (ff in faFiles){
      message(ff)
      dss = readDNAStringSet(ff)
      chromSizes[names(dss)] = width(dss)
    }
    ezWrite.table(chromSizes, file=param$ezRef@refChromSizesFile, col.names = FALSE)    
  }
  return(param$ezRef@refChromSizesFile)
}

