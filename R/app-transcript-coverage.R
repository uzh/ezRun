###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


### computes per transcript read counts and per-transcript-base coverage
ezMethodTranscriptCoverage = function(input=NA, output=NA, param=NA){
  
  require(GenomicFeatures)
  
  seqAnno = ezRead.table(param$ezRef@refAnnotationFile)
  bamFile = input$getFullPaths("trBAM")
  strand = switch(param$strandMode,
                  sense="+",
                  antisense="-",
                  both="*")
  countVec = countTranscriptBam(bamFile, strand=strand)
  stopifnot(names(countVec) %in% rownames(seqAnno))
  counts = ezMatrix(countVec, rows=names(countVec), cols="matchCounts")
  ezWrite.table(counts, head="transcript_id", file=basename(output$getColumn("Count")))
  
  stopifnot(param$paired == FALSE)

  ### compute the coverage
  txdb = makeTxDbFromGFF(param$ezRef@refFeatureFile,
                         dataSource="FGCZ", taxonomyId = 2759)# organism=organism, chrominfo=NULL)
  tisPos = getTranslationInitiationSite(txdb)
  
  profileList = getTranscriptProfiles(bamFile, param)
  save(profileList, file=basename(output$getColumn("Profiles")))
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodTranscriptCoverage(input=NA, output=NA, param=NA)
##' @description Process alignments to transcripts in trBAM files. See also Bowtie2TranscriptomeApp.
EzAppTranscriptCoverage <-
  setRefClass("EzAppTranscriptCoverage",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTranscriptCoverage
                  name <<- "EzAppTranscriptCoverage"
                  appDefaults <<- rbind(minReadLength=ezFrame(Type="integer", DefaultValue="28", Description="minimum read length to consider"),
                                        maxReadLength=ezFrame(Type="integer", DefaultValue="32", Description="maximum read length to consider"),
                                        getCoverageByReadLength=ezFrame(Type="logical", DefaultValue="TRUE", Description="should the coverage also computed for each read length separatedly"),
                                        readStartOnlyCoverage=ezFrame(Type="logical", DefaultValue="TRUE", Description="should the coverage profile show the read start or the entire read length"))
                }
              )
  )
