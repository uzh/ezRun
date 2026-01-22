###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpaceRangerSegQC <-
  setRefClass(
    "EzAppSpaceRangerSeqQC",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpaceRangerSegQC
        name <<- "EzAppSpaceRangerSegQC"
      }
    )
  )

ezMethodSpaceRangerSegQC <- function(input = NA, output = NA, param = NA) {
  sampleName <- input$getNames()

  #Fix image because of a bug in spaceranger 4.0.1
  myImage <- input$getFullPaths("Image")
  imageSize <- file.size(myImage) / 1024^3
  if (imageSize < 4) {
    cmd_tiffSplit <- paste(
      '/usr/local/ngseq/src/tiff-4.7.0/bin/bin/tiffsplit',
      myImage,
      'output_'
    )
    ezSystem(cmd_tiffSplit)
    highResName <- sub('.tif$', '_highRes.tif', basename(myImage))
    highresImage <- ezSystem('ls -S output_*.tif | head -n 1', intern = TRUE)
    ezSystem(paste('mv', highresImage, highResName))
  } else {
    highResName <- myImage
  }

  cmd <- paste(
    "spaceranger segment",
    paste0("--id=", sampleName),
    paste0("--localmem=", param$ram),
    paste0("--localcores=", param$cores)
  )

  cmd <- paste(cmd, paste0("--tissue-image=", highResName))

  if (ezIsSpecified(param$cmdOptions)) {
    cmd <- paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  ezSystem(paste("mv", file.path(sampleName, "outs"), "result"))
  unlink(sampleName, recursive = TRUE)
  file.rename('result', sampleName)
  ezSystem(paste(
    "mv",
    file.path(sampleName, "websummary.html"),
    file.path(sampleName, "web_summary.html")
  ))

  return("Success")
}
