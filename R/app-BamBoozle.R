###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
###################################################################

##' @template app-template
##' @templateVar method ezMethodBamBoozle(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppBamBoozle <-
  setRefClass(
    "EzAppBamBoozle",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBamBoozle
        name <<- "EzAppBamBoozle"
      }
    )
  )

ezMethodBamBoozle = function(input = NA, output = NA, param = NA) {
  bamFile = input$getFullPaths("BAM")
  file.copy(bamFile, basename(bamFile))
  genomeFa = param$ezRef["refFastaFile"]
  outFile <- sub('.bam', '_cleaned.bam', basename(bamFile))

  cmd <- paste(
    'BAMboozle --bam',
    basename(bamFile),
    '--fa',
    genomeFa,
    '--out',
    outFile,
    '--p',
    param$cores,
    param$cmdOptions
  )
  ezSystem(cmd)

  fq1 <- sub('.bam', '_R1.fastq', outFile)
  if (param$paired) {
    system(paste(
      'samtools sort -n',
      outFile,
      '-@',
      param$cores,
      '-o sortedBam.bam'
    ))
    fq2 <- sub('.bam', '_R2.fastq', outFile)
    cmd <- paste(
      'bedtools bamtofastq -i sortedBam.bam',
      '-fq',
      fq1,
      '-fq2',
      fq2
    )
  } else {
    fq2 <- ''
    cmd <- paste('bedtools bamtofastq -i', outFile, '-fq', fq1)
  }
  system(cmd)
  system(paste('pigz --best -p', param$cores, fq1, fq2))
  return('suceess')
}
