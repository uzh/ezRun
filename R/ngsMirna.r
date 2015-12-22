###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


.getLeftClipping = function(cigar, batchSize=100000){
  
  hRaw = as.raw(72)
  result = integer(length(cigar))
  for (i in seq(from=1, to=length(cigar), by=batchSize)){
    j = min(i+batchSize-1, length(cigar))
    cigarList = splitCigar(cigar[i:j]) ## TODOMF: function splitCigar doesn't exist
    gc()
    result[i:j] = sapply(cigarList, function(x){if (x[[1]][1] == hRaw) x[[2]][1] else 0})
    gc()
  }
  return(result)
}

## TODO: refactor: use the general ezMethodTrim
trimMirna = function(input, output, adapter=NA, param){
  qtOut = sub(".fastq.gz", "_qualtrim.fastq", basename(input))
  qtBad = sub(".fastq.gz", "_qualtrimBad.fastq", basename(input))
  catCmd = paste("gunzip -c", input)
  filtCmd = paste("|", PRINSEQ_LITE, 
                  "-no_qual_header",
                  "-trim_qual_right", 20,
                  "-trim_qual_type", "mean",
                  "-trim_qual_window", 4,
                  "-fastq", "stdin",
                  "-out_bad", sub(".fastq$", "", qtBad),
                  "-out_good", sub(".fastq$", "", qtOut), 
                  ">", sub(".fastq", ".prinseq.out", output))
  cmd = paste(catCmd, filtCmd)
  ezSystem(cmd)
  atOut = sub("_qualtrim.fastq", "_allTrimmed.fastq", qtOut)
  stopifnot(atOut != qtOut)
  cmd = paste("/usr/local/ngseq/bin/flexbar",
              "--adapter-seq", adapter,
              "--adapter-trim-end", "RIGHT",
              "--adapter-min-overlap", 10,
              "--adapter-threshold", 1.5,
              "--min-read-length", 18,
              "--max-uncalled", 4,
              "--format", "i1.8",
              "--reads", qtOut,
              "--target",  sub(".fastq$", "", atOut),
              ">", sub(".fastq", ".flexbar.out", output))
  ezSystem(cmd)
  ezSystem(paste("rm -f ", qtOut, qtBad))
  if (atOut != output){
    ezSystem(paste("mv", atOut, output))
  } else {
    ezSystem(paste("rm", atOut))    
  }
  return(output)
}
