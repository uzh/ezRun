


## if the outFile is not specified, the inFile will be overwritten
trimFastq <- function(inFile, outFile=NULL, length=NULL){
  stopifnot(!is.null(length))
  fqTmp <- tempfile(pattern = "trimmed", tmpdir = ".", fileext = ".fastq.gz")
  cmd = paste("seqtk seq -L", length, inFile, "| pigz --fast -p1 >", fqTmp)
  ezSystem(cmd)
  if (is.null(outFile)){
    outFile <- inFile
  }
  ezSystem(paste("mv", fqTmp, outFile))
  return(invisible(outFile))
}