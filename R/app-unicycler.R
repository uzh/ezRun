ezMethodUnicycler = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
  sampleName = input$getNames()
  ##stopifnot((param$paired))
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  basicOpts <- param$unicyclerOpt
  if (param$long == "YES" & param$paired == "true"){
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    long=param$pathToLong
    readOpt = paste("-1", read1, "-2", read2, "-l", long)
    cmd = paste("unicycler", readOpt, basicOpts, "--keep 1",
                "--mode", param$mode, "-o", "unicycler", '-t', ezThreads())
    ezSystem(cmd) 
    } else if (param$long == "NO" & param$paired == "true") {
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste("-1", read1, "-2", read2)
    cmd = paste("unicycler", readOpt, basicOpts, "--keep 1",
                "--mode", param$mode, "-o", "unicycler", '-t', ezThreads())
    ezSystem(cmd)
    } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste("-s", read1)
    cmd = paste("unicycler", readOpt, basicOpts, "--keep 1",
                "--mode", param$mode, "-o", "unicycler", '-t', ezThreads())
    ezSystem(cmd)
    }
  wddir <- "."
  gfafile <- file.path(wddir, "unicycler/assembly.gfa")
  ezSystem(paste("cp", gfafile, basename(output$getColumn("Graph"))))
  afile <- file.path(wddir, "unicycler/assembly.fasta")
  ezSystem(paste("cp", afile, basename(output$getColumn("Draft"))))
  logfile <- file.path(wddir, "unicycler/unicycler.log")
  ezSystem(paste("cp", logfile, basename(output$getColumn("Log"))))
  ezSystem(paste("mv", "unicycler", sampleName))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodUnicycler()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppUnicycler <-
  setRefClass("EzAppUnicycler",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodUnicycler
                  name <<- "EzAppUnicycler"
                  appDefaults <<- rbind(
                    pathToLong = ezFrame(Type="character", DefaultValue="", Description="Specify path to the FASTA or FASTQ file for your long reads."),
                    unicyclerOpt = ezFrame(Type="character",  DefaultValue="",  Description="Predefine any options not already specified by default. By default is empty for genome assembly"))
                }
              )
  )