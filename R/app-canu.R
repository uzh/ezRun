###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCanu = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  SMRT_Path = input$getFullPaths("Read1")
  SMRT_File = basename(input$getColumn("Read1"))
  #inputFileType = param$inputType
  #if (inputFileType == "pacbioSmrtCell"){
  #ezSystem(paste("cp -r", SMRT_Path, "."))
  #ezSystem(paste("mkdir", "smrt_input"))
  #ezSystem(paste("tar -zxf", SMRT_File, "--strip-components=4 -C smrt_input"))
  #readFile = file.path(getwd(), "smrt_input", "Analysis_Results", "*.subreads.fastq")
  #ezSystem(paste("cat", readFile, ">", paste0(sampleName,".fastq")))
  #cmd = paste("canu", "-p", "canu", "-d", "canu", paste0("genomeSize=", param$canuGenomeSize, "k"), paste0("maxMemory=", param$ram), paste0("maxThreads=", ezThreads()), opt, param$canuReadOpt, paste0(sampleName,".fastq"), "1> ", paste0(sampleName,"_canu.log"))
  #ezSystem(cmd)
  #ezSystem(paste("cp", "canu/canu.contigs.fasta", basename(output$getColumn("Draft"))))
  #ezSystem(paste("mv", "canu", sampleName))
  #return("Success")
  #}else{
  cmd = paste(
    "canu",
    "-p",
    sampleName,
    "-d",
    sampleName,
    paste0("genomeSize=", param$canuGenomeSize, "k"),
    paste0("maxMemory=", param$ram),
    paste0("maxThreads=", ezThreads()),
    opt,
    param$canuReadOpt,
    input$getFullPaths("Read1"),
    "1> ",
    paste0(sampleName, "_canu.log")
  )
  ezSystem(cmd)
  contigOutputFile <- paste0(sampleName, "/", sampleName, ".contigs.fasta")
  ezSystem(paste("cp", contigOutputFile, basename(output$getColumn("Draft"))))
  return("Success")
  #}
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppCanu <-
  setRefClass(
    "EzAppCanu",
    contains = "EzApp",
    methods = list(
      ## Canu: always runs (ezMethodCanu).
      ## HiCanu: Canu's mode for PacBio HiFi reads, when canuReadOpt is -pacbio-hifi
      ## (canuReadOpt is passed as the read-type flag in ezMethodCanu).
      citation = function(param = list()) {
        c(
          "Koren, S., Walenz, B.P., Berlin, K., Miller, J.R., Bergman, N.H. & Phillippy, A.M. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Research 27(5), 722-736 (2017). https://doi.org/10.1101/gr.215087.116",
          if (isTRUE(grepl("pacbio-hifi", param$canuReadOpt))) "Nurk, S. et al. HiCanu: accurate assembly of segmental duplications, satellites, and allelic variants from high-fidelity long reads. Genome Research 30(9), 1291-1305 (2020). https://doi.org/10.1101/gr.263566.120"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCanu
        name <<- "EzAppCanu"
        appDefaults <<- rbind(
          canuReadOpt = ezFrame(
            Type = "character",
            DefaultValue = "-pacbio-raw",
            Description = "input read types: -pacbio-raw, -pacbio-corrected, -nanopore-raw, -nanopore-corrected. Default is pacbio raw data"
          ),
          canuGenomeSize = ezFrame(
            Type = "integer",
            DefaultValue = "5000",
            Description = "estimated genome size in Kbp"
          )
        )
      }
    )
  )
