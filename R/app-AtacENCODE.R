###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### ATAC-seq ENCODE pipeline
### https://www.encodeproject.org/atac-seq/

EzAppAtacENCODE <-
  setRefClass("EzAppAtacENCODE",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodAtacENCODE
                  name <<- "EzAppAtacENCODE"
                }
              )
  )

ezMethodAtacENCODE <- function(input=NA, output=NA, param=NA){
    require(rjson)
    #    require(Herper)
   #local_CondaEnv("encode-atac-seq-pipeline", pathToMiniConda = "/usr/local/ngseq/miniconda3")
  #setwdNew(basename(output$getColumn("Report")))
  sampleName <- input$getNames()
  createJson(input, param)
  
  cmd <- paste('caper run /usr/local/ngseq/opt/atac-seq-pipeline/atac.wdl -i atac.json')
  system(cmd)
  
  ##find qc.html
  system(paste('find . -name qc.html -type f -print0| xargs -r0 cp -t .;'))
  file.rename('qc.html', paste0(sampleName, '_qc.html'))
  system(paste('find . -name qc.json -type f -print0| xargs -r0 cp -t .;'))
  file.rename('qc.json', paste0(sampleName, '_qc.json'))
  return("Success")
}


createJson <- function(input, param){
    atac_qc = list()
    sampleName <- input$getNames()
    atac_qc[["atac.title"]] <- sampleName
    atac_qc[["atac.description"]] <- paste('ATAC-seq on', sampleName)
    atac_qc[["atac.pipeline_type"]] <- "atac"
    atac_qc[["atac.align_only"]] <- FALSE
    atac_qc[["atac.true_rep_only"]] <- FALSE
    
    ## This pipeline only supports human and mouse
    species <- input$getColumn("Species")
    stopifnot(length(species)==1L)
    
    if(grepl("(Homo|human)", species, ignore.case = TRUE)){
        refBuild <- "hg38"
    } else if(grepl("(Mus|mouse)", species, ignore.case = TRUE)){
        refBuild <- "mm10"
    } else{
        stop("Only human and mouse are supported by this pipeline!")
    }
    atac_qc[["atac.genome_tsv"]] <- paste0("https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/",refBuild, '.tsv')
    ## We only work on paired-end data for now; although the pipeline supporst SE
    atac_qc[["atac.paired_end"]] <- TRUE
    atac_qc[["atac.fastqs_rep1_R1"]] <- as.vector(input$getFullPaths("Read1"))
    atac_qc[["atac.fastqs_rep1_R2"]] <- as.vector(input$getFullPaths("Read2"))
    if(param$nReads>0){
        atac_qc[["atac.subsample_reads"]] <- param$nReads
        atac_qc[["atac.xcor_subsample_reads"]] <- param$nReads
    }
    atac_qc[["atac.auto_detect_adapter"]] <- TRUE
    atac_qc[["atac.enable_xcor"]] <- TRUE
    atac_qc[["atac.multimapping"]] <- 4
    atac_qc[["atac.cap_num_peak"]] <- 300000
    atac_qc[["atac.pval_thresh"]] <- 0.01
    atac_qc[["atac.smooth_win"]] <- 150
    atac_qc[["atac.align_cpu"]] <- param$cores
    atac_qc[["atac.align_mem_factor"]] <- 0.15
    atac_qc[["atac.align_time_hr"]] <- 48
    atac_qc[["atac.align_disk_factor"]] <- 8.0
    
    myfile <- toJSON(atac_qc)
    myfile <- gsub(',', ',\n', myfile)
    myfile <- sub('R1\\":', 'R1\\":[', myfile)
    myfile <- sub('R2\\":', 'R2\\":[', myfile)
    myfile <- gsub('fastq.gz\\"', 'fastq.gz\\"]', myfile)
    write(myfile, "atac.json")
}
