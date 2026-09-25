###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodSamsa2 = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  ### metatranscriptomics assemby with Samsa2, annotation with RefSeq

  library(plyr)

  sampleName = input$getNames()
  file1PathInDatset <- input$getFullPaths("Read1")
  #fastqName1 <- paste0(sampleName,".R1.fastq")
  #cpCmd1 <- paste0("gunzip -c ", file1PathInDatset, "  > ", fastqName1)
  #ezSystem(cpCmd1)
  if (param$paired) {
    file2PathInDatset <- input$getFullPaths("Read2")
    #fastqName2 <- paste0(sampleName,".R2.fastq")
    #cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", fastqName2)
    #ezSystem(cpCmd2)
  }
  ##make input directory
  make_inputDir <- paste(
    "if [ -d input_files ]; then echo dir_exists; else mkdir input_files; fi"
  )
  ezSystem(make_inputDir)
  move_input1 <- paste("cp", file1PathInDatset, "input_files")
  ezSystem(move_input1)
  if (param$paired) {
    move_input2 <- paste("cp", file2PathInDatset, "input_files")
    ezSystem(move_input2)
  }
  ##copy samsa2
  cpSamsa2 <- paste("cp -r /usr/local/ngseq/src/samsa2/ .")
  ezSystem(cpSamsa2)
  #updateTemplateScriptCmd <- paste("cp",
  #                                 file.path(METAGENOMICS_ROOT,SAMSA2_TEMPLATE_SCRIPT),
  #                                  SAMSA2_TEMPLATE_SCRIPT)
  #ezSystem(updateTemplateScriptCmd)
  ## run bash script
  samsa2ToBeExecCmd <- paste(
    "bash ./samsa2/bash_scripts/master_script.sh ./input_files/ ."
  )
  ezSystem(samsa2ToBeExecCmd)

  ## place output files

  oldAnnFile <- list.files(
    "step_4_output",
    pattern = "RefSeq_annotated",
    full.names = T
  )
  newAnnFile <- basename(output$getColumn("annotationFileRefSeq"))
  ezSystem(paste("mv", oldAnnFile, newAnnFile))
  oldAnnFile_org <- list.files(
    "step_5_output/RefSeq_results/org_results",
    pattern = "RefSeq_annot",
    full.names = T
  )
  newAnnFile_org <- basename(output$getColumn("annotationORGFileRefSeq"))
  ezSystem(paste("mv", oldAnnFile_org, newAnnFile_org))
  oldAnnFile_func <- list.files(
    "step_5_output/RefSeq_results/func_results",
    pattern = "RefSeq_annot",
    full.names = T
  )
  newAnnFile_func <- basename(output$getColumn("annotationFUNCFileRefSeq"))
  ezSystem(paste("mv", oldAnnFile_func, newAnnFile_func))
}

##' @template app-template
##' @templateVar method ezMethodSamsa2()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppSamsa2 <-
  setRefClass(
    "EzAppSamsa2",
    contains = "EzApp",
    methods = list(
      ## The pipeline is the installed SAMSA2 master script, run unmodified; its steps
      ## are not in the job log, only the one `bash master_script.sh` call.
      ## Lines below refer to /usr/local/ngseq/src/samsa2/bash_scripts/master_script.sh.
      ## SAMSA2 master_script.sh: Trimmomatic, SortMeRNA, DIAMOND vs RefSeq unconditional; PEAR gated
      ## on paired input.
      citation = function() {
        c(
          "Westreich, S.T. et al. SAMSA2: a standalone metatranscriptome analysis pipeline. BMC Bioinformatics 19, 175 (2018). https://doi.org/10.1186/s12859-018-2189-z",
          "Bolger, A.M. et al. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30, 2114-2120 (2014). https://doi.org/10.1093/bioinformatics/btu170",
          "Zhang, J. et al. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30, 614-620 (2014). https://doi.org/10.1093/bioinformatics/btt593",
          "Kopylova, E. et al. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics 28, 3211-3217 (2012). https://doi.org/10.1093/bioinformatics/bts611",
          "Buchfink, B. et al. Fast and sensitive protein alignment using DIAMOND. Nature Methods 12, 59-60 (2015). https://doi.org/10.1038/nmeth.3176",
          "O'Leary, N.A. et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Research 44, D733-D745 (2016). https://doi.org/10.1093/nar/gkv1189"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## master_script.sh:113-115
          "Reads were quality-trimmed with Trimmomatic (SLIDINGWINDOW:4:15, MINLEN:70, phred33), in paired-end mode when an R2 file matching the R1 file name was present.",
          ## master_script.sh:149
          if (isTRUE(as.logical(param$paired))) "Paired reads were merged with PEAR and only the merged (assembled) reads were carried forward; unmerged pairs were not analysed further.",
          ## master_script.sh:209
          "Ribosomal RNA reads were removed with SortMeRNA against the SILVA bacterial 16S database (silva-bac-16s-id90) only.",
          ## master_script.sh:239
          "The remaining reads were aligned with DIAMOND blastx against the SAMSA2 RefSeq bacterial protein database, keeping only the best hit per read (-k 1).",
          ## master_script.sh:265-266
          "Hits were aggregated into organism-level and function-level read counts with the SAMSA2 script DIAMOND_analysis_counter.py.",
          ## master_script.sh:280 exits before the Subsystems and DESeq2 steps
          "The SEED Subsystems annotation and the SAMSA2 DESeq2 statistics were not run: the installed script stops after the RefSeq annotation, whatever useSubsystemDB was set to."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSamsa2
        name <<- "EzAppSamsa2"
        appDefaults <<- rbind(
          useSubsystemDB = ezFrame(
            Type = "logical",
            DefaultValue = "RefSeq",
            Description = "database"
          )
        )
      }
    )
  )
