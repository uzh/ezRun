###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastQC_10x <- function(input=NA, output=NA, param=NA, 
                               htmlFile="00index.html"){
  require(rmarkdown)
  require(tidyverse)
  setwdNew(basename(output$getColumn("Report")))
  
  dataset <- input$meta
  samples <- rownames(dataset)
  sampleDirs <- input$getFullPaths("RawDataDir")
  stopifnot(all(grepl("\\.tar$", sampleDirs)))
  
  taredfiles <- lapply(sampleDirs, untar, list=TRUE)
  taredfiles_R1 <- sapply(taredfiles, function(x){grep("_R1_", x, value=TRUE) %>% head(1)})
  taredfiles_R2 <- sapply(taredfiles, function(x){grep("_R2_", x, value=TRUE) %>% head(1)})
  for(i in 1:length(sampleDirs)){
    untar(sampleDirs[i], files=c(taredfiles_R1[i], taredfiles_R2[i]))
  }
  taredfiles_R1 <- normalizePath(taredfiles_R1)
  taredfiles_R2 <- normalizePath(taredfiles_R2)
  
  dataset$`Read1` <- taredfiles_R1
  dataset$`Read2` <- taredfiles_R2
  input <- EzDataset(meta=dataset, dataRoot=NULL)
  if(sum(input$getColumn("Read Count")) > 1e9){
    doSubsample <- TRUE # subsample to 1Mio reads
    input <- ezMethodSubsampleFastq(input=input, param=param)
    dataset = input$meta
  }else{
    doSubsample <- FALSE
  }
  
  files = c()
  for (sm in samples){
    files[paste0(sm, "_R1")] = input$getFullPaths("Read1")[sm]
    if(isTRUE(param$paired)){
      files[paste0(sm, "_R2")] = input$getFullPaths("Read2")[sm]
    }
  }
  nFiles = length(files)
  
  ## guess the names of the report directories that will be creatd by fastqc
  reportDirs = sub("\\.fastq(\\.gz)*$", "_fastqc", basename(files))
  reportDirs = sub("\\.fq(\\.gz)*$", "_fastqc", reportDirs)
  reportDirs = sub("\\.bam$", "_fastqc", reportDirs)
  stopifnot(!duplicated(reportDirs))
  filesUse = files[!file.exists(reportDirs)]
  if (length(filesUse) > 0){
    cmd = paste("fastqc", "--extract -o . -t", min(param$cores, 8),
                "-a", FASTQC_ADAPTERS, "--kmers 7",
                param$cmdOptions, paste(filesUse, collapse=" "),
                "> fastqc.out", "2> fastqc.err")
    result = ezSystem(cmd)
  }
  statusToPng = c(PASS="tick.png", WARN="warning.png", FAIL="error.png")

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "FastQC.Rmd", "FastQC_overview.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  
  ## collect the overview table
  plots = c("Per base sequence quality"="per_base_quality.png",
            "Per sequence quality scores"="per_sequence_quality.png",
            "Per tile sequence quality"="per_tile_quality.png",
            "Per base sequence content"="per_base_sequence_content.png",
            "Per sequence GC content"="per_sequence_gc_content.png",
            "Per base N content"="per_base_n_content.png",
            "Sequence Length Distribution"="sequence_length_distribution.png",
            "Sequence Duplication Levels"="duplication_levels.png",
            "Adapter Content"="adapter_content.png",
            "Kmer Content"="kmer_profiles.png")
  
  ## make for each plot type an html report with all samples
  plotPages = sub(".png", ".html", plots)
  for(i in 1:length(plots)){
    plotPage <- plotPages[i]
    pngs <- file.path(reportDirs, "Images", plots[i])
    rmarkdown::render(input="FastQC_overview.Rmd", envir = new.env(),
                      output_dir=".", output_file=plotPage, quiet=TRUE)
  }
  
  ## establish the main report
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  if (!is.null(dataset$"Read Count")){
    readCount = signif(dataset$"Read Count" / 1e6, digits=3)
    names(readCount) = rownames(dataset)
  } else {
    readCount = integer()
    for (i in 1:nFiles){
      x = ezRead.table(file.path(reportDirs[i], "fastqc_data.txt"), 
                       header=FALSE, nrows=7, fill=TRUE)
      readCount[names(files)[i]] = signif(as.integer(x["Total Sequences", 1]) / 
                                            1e6, digits=3)
    }
  }
  ans4Report[["Read Counts"]] <- readCount
  
  ## Each sample can have different number of reports.
  ## Especially per tile sequence quality
  nrReports <- sapply(reportDirs, 
                      function(x){smy = ezRead.table(file.path(x, "summary.txt"), 
                                                     row.names=NULL, header=FALSE);
                      nrow(smy)})
  i <- which.max(nrReports)
  smy = ezRead.table(file.path(reportDirs[i], "summary.txt"), 
                     row.names=NULL, header=FALSE)
  rowNames = paste0("<a href=", reportDirs, "/fastqc_report.html>", 
                    names(files), "</a>")
  tbl = ezMatrix("", rows=rowNames, cols=smy[[2]])
  for (i in 1:nFiles){
    smy = ezRead.table(file.path(reportDirs[i], "summary.txt"), 
                       row.names=NULL, header=FALSE)
    href = paste0(reportDirs[i], "/fastqc_report.html#M", 
                  0:(ncol(tbl)-1))[colnames(tbl) %in% smy[[2]]]
    img = paste0(reportDirs[1], 	"/Icons/", statusToPng[smy[[1]]]) ##
    tbl[i, colnames(tbl) %in% smy[[2]]] = paste0("<a href=", href, 
                                                 "><img src=", img, "></a>")
  }
  colnames(tbl) <- ifelse(colnames(tbl) %in% names(plotPages),
                          paste0("<a href=", plotPages[colnames(tbl)], 
                                 ">", colnames(tbl),
                                 "</a>"),
                          colnames(tbl))
  
  ans4Report[["Fastqc quality measures"]] <- tbl
  
  qualMatrixList = ezMclapply(files, getQualityMatrix, mc.cores=param$cores)
  ans4Report[["Per Base Read Quality"]] <- qualMatrixList
  
  ## debug
  ## save(ans4Report, file="ans4Report.rda")
  
  ## generate the main reports
  rmarkdown::render(input="FastQC.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  prepareRmdLib()
  
  ## generate multiQC report
  ezSystem("multiqc .")
  
  unlink(dirname(taredfiles_R1), recursive=TRUE)
  unlink(paste0(reportDirs, ".zip"), recursive = TRUE)
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }
EzAppFastqc_10x <-
  setRefClass("EzAppFastqc_10x",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFastQC_10x
                  name <<- "EzAppFastqc_10x"
                }
              )
  )
