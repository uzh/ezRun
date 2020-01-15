###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title exceRpt_smallRNA report app
##' @description Use this reference class to process exceRpt_smallRNA's output after execution and perform a QC.
##' @author Miquel Anglada Girotto
EzAppExceRptReport =
  setRefClass( "EzAppExceRptReport",
               contains = "EzApp",
               methods = list(
                 initialize = function()
                 {
                   "Initializes the application using its specific defaults."
                   runMethod <<- ezMethodExceRptReport
                   name <<- "EzAppExceRptReport"
                   appDefaults <<- rbind(name=ezFrame(Type="character",  DefaultValue="Excerpt_Report",  Description=""))
                 }
               )
               
  )

##' @title execute exceRpt_smallRNA's output processing.
##' @description create QC plots and count smallRNAs detected after executing ezMethodExceRpt().
##' @param input ezDataFrame()
##' @param output ezDataFrame()
##' @param param list() environment parameters.
##' @author Miquel Anglada Girotto
ezMethodExceRptReport = function(input=NA, output=NA, param=NA){
  ## set report's working directory
  setwdNew(basename(output$getColumn("Report")))
  
  ## Process output
  processedOutputDir = "processed_output"
  selectedSamples = param[['samples']]
  
  samplePaths = input$getFullPaths("excerpt")
  idx = basename(samplePaths) %in% selectedSamples
  samplePathsFilt = samplePaths[idx]
  
  plots = processSamples(samplePaths = samplePathsFilt, outputDir = processedOutputDir, getPlotsObjects=TRUE)

  ## list files generated
  dataFiles = list.files(processedOutputDir, pattern = '*.txt')
  
  ## Create report
  makeRmdReport(plots=plots, dataFiles=dataFiles, rmdFile = "excerpt.Rmd")
}


##' @title generate smallRNA counts and QC.
##' @description Obtain counts and QC plots across all samples.
##' @param samplePaths <string> vector of all the files that need to be processed together.
##' @param outputDir <string> directory name where to output generated objects and report after processing.
##' @param getPlotsObjects <bool> TRUE to output plots as an object.
##' @author (Rob Kitchen) Miquel Anglada Girotto
processSamples = function(samplePaths, outputDir, getPlotsObjects=FALSE){
  ## Load required dependencies
  loadDependencies()
  
  ## delete -e from .stats (they do not appear in the ExampleData, but they do when I ran the program)
  #delete_e(samplePaths)
  
  ## Get directories containing counts
  samplePathList = list.dirs(samplePaths,recursive=F)
  
  ## create output dir
  dir.create(outputDir)
  
  ## reads, normalises, and saves individual sample results
  sampleIDs = readData(samplePathList = samplePathList, output.dir = outputDir)
  
  ## plot the data
  plotsList = PlotData(sampleIDs, outputDir)
  if(getPlotsObjects==TRUE){return(plotsList)}
}


##' @title load dependencies
##' @description load dependencies required within processSamples()
##' @author Miquel Anglada Girotto
loadDependencies = function(){
  ## load
  require(plyr)
  require(gplots)
  require(marray)
  require(reshape2)
  require(ggplot2)
  require(tools)
  require(Rgraphviz)
  require(scales)
}


##' @title Get smallRNA counts and stats
##' @description Get smallRNA counts and stats to return as final output and make QC plots.
##' @author (Rob Kitchen) Miquel Anglada Girotto
readData = function(samplePathList, output.dir){
  
  ##
  ## Create objects to contain the data
  ##
  sample.data = vector(mode="list",length=length(samplePathList))
  allIDs.calibrator = NULL
  allIDs.miRNA = NULL
  allIDs.tRNA = NULL
  allIDs.piRNA = NULL
  allIDs.gencode = NULL
  allIDs.circularRNA = NULL
  allIDs.exogenous_miRNA = NULL
  allIDs.exogenous_rRNA = NULL
  allIDs.exogenous_genomes = NULL
  taxonomyInfo.exogenous_rRNA = NULL
  taxonomyInfo.exogenous_genomes = NULL
  
  mapping.stats = matrix(0,nrow=length(samplePathList),ncol=30, dimnames=list(1:length(samplePathList), c("input","successfully_clipped","failed_quality_filter","failed_homopolymer_filter","calibrator","UniVec_contaminants","rRNA","reads_used_for_alignment","genome","miRNA_sense","miRNA_antisense","miRNAprecursor_sense","miRNAprecursor_antisense","tRNA_sense","tRNA_antisense","piRNA_sense","piRNA_antisense","gencode_sense","gencode_antisense","circularRNA_sense","circularRNA_antisense","not_mapped_to_genome_or_libs","repetitiveElements","endogenous_gapped","input_to_exogenous_miRNA","exogenous_miRNA","input_to_exogenous_rRNA","exogenous_rRNA","input_to_exogenous_genomes","exogenous_genomes")))
  qc.results = matrix(0,nrow=length(samplePathList),ncol=5, dimnames=list(1:length(samplePathList), c("InputReads","GenomeReads","TranscriptomeReads","TranscriptomeGenomeRatio","TranscriptomeComplexity")))
  maxReadLength = 10000
  read.lengths = matrix(0,nrow=length(samplePathList),ncol=maxReadLength+1,dimnames=list(1:length(samplePathList), 0:maxReadLength))
  
  
  ##
  ## Loop through all samples and read the pipeline output
  ##
  printMessage(c("Reading sample data..."))
  removeSamples = NULL
  for(i in 1:length(samplePathList)){
    ## Parse the sampleID from the path:
    tmp = unlist(strsplit(samplePathList[i], "/"))
    thisSampleID = tmp[length(tmp)]
    
    ## Get timings and check this sample finished successfully
    tmp.stats = readLines(paste(samplePathList[i],".stats",sep=""))
    tmp.stats = delete_e(tmp.stats)
    tmp.stats = do.call(rbind,strsplit(tmp.stats,'\t'))
    tmp.stats = data.frame(tmp.stats,stringsAsFactors = FALSE)
    tmp.stats[c(1,nrow(tmp.stats)),2] = ""
    x.start = grep("#STATS",tmp.stats[,1])
    x.end = grep("#END OF STATS",tmp.stats[,1])
    if(length(x.start) > 0  &&  length(x.end) > 0){
      tmp.start = strptime(unlist(strsplit(tmp.stats[x.start[1],1],"Run started at "))[2],"%Y-%m-%d--%H:%M:%S")
      tmp.end = strptime(unlist(strsplit(tmp.stats[x.end[1],1],"Run completed at "))[2],"%Y-%m-%d--%H:%M:%S")
      runTiming = data.frame(start=tmp.start, completed=tmp.end, duration=difftime(tmp.end,tmp.start), duration_secs=as.numeric(difftime(tmp.end,tmp.start,units="secs")))
      continue = T
    }else{
      continue = F
      removeSamples = c(removeSamples, i)
      printMessage(c("[",i,"/",length(samplePathList),"] WARNING: Incomplete run for sample \'",thisSampleID,"\', ignoring"))
    }
    
    if(continue == T){
      ##
      ## Read sample mapping stats
      ##
      tmp.stats = tmp.stats[-c(1,2,nrow(tmp.stats)),]
      tmp.stats[tmp.stats[,1] %in% "clipped", 1] = "successfully_clipped"
      mapping.stats[i, match(tmp.stats[,1], colnames(mapping.stats))] = as.numeric(tmp.stats[,2])
      rownames(mapping.stats)[i] = thisSampleID
      
      ##
      ## Read the QC result
      ##
      adapterConfidence = NA
      qcOutcome = NA
      if(file.exists(paste(samplePathList[i],".qcResult",sep=""))){
        tmp.qc = read.table(paste(samplePathList[i],".qcResult",sep=""), stringsAsFactors=F, fill=T, header=F, sep=" ",skip=0)
        if(tmp.qc[1,1] == "Adapter_confidence:"){
          adapterConfidence = tmp.qc[1,2]
          tmp.qc = tmp.qc[-1,]
        }
        qcOutcome = tmp.qc[1,2]
        
        qc.results[i, match(gsub(":$","",tmp.qc[-1,1]), colnames(qc.results))] = as.numeric(tmp.qc[-1,2])
        #qc.results[i, ] = as.numeric(tmp.qc[-1,2])
        rownames(qc.results)[i] = thisSampleID
      }
      
      ##
      ## Read the adapter sequence
      ##
      if(paste(thisSampleID,".adapterSeq",sep="") %in% dir(samplePathList[i])){
        tmp.seq = try(read.table(paste(samplePathList[i],"/",thisSampleID,".adapterSeq",sep="")), silent=T)
        if(class(tmp.seq) == "try-error"){
          adapterSeq = NA
        }else{
          adapterSeq = as.character(tmp.seq[1,1])
        }
      }
      
      
      ##
      ## Read the calibrator counts, if available
      ##
      calibratorCounts = NULL
      if(paste(thisSampleID,".clipped.trimmed.filtered.calibratormapped.counts",sep="") %in% dir(samplePathList[i])){
        calibratorCounts = try(read.table(paste(samplePathList[i],"/",thisSampleID,".clipped.trimmed.filtered.calibratormapped.counts",sep=""), stringsAsFactors=F)[,2:1], silent=T)
        if(class(calibratorCounts) == "try-error"){
          calibratorCounts = NULL
        }else{
          colnames(calibratorCounts) = c("calibratorID","readCount")
        }
      }
      
      
      ##
      ## Read the clipped read lengths  
      ##
      if(length(grep(".readLengths.txt$", dir(samplePathList[i]))) == 1){
        tmp = read.table(paste(samplePathList[i], dir(samplePathList[i])[grep(".readLengths.txt$", dir(samplePathList[i]))], sep="/"))
        read.lengths[i, 1:ncol(tmp)] = as.numeric(tmp[2,])
        rownames(read.lengths)[i] = thisSampleID
      }
      
      
      ##
      ## Read sample data
      ##
      availableFiles = dir(samplePathList[i])
      miRNA_sense=miRNA_antisense = tRNA_sense=tRNA_antisense = piRNA_sense=piRNA_antisense = gencode_sense=gencode_antisense = circRNA_sense=circRNA_antisense = NULL
      
      if("readCounts_miRNAmature_sense.txt" %in% availableFiles){
        miRNA_sense = read.table(paste(samplePathList[i],"readCounts_miRNAmature_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        miRNA_sense = cbind(miRNA_sense, ID=sapply(rownames(miRNA_sense), function(id){ multiID = unlist(strsplit(id,"\\|")); multiIDs = sapply(multiID, function(idPart){unlist(strsplit(idPart,":"))[1]});   if(length(multiIDs) == 1){ multiIDs }else{ paste(sort(multiIDs),collapse="|") }  }))
      }
      if("readCounts_miRNAmature_antisense.txt" %in% availableFiles){
        miRNA_antisense = read.table(paste(samplePathList[i],"readCounts_miRNAmature_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        miRNA_antisense = cbind(miRNA_antisense, ID=sapply(rownames(miRNA_antisense), function(id){ multiID = unlist(strsplit(id,"\\|")); multiIDs = sapply(multiID, function(idPart){unlist(strsplit(idPart,":"))[1]});   if(length(multiIDs) == 1){ multiIDs }else{ paste(sort(multiIDs),collapse="|") }  }))
      }
      
      if("readCounts_tRNA_sense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_tRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), function(id){ unlist(strsplit(id,"-"))[2]  }))
        tRNA_sense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(tRNA_sense)[-1] = colnames(tmp)[1:4]
        tRNA_sense = tRNA_sense[order(tRNA_sense$multimapAdjustedReadCount,decreasing=T), ]
      }
      if("readCounts_tRNA_antisense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_tRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), function(id){ unlist(strsplit(id,"-"))[2]  }))
        tRNA_antisense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(tRNA_antisense)[-1] = colnames(tmp)[1:4]
        tRNA_antisense = tRNA_antisense[order(tRNA_antisense$multimapAdjustedReadCount,decreasing=T), ]
      }
      
      if("readCounts_piRNA_sense.txt" %in% availableFiles){
        piRNA_sense = read.table(paste(samplePathList[i],"readCounts_piRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      if("readCounts_piRNA_antisense.txt" %in% availableFiles){
        piRNA_antisense = read.table(paste(samplePathList[i],"readCounts_piRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      
      
      makeGeneID = function(id){ 
        bits=unlist(strsplit(id,":")); geneNameBits=unlist(strsplit(bits[3],"-")); 
        geneName = geneNameBits[1]
        if(length(geneNameBits) > 2){ geneName=paste(geneNameBits[-length(geneNameBits)],collapse="-") }
        paste(geneName,bits[2],sep=":") 
      }
      if("readCounts_gencode_sense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_gencode_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), makeGeneID))
        gencode_sense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(gencode_sense)[-1] = colnames(tmp)[1:4]
        gencode_sense = gencode_sense[order(gencode_sense$multimapAdjustedReadCount,decreasing=T), ]
      }
      if("readCounts_gencode_antisense.txt" %in% availableFiles){
        tmp = read.table(paste(samplePathList[i],"readCounts_gencode_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        tmp = cbind(tmp, ID=sapply(rownames(tmp), makeGeneID))
        gencode_antisense = ddply(tmp, "ID", function(mat){ c(as.numeric(mat[1,1:2]),sum(mat$multimapAdjustedReadCount),sum(mat$multimapAdjustedBarcodeCount)) })
        colnames(gencode_antisense)[-1] = colnames(tmp)[1:4]
        gencode_antisense = gencode_antisense[order(gencode_antisense$multimapAdjustedReadCount,decreasing=T), ]
      }
      
      if("readCounts_circRNA_sense.txt" %in% availableFiles){
        circRNA_sense = read.table(paste(samplePathList[i],"readCounts_circRNA_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      if("readCounts_circRNA_antisense.txt" %in% availableFiles){
        circRNA_antisense = read.table(paste(samplePathList[i],"readCounts_circRNA_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
      }
      
      
      ##
      ## Read exogenous miRNA alignments (if applicable)
      ##
      exogenous_miRNA_sense = NA
      exogenous_miRNA_IDs = NULL
      if("EXOGENOUS_miRNA" %in% availableFiles){
        tmp.dir = paste(samplePathList[i],"EXOGENOUS_miRNA",sep="/")
        if("readCounts_miRNAmature_sense.txt" %in% dir(tmp.dir)){
          exogenous_miRNA_sense = read.table(paste(tmp.dir,"readCounts_miRNAmature_sense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
          exogenous_miRNA_IDs = rownames(exogenous_miRNA_sense)
        }
        if("readCounts_miRNAmature_antisense.txt" %in% dir(tmp.dir)){
          exogenous_miRNA_antisense = read.table(paste(tmp.dir,"readCounts_miRNAmature_antisense.txt",sep="/"), header=T, sep="\t", comment.char="", stringsAsFactors=F,colClasses=c("character","numeric","numeric","numeric","numeric"), row.names=1)
        }
      }
      
      
      ##
      ## Read exogenous rRNA alignments (if applicable)
      ##
      exogenous_rRNA = NA
      exogenous_rRNA_IDs = NULL
      if("EXOGENOUS_rRNA" %in% availableFiles){
        tmp.dir = paste(samplePathList[i],"EXOGENOUS_rRNA",sep="/")
        if("ExogenousRibosomalAlignments.result.taxaAnnotated.txt" %in% dir(tmp.dir)){
          exogenous_rRNA = try(read.table(paste(tmp.dir,"/ExogenousRibosomalAlignments.result.taxaAnnotated.txt",sep=""), sep="\t", stringsAsFactors = F, quote="", comment.char="",header=T), silent=T) 
          if(class(exogenous_rRNA) == "try-error"){
            exogenous_rRNA = NULL
          }#else{
          #  colnames(exogenous_rRNA) = c("indent","distFromRoot","level","name","uniqueReads","allSumReads")
          #}
          ## remove the 'Bacteria' stick insect!
          i.toRemove = which(exogenous_rRNA$name == "Bacteria"  &  exogenous_rRNA$level == "genus")
          if(length(i.toRemove) > 0)
            exogenous_rRNA = exogenous_rRNA[-i.toRemove, ]
          exogenous_rRNA_IDs = exogenous_rRNA$name
          taxonomyInfo.exogenous_rRNA = unique(rbind(taxonomyInfo.exogenous_rRNA, exogenous_rRNA[,1:5]))
        }
      }
      
      
      ##
      ## Read exogenous genome alignments (if applicable)
      ##
      exogenous_genomes = NA
      exogenous_genomes_IDs = NULL
      if("EXOGENOUS_genomes" %in% availableFiles){
        tmp.dir = paste(samplePathList[i],"EXOGENOUS_genomes",sep="/")
        if("ExogenousGenomicAlignments.result.taxaAnnotated.txt" %in% dir(tmp.dir)){
          exogenous_genomes = read.table(paste(tmp.dir,"/ExogenousGenomicAlignments.result.taxaAnnotated.txt",sep=""), sep="\t", stringsAsFactors = F, quote="", comment.char="", header=T)
          if(class(exogenous_genomes) == "try-error"){
            exogenous_genomes = NULL
          }#else{
          #  colnames(exogenous_genomes) = c("indent","distFromRoot","level","name","uniqueReads","allSumReads")
          #}
          
          i.toRemove = which(exogenous_genomes$name == "Bacteria"  &  exogenous_genomes$level == "genus")
          if(length(i.toRemove) > 0)
            exogenous_genomes = exogenous_genomes[-i.toRemove, ]
          exogenous_genomes_IDs = exogenous_genomes$name
          taxonomyInfo.exogenous_genomes = unique(rbind(taxonomyInfo.exogenous_genomes, exogenous_genomes[,1:5]))
        }
      }
      
      # Update list of detected smallRNA IDs
      allIDs.calibrator = unique(c(allIDs.calibrator, as.character(calibratorCounts$calibratorID)))
      allIDs.miRNA = unique(c(allIDs.miRNA, as.character(miRNA_sense$ID)))
      allIDs.tRNA = unique(c(allIDs.tRNA, as.character(tRNA_sense$ID)))
      allIDs.piRNA = unique(c(allIDs.piRNA, rownames(piRNA_sense)))
      allIDs.gencode = unique(c(allIDs.gencode, as.character(gencode_sense$ID)))
      allIDs.circularRNA = unique(c(allIDs.circularRNA, rownames(circRNA_sense)))
      allIDs.exogenous_miRNA = unique(c(allIDs.exogenous_miRNA, exogenous_miRNA_IDs))
      allIDs.exogenous_rRNA = unique(c(allIDs.exogenous_rRNA, exogenous_rRNA_IDs))
      allIDs.exogenous_genomes = unique(c(allIDs.exogenous_genomes, exogenous_genomes_IDs))
      
      sample.data[[i]] = list("miRNA_sense"=miRNA_sense,"miRNA_antisense"=miRNA_antisense, "tRNA_sense"=tRNA_sense,"tRNA_antisense"=tRNA_antisense, "piRNA_sense"=piRNA_sense,"piRNA_antisense"=piRNA_antisense, "gencode_sense"=gencode_sense,"gencode_antisense"=gencode_antisense, "circRNA_sense"=circRNA_sense,"circRNA_antisense"=circRNA_antisense, "exogenous_miRNA_sense"=exogenous_miRNA_sense, "exogenous_rRNA"=exogenous_rRNA, "exogenous_genomes"=exogenous_genomes, "adapterSeq"=adapterSeq, "adapterConfidence"=adapterConfidence, "qcOutcome"=qcOutcome, "runTiming"=runTiming, "calibratorCounts"=calibratorCounts)
      names(sample.data)[i] = thisSampleID
      
      
      printMessage(c("[",i,"/",length(samplePathList),"] Added sample \'",thisSampleID,"\'"))
    }
  }
  
  
  
  ##
  ## Remove failed/incomplete samples
  ##
  stopifnot(length(removeSamples) < length(sample.data))
  if(length(removeSamples) > 0){
    read.lengths = read.lengths[-removeSamples, ]
    sample.data = sample.data[-removeSamples]
    mapping.stats = mapping.stats[-removeSamples,]
    qc.results = qc.results[-removeSamples,]
  }
  
  
  
  ##
  ## Trim read-length matrix
  ##
  read.lengths = read.lengths[,0:(max(as.numeric(colnames(read.lengths[, colSums(read.lengths) > 0, drop=F])))+1), drop=F]
  #read.lengths = read.lengths[,colSums(read.lengths) > 0, drop=F]
  
  
  ##
  ## Collect IDs
  ##
  allIDs = list("calibrator"=allIDs.calibrator, "miRNA_sense"=allIDs.miRNA, "tRNA_sense"=allIDs.tRNA, "piRNA_sense"=allIDs.piRNA, "gencode_sense"=allIDs.gencode, "circRNA_sense"=allIDs.circularRNA, "exogenous_miRNA"=allIDs.exogenous_miRNA, "exogenous_rRNA"=allIDs.exogenous_rRNA, "exogenous_genomes"=allIDs.exogenous_genomes)
  
  
  ##
  ## Convert sample data to large per-smallRNA expression matrices
  ##
  printMessage("Creating raw read-count matrices for available libraries")
  #run.duration = data.frame(runDuration_string=rep("",length(sample.data)), runDuration_secs=rep(0,length(sample.data)),stringsAsFactors = F)
  run.duration = data.frame(runDuration_secs=rep(0,length(sample.data)),stringsAsFactors = F)
  rownames(run.duration) = names(sample.data)
  
  exprs.calibrator = matrix(0,ncol=length(sample.data),nrow=length(allIDs$calibrator), dimnames=list(allIDs$calibrator, names(sample.data)))
  
  exprs.miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$miRNA_sense), dimnames=list(allIDs$miRNA_sense, names(sample.data)))
  exprs.tRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$tRNA_sense), dimnames=list(allIDs$tRNA_sense, names(sample.data)))
  exprs.piRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$piRNA_sense), dimnames=list(allIDs$piRNA_sense, names(sample.data)))
  exprs.gencode = matrix(0,ncol=length(sample.data),nrow=length(allIDs$gencode_sense), dimnames=list(allIDs$gencode_sense, names(sample.data)))
  exprs.circRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$circRNA_sense), dimnames=list(allIDs$circRNA_sense, names(sample.data)))
  exprs.exogenous_miRNA = matrix(0,ncol=length(sample.data),nrow=length(allIDs$exogenous_miRNA), dimnames=list(allIDs$exogenous_miRNA, names(sample.data)))
  
  if(is.null(taxonomyInfo.exogenous_rRNA))
    tmp.nrow = 0
  else
    tmp.nrow = nrow(taxonomyInfo.exogenous_rRNA)
  exprs.exogenousRibosomal_specific = matrix(0,ncol=length(sample.data),nrow=tmp.nrow, dimnames=list(taxonomyInfo.exogenous_rRNA$ID, names(sample.data)))
  exprs.exogenousRibosomal_cumulative = matrix(0,ncol=length(sample.data),nrow=tmp.nrow, dimnames=list(taxonomyInfo.exogenous_rRNA$ID, names(sample.data)))
  
  if(is.null(taxonomyInfo.exogenous_genomes))
    tmp.nrow = 0
  else
    tmp.nrow = nrow(taxonomyInfo.exogenous_genomes)
  exprs.exogenousGenomes_specific = matrix(0,ncol=length(sample.data),nrow=tmp.nrow, dimnames=list(taxonomyInfo.exogenous_genomes$ID, names(sample.data)))
  exprs.exogenousGenomes_cumulative = matrix(0,ncol=length(sample.data),nrow=tmp.nrow, dimnames=list(taxonomyInfo.exogenous_genomes$ID, names(sample.data)))
  
  for(i in 1:length(sample.data)){
    run.duration[i,] = sample.data[[i]]$runTiming[1,4,drop=F]
    
    if(!is.null(nrow(sample.data[[i]]$calibratorCounts)))
      exprs.calibrator[match(sample.data[[i]]$calibratorCounts$calibratorID, rownames(exprs.calibrator)),i] = as.numeric(sample.data[[i]]$calibratorCounts$readCount)
    
    exprs.miRNA[match(sample.data[[i]]$miRNA_sense$ID, rownames(exprs.miRNA)),i] = as.numeric(sample.data[[i]]$miRNA_sense$multimapAdjustedReadCount)
    exprs.tRNA[match(sample.data[[i]]$tRNA_sense$ID, rownames(exprs.tRNA)),i] = as.numeric(sample.data[[i]]$tRNA_sense$multimapAdjustedReadCount)
    exprs.piRNA[match(rownames(sample.data[[i]]$piRNA_sense), rownames(exprs.piRNA)),i] = as.numeric(sample.data[[i]]$piRNA_sense$multimapAdjustedReadCount)
    exprs.gencode[match(sample.data[[i]]$gencode_sense$ID, rownames(exprs.gencode)),i] = as.numeric(sample.data[[i]]$gencode_sense$multimapAdjustedReadCount)
    exprs.circRNA[match(rownames(sample.data[[i]]$circRNA_sense), rownames(exprs.circRNA)),i] = as.numeric(sample.data[[i]]$circRNA_sense$multimapAdjustedReadCount)
    ## Exogenous miRNA
    if(!is.null(nrow(sample.data[[i]]$exogenous_miRNA)))
      exprs.exogenous_miRNA[match(rownames(sample.data[[i]]$exogenous_miRNA), rownames(exprs.exogenous_miRNA)),i] = as.numeric(sample.data[[i]]$exogenous_miRNA$multimapAdjustedReadCount)
    ## Exogenous rRNA
    if(!is.null(nrow(sample.data[[i]]$exogenous_rRNA))){
      exprs.exogenousRibosomal_specific[match(sample.data[[i]]$exogenous_rRNA$ID, rownames(exprs.exogenousRibosomal_specific)),i] = as.numeric(sample.data[[i]]$exogenous_rRNA$readCount_direct)
      exprs.exogenousRibosomal_cumulative[match(sample.data[[i]]$exogenous_rRNA$ID, rownames(exprs.exogenousRibosomal_cumulative)),i] = as.numeric(sample.data[[i]]$exogenous_rRNA$readCount_inherited)
    }
    ## Exogenous Genomes
    if(!is.null(nrow(sample.data[[i]]$exogenous_genomes))){
      exprs.exogenousGenomes_specific[match(sample.data[[i]]$exogenous_genomes$ID, rownames(exprs.exogenousGenomes_specific)),i] = as.numeric(sample.data[[i]]$exogenous_genomes$readCount_direct)
      exprs.exogenousGenomes_cumulative[match(sample.data[[i]]$exogenous_genomes$ID, rownames(exprs.exogenousGenomes_cumulative)),i] = as.numeric(sample.data[[i]]$exogenous_genomes$readCount_inherited)
    }
  }
  
  
  ##
  ## Calculate the total number of mapped reads to the rRNA, genome, and exogenous sequences
  ##
  mapping.stats[is.na(mapping.stats)] = 0
  mapping.stats = as.data.frame(mapping.stats)
  libSizes = list()
  libSizes$input = mapping.stats[,colnames(mapping.stats) %in% c("input")]
  libSizes$successfully_clipped = mapping.stats[,colnames(mapping.stats) %in% c("successfully_clipped")]
  libSizes$reads_used_for_alignment = mapping.stats[,colnames(mapping.stats) %in% c("reads_used_for_alignment")]
  libSizes$all = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome","miRNA_exogenous_sense")])
  libSizes$endogenous = rowSums(mapping.stats[,colnames(mapping.stats) %in% c("rRNA","genome")])
  libSizes$genome = mapping.stats[,colnames(mapping.stats) %in% "genome"]
  libSizes$smRNA = mapping.stats[,grep("sense",colnames(mapping.stats))]
  libSizes$miRNA = colSums(exprs.miRNA)
  libSizes$exogenous_miRNA = colSums(exprs.exogenous_miRNA)
  libSizes$exogenous_rRNA = exprs.exogenousRibosomal_cumulative[rownames(exprs.exogenousRibosomal_cumulative)=="1",]
  libSizes$exogenous_genomes = exprs.exogenousGenomes_cumulative[rownames(exprs.exogenousGenomes_cumulative)=="1",]
  
  
  ##
  ## Save the raw count data
  ##
  printMessage("Saving raw data to disk")
  #save(exprs.miRNA, exprs.tRNA, exprs.piRNA, exprs.gencode, exprs.circRNA, exprs.exogenous_miRNA, exprs.exogenous_genomes, mapping.stats, libSizes, read.lengths, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadCounts.RData", sep="/"))
  save(exprs.miRNA, exprs.tRNA, exprs.piRNA, exprs.gencode, exprs.circRNA, exprs.exogenous_miRNA, exprs.exogenousRibosomal_specific, exprs.exogenousRibosomal_cumulative, taxonomyInfo.exogenous_rRNA, exprs.exogenousGenomes_specific, exprs.exogenousGenomes_cumulative, taxonomyInfo.exogenous_genomes, mapping.stats, qc.results, libSizes, read.lengths, run.duration, exprs.calibrator, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadCounts.RData", sep="/"))
  
  if(nrow(exprs.calibrator) > 0)
    write.table(exprs.calibrator, file=paste(output.dir, "exceRpt_CALIBRATOR_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  if(nrow(exprs.miRNA) > 0)
    write.table(exprs.miRNA, file=paste(output.dir, "exceRpt_miRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.tRNA) > 0)
    write.table(exprs.tRNA, file=paste(output.dir, "exceRpt_tRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.piRNA) > 0)
    write.table(exprs.piRNA, file=paste(output.dir, "exceRpt_piRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.gencode) > 0)
    write.table(exprs.gencode, file=paste(output.dir, "exceRpt_gencode_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.circRNA) > 0)
    write.table(exprs.circRNA, file=paste(output.dir, "exceRpt_circularRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.exogenous_miRNA) > 0)
    write.table(exprs.exogenous_miRNA, file=paste(output.dir, "exceRpt_exogenous_miRNA_ReadCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.exogenousRibosomal_specific) > 0){
    tmp = cbind(taxonomyInfo.exogenous_rRNA[match(rownames(exprs.exogenousRibosomal_specific), taxonomyInfo.exogenous_rRNA$ID), ], exprs.exogenousRibosomal_specific)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousRibosomal_taxonomySpecific_ReadCounts.txt", sep="/"), sep="\t", row.names=F, quote=F)
    
    tmp = cbind(taxonomyInfo.exogenous_rRNA[match(rownames(exprs.exogenousRibosomal_cumulative), taxonomyInfo.exogenous_rRNA$ID), ], exprs.exogenousRibosomal_cumulative)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousRibosomal_taxonomyCumulative_ReadCounts.txt", sep="/"), sep="\t", row.names=F, quote=F)
  }
  
  if(nrow(exprs.exogenousGenomes_specific) > 0){
    tmp = cbind(taxonomyInfo.exogenous_genomes[match(rownames(exprs.exogenousGenomes_specific), taxonomyInfo.exogenous_genomes$ID), ], exprs.exogenousGenomes_specific)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousGenomes_taxonomySpecific_ReadCounts.txt", sep="/"), sep="\t", row.names=F, quote=F)
    
    tmp = cbind(taxonomyInfo.exogenous_genomes[match(rownames(exprs.exogenousGenomes_cumulative), taxonomyInfo.exogenous_genomes$ID), ], exprs.exogenousGenomes_cumulative)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousGenomes_taxonomyCumulative_ReadCounts.txt", sep="/"), sep="\t", row.names=F, quote=F)
  }
  
  write.table(read.lengths, file=paste(output.dir, "exceRpt_ReadLengths.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ##
  ## Keep a record of the adapter sequences for QC
  ##
  adapterSeq = unlist(lapply(sample.data, function(l){ l$adapterSeq }))
  write.table(as.data.frame(adapterSeq), file=paste(output.dir, "exceRpt_adapterSequences.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ##
  ## Write the numbers of reads mapping at each stage and the QC results
  ##
  write.table(mapping.stats, file=paste(output.dir,"exceRpt_readMappingSummary.txt",sep="/"), sep="\t", col.names=NA, quote=F)
  write.table(qc.results, file=paste(output.dir,"exceRpt_QCresults.txt",sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ##
  ## Calculate reads per million (RPM)
  ##
  printMessage("Normalising to RPM")
  #libSize.use = libSizes$all
  #libSize.use = libSizes$miRNA
  libSize.use = libSizes$genome
  exprs.miRNA.rpm = t(10^6 * t(exprs.miRNA) / libSize.use)
  exprs.tRNA.rpm = t(10^6 * t(exprs.tRNA) / libSize.use)
  exprs.piRNA.rpm = t(10^6 * t(exprs.piRNA) / libSize.use)
  exprs.gencode.rpm = t(10^6 * t(exprs.gencode) / libSize.use)
  exprs.circRNA.rpm = t(10^6 * t(exprs.circRNA) / libSize.use)
  exprs.exogenous_miRNA.rpm = t(10^6 * t(exprs.exogenous_miRNA) / libSizes$exogenous_miRNA)
  
  exprs.exogenousRibosomal_specific.rpm = exprs.exogenousRibosomal_specific
  exprs.exogenousRibosomal_cumulative.rpm = exprs.exogenousRibosomal_cumulative
  if(nrow(exprs.exogenousRibosomal_specific) > 0){
    exprs.exogenousRibosomal_specific.rpm = t(10^6 * t(exprs.exogenousRibosomal_specific) / libSizes$exogenous_rRNA)
    exprs.exogenousRibosomal_cumulative.rpm = t(10^6 * t(exprs.exogenousRibosomal_cumulative) / libSizes$exogenous_rRNA)
  }
  
  exprs.exogenousGenomes_specific.rpm = exprs.exogenousGenomes_specific
  exprs.exogenousGenomes_cumulative.rpm = exprs.exogenousGenomes_cumulative
  if(nrow(exprs.exogenousGenomes_specific) > 0){
    exprs.exogenousGenomes_specific.rpm = t(10^6 * t(exprs.exogenousGenomes_specific) / libSizes$exogenous_genomes)
    exprs.exogenousGenomes_cumulative.rpm = t(10^6 * t(exprs.exogenousGenomes_cumulative) / libSizes$exogenous_genomes)
  }
  
  
  ##
  ## Save the RPM normalised data
  ##
  printMessage("Saving normalised data to disk")
  save(exprs.miRNA.rpm, exprs.tRNA.rpm, exprs.piRNA.rpm, exprs.gencode.rpm, exprs.circRNA.rpm, exprs.exogenous_miRNA.rpm, exprs.exogenousRibosomal_specific.rpm, exprs.exogenousRibosomal_cumulative.rpm, exprs.exogenousGenomes_specific.rpm, exprs.exogenousGenomes_cumulative.rpm, file=paste(output.dir, "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep="/"))
  
  if(nrow(exprs.miRNA.rpm) > 0)
    write.table(exprs.miRNA.rpm, file=paste(output.dir, "exceRpt_miRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  if(nrow(exprs.tRNA.rpm) > 0)  
    write.table(exprs.tRNA.rpm, file=paste(output.dir, "exceRpt_tRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  if(nrow(exprs.piRNA.rpm) > 0)
    write.table(exprs.piRNA.rpm, file=paste(output.dir, "exceRpt_piRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  if(nrow(exprs.gencode.rpm) > 0)
    write.table(exprs.gencode.rpm, file=paste(output.dir, "exceRpt_gencode_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  if(nrow(exprs.circRNA.rpm) > 0)
    write.table(exprs.circRNA.rpm, file=paste(output.dir, "exceRpt_circularRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.exogenous_miRNA.rpm) > 0)
    write.table(exprs.exogenous_miRNA.rpm, file=paste(output.dir, "exceRpt_exogenous_miRNA_ReadsPerMillion.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  if(nrow(exprs.exogenousRibosomal_specific) > 0){
    tmp = cbind(taxonomyInfo.exogenous_rRNA[match(rownames(exprs.exogenousRibosomal_specific.rpm), taxonomyInfo.exogenous_rRNA$ID), ], exprs.exogenousRibosomal_specific.rpm)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousRibosomal_taxonomySpecific_ReadsPerMillion.txt", sep="/"), sep="\t", row.names=F, quote=F)
    
    tmp = cbind(taxonomyInfo.exogenous_rRNA[match(rownames(exprs.exogenousRibosomal_cumulative.rpm), taxonomyInfo.exogenous_rRNA$ID), ], exprs.exogenousRibosomal_cumulative.rpm)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousRibosomal_taxonomyCumulative_ReadsPerMillion.txt", sep="/"), sep="\t", row.names=F, quote=F)
  }
  
  if(nrow(exprs.exogenousGenomes_specific) > 0){
    tmp = cbind(taxonomyInfo.exogenous_genomes[match(rownames(exprs.exogenousGenomes_specific.rpm), taxonomyInfo.exogenous_genomes$ID), ], exprs.exogenousGenomes_specific.rpm)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousGenomes_taxonomySpecific_ReadsPerMillion.txt", sep="/"), sep="\t", row.names=F, quote=F)
    tmp = cbind(taxonomyInfo.exogenous_genomes[match(rownames(exprs.exogenousGenomes_cumulative.rpm), taxonomyInfo.exogenous_genomes$ID), ], exprs.exogenousGenomes_cumulative.rpm)
    write.table(tmp, file=paste(output.dir, "exceRpt_exogenousGenomes_taxonomyCumulative_ReadsPerMillion.txt", sep="/"), sep="\t", row.names=F, quote=F)
  }
  
  return(rownames(mapping.stats))
}


##' @title delete "-e " strings
##' @description delete "-e " strings in *.stats processed output document to create the required dataframe
##' @description for downstream processing
##' @param linesRead <string> vector output from readLines()
##' @author Miquel Anglada Girotto
delete_e = function(linesRead){
  f = gsub('-e ','',linesRead)
  return(f)
}


##' @title print info messages
##' @description to print formatted info messages.
##' @param message <string>
##' @author (Rob Kitchen)
printMessage = function(message=""){
  cat(as.character(Sys.time()),":  ",paste(message,sep=""),"\n",sep="")
}


##' @title Make QC plots
##' @description Create QC plots from data obtained through readData()
##' @param sampleIDs <string> names of the samples obtained through readData()
##' @param output.dir <string> directory to which output the files to create the report.
##' @param sampleGroups <string> if any sample groups were defined. Deprecated.
##' @param minPercent_exogenousRibosomal <double> Not used.
##' @param minPercent_exogenousGenomes <double> Not used.
##' @author (Rob Kitchen) Miquel Anglada Girotto
PlotData = function(sampleIDs, output.dir, sampleGroups=NA, minPercent_exogenousRibosomal=0.5, minPercent_exogenousGenomes=0.5){
  
  load(paste(output.dir, "exceRpt_smallRNAQuants_ReadCounts.RData", sep="/"))
  
  load(paste(output.dir, "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep="/"))
  
  ##
  ## Order samples based on similarity of mapping statistics
  ##
  sampleOrder = 1
  if(nrow(mapping.stats) > 1){
    h = hclust(dist(1-cor(t(mapping.stats))))
    sampleOrder = h$order
  }
  
  ##
  ## Create list where to save plots that can be returned
  ##
  plotsList = list()
  
  ##
  ## Open PDF for diagnostic plots
  ##
  printMessage("Creating QC plots")
  pdf(paste(output.dir,"exceRpt_DiagnosticPlots.pdf",sep="/"), height=10, width=20)
  
  
  if(ncol(read.lengths) > 1){
    printMessage("Plotting read-length distributions")
    ##
    ## plot distribution of clipped read lengths - read count
    ##
    tmp = melt(read.lengths); colnames(tmp) = c("sample","length","count")
    if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sample, sampleGroups$sampleID), 2] }
    maxX = min(c(100,max(tmp$length)))
    p = ggplot(tmp, aes(x=length, y=count, color=sample)) + geom_line(alpha=1) +
      xlab("read length (nt)") + ylab("# reads") + 
      ggtitle("read-length distributions:\nraw read count")+
      scale_x_continuous(limits=c(15,maxX), minor_breaks=1:maxX, breaks=seq(15,maxX,by=5))+
      scale_colour_manual(values=getSampleColors(tmp$sample))
    if(nrow(read.lengths) > 30){ p = p +guides(colour=FALSE) }
    if(is.data.frame(sampleGroups)){ p = p +facet_wrap(~sampleGroup,ncol=1)}
    print(p)
    
    
    # save
    plotsList[["read-length distributions: raw read count"]] = p
    
    ##
    ## plot distribution of clipped read lengths - fraction
    ##
    tmp = melt(t(apply(read.lengths, 1, function(row){ row/sum(row) }))); colnames(tmp) = c("sample","length","fraction")
    if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sample, sampleGroups$sampleID), 2] }
    maxX = min(c(100,max(tmp$length)))
    p = ggplot(tmp, aes(x=length, y=fraction, colour=sample)) +geom_line(alpha=1) +xlab("read length (nt)") +ylab("fraction of reads") +ggtitle("read-length distributions:\nnormalised read fraction") +scale_x_continuous(limits=c(15,maxX), minor_breaks=1:maxX, breaks=seq(15,maxX,by=5))+
      scale_colour_manual(values=getSampleColors(tmp$sample))
    if(nrow(read.lengths) > 30){ p = p +guides(colour=FALSE) }
    if(is.data.frame(sampleGroups)){ p = p +facet_wrap(~sampleGroup,ncol=1)}
    print(p)
    
    # save
    plotsList[["read-length distributions: normalised read fraction"]] = p
  }
  
  
  
  ##
  ## Plot run duration of each sample
  ##
  printMessage("Plotting run-duration")
  tmp=reshape2::melt(as.matrix(run.duration))
  colnames(tmp) = c("sampleID","stuff","runDuration_seconds")
  tmp = cbind(tmp, category=.bincode(tmp[,3], breaks=c(0,as.numeric(quantile(tmp[,3],probs=c(0.10,0.90,1))))))
  tmp = cbind(tmp, colour=tmp$category)
  tmp = cbind(tmp, inputReadCount=mapping.stats$input)
  tmp$category[tmp$category == 1] = "fast"
  tmp$category[tmp$category == 2] = "normal"
  tmp$category[tmp$category == 3] = "slow"
  tmp$colour = getSampleColors(tmp$colour)
  tmp$runDuration_minutes = tmp$runDuration_seconds/60
  tmp$runDuration_hours = tmp$runDuration_minutes/60
  p = ggplot(tmp, aes(x=sampleID,y=runDuration_hours)) + geom_bar(stat="identity",fill=tmp$colour) +
    facet_grid(~category,scales="free_x",space="free_x") +guides(fill=FALSE) +
    theme(axis.text.x=element_text(angle=60, hjust=1.0, vjust=1)) +
    ggtitle("Duration of exceRpt run for each sample") +ylab("Run duration (hours)")
  print(p)
  # save
  plotsList[["Duration of exceRpt run for each sample"]] = p
  
  if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sampleID, sampleGroups$sampleID), 2] }
  p = ggplot(tmp, aes(x=inputReadCount,y=runDuration_hours,colour=colour)) + geom_point(size=5) + guides(colour=FALSE) +
    scale_x_log10(limits=c(min(c(100000,10^floor(log10(min(tmp$inputReadCount+1))))),10^ceiling(log10(max(tmp$inputReadCount)))), breaks=10^seq(min(c(100000,floor(log10(min(tmp$inputReadCount+1))))),ceiling(log10(max(tmp$inputReadCount))))) +
    ggtitle("Duration of exceRpt run per sequencing yield") + ylab("Run duration (hours)") + xlab("Total number of reads input") + 
    scale_colour_manual(values=tmp$colour)
  print(p)
  
  # save
  plotsList[["Duration of exceRpt run per sequencing yield"]] = p
  
  ##
  ## plot distribution of # mapped reads per sample
  ##
  printMessage("Plotting # mapped reads")
  tmp = log10(libSizes$all)
  ## ORIGINAL: hist(tmp, breaks=seq(0,ceiling(max(tmp)), by=0.1), col="grey", border="white", xlab="log10(# mapped reads)", main="Library size (all mapped reads)", ylab="# samples")
  p = ggplot(as.data.frame(tmp),aes(x=tmp))+geom_histogram(breaks=seq(0,ceiling(max(tmp)), by=0.1),fill='lightblue')+xlab("log10(# mapped reads)")+ylab("# samples")+ggtitle("Library size (all mapped reads)")
  
  # save
  plotsList[["Library size (all mapped reads)"]] = p
  
  ##
  ## Plot the rRNA contamination
  ##
  mapping.stats.orig = mapping.stats
  mapping.stats = mapping.stats[,-grep("input_to_",colnames(mapping.stats))]
  
  ## remove the exogenous stuff from the stats if this wasn't used in the run
  if(sum(mapping.stats[,23:27]) == 0)
    mapping.stats = mapping.stats[, -c(23:27)]
  
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  printMessage("Plotting mapping stats heatmap (1/3)")
  toplot = melt(as.matrix(mapping.stats / mapping.stats$input)); colnames(toplot) = c("Sample","Stage","ReadFraction")
  toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
  toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
  if(is.data.frame(sampleGroups)){ toplot$sampleGroup = sampleGroups[match(toplot$Sample, sampleGroups$sampleID), 2] }
  p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=40, hjust=1.0, vjust=1)) +ggtitle("fraction aligned reads\n(normalised by # input reads)")
  if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
  print(p)
  
  # save
  plotsList[["fraction aligned reads (normalised by # input reads)"]] = p
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  printMessage("Plotting mapping stats heatmap (2/3)")
  if(max(mapping.stats$successfully_clipped) > 0){
    tmp = mapping.stats
    i.toFix = which(tmp$successfully_clipped == 0)
    if(length(i.toFix) > 0)
      tmp$successfully_clipped[i.toFix] = tmp$input[i.toFix]
    toplot = melt(as.matrix(tmp / tmp$successfully_clipped)[,-1,drop=F]); colnames(toplot) = c("Sample","Stage","ReadFraction")
    toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
    toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
    if(is.data.frame(sampleGroups)){ toplot$sampleGroup = sampleGroups[match(toplot$Sample, sampleGroups$sampleID), 2] }
    p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=40, hjust=1.0, vjust=1)) +ggtitle("fraction aligned reads\n(normalised by # adapter-clipped reads)")
    if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
    if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
    print(p)
    
    # save
    plotsList[["fraction aligned reads (normalised by # adapter-clipped reads)"]] = p
  }
  
  ##
  ## Plot heatmap of mapping percentages through the pipeline
  ##
  printMessage("Plotting mapping stats heatmap (3/3)")
  toplot = melt(as.matrix(mapping.stats / mapping.stats$reads_used_for_alignment)[,-c(1:7),drop=F]); colnames(toplot) = c("Sample","Stage","ReadFraction")
  toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
  toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
  if(is.data.frame(sampleGroups)){ toplot$sampleGroup = sampleGroups[match(toplot$Sample, sampleGroups$sampleID), 2] }
  p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=40, hjust=1.0, vjust=1)) +ggtitle("fraction aligned reads\n(normalised by # non-contaminant reads)")
  if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
  print(p)
  
  # save
  plotsList[["fraction aligned reads (normalised by # non-contaminant reads)"]] = p
  
  ##
  ## Plot QC results
  ##
  printMessage("Plotting QC result")
  toplot = as.data.frame(qc.results)
  toplot$Sample = factor(rownames(toplot), levels=rownames(mapping.stats)[sampleOrder])
  p = ggplot(toplot, aes(x=TranscriptomeReads, y=TranscriptomeGenomeRatio)) +
    geom_text(label=toplot$Sample,hjust = 0,nudge_x = 0.02, nudge_y = 0.02, colour='black') +
    scale_colour_manual(values=getSampleColors(toplot$Sample))
  if(is.data.frame(sampleGroups)){ 
    toplot$sampleGroup = sampleGroups[match(toplot$Sample, sampleGroups$sampleID), 2] 
    p = ggplot(toplot, aes(x=TranscriptomeReads, y=TranscriptomeGenomeRatio, colour=sampleGroup)) +
      geom_text(label=toplot$Sample,hjust = 0, nudge_x = 0.02, nudge_y = 0.02, colour='black') +
      scale_colour_manual(values=getSampleColors(toplot$sampleGroup))
  }
  minX = floor(log10(min(toplot$TranscriptomeReads)+0.001))
  maxX = ceiling(log10(max(toplot$TranscriptomeReads)+0.001))
  p = p + scale_x_log10(breaks=10^c(minX:maxX)) +coord_cartesian(xlim=c(10^(minX),10^(maxX)),ylim=c(0,1)) +geom_vline(xintercept=100000,col="red",alpha=0.5) +geom_hline(yintercept=0.5,col="red",alpha=0.5) +annotate("rect",xmin=0,xmax=Inf,ymin=-1,ymax=0.5,alpha=0.2,fill="red") +annotate("rect",xmin=0,xmax=100000,ymin=-1,ymax=1.1,alpha=0.2,fill="red") +ylab("# transcriptome reads / # genome reads") +xlab("# transcriptome reads (log10)") +ggtitle("QC result: overall")
  
  print(p+geom_point(size=4))
  
  # save
  plotsList[["QC result: overall"]] = p + geom_point(size=4)
  
  if(is.data.frame(sampleGroups)){
    print(p + facet_wrap(~sampleGroup) + theme(legend.position="none") +geom_point(size=2)) 
    
    # save
    plotsList[["QC result: overall"]] = p +facet_wrap(~sampleGroup) +theme(legend.position="none") +
      geom_point(size=2)
    
  }
  
  
  
  ##
  ## Heatmap
  ##
  qc.results[,4] = round(qc.results[,4]*100)/100
  qc.results[,5] = round(qc.results[,5]*1000)/1000
  tmp.mat = qc.results
  tmp.mat[,1] = rep(1,nrow(tmp.mat))
  tmp.mat[,2] = rep(1,nrow(tmp.mat))
  tmp.mat[,5] = rep(1,nrow(tmp.mat))
  tmp.pass=tmp.mat[,3] >= 100000;  tmp.mat[tmp.pass,3] = "pass"; tmp.mat[!tmp.pass,3] = "fail"
  tmp.pass=tmp.mat[,4] >= 0.5;  tmp.mat[tmp.pass,4] = "pass"; tmp.mat[!tmp.pass,4] = "fail"
  
  toplot=cbind(melt(tmp.mat), Actual=melt(qc.results)[,3]); colnames(toplot)[1:3]=c("Sample","Stage","Value")
  toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
  if(is.data.frame(sampleGroups)){ toplot$sampleGroup = sampleGroups[match(toplot$Sample, sampleGroups$sampleID), 2] }
  p = ggplot(toplot, aes(y=Sample, x=Stage, fill=Value, label=Actual)) +
    scale_fill_manual(values=c("fail"="red","pass"="palegreen","1"="lightgrey")) +
    geom_label() +theme(plot.background=element_rect(fill="white"),panel.background=element_rect(fill=rgb(0.97,0.97,0.97)), axis.text.x=element_text(angle=20, hjust=1, vjust=1)) +
    ggtitle("QC result: per-sample results") +xlab("") +ylab("")
  print(p)
  
  # save
  plotsList[["QC result: per-sample results"]] = p
  
  ##
  ## Plot breakdown of counts in each biotype 
  ##
  printMessage("Plotting biotype counts")
  require(plyr)
  sampleTotals=matrix(NA,ncol=nrow(mapping.stats),nrow=0); colnames(sampleTotals) = rownames(mapping.stats)
  if(nrow(exprs.miRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.miRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "miRNA"
  }
  if(nrow(exprs.tRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.tRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "tRNA"
  }
  if(nrow(exprs.piRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.piRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "piRNA"
  }
  if(nrow(exprs.gencode) > 0){
    tmp = data.frame(biotype=sapply(rownames(exprs.gencode), function(id){bits=unlist(strsplit(id,":")); bits[length(bits)]}), exprs.gencode)
    tmp = ddply(tmp, "biotype", function(mat){ colSums(mat[,-1,drop=F]) })
    rownames(tmp) = tmp[,1]; tmp = tmp[,-1,drop=F]
    colnames(tmp) = colnames(sampleTotals)
    ## add gencode miRNA to the existing count...
    if("miRNA" %in% rownames(tmp)  &&  "miRNA" %in% rownames(sampleTotals)){
      i = which(rownames(tmp) %in% "miRNA")
      j = which(rownames(sampleTotals) %in% "miRNA")
      for(x in 1:ncol(sampleTotals)){
        sampleTotals[j,x] = sampleTotals[j,x] + tmp[i,x]
      }
      tmp = tmp[-i,,drop=F]
    }
    sampleTotals = rbind(sampleTotals, tmp)
  }
  if(nrow(exprs.circRNA) > 0){
    sampleTotals = rbind(sampleTotals, colSums(exprs.circRNA))
    rownames(sampleTotals)[nrow(sampleTotals)] = "circularRNA"
  }
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_miRNA)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_miRNA"
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_rRNA)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_rRNA"
  sampleTotals = rbind(sampleTotals, mapping.stats$exogenous_genomes)
  rownames(sampleTotals)[nrow(sampleTotals)] = "exogenous_genomes"
  
  sampleTotals = sampleTotals[order(apply(sampleTotals, 1, median, na.rm=T), decreasing=F), ,drop=F]
  tmp = melt(as.matrix(sampleTotals))
  colnames(tmp) = c("biotype","sampleID","readCount")
  if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sampleID, sampleGroups$sampleID), 2] }
  p = ggplot(na.omit(tmp), aes(y=readCount,x=biotype, colour=biotype)) +geom_hline(aes(yintercept=1),linetype="dashed") +
    geom_boxplot() +scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000)) +
    guides(colour=FALSE) +coord_flip() +ggtitle("Biotypes: distributions, raw read-counts")+
    scale_color_manual(values = getSampleColors(tmp$biotype))
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x")}
  print(p)
  
  # save
  plotsList[["Biotypes: distributions, raw read-counts"]] = p
  
  ## save the biotype counts
  write.table(sampleTotals[order(apply(sampleTotals, 1, median, na.rm=T), decreasing=T), ,drop=F], file=paste(output.dir, "exceRpt_biotypeCounts.txt", sep="/"), sep="\t", col.names=NA, quote=F)
  
  
  ## plot biotype breakdown as RPM:
  tmp = melt(as.matrix(apply(sampleTotals, 2, function(col){ col*1000000/sum(col) })))
  colnames(tmp) = c("biotype","sampleID","readPerMillion")
  if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sampleID, sampleGroups$sampleID), 2] }
  p = ggplot(na.omit(tmp), aes(y=readPerMillion,x=biotype, colour=biotype)) +geom_hline(aes(yintercept=1),linetype="dashed") +geom_boxplot() +scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000)) +guides(colour=FALSE) +coord_flip() +ggtitle("Biotypes: distributions, normalised")+
    scale_color_manual(values = getSampleColors(tmp$biotype))
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x")}
  print(p)
  
  # save
  plotsList[["Biotypes: distributions, normalised"]] = p
  
  ## plot top N biotypes for each sample as a barplot - normalised to INPUT reads
  tmp = sampleTotals
  for(i in 1:ncol(sampleTotals))
    tmp[,i] = sampleTotals[,i]*1000000/libSizes$reads_used_for_alignment[i]
  biotypeOrder = order(apply(tmp, 1, mean), decreasing=T)
  tmp = tmp[biotypeOrder, , drop=F]
  
  N = min(7, length(biotypeOrder)) # number of biotypes detected
  
  tmp = as.matrix(rbind(tmp[1:N, , drop=F], other=colSums(tmp[-c(1:N), , drop=F])))
  tmp = melt(rbind(tmp, unmapped=1000000-colSums(tmp)))
  colnames(tmp) = c("biotype","sampleID","readsPerMillion")
  if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sampleID, sampleGroups$sampleID), 2] }
  p = ggplot(na.omit(tmp), aes(y=readsPerMillion,x=sampleID,fill=biotype)) +geom_bar(stat="identity") +scale_fill_brewer(palette = "Paired") +theme(axis.text.x=element_text(angle=50, hjust=1.0, vjust=1)) +ggtitle("Biotypes: per-sample, normalised") +ylab("reads per million reads used for alignment") +xlab("") +ylim(limits=c(0,1E6)) #+scale_fill_discrete(rich.colors(10))
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
  suppressWarnings(print(p))
  
  # save
  plotsList[["Biotypes: per-sample, normalised"]] = p
  
  ## plot top N biotypes for each sample as a barplot - normalised to MAPPED reads
  tmp = as.matrix(apply(sampleTotals, 2, function(col){ col*1000000/sum(col) }))
  tmp = tmp[biotypeOrder, , drop=F]
  tmp = as.matrix(rbind(tmp[1:N, , drop=F], other=colSums(tmp[-c(1:N), , drop=F])))
  tmp = melt(tmp)
  
  colnames(tmp) = c("biotype","sampleID","readsPerMillion")
  if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sampleID, sampleGroups$sampleID), 2] }
  p = ggplot(na.omit(tmp), aes(y=readsPerMillion,x=sampleID,fill=biotype)) +geom_bar(stat="identity") +theme(axis.text.x=element_text(angle=50, hjust=1.0, vjust=1)) +ggtitle("Biotypes: per-sample, normalised") +ylab("reads per million mapped reads") +xlab("") +scale_fill_brewer(palette = "Paired") #+scale_fill_manual( values = c(colorRampPalette( brewer.pal( 6 , "Paired" ) )(8), "grey") )
  if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
  suppressWarnings(print(p))    
  
  # save
  plotsList[["Biotypes: per-sample, normalised"]] = p
  
  ## Plot miRNA expression distributions
  if(nrow(exprs.miRNA) > 0){
    printMessage("Plotting miRNA expression distributions")
    tmp = melt(exprs.miRNA)
    colnames(tmp) = c("miRNA","sample","abundance")
    if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sample, sampleGroups$sampleID), 2] }
    p = ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +
      geom_violin() +
      geom_boxplot(alpha=0.2) +ylab("Read count") +ggtitle("miRNA abundance distributions (raw counts)") +
      scale_y_log10() +guides(colour=FALSE) +
      scale_color_manual(values = setNames(getSampleColors(tmp$sample),tmp$sample))
    if(ncol(exprs.miRNA) < 30){ 
      p = p +theme(axis.text.x=element_text(angle=50, hjust=1.0, vjust=1))
    }else{
      p = p+theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    }
    if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
    print(p)
    
    # save
    plotsList[["miRNA abundance distributions (raw counts)"]] = p
    
    p = ggplot(tmp, aes(x=abundance, colour=sample)) +geom_density() +xlab("Read count") +ggtitle("miRNA abundance distributions (raw counts)") +scale_x_log10()+
      scale_color_manual(values = setNames(getSampleColors(tmp$sample),tmp$sample))
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup)}
    print(p)
    
    # save
    plotsList[["miRNA abundance distributions (raw counts)"]] = p
    
    tmp = melt(exprs.miRNA.rpm)
    colnames(tmp) = c("miRNA","sample","abundance")
    if(is.data.frame(sampleGroups)){ tmp$sampleGroup = sampleGroups[match(tmp$sample, sampleGroups$sampleID), 2] }
    p = ggplot(tmp, aes(y=abundance, x=sample, colour=sample)) +geom_violin() +
      geom_boxplot(alpha=0.2) +ylab("Reads per million (RPM)") +
      ggtitle("miRNA abundance distributions (RPM)") +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
      scale_y_log10() +guides(colour=FALSE)+
      scale_color_manual(values = setNames(getSampleColors(tmp$sample),tmp$sample))
    if(ncol(exprs.miRNA.rpm) < 30){ 
      p = p +theme(axis.text.x=element_text(angle=50, hjust=1.0, vjust=1))
    }else{
      p = p+theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    }
    if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup, scales="free_x",space="free_x")}
    print(p)
    
    # save
    plotsList[["miRNA abundance distributions (RPM)"]] = p
    
    p = ggplot(tmp, aes(x=abundance, colour=sample)) +geom_density() +
      xlab("Reads per million (RPM)") +ggtitle("miRNA abundance distributions (RPM)") +scale_x_log10()+
      scale_color_manual(values = setNames(getSampleColors(tmp$sample),tmp$sample))
    if(ncol(exprs.miRNA.rpm) > 30){ p = p +guides(colour=FALSE) }
    if(is.data.frame(sampleGroups)){ p = p +facet_grid(~sampleGroup)}
    print(p)
    
    # save
    plotsList[["miRNA abundance distributions (RPM)"]] = p
  }
  
  
  
  ##
  ## Finally, plot exogenous if there are any (not used by default)
  ##
  if(nrow(exprs.exogenousGenomes_specific) > 0){
    printMessage("Plotting exogenous counts")
    par(oma=c(20,2,0,0))
    barplot(exprs.exogenousGenomes_cumulative[1,,drop=F], las=2, main="Total # reads mapped to NCBI taxonomy")
    
    ## if we have more than one sample, plot some heatmaps
    if(ncol(exprs.exogenousGenomes_specific) > 1){
      
      par(oma=c(8,0,0,20))
      maxRow = 50; if(nrow(exprs.exogenousGenomes_specific) < maxRow){ maxRow = nrow(exprs.exogenousGenomes_specific) }
      tmp.order = order(apply(t(t(exprs.exogenousGenomes_specific)/colSums(exprs.exogenousGenomes_specific)), 1, median), decreasing=T)
      tmp = t(log10(t(t(exprs.exogenousGenomes_specific)*1000000/colSums(exprs.exogenousGenomes_specific))[tmp.order, ][1:maxRow,]+0.1))
      colnames(tmp) = taxonomyInfo.exogenous_genomes[match(colnames(tmp), taxonomyInfo.exogenous_genomes$ID), ]$name
      heatmap.2(tmp,trace="none",main="top taxa nodes: specific normalised read count", symbreaks=F,col=rich.colors(50))
      
      tmp.order = order(apply(t(t(exprs.exogenousGenomes_specific)), 1, median), decreasing=T)
      tmp = t(log10(exprs.exogenousGenomes_specific[tmp.order, ][1:maxRow,]+0.1))
      colnames(tmp) = taxonomyInfo.exogenous_genomes[match(colnames(tmp), taxonomyInfo.exogenous_genomes$ID), ]$name
      heatmap.2(tmp,trace="none",main="top taxa nodes: specific absolute read count", symbreaks=F,col=rich.colors(50))
      
      maxRow = 50; if(nrow(exprs.exogenousGenomes_cumulative) < maxRow){ maxRow = nrow(exprs.exogenousGenomes_cumulative) }
      tmp.order = order(apply(t(t(exprs.exogenousGenomes_cumulative)/libSizes$exogenous_genomes), 1, median), decreasing=T)
      tmp = t(log10(t(t(exprs.exogenousGenomes_cumulative)*1000000/libSizes$exogenous_genomes)[tmp.order, ][1:maxRow,]+0.1))
      colnames(tmp) = taxonomyInfo.exogenous_genomes[match(colnames(tmp), taxonomyInfo.exogenous_genomes$ID), ]$name
      heatmap.2(tmp,trace="none",main="top taxa nodes: cumulative normalised read count", symbreaks=F,col=rich.colors(50))
      
      tmp.order = order(apply(t(t(exprs.exogenousGenomes_cumulative)), 1, median), decreasing=T)
      tmp = t(log10(exprs.exogenousGenomes_cumulative[tmp.order, ][1:maxRow,]+0.1))
      colnames(tmp) = taxonomyInfo.exogenous_genomes[match(colnames(tmp), taxonomyInfo.exogenous_genomes$ID), ]$name
      heatmap.2(tmp,trace="none",main="top taxa nodes: cumulative absolute read count", symbreaks=F,col=rich.colors(50))
    }
  }
  dev.off()
  
  printMessage("All done!")
  return(plotsList)
}


## Deprecated below:
##
## Plots a taxonomy tree with a given set of weights
## (Rob Kitchen)
plotTree = function(rEG, taxonomyInfo, counts_uniq, counts_cum, title="", what){

## node parameters
nNodes = length(nodes(rEG))
nA <- list()
nA$shape = rep("circle",nNodes)
nA$fixedSize<-rep(FALSE, nNodes)
nA$height <- nA$width <- rescale(sqrt(counts_cum/10), to=c(0.25,7))
nA$color <- rep(rgb(0,0,0,0.25),nNodes)
nA$style <- rep("bold", nNodes)
if(what == "exogenousRibosomal"){
  nA$fillcolor <- sapply(counts_uniq*10, function(val){ if(val>100){val=100}; rgb(100-val,100,100-val,maxColorValue=100)})
}else{
  nA$fillcolor <- sapply(counts_uniq*10, function(val){ if(val>100){val=100}; rgb(100-val,100-val,100,maxColorValue=100)})
}

newNodeIDs = sapply(taxonomyInfo[match(as.numeric(nodes(rEG)), taxonomyInfo$ID), ]$name, function(id){ newID=unlist(strsplit(id," ")); if(length(newID) == 1){id}else{paste(newID[1], "\n", paste(newID[-1],collapse=" "), sep="") }})
nA$label <- paste(newNodeIDs,"\n",round(counts_cum*10)/10,"%",sep="")
nA <- lapply(nA, function(x) { names(x) <- nodes(rEG); x})

## edge parameters
eA <- list(arrowsize=rep(0.1,length(names(rEG@edgeData))), arrowhead=rep("none",length(names(rEG@edgeData))))
eA <- lapply(eA, function(x) { names(x) <- names(rEG@edgeData); x})

## layout the graph
tmp = layoutGraph(rEG, nodeAttrs=nA, edgeAttrs=eA)

## hack to make sure the node labels are visible!
sizes = rescale(tmp@renderInfo@nodes$rWidth, to=c(0.2,1.5))
names(sizes) = nodes(rEG)
nodeRenderInfo(tmp) <- list(cex=sizes)

graphRenderInfo(tmp) <- list(main=title)

## plot the graph
renderGraph(tmp)
}


##
## Plot exogenous genomes
## (Rob Kitchen)
plotExogenousTaxonomyTrees = function(counts, cumcounts, what, output.dir, taxonomyInfo, fontScale=2, sampleGroups=NA, minPercent=0.5){
  
  
  # counts = exprs.exogenousGenomes_specific
  # cumcounts = exprs.exogenousGenomes_cumulative
  # taxonomyInfo = taxonomyInfo.exogenous_genomes
  # 
  # counts = exprs.exogenousRibosomal_specific
  # cumcounts = exprs.exogenousRibosomal_cumulative
  # taxonomyInfo = taxonomyInfo.exogenous_rRNA
  
  ## add direct count to the cumulative counts matrix
  cumcounts = cumcounts+counts
  
  #counts.norm = t(t(counts*100)/colSums(counts))
  counts.norm = apply(counts, 2, function(col){ col*100/sum(col) })
  cumcounts.norm = apply(cumcounts, 2, function(col){ col*100/col[1] })
  dim(counts)
  
  ## remove nodes with < 0.1% of all reads
  #minPercent = 1
  keepRows = which(apply(counts.norm, 1, max) >= minPercent)
  keepRows = sort(unique(c(keepRows, which(apply(cumcounts.norm, 1, max) >= minPercent))))
  
  # use only paths through the tree that capture above a certain fraction of reads
  counts = counts[keepRows, , drop=F]
  cumcounts = cumcounts[keepRows, , drop=F]
  nrow(counts)
  #data_uniq = counts.norm[keepRows, , drop=F]
  #data_cum = cumcounts.norm[keepRows, , drop=F]
  #nrow(data_cum)
  
  ## Re-scale the node percentages after trimming branches to make the numbers make more sense - shouldn't make much diff to the cumcounts
  data_uniq = apply(counts, 2, function(col){ col*100/sum(col) })
  data_cum = apply(cumcounts, 2, function(col){ col*100/col[1] })
  
  #if("significantDEX" %in% names(combinedSamples)){
  #  significant = combinedSamples$significantDEX[keepRows]
  #  foldChange = combinedSamples$foldChange[keepRows]
  #}
  
  ## remove edges with no useable counts (based on minPercent threshold)
  taxonomyInfo = taxonomyInfo[taxonomyInfo$ID %in% rownames(data_cum), ]
  
  ## Build the graph object
  rEG <<- new("graphNEL", nodes=as.character(taxonomyInfo$ID), edgemode="directed")
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  apply(taxonomyInfo[-1,], 1, function(row){ 
    from = trim(as.character(row[4]));
    if(from %in% taxonomyInfo$ID){ rEG <<- addEdge(trim(as.character(row[4])), trim(as.character(row[3])), rEG, 1) }
    NULL })
  
  
  data_uniq = data_uniq[match(taxonomyInfo$ID, rownames(data_uniq)), , drop=F]
  data_cum = data_cum[match(taxonomyInfo$ID, rownames(data_cum)), , drop=F]
  data_uniq[is.na(data_uniq)] = 0
  data_cum[is.na(data_cum)] = 0
  
  
  ##
  ## Write to PDF
  ##
  ## plot an average tree over all samples
  printMessage(c("Plotting a taxonomy tree based on the average of all samples "))
  pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_aggregateSamples.pdf",sep=""),height=7,width=15)
  plotTree(rEG, taxonomyInfo, apply(data_uniq, 1, max), rowMeans(data_cum), what=what)
  dev.off()
  
  ## plot samples individually
  printMessage(c("Plotting a separate taxonomy tree for each sample"))
  pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_perSample.pdf",sep=""), height=7, width=15)
  for(i in 1:ncol(data_uniq))
    plotTree(rEG, taxonomyInfo, data_uniq[,i], data_cum[,i], title=paste(colnames(data_uniq)[i]," (total reads: ",cumcounts[1,i],")", sep=""), what=what)
  dev.off()
  
  ## if there are groups of samples
  if(is.data.frame(sampleGroups)){
    printMessage(c("Plotting a separate taxonomy tree for each sample-group"))
    pdf(file=paste(output.dir,"/exceRpt_",what,"_TaxonomyTrees_perGroup.pdf",sep=""), height=7, width=15)
    for(thisgroup in levels(as.factor(sampleGroups$sampleGroup))){
      tmpDat_uniq = rowMeans(data_uniq[, match(sampleGroups[sampleGroups$sampleGroup %in% thisgroup, ]$sampleID, colnames(data_uniq)), drop=F])
      tmpDat_cum = rowMeans(data_cum[, match(sampleGroups[sampleGroups$sampleGroup %in% thisgroup, ]$sampleID, colnames(data_cum)), drop=F])
      plotTree(rEG, taxonomyInfo, tmpDat_uniq, tmpDat_cum, title=paste(thisgroup,sep=""), what=what)
    }
    dev.off()
  }
  
}

