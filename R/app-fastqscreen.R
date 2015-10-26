###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Fastq Screen
##' @seealso \code{\link{EzAppFastqScreen}}
ezMethodFastqScreen = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  dataset = input$meta
  files = input$getFullPaths(param, "Read1")
  resultFiles = executeFastqscreenCMD(param, files)
  countFiles = executeBowtie2CMD(param, files)
  fastqData = collectFastqscreenOutput(dataset, files, resultFiles)
  speciesPercentageTop = collectBowtie2Output(param, dataset, countFiles)
  fastqscreenReport(dataset, fastqData, param, htmlFile, resultFiles, speciesPercentageTop)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastqScreen()
##' @seealso \code{\link{ezMethodFastqScreen}}
EzAppFastqScreen <-
  setRefClass("EzAppFastqScreen",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFastqScreen
                  name <<- "EzAppFastqScreen"
                  appDefaults <<- rbind(nTopSpecies=ezFrame(Type="integer",  DefaultValue=10,  Description="number of species to show in the plots"),
                                        confFile=ezFrame(Type="character",  DefaultValue="",  Description="the configuration file for fastq screen"),
                                        minAlignmentScore=ezFrame(Type="integer",  DefaultValue="-20",  Description="the min alignment score for bowtie2"))
                }
              )
  )


## NOTEP: all 5 functions below get only called once each in ezMethodFastqScreen()
executeFastqscreenCMD = function(param, files){
  confFile = paste(FASTQSCREEN_CONF_DIR, param$confFile, sep="")
  opt = ""
  if (param$nReads > 0){
    opt = paste(opt, "--subset", ezIntString(param$nReads))
  }
  cmd = paste(FASTQSCREEN, opt, " --threads", param$cores, " --conf ", confFile,
              paste(files, collapse=" "), "--outdir . --aligner bowtie2",
              "> fastqscreen.out", "2> fastqscreen.err")
  ezSystem(cmd)
  resultFiles = paste(sub(".fastq.gz", "", basename(files)), '_screen.txt',sep='')
  names(resultFiles) = names(files)
  return(resultFiles)
}

# TODO merge with 
executeBowtie2CMD = function(param, files){
  countFiles = character()
  for (nm in names(files)){
    countFiles[nm] = paste(nm, "-counts.txt", sep="")
    bowtie2options = param$cmdOptions
    if(!param$paired){
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                  " -U ",files[i], bowtie2options ,"-p",param$cores,"-u",param$nReads,
                  "--no-unal", "2> ", paste(nm, "_bowtie2.err",sep=''),
                  "|", SAMTOOLS, "view -S -b -", ">", "bowtie.bam")
    } else {
      R2_file = sub('R1','R2',files[i])  ## TODO: this is a hack, R2 files should be passed to the function
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                  " -1 ",files[i]," -2 ", R2_file, bowtie2options, "-p",param$cores,"-u",param$nReads,
                  "--no-discordant --no-mixed --no-unal",
                  "2> ",paste(nm, "_bowtie2.err",sep=''),
                  "|", SAMTOOLS, "view -S -b -", ">", "bowtie.bam")
    }
    ezSystem(cmd)
    ezSystem(paste(SAMTOOLS,"view", "bowtie.bam", "|cut -f1,3,12 |sort|sed 's/AS:i://g' >",countFiles[nm]))
    #system(paste(SAMTOOLS,"view -F 256",sbamFile,"|cut -f1,12 |sort|sed 's/AS:i://g' >",bestScoreFile))
  }
  system('rm bowtie.bam *.err')
  return(countFiles)
}

collectFastqscreenOutput = function(dataset, files, resultFiles){
  fastqData = list()
  fastqData$MappingRate = vector(length=length(resultFiles),mode='double')
  names(fastqData$MappingRate) = rownames(dataset)
  fastqData$Reads = vector(length=length(resultFiles),mode='integer')
  names(fastqData$Reads) = rownames(dataset)
  fastqData$CommonResults = list()
  for(i in 1:length(resultFiles)){
    cat('Process ',files[i],':')
    x = ezRead.table(resultFiles[i], skip=1, stringsAsFactors=F, blank.lines.skip=T, fill=T, row.names=NULL)
    fastqData$Reads[i] = x$"#Reads_processed"[1]
    UnmappedReads = as.numeric(unlist(strsplit(x$Genome[nrow(x)], split = " "))[2])
    fastqData$MappingRate[i] = round((100 - UnmappedReads), digits = 2)
    x = x[-nrow(x), grep('%.*hit',colnames(x))]
    rownames(x) = x$Genome
    fastqData$CommonResults[[i]] = x
  }
  return(fastqData)
}

# TODO: merge with executBowtie2Cmd; return data; don't make plots
collectBowtie2Output = function(param, dataset, countFiles){
  require(limma)
  tax2name = read.table('/srv/GT/reference/RefSeq/mRNA/20150301/Annotation/tax2name.txt',header=F,stringsAsFactors=F,sep='|', 
                        colClasses="character",quote='', comment.char="")
  colnames(tax2name) = c('TAX_ID','Name')
  tax2name$TAX_ID = sub("\t", "", tax2name$TAX_ID)
  tax2name$Name = gsub('\t','',tax2name$Name)
  rownames(tax2name) = tax2name$TAX_ID
  speciesPercentageTop = list()
  for (nm in names(countFiles)){
    countData = read.table(countFiles[nm], header=F, sep='', stringsAsFactors=F, comment.char='')
    colnames(countData) = c('readId', 'hit', 'aScore')
    bestScores = tapply(countData$aScore, countData$readId, max)
    countData = countData[ countData$aScore == bestScores[countData$readId], ]
    countData = countData[countData$aScore >= param$minAlignmentScore, ]
    ## TODO: strsplit2 fails if no hits satsify the criteria and countData has zero rows
    countData$species = strsplit2(countData$hit,'_')[,1]
    speciesHitsPerRead = tapply(countData$species, by=countData$readId, unique)
    uniqSpeciesHitsPerRead = names(speciesHitsPerRead)[which(sapply(speciesHitsPerRead, length)==1)]
    ###Result UniqHits:
    uniqSpeciesHits = sort(table(unlist(speciesHitsPerRead[uniqSpeciesHitsPerRead])), decreasing = T)
    ###Results MultipleHits:
    multipleSpeciesHitsPerRead = countData[-which(countData$readId %in% uniqSpeciesHitsPerRead),]
    by = paste(multipleSpeciesHitsPerRead$readId, multipleSpeciesHitsPerRead$species,sep='_')
    ##multipleSpeciesHits = sort(table(multipleSpeciesHitsPerRead$species[!duplicated(by)]), decreasing=T) # is equivalent to the row below
    multipleSpeciesHits = sort(table(tapply(multipleSpeciesHitsPerRead$species,by,unique)),decreasing = T)
    
    topSpeciesUniq = uniqSpeciesHits[1:param$nTopSpecies]
    topSpeciesMultiple = multipleSpeciesHits[names(topSpeciesUniq)]
    
    if(param$nReads > 0){
      totalCount = param$nReads
    } else {
      totalCount = dataset[nm, 'Read Count']
    }
    taxIds = names(topSpeciesMultiple)
    taxNames = tax2name[taxIds, 'Name']
    hasNoName = is.na(taxNames)
    taxNames[hasNoName] = taxIds[hasNoName]
    x = ezFrame(UniqueSpeciesHits=signif(100 * topSpeciesUniq/totalCount, digits=4),
                MultipleSpeciesHits=signif(100 * topSpeciesMultiple/totalCount, digits=4),
                row.names = taxNames)
    speciesPercentageTop[[nm]] = x
  }
  return(speciesPercentageTop)
}

## EzPlotter$new(expr="plotMappingRate(data$MappingRate)", mouseOvertext="ädfasfd", )
fastqscreenReport = function(dataset, fastqData=fastqData, param, htmlFile="00index.html", resultFiles, speciesPercentageTop){
  require(ReporteRs, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  html = openBsdocReport(title=paste("FastQ Screen:", param$name), dataset=dataset)
  html = addTitle(html, "Settings", level=2)
  settings = ezFrame(Value=character(0))
  settings["Configuration File:", "Value"] = param$confFile
  settings["RefSeq mRNA Reference:", "Value"] = REFSEQ_mRNA_REF
  settings["FastqScreen Version:", "Value"] = basename(dirname(FASTQSCREEN))
  settings["Bowtie2 Version:", "Value"] = basename(BOWTIE2_DIR)
  settings["Bowtie2 Parameters:", "Value"] = param$cmdOptions
  settings["Minimum AlignmentScore:", "Value"] = param$minAlignmentScore
  settings["TopSpecies:", "Value"] = param$nTopSpecies
  settings["Subset:", "Value"] = param$subset
  html = addFlexTable(html, ezFlexTable(settings, add.rownames=TRUE))
  html = addTitle(html, "rRNA-Check", level=2)
  html = addTitle(html, "Per Dataset", level=3)
  plotter = EzPlotterFastqScreen$new(x=fastqData$MappingRate)
  mappingRateLink = ezImageFileLink(plotter, file="MappingRate.png", width=600, height=450, las=2, ylim=c(0,100),
                                    ylab='MappedReads in %', main="MappingRate", col="blue",
                                    addText=expression(text(y=data$x+5, x=data$bplot,
                                                            labels=paste(as.character(data$x), '%', sep=''), cex=0.7, xpd=TRUE)))
  plotter = EzPlotterFastqScreen$new(x=fastqData$Reads)
  readsLink = ezImageFileLink(plotter, file="Reads.png", width=600, height=450, las=2,    #seq(0.7,0.7+1.2*(length(x)-1),by=1.2)
                              ylab="#Reads", main="ProcessedReads", col="lightblue")
  html = addFlexTable(html, ezFlexTable(cbind(mappingRateLink, readsLink)))
  html = addTitle(html, "Per Sample", level=3)
  screenLinks = list()
  for (i in 1:length(fastqData$CommonResults)){
    plotter = EzPlotterFastqScreen$new(x=t(fastqData$CommonResults[[i]]))
    link = ezImageFileLink(plotter, file=gsub('.txt', '.png', resultFiles[i], '.png'), las=2, ylim=c(0,100),
                           legend.text=T, ylab='MappedReads in %', main=rownames(dataset)[i])
    if (grepl("screen", resultFiles[i])){
      screenLinks = append(screenLinks, link)
    }
  }
  detectedSpeciesLinks = list()
  for (i in 1:length(speciesPercentageTop)){
    plotter = EzPlotterFastqScreen$new(x=t(speciesPercentageTop[[i]]))
    link = ezImageFileLink(plotter, file=paste('DetectedSpecies_', rownames(dataset)[i], '.png', sep=''),
                           col=c('blue', 'lightblue'), las=2, ylim=c(0,100), legend.text=T, ylab='MappedReads in %',
                           main=rownames(dataset)[i], addText=expression(text(y=3, x=data$bplot,
                                                                              labels=paste(data$x[1,], '%', sep=''), xpd=TRUE)))
    detectedSpeciesLinks = append(detectedSpeciesLinks, link)
  }
  IMAGESperROW = 4
  if (ezIsSpecified(screenLinks)){
    if(length(screenLinks) <= IMAGESperROW){
      html = addFlexTable(html, ezFlexTable(rbind(screenLinks)))
    } else {
      html = addFlexTable(html, ezFlexTable(rbind(screenLinks[1:IMAGESperROW])))
      for (i in 1:(ceiling(length(screenLinks)/IMAGESperROW)-1)){
        html = addFlexTable(html, ezFlexTable(rbind(screenLinks[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(screenLinks))])))
      }
    }
  }
  if (ezIsSpecified(detectedSpeciesLinks)){
    html = addTitle(html, "Mapping to RefSeq mRNA", level=2)
    if(length(detectedSpeciesLinks) <= IMAGESperROW){
      html = addFlexTable(html, ezFlexTable(rbind(detectedSpeciesLinks)))
    } else {
      html = addFlexTable(html, ezFlexTable(rbind(detectedSpeciesLinks[1:IMAGESperROW])))
      for (i in 1:(ceiling(length(detectedSpeciesLinks)/IMAGESperROW)-1)){
        html = addFlexTable(html, ezFlexTable(rbind(detectedSpeciesLinks[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(detectedSpeciesLinks))])))
      }
    }
  }
  html = addTitle(html, "Misc", level=2)
  txts = list.files(".",pattern="screen\\.txt")
  for (each in txts){
    html = addParagraph(html, pot(each, hyperlink = each))
  }
  ezSessionInfo()
  html = addParagraph(html, pot("sessionInfo.txt", hyperlink = "sessionInfo.txt"))
  closeBsdocReport(doc=html, file=htmlFile)
}
