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
  samples = rownames(dataset)
  files = input$getFullPaths(param,"Read1")
  names(files) = samples
  executeFastqscreenCMD(param,files)
  executeBowtie2CMD(param,files)
  data = collectFastqscreenOutput(dataset,files)
  collectBowtie2Output(param,dataset)
  generatePlots(dataset,data)
  generateFastqscreenHtmlReport(param=param, dataset=dataset, htmlFile=htmlFile)
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


## NOTEP: all 6 functions below get only called once each in ezMethodFastqScreen()
executeFastqscreenCMD = function(param,files){
  confFile = paste(FASTQSCREEN_CONF_DIR,param$confFile,sep="")
  opt = ""
  if (param$nReads > 0){
    opt = paste(opt, "--subset", ezIntString(param$nReads))
  }
  cmd = paste(FASTQSCREEN, opt, " --threads", param$cores," --conf ",confFile,
                paste(files, collapse=" "),"--outdir . --aligner bowtie2",
                "> fastqscreen.out", "2> fastqscreen.err")
  ezSystem(cmd)
}

# TODO merge with 
executeBowtie2CMD = function(param,files){
  # TODO build filenames from samplenames; return result files; use ezSortIndexBam
  for (i in 1:length(files)){
    countFile = basename(files[i])
    countFile = sub('f.*q.gz$','counts.txt',countFile)
    samFile = sub('counts.txt','sam',countFile)
    bamFile = sub('counts.txt','bam',countFile)
    sbamFile = sub('bam','sorted.bam',bamFile)
    bestScoreFile = sub('counts.txt','bestScore.txt',countFile)
    bowtie2options = param$cmdOptions
    if(!param$paired){
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                " -U ",files[i], bowtie2options ,"-p",param$cores,"-u",param$nReads,
                "--no-unal -S",samFile,"2> ",paste(sub('.sam','',samFile),"_bowtie2.err",sep=''))
    } else {
      R2_file = sub('R1','R2',files[i])
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                    " -1 ",files[i]," -2 ", R2_file, bowtie2options, "-p",param$cores,"-u",param$nReads,
                    "--no-discordant --no-mixed --no-unal -S",samFile,"2> ",paste(sub('.sam','',samFile),"_bowtie2.err",sep=''))
    }
    ezSystem(cmd)
    system(paste(SAMTOOLS,"view -bS ",samFile,"-o", bamFile))
    system(paste(SAMTOOLS,"sort",bamFile, sub('.bam','',sbamFile)))       
    system(paste(SAMTOOLS,"view",sbamFile,"|cut -f1,3,12 |sort|sed 's/AS:i://g' >",countFile))
    system(paste(SAMTOOLS,"view -F 256",sbamFile,"|cut -f1,12 |sort|sed 's/AS:i://g' >",bestScoreFile))
     }
}

collectFastqscreenOutput = function(dataset,files){
  resultFiles = paste(sub(".fastq.gz", "", basename(dataset$"Read1 [File]")),'_screen.txt',sep='') 
  data = list()
  data$MappingRate = vector(length=length(resultFiles),mode='double')
  names(data$MappingRate) = rownames(dataset)
  data$Reads = vector(length=length(resultFiles),mode='integer')
  names(data$Reads) = rownames(dataset)
  data$CommonResults = list()
  for(i in 1:length(resultFiles)){
    cat('Process ',files[i],':')
    data$CommonResults[[i]] = ezRead.table(resultFiles[i], skip=1, stringsAsFactors=F, blank.lines.skip=T, fill=T, row.names=NULL)
    data$Reads[i] = data$CommonResults[[i]]$"#Reads_processed"[1]
    UnmappedReads = as.numeric(unlist(strsplit(data$CommonResults[[i]]$Genome[nrow(data$CommonResults[[i]])], split = " "))[2])
    data$MappingRate[i] = round((100 - UnmappedReads), digits = 2)
    lastLine = nrow(data$CommonResults[[i]])
    LibNames = data$CommonResults[[i]]$Genome[-c(lastLine)]
    data$CommonResults[[i]] = data$CommonResults[[i]][c(1:(lastLine-1)),grep('%.*hit',colnames(data$CommonResults[[i]]))]
    rownames(data$CommonResults[[i]]) = LibNames 
  }
  return(data)
}

# TODO: merge with executBowtie2Cmd; return data; don't make plots
collectBowtie2Output = function(param,dataset){
  require(limma)
  tax2name = read.table('/srv/GT/reference/RefSeq/mRNA/20150301/Annotation/tax2name.txt',header=F,stringsAsFactors=F,sep='|', 
                        colClasses="character",quote='', comment.char="")
  colnames(tax2name) = c('TAX_ID','Name')
  tax2name$TAX_ID = sub("\t", "", tax2name$TAX_ID)
  tax2name$Name = gsub('\t','',tax2name$Name)
  rownames(tax2name) = tax2name$TAX_ID
  files = dataset[['Read1 [File]']]
  for (i in 1:length(files)){
    countFile = basename(files[i])
    countFile = sub('f.*q.gz$','counts.txt',countFile)
    bestScoreFile = sub('counts.txt','bestScore.txt',countFile)
    SampleName = rownames(dataset)[i]
    countData = read.table(countFile,header=F,sep='',stringsAsFactors=F,comment.char = '')
    colnames(countData)=c('readId','hit','aScore')
    bestScoreData = read.table(bestScoreFile,header=F,sep='',stringsAsFactors=F,comment.char = '')
    colnames(bestScoreData) = c('readId','bestAScore')
    combinedCountData = merge(countData,bestScoreData,by.x='readId',by.y='readId')
    combinedCountData = combinedCountData[which(combinedCountData$aScore == combinedCountData$bestAScore),]
    combinedCountData = combinedCountData[which(combinedCountData$aScore >= param$minAlignmentScore),]
    ## TODO: strsplit2 fails if no hits satsify the criteria and combinedCountData has zero rows
    combinedCountData$hit = strsplit2(combinedCountData$hit,'_')[,1]
    by = combinedCountData$readId
    speciesHitsPerRead = tapply(combinedCountData$hit,by,unique)
    uniqSpeciesHitsPerRead = names(speciesHitsPerRead)[which(sapply(speciesHitsPerRead,length)==1)]
    ###Result UniqHits:
    uniqSpeciesHits = sort(table(unlist(speciesHitsPerRead[uniqSpeciesHitsPerRead])),decreasing = T)
    ###Results MultipleHits:
    multipleSpeciesHitsPerRead = combinedCountData[-which(combinedCountData$readId %in% uniqSpeciesHitsPerRead),]
    by = paste(multipleSpeciesHitsPerRead$readId,multipleSpeciesHitsPerRead$hit,sep='_')
    multipleSpeciesHits = sort(table(tapply(multipleSpeciesHitsPerRead$hit,by,unique)),decreasing = T)
    
    topSpeciesUniq = uniqSpeciesHits[1:param$nTopSpecies]
    topSpeciesMultiple = multipleSpeciesHits[names(topSpeciesUniq)]
    
    if(param$nReads > 0){
      totalCount = param$nReads
    } else {
      totalCount = dataset[['Read Count']][i]
    }
    speciesPercentageTop = signif(100 * cbind(topSpeciesUniq,topSpeciesMultiple)/totalCount, digits=4)
    taxIdsTopSpecies = rownames(speciesPercentageTop)
    namesTopSpecies = tax2name[taxIdsTopSpecies,c('Name')]
    rownames(speciesPercentageTop)[which(!is.na(namesTopSpecies))] = namesTopSpecies[which(!is.na(namesTopSpecies))]
    colnames(speciesPercentageTop) = c('UniqueSpeciesHits','MultipleSpeciesHits')
    png(paste('DetectedSpecies_',SampleName,'.png',sep=''))
    par(mar=c(18.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    bp = barplot(t(speciesPercentageTop), las=2, ylab='MappedReads in %',main=SampleName, ylim=c(0, 100),
                 col=c('blue','lightblue'),legend.text =T)
    text(y= 3,x= bp, labels=paste(speciesPercentageTop[,1],'%',sep=''), xpd=TRUE)
    dev.off()
    }
  system('rm *.sam')
  system('rm *.bam')
  system('rm *.bestScore.txt')
}
 

generatePlots = function(dataset, data){
  resultFiles = paste(basename(dataset$"Read1 [File]"),'_screen.txt',sep='')
  resultFiles = sub('\\.fastq.gz','',resultFiles)
  for(i in 1:length(data$CommonResults)){
    png(gsub('.txt','.png',resultFiles[i],'.png'))
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    barplot(t(data$CommonResults[[i]]), las=2, ylim=c(0,100),legend.text=T, ylab='MappedReads in %', main=rownames(dataset)[i])
    dev.off()
  }
  png('MappingRate.png',width=800,height=600)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  bplt = barplot(data$MappingRate,las=2,ylim=c(0,100),ylab='MappedReads in %',main="MappingRate",col="blue")
  text(y= data$MappingRate+5, x= bplt, labels=paste(as.character(data$MappingRate),'%',sep=''),cex=0.7, xpd=TRUE)
  dev.off()
  
  png('Reads.png',width=800,height=600)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  bplt = barplot(data$Reads,las=2,ylab="#Reads",main="ProcessedReads",col="lightblue")
  dev.off()
}

generateFastqscreenHtmlReport = function(param,dataset, htmlFile="00index.html"){
  html = openHtmlReport(htmlFile, param=param, title=paste("FastQ Screen:", param$name),
                        dataset=dataset)
  ezWrite("<h2>Settings</h2>", con=html)
  ezWrite("<table border='0'>", con=html)
  ezWrite("<tr><td>Configuration File:</td><td>",param$confFile, "</td></tr>", con=html)
  ezWrite("<tr><td>RefSeq mRNA Reference:</td><td>",REFSEQ_mRNA_REF,"</td></tr>", con=html)
  ezWrite("<tr><td>FastqScreen Version:</td><td>",basename(dirname(FASTQSCREEN)), "</td></tr>", con=html)
  ezWrite("<tr><td>Bowtie2 Version:</td><td>",basename(BOWTIE2_DIR), "</td></tr>", con=html)
  ezWrite("<tr><td>Bowtie2 Parameters:</td><td>", param$cmdOptions, "</td></tr>", con=html)
  ezWrite("<tr><td>Minimum AlignmentScore:</td><td>",param$minAlignmentScore, "</td></tr>", con=html)
  ezWrite("<tr><td>TopSpecies:</td><td>",param$nTopSpecies, "</td></tr>", con=html)
  ezWrite("<tr><td>Subset:</td><td>", param$subset, "</td></tr>", con=html)
  ezWrite("</table>", con=html)
  ezWrite("<h2>rRNA-Check</h2>", con=html)
  ezWrite("<h3>Per Dataset</h3>", con=html)
  pngs = list.files(".",pattern="^Reads.png|^MappingRate.png")
  ezWrite("<br>", con=html)
  writeImageRowToHtml(pngs, con=html)
  ezWrite("</br>", con=html)
  ezWrite("<h3>Per Sample</h3>", con=html)
  IMAGESperROW = 4
  pngs = list.files(".",pattern="screen\\.png")
  ezWrite("<br>", con=html)
  if(length(pngs) <= IMAGESperROW){
    writeImageRowToHtml(pngs, con=html)
  } else {
    writeImageRowToHtml(pngs[1:IMAGESperROW], con=html)
    for(i in 1:(ceiling(length(pngs)/IMAGESperROW)-1)){
      writeImageRowToHtml(pngs[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(pngs))], con=html)  
    } 
  }
  ezWrite("<h2>Mapping to RefSeq mRNA</h2>", con=html)
  pngs = list.files(".",pattern="DetectedSpecies.*.png")
  ezWrite("<br>", con=html)
  if(length(pngs) <= IMAGESperROW){
    writeImageRowToHtml(pngs, con=html)
  } else {
    writeImageRowToHtml(pngs[1:IMAGESperROW], con=html)
    for(i in 1:(ceiling(length(pngs)/IMAGESperROW)-1)){
      writeImageRowToHtml(pngs[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(pngs))], con=html)  
    } 
  }
  
  ezWrite("<h2>Misc</h2>", con=html)
  txts = list.files(".",pattern="screen\\.txt")
  writeTxtLinksToHtml(txts, con=html)
  ezSessionInfo()
  writeTxtLinksToHtml('sessionInfo.txt',con=html)
  flush(html)
  closeHTML(html)
}
