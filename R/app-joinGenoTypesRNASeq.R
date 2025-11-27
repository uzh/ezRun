###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodJoinGenoTypesRNASeq <- function(input=NA, output=NA, param=NA){
    setwdNew(param[['name']])
    param[['genomeSeq']] = param$ezRef["refFastaFile"]
    param[['species']] = limma::strsplit2(param$ezRef['refBuild'],'/')[1]
    param[['javaCall']] = paste("java", "-Djava.io.tmpdir=.")
    param[['gatk']] = file.path(Sys.getenv("GATK"),'gatk')
    
    dataset = input$meta
    dataset[['GVCF [File]']] = input$getFullPaths("GVCF")
    datasetCaseList = split(dataset,input$getColumn(param$grouping))
    results <- ezMclapply(names(datasetCaseList),runGatkPipelineRNASeq, param=param, datasetCaseList=datasetCaseList, mc.cores = param$cores)
    if(length(results) > 1){
        results <- runGatkPipelineRNASeq('allSamples', param=param, datasetCaseList = list(allSamples = dataset))
        
    }
    vcfOutputFile = results[[1]]
    chromSizes = ezChromSizesFromVcf(vcfOutputFile)
    genotype = geno(readVcf(vcfOutputFile, genome="genomeDummy"))
    gt = genotype$GT
    gt[genotype$DP < param$minReadDepth] = "lowCov" ## those calls will become NA in subsequent analyses
    
    makeRmdReport(
        input=input,
        output = output, param = param, chromSizes=chromSizes, gt=gt,
        rmdFile = "Mpileup.Rmd", reportTitle = 'GATK Joint Genotyping RNA-Seq Report'
    )
    return("Success")
}

runGatkPipelineRNASeq <- function(caseName, param=NA, datasetCaseList=NULL){
    gatk = param[['gatk']]
    datasetCase <- datasetCaseList[[caseName]]
    myLog = paste0('log_',caseName,'.txt')

    ##Create CombinedGVCF:
    if(nrow(datasetCase) >1){
        GenotypeGVCF = paste(gatk,'CombineGVCFs')
        myGVCF = paste0(caseName, ".g.vcf")
        cmd = paste(GenotypeGVCF,
                    "-R", param$genomeSeq,
                    paste('--variant', datasetCase[['GVCF [File]']], collapse = ' '),
                    "--output", myGVCF)
        ezSystem(paste(cmd,'2>',myLog))
        fileCmd = paste("--variant", myGVCF)
    } else {
        fileCmd = paste("--variant", datasetCase[['GVCF [File]']])
    }
    
    
    GenotypeGVCF = paste(gatk,'GenotypeGVCFs')
    gvcfFile = paste0(caseName,'.g.vcf')
    tmpGvcf = paste0(caseName,'_temp.vcf')
    cmd = paste(GenotypeGVCF, "-R", param$genomeSeq, fileCmd, "--output", tmpGvcf)
    ezSystem(paste(cmd,'2>',myLog))
    ezSystem(paste('mv', tmpGvcf, gvcfFile))
    ezSystem(paste('mv', paste0(tmpGvcf, ".idx"), paste0(gvcfFile, ".idx")))
    ezSystem(paste("bgzip", gvcfFile))
    ezSystem(paste0("tabix -p vcf ", gvcfFile, ".gz"))
    return(paste0(gvcfFile,'.gz'))
}


##' @template app-template
##' @templateVar method ezMethodJoinGenoTypesRNASeq(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppJoinGenoTypesRNASeq <-
    setRefClass("EzAppJoinGenoTypesRNASeq",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodJoinGenoTypesRNASeq
                        name <<- "JoinGenoTypesRNASeq"
                        appDefaults <<- rbind(
                        minReadDepth=ezFrame(Type="integer", DefaultValue="20", Description="use for clustering only SNV with coverage higher than")
                        )
                    }
                )
    )
