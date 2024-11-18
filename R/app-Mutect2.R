###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMutect2 = function(input=NA, output=NA, param=NA){
    sampleBamFile <- input$getFullPaths("BAM")
    genomeSeq = param$ezRef["refFastaFile"]
    ctrlBamFile <- input$getFullPaths("CtrlBam")
    sampleName = input$getNames()
    sampleNameNormal <- sub('.bam', '', basename(ctrlBamFile))
    if(param$snpEffDB == 'mm39'){
        param$snpEffConfig <- '/srv/GT/reference/Mus_musculus/UCSC/mm39/Annotation/Genes/snpEff/snpEff.config'
        param$dataDir <- paste0(dirname(param$snpEffConfig),'/data')
    }
    param$javaCall <- paste("java", "-Djava.io.tmpdir=.")
    outFile <- paste0(sampleName,".somatic.vcf.gz")
    cmd <- paste("gatk Mutect2 -R" ,genomeSeq, "-I", sampleBamFile, "-I", ctrlBamFile, "-normal", sampleNameNormal, "-O", outFile, 
                 "--native-pair-hmm-threads", param$cores, param$cmdOptions)
    ezSystem(cmd)
    filteredOutFile <- paste0(sampleName, '.filtered.somatic.vcf.gz')
    cmd <- paste("gatk FilterMutectCalls -V",outFile, "-O", filteredOutFile, "-R", genomeSeq)
    ezSystem(cmd)
    annotatedOutFile <- paste0(sampleName,'.somatic.ann.vcf')
    if(param$snpEffDB == 'mm39'){
    cmd <- paste(param$javaCall, "-jar", "$SnpEff/snpEff.jar", "ann -c",param$snpEffConfig, "-dataDir", param$dataDir,
                 "-csvStats", paste0(sampleName,".snpeff.csv"), "-s", paste0(sampleName,".snpeff.html"), param$snpEffDB, filteredOutFile, ">", annotatedOutFile)
    ezSystem(cmd)
    cmd <- paste("bgzip", annotatedOutFile)
    ezSystem(cmd)
    cmd <- paste("tabix -p vcf", paste0(annotatedOutFile,".gz"))
    ezSystem(cmd)
    ezSystem(paste('zip', paste0(sampleName, '_misc.zip'), '*snpeff*'))
    ezSystem(paste('rm *snpeff*'))
    } else {
       system(paste('mv', filteredOutFile, paste0(annotatedOutFile,".gz")))
        ezSystem(paste('touch', paste0(sampleName, '_misc.zip')))
    }
    ezSystem(paste('rm', paste0(outFile,'*'), paste0(filteredOutFile,'*')))
}

##' @template app-template
##' @templateVar method ezMethodMutect2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMutect2 <-
    setRefClass("EzAppMutect2",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodMutect2
                        name <<- "ezMethodMutect2"
                        appDefaults <<- rbind(snpEffDB=ezFrame(Type="character", DefaultValue="",	Description="snpEffDB Name")
                            )
                    }
                )
    )
