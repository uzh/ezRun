###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodJoinGenoTypes = function(input=NA, output=NA, param=NA){
    setwdNew(param[['name']])
    dataset = input$meta
    dataset[['GVCF [File]']] = input$getFullPaths("GVCF")
    datasetCaseList = split(dataset,input$getColumn(param$grouping))
    if(is.null(param$mc.cores)){
        param[['mc.cores']] = max(1L, round(param$cores/length(datasetCaseList)))
    }
    param[['genomeSeq']] = param$ezRef["refFastaFile"]
    param[['species']] = limma::strsplit2(param$ezRef['refBuild'],'/')[1]
    param[['knownSites']] = list.files(param$ezRef["refVariantsDir"],pattern='vcf.gz$',full.names = T) 
    param[['dbsnpFile']] = param$knownSites[grep('dbsnp.*vcf.gz$', param$knownSites)]
    param[['snpEffConfig']] = file.path(dirname(param$ezRef["refFeatureFile"]), "snpEff/snpEff.config")
    snpEffDB = basename(list.dirs(dirname(param$snpEffConfig))) ####bad Hack
    param[['snpEffDB']] = snpEffDB[grep(param$species, snpEffDB)]
    param[['javaCall']] = paste0("java", " -Djava.io.tmpdir=. -Xmx", param$ram/param$mc.cores, "g")
    
    require(parallel)
    ezMclapply(datasetCaseList,runGatkPipeline, param, mc.cores = param$cores)
    if(param$recalibrateVariants){
        ezSystem('rm *.tranches *recal* *_annotated.vcf*')
    } else {
        ezSystem('rm *_annotated.vcf*')
    }
    #Further processing for Kispi-Samples ->make Report per Family (per flag)
    return("Success")
}

runGatkPipeline = function(datasetCase, param=NA){
    caseName = unique(datasetCase[[paste(param$grouping,'[Factor]')]])
    myLog = paste0('log_',caseName,'.txt')
    fileCmd = vector(mode = 'character',length = 1)
    for (j in 1:nrow(datasetCase)){
        fileCmd = paste(fileCmd,paste("--variant", datasetCase[['GVCF [File]']][j], collapse=','))
    }  
    
    GenotypeGVCF = paste(param$javaCall,"-jar", Sys.getenv("GATK_jar"), " -T GenotypeGVCFs")
    outputFile = paste0(caseName,'.vcf')
    cmd = paste(GenotypeGVCF, "-R", param$genomeSeq,
                fileCmd,
                "--dbsnp", param$dbsnpFile,
                "-o", outputFile,
                "-nt", 1) ####It's often crashing with more than one core
    if(param$targetFile != ''){
        param$targetFile = file.path(TARGET_ENRICHMENT_DESIGN_DIR, param$targetFile)
        cmd = paste(cmd,
                    "-L", param$targetFile)
    }
    ezSystem(paste(cmd,'2>',myLog))
    
    if(param$species == 'Homo_sapiens'){
        VariantRecalibrator1 = paste(param$javaCall,"-jar", Sys.getenv("GATK_jar"), " -T VariantRecalibrator")
        hapmapFile = param$knownSites[grep('hapmap_.*vcf.gz$', param$knownSites)]
        h1000G_omniFile = param$knownSites[grep('1000G_omni.*vcf.gz$', param$knownSites)]
        h1000G_phase1File = param$knownSites[grep('1000G_phase1.snps.*vcf.gz$', param$knownSites)]
        
        if(param$recalibrateVariants){
            cmd = paste(VariantRecalibrator1, "-R", param$genomeSeq,
                        "-input", outputFile,
                        "-resource:hapmap,known=false,training=true,truth=true,prior=15.0", hapmapFile,
                        "-resource:omni,known=false,training=true,truth=false,prior=12.0",  h1000G_omniFile,
                        "-resource:1000G,known=false,training=true,truth=false,prior=10.0", h1000G_phase1File,
                        "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0",  param$dbsnpFile,
                        "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ",
                        "-mode SNP",
                        "-recalFile", paste0(caseName,'_raw.SNPs.recal'),
                        "-tranchesFile",paste0(caseName,'_raw.SNPs.tranches'),
                        "-rscriptFile", paste0(caseName,'_recal.SNPs.plots.R'))
            
            if(param$targetFile != ''){
                cmd = paste(cmd,
                            "-L", param$targetFile, "--max_attempts 3")
            } else {
                cmd = paste(cmd, "-an DP")
            }
            ezSystem(paste(cmd,'2>>',myLog))
            
            #Apply RecalibrationOutput for SNPs:
            VariantRecalibrator2 = paste(param$javaCall,"-jar", Sys.getenv("GATK_jar"), "-T ApplyRecalibration")
            cmd = paste(VariantRecalibrator2, "-R", param$genomeSeq,
                        "-input", outputFile,
                        "-mode SNP",
                        "-recalFile", paste0(caseName,'_raw.SNPs.recal'),
                        "-tranchesFile",paste0(caseName,'_raw.SNPs.tranches'),
                        "-o",paste0(caseName,'_recal.SNPs.vcf'),
                        "-ts_filter_level 99.0",
                        "-nt", round(param$cores/param$mc.cores))
            if(param$targetFile != ''){
                cmd = paste(cmd,
                            "-L", param$targetFile)
            } 
            ezSystem(paste(cmd,'2>>',myLog))
            
            #2.Run VariantRecalibration for InDels:
            if(param$recalibrateInDels){
                millsFile = param$knownSites[grep('Mills.*vcf.gz$', param$knownSites)]
                cmd = paste(VariantRecalibrator1, "-R", param$genomeSeq,
                            "-input", paste0(caseName,"_recal.SNPs.vcf"),
                            "-resource:mills,known=false,training=true,truth=true,prior=12.0", millsFile,
                            "--maxGaussians 4 --minNumBadVariants 500 --max_attempts 3", #might be suboptimal
                            "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ",
                            "-mode INDEL",
                            "-recalFile", paste0(caseName,"_raw.InDels.recal"),
                            "-tranchesFile", paste0(caseName,"_raw.InDels.tranches"),
                            "-rscriptFile", paste0(caseName,"_recal.InDels.plots.R"))
                if(param$targetFile != ''){
                    cmd = paste(cmd,
                                "-L", param$targetFile)
                } else {
                    cmd = paste(cmd, "-an DP")
                }
                ezSystem(paste(cmd,'2>>',myLog))
                
                #Apply RecalibrationOutput for InDels:
                VariantRecalibrator2 = paste(param$javaCall,"-jar", Sys.getenv("GATK_jar"), "-T ApplyRecalibration")
                cmd = paste(VariantRecalibrator2, "-R", param$genomeSeq,
                            "-input", paste0(caseName,"_recal.SNPs.vcf"),
                            "-mode INDEL",
                            "-recalFile",paste0(caseName,"_raw.InDels.recal"),
                            "-tranchesFile",paste0(caseName,"_raw.InDels.tranches"),
                            "-o ",outputFile,
                            "-ts_filter_level 99.0",
                            "-nt", round(param$cores/param$mc.cores))
                if(param$targetFile != ''){
                    cmd = paste(cmd,
                                "-L", param$targetFile)
                } 
                ezSystem(paste(cmd,'2>>',myLog))
            } else {
                ezSystem(paste('mv', paste0(caseName,"_recal.SNPs.vcf"), outputFile))
            }
        }
        #Add ExAc-Annotation:
        ExAcFile = param$knownSites[grep('ExAC.*vcf.gz$', param$knownSites)]
        VariantAnnotation = paste(param$javaCall,"-jar", Sys.getenv("GATK_jar"), "-T VariantAnnotator")
        cmd = paste(VariantAnnotation, "-R", param$genomeSeq,
                    "--resource:ExAC",  ExAcFile,
                    "-o ",paste0(outputFile, "_annotated.vcf"),
                    "--expression ExAC.AF",
                    "--expression ExAC.AC",
                    "--expression ExAC.AN",
                    "--expression ExAC.AC_Het",
                    "--expression ExAC.AC_Hom",
                    "-V", outputFile,
                    "-nt", 1) 
        if(param$targetFile != ''){
            cmd = paste(cmd,
                        "-L", param$targetFile)
        } 
        ezSystem(paste(cmd,'2>>',myLog))
        cmd = paste(param$javaCall, "-jar", "$SnpEff/SnpSift.jar", "dbnsfp -f 1000Gp1_EUR_AF,Uniprot_acc,Interpro_domain,phastCons100way_vertebrate,CADD_phred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_score,SIFT_pred -v -db /srv/GT/databases/dbNSFP/dbNSFP2.9.txt.gz",paste0(outputFile, "_annotated.vcf"), ">", outputFile)
        ezSystem(paste(cmd,'2>>',myLog))
    } 
    
    ezSystem(paste("mv",outputFile,paste0(outputFile, "_annotated.vcf")))
    #SnpEff:
    if(param$annotateVariants){
        htmlOutputFile = paste0(caseName,'.html')
        if(param$proteinCodingTranscriptsOnly){
            gtfFile = param$ezRef@refFeatureFile
            gtf = rtracklayer::import(gtfFile)
            idx = gtf$type == 'transcript'
            gtf = data.frame(gtf[idx])
            protCodingTranscripts = gtf$transcript_id[gtf$source=='protein_coding']
            write.table(protCodingTranscripts, 'protCodingTranscripts.txt',col.names = F, row.names = F, quote =F)
            cmd = paste(param$javaCall, "-jar", "$SnpEff/snpEff.jar","-onlyTr protCodingTranscripts.txt", "-s", htmlOutputFile, "-c", param$snpEffConfig, param$snpEffDB, "-v", paste0(outputFile, "_annotated.vcf"),">", 
                        outputFile)
        } else {
            cmd = paste(param$javaCall, "-jar", "$SnpEff/snpEff.jar", "-s", htmlOutputFile, "-c", param$snpEffConfig, param$snpEffDB, "-v", paste0(outputFile, "_annotated.vcf"),">", 
                        outputFile)
        }
        ezSystem(paste(cmd,'2>>',myLog))
    } else {
        ezSystem(paste("mv",paste0(outputFile, "_annotated.vcf"), outputFile))
    }
}


##' @template app-template
##' @templateVar method ezMethodJoinGenoTypes(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppJoinGenoTypes <-
    setRefClass("EzAppJoinGenoTypes",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodJoinGenoTypes
                        name <<- "EzAppJoinGenoTypes"
                        appDefaults <<- rbind(targetFile = ezFrame(Type="character",  DefaultValue="", Description="restrict to targeted genomic regions"),
                                              vcfFilt.minAltCount = ezFrame(Type="integer",  DefaultValue=3,  Description="minimum coverage for the alternative variant"),
                                              vcfFilt.minReadDepth = ezFrame(Type="integer",  DefaultValue=20,  Description="minimum read depth"),
                                              recalibrateVariants = ezFrame(Type="logical",  DefaultValue=TRUE,  Description="recalibrateVariants"),
                                              recalibrateInDels = ezFrame(Type="logical",  DefaultValue=TRUE,  Description="recalibrate InDels"),
                                              annotateVariants = ezFrame(Type="logical",  DefaultValue=TRUE,  Description="annotate Variants with SnpEff"),
                                              proteinCodingTranscriptsOnly = ezFrame(Type="logical",  DefaultValue=TRUE,  Description="annotate variants only with protCod variants"))
                    }
                )
    )
