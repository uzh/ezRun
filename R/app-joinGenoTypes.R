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
    param[['javaCall']] = paste("java", "-Djava.io.tmpdir=.")
    param[['gatk']] = file.path(Sys.getenv("GATK"),'gatk')
    
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
    gatk = param[['gatk']]
    caseName = unique(datasetCase[[paste(param$grouping,'[Factor]')]])
    myLog = paste0('log_',caseName,'.txt')
    if(param$targetFile != ''){
        param$targetFile = file.path(TARGET_ENRICHMENT_DESIGN_DIR, param$targetFile)
    }
    
    
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
        fileCmd = paste("--variant", datasetCase[['GVCF [File]']][j])
    }
    

    GenotypeGVCF = paste(gatk,'GenotypeGVCFs')
    outputFile = paste0(caseName,'.vcf')
    cmd = paste(GenotypeGVCF, "-R", param$genomeSeq,
                fileCmd,
                "--dbsnp", param$dbsnpFile,
                "--output", outputFile)
    if(param$targetFile != ''){
        cmd = paste(cmd, "-L", param$targetFile) }
    ezSystem(paste(cmd,'2>',myLog))
    
    if(param$species == 'Homo_sapiens'){
        VariantRecalibrator1 = paste(gatk,'VariantRecalibrator')
        hapmapFile = param$knownSites[grep('hapmap_.*vcf.gz$', param$knownSites)]
        h1000G_omniFile = param$knownSites[grep('1000G_omni.*vcf.gz$', param$knownSites)]
        h1000G_phase1File = param$knownSites[grep('1000G_phase1.snps.*vcf.gz$', param$knownSites)]
        
        if(param$recalibrateVariants){
            cmd = paste(VariantRecalibrator1, "-R", param$genomeSeq,
                        "--variant", outputFile,
                        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0", hapmapFile,
                        "--resource:omni,known=false,training=true,truth=false,prior=12.0",  h1000G_omniFile,
                        "--resource:1000G,known=false,training=true,truth=false,prior=10.0", h1000G_phase1File,
                        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0",  param$dbsnpFile,
                        "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ -an DP",
                        "-mode SNP",
                        "--output", paste0(caseName,'_raw.SNPs.recal'),
                        "--tranches-file", paste0(caseName,'_raw.SNPs.tranches'),
                        "--rscript-file", paste0(caseName,'_recal.SNPs.plots.R'),
                        "--max-attempts 3")
            if(param$targetFile != ''){
                cmd = paste(cmd, "-L", param$targetFile) }
            ezSystem(paste(cmd,'2>>',myLog))
            
            #Apply RecalibrationOutput for SNPs:
            VariantRecalibrator2 = paste(gatk,'ApplyVQSR')
            cmd = paste(VariantRecalibrator2, "-R", param$genomeSeq,
                        "--variant", outputFile,
                        "-mode SNP",
                        "--recal-file", paste0(caseName,'_raw.SNPs.recal'),
                        "--tranches-file",paste0(caseName,'_raw.SNPs.tranches'),
                        "--output", paste0(caseName,'_recal.SNPs.vcf'),
                        "--truth-sensitivity-filter-level 99.0")
            if(param$targetFile != ''){
                cmd = paste(cmd, "-L", param$targetFile) }
            }
            ezSystem(paste(cmd,'2>>',myLog))
            
            #2.Run VariantRecalibration for InDels:
            if(param$recalibrateInDels){
                millsFile = param$knownSites[grep('Mills.*vcf.gz$', param$knownSites)]
                cmd = paste(VariantRecalibrator1, "-R", param$genomeSeq,
                            "--variant", paste0(caseName,"_recal.SNPs.vcf"),
                            "-resource:mills,known=false,training=true,truth=true,prior=12.0", millsFile,
                            "--max-gaussians 4 --minimum-bad-variants 500 --max-attempts 3", #might be suboptimal
                            "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ -an DP",
                            "-mode INDEL",
                            "--output", paste0(caseName,"_raw.InDels.recal"),
                            "--tranches-file", paste0(caseName,"_raw.InDels.tranches"),
                            "--rscript-file", paste0(caseName,"_recal.InDels.plots.R"))
                if(param$targetFile != ''){
                    cmd = paste(cmd, "-L", param$targetFile) }
                ezSystem(paste(cmd,'2>>',myLog))
                
                #Apply RecalibrationOutput for InDels:
                cmd = paste(VariantRecalibrator2, "-R", param$genomeSeq,
                            "--variant", paste0(caseName,"_recal.SNPs.vcf"),
                            "-mode INDEL",
                            "--recal-file", paste0(caseName,'_raw.InDels.recal'),
                            "--tranches-file",paste0(caseName,'_raw.InDels.tranches'),
                            "--output", paste0(caseName,'_recal.InDels.vcf'),
                            "--truth-sensitivity-filter-level 99.0")
                if(param$targetFile != ''){
                    cmd = paste(cmd, "-L", param$targetFile) }
                ezSystem(paste(cmd,'2>>',myLog))
                ezSystem(paste('mv', paste0(caseName,"_recal.InDels.vcf"), outputFile))
            } else {
                ezSystem(paste('mv', paste0(caseName,"_recal.SNPs.vcf"), outputFile))
            }
        #Add ExAc-Annotation:
            runExAc = FALSE
            if(runExAC){
        ExAcFile = param$knownSites[grep('ExAC.*vcf.gz$', param$knownSites)]
        VariantAnnotation = paste(gatk,'VariantAnnotator')
        cmd = paste(VariantAnnotation, "-R", param$genomeSeq,
                    "--variant", outputFile,
                    "--output ",paste0(outputFile, "_annotated.vcf"),
                    "--resource:ExAC",  ExAcFile,
                    "--expression ExAC.AF",
                    "--expression ExAC.AC",
                    "--expression ExAC.AN",
                    "--expression ExAC.AC_Het",
                    "--expression ExAC.AC_Hom")
        if(param$targetFile != ''){
            cmd = paste(cmd, "-L", param$targetFile) }
        
        #cmd = '/usr/local/ngseq/packages/Variants/GATK/4.1.8.0/gatk SelectVariants --reference /srv/GT/reference/Homo_sapiens/Ensembl/GRCh37.p13/Sequence/WholeGenomeFasta/genome.fa  -V S5008.vcf --keep-original-ac true --keep-original-dp true --set-filtered-gt-to-nocall true --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --output test.vcf'
        ezSystem(paste(cmd,'2>>',myLog))
            } else {
                ezSystem(paste('mv' outputFile, paste0(outputFile, "_annotated.vcf")))
            }
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
