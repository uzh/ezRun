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
  genomeSeq = param$ezRef["refFastaFile"]
  species = limma::strsplit2(param$ezRef['refBuild'],'/')[1]
  knownSites = list.files(param$ezRef["refVariantsDir"],pattern='vcf$',full.names = T) 
  dbsnpFile = knownSites[grep('dbsnp.*vcf$', knownSites)]
  ExAcFile = knownSites[grep('ExAC.*vcf$', knownSites)]
  snpEffConfig = file.path(dirname(param$ezRef["refFeatureFile"]),"snpEff/snpEff.config")
  snpEffDB = basename(list.dirs(dirname(snpEffConfig))) ####bad Hack
  snpEffDB = snpEffDB[grep(species, snpEffDB)]
  javaCall = paste0(JAVA, " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  datasetCaseList = split(dataset,input$getColumn(param$grouping))
  for (i in 1:length(datasetCaseList)) {
    datasetCase = datasetCaseList[[i]]
    fileCmd = vector(mode = 'character',length = 1)
    for (j in 1:nrow(datasetCase)){
      fileCmd = paste(fileCmd,paste("--variant", datasetCase[['GVCF [File]']][j], collapse=','))
    }  
    
    GenotypeGVCF = paste(javaCall,"-jar", GATK_JAR, " -T GenotypeGVCFs")
    outputFile = paste0(names(datasetCaseList)[i],'.vcf')
    cmd = paste(GenotypeGVCF, "-R", genomeSeq,
                fileCmd,
                "-o", outputFile,
                "-nt", param$cores)
    if(!is.null(param$targetFile)){
      cmd = paste(cmd,
                  "-L", param$targetFile)
    }
    ezSystem(cmd)
    
    if(species == 'Homo_sapiens'){
      VariantRecalibrator1 = paste(javaCall,"-jar", GATK_JAR, " -T VariantRecalibrator")
      hapmapFile = knownSites[grep('hapmap_.*vcf$', knownSites)]
      h1000G_omniFile = knownSites[grep('1000G_omni.*vcf$', knownSites)]
      h1000G_phase1File = knownSites[grep('1000G_phase1.snps.*vcf$', knownSites)]
      
      cmd = paste(VariantRecalibrator1, "-R", genomeSeq,
                  "-input", outputFile,
                  "-resource:hapmap,known=false,training=true,truth=true,prior=15.0", hapmapFile,
                  "-resource:omni,known=false,training=true,truth=false,prior=12.0",  h1000G_omniFile,
                  "-resource:1000G,known=false,training=true,truth=false,prior=10.0", h1000G_phase1File,
                  "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0",  dbsnpFile,
                  "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ",
                  "-mode SNP",
                  "-recalFile raw.SNPs.recal",
                  "-tranchesFile raw.SNPs.tranches",
                  "-rscriptFile recal.SNPs.plots.R")
      
      if(!is.null(param$targetFile)){
        cmd = paste(cmd,
                    "-L", param$targetFile, "--max_attempts 3")
      } else {
        cmd = paste(cmd, "-an DP")
      }
      #ezSystem(cmd)
      
      #Apply RecalibrationOutput for SNPs:
      VariantRecalibrator2 = paste(javaCall,"-jar", GATK_JAR, "-T ApplyRecalibration")
      cmd = paste(VariantRecalibrator2, "-R", genomeSeq,
                  "-input", outputFile,
                  "-mode SNP",
                  "-recalFile raw.SNPs.recal",
                  "-tranchesFile raw.SNPs.tranches",
                  "-o recal.SNPs.vcf",
                  "-ts_filter_level 99.0",
                  "-nt", param$cores)
      if(!is.null(param$targetFile)){
        cmd = paste(cmd,
                    "-L", param$targetFile)
      } 
      #ezSystem(cmd)
      #ezSystem(paste("mv raw.SNPs.tranches.pdf", paste0(names(datasetCaseList)[i],"_raw.SNPs.tranches.pdf")))
      
      #2.Run VariantRecalibration for InDels:
      millsFile = knownSites[grep('Mills.*vcf$', knownSites)]
      cmd = paste(VariantRecalibrator1, "-R", genomeSeq,
                  "-input recal.SNPs.vcf",
                  "-resource:mills,known=false,training=true,truth=true,prior=12.0", millsFile,
                  "--maxGaussians 4 --minNumBadVariants 100 --max_attempts 3", #might be suboptimal
                  "-an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ",
                  "-mode INDEL",
                  "-recalFile raw.InDels.recal",
                  "-tranchesFile raw.InDels.tranches",
                  "-rscriptFile recal.InDels.plots.R")
      if(!is.null(param$targetFile)){
        cmd = paste(cmd,
                    "-L", param$targetFile)
      } else {
        cmd = paste(cmd, "-an DP")
      }
      #ezSystem(cmd)
      
      #Apply RecalibrationOutput for InDels:
      VariantRecalibrator2 = paste(javaCall,"-jar", GATK_JAR, "-T ApplyRecalibration")
      cmd = paste(VariantRecalibrator2, "-R", genomeSeq,
                  "-input recal.SNPs.vcf",
                  "-mode INDEL",
                  "-recalFile raw.InDels.recal",
                  "-tranchesFile raw.InDels.tranches",
                  "-o ",outputFile,
                  "-ts_filter_level 99.0",
                  "-nt", param$cores)
      if(!is.null(param$targetFile)){
        cmd = paste(cmd,
                    "-L", param$targetFile)
      } 
     # ezSystem(cmd)
    }
    #1.dbSnp-Annotation/ExAc-Annotation:
    cmd = paste(javaCall, "-jar", file.path(SNPEFF_DIR, "SnpSift.jar"), "annotate -id", dbsnpFile, outputFile, ">", paste0(outputFile, "_annotated.vcf"))
    ezSystem(cmd)
    if(species == 'Homo_sapiens'){
      ezSystem(paste("sed -i 's/;AF=/;Case_AF=/g' ", paste0(outputFile, "_annotated.vcf")))
      ezSystem(paste("sed -i 's/\tAC=/\tCase_AC=/g' ",paste0(outputFile, "_annotated.vcf")))
      ezSystem(paste("sed -i 's/;AN=/;Case_AN=/g' ",paste0(outputFile, "_annotated.vcf")))
      cmd = paste(javaCall, "-jar", file.path(SNPEFF_DIR, "SnpSift.jar"), "annotate -noId -info AF,AC,AN", ExAcFile, paste0(outputFile, "_annotated.vcf"), ">", outputFile)
      ezSystem(cmd)
    } else {
      ezSystem(paste("mv",paste0(outputFile, "_annotated.vcf"), outputFile))
    }
    #2.SnpEff:
    htmlOutputFile = paste0(names(datasetCaseList)[i],'.html')
    cmd = paste(javaCall, "-jar", file.path(SNPEFF_DIR, "snpEff.jar"), "-s", htmlOutputFile, "-c", snpEffConfig, snpEffDB, "-v", outputFile,">", 
                paste0(outputFile, "_annotated.vcf"))
    ezSystem(cmd)
    
    #3.Add SIFT & CO Annotation if Human (dbnsfp):
    if(species == 'Homo_sapiens'){
      cmd = paste(javaCall, "-jar", file.path(SNPEFF_DIR, "SnpSift.jar"), "dbnsfp -f SIFT_score,SIFT_pred,Polyphen2_HDIV_pred,1000Gp1_EUR_AF,CADD_phred -v -db /srv/GT/databases/dbNSFP/dbNSFP2.9.txt.gz",paste0(outputFile, "_annotated.vcf"), ">", outputFile)
      ezSystem(cmd)
    } else {
      ezSystem(paste("mv",paste0(outputFile, "_annotated.vcf"), outputFile))
    }
  }
  ezSystem('rm *.tranches *recal* *_annotated.vcf')
  
  
  
  #Further processing for Kispi-Samples ->make Report per Family (per flag)
  return("Success")
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
                                        vcfFilt.minAltCount = ezFrame(Type="integer",  DefaultValue=8,  Description="minimum coverage for the alternative variant"),
                                        vcfFilt.minReadDepth = ezFrame(Type="integer",  DefaultValue=20,  Description="minimum read depth"))
                }
              )
  )
