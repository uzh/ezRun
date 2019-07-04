###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppGatkRnaHaplotyper <-
  setRefClass("EzAppGatkRnaHaplotyper",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGatkRnaHaplotyper
                  name <<- "EzAppGatkRnaHaplotyper"
                  appDefaults <<- rbind(vcfFilt.minAltCount = ezFrame(Type="integer",
                                                                      DefaultValue=10,
                                                                      Description="minimum coverage for the alternative variant"),
                                        vcfCall.minReadDepth = ezFrame(Type="integer",
                                                                       DefaultValue=10,
                                                                       Description="minimum read deapth"))
                }
              )
  )


ezMethodGatkRnaHaplotyper = function(input=NA, output=NA, param=NA,
                                     htmlFile="00index.html"){
  require(Rsamtools)
  require(VariantAnnotation)
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  bamDataset = input$meta
  
  reportDir = basename(output$getColumn("Report"))
  htmlFile = basename(output$getColumn("Html"))
  vcfOutputFile = output$getColumn("VCF")
  
  bamFiles = input$getFullPaths("BAM")
  genomeSeq = param$ezRef["refFastaFile"]
  nBamsInParallel = min(8, param$cores)
  bamFilesClean = ezMclapply(names(bamFiles), function(sampleName){
    setwdNew(paste(sampleName, "proc", sep="-"))
    bf = bamFiles[sampleName]
    obf = file.path(getwd(), basename(bf))
    if(!file.exists(obf)){
      ## Easier for debug when having obf already
      if (ezIsSpecified(param$seqNames)){
        bamParam = ScanBamParam(what=scanBamWhat())
        seqLengths = ezBamSeqLengths(bf)
        bamWhich(bamParam) = GRanges(seqnames=param$seqNames,
                                     ranges=IRanges(start=1, 
                                                    end=seqLengths[param$seqNames]))
        filterBam(bf, "local.bam", param=bamParam)
        indexBam("local.bam")
      } else {
        file.copy(bf, "local.bam")
        file.copy(paste0(bf, ".bai"), "local.bam.bai")
      }
      
      cmd = paste(prepareGATK(), "SplitNCigarReads", "-R", genomeSeq,
                  "-I local.bam",
                  "-O", obf,
                  "> splitncigars.stdout 2> splitncigars.stderr")
      ezSystem(cmd)
      file.remove(c("local.bam", "local.bam.bai"))
      indexBam(obf)
    }
    setwd("..")
    return(obf)
  }, mc.cores=nBamsInParallel, mc.preschedule=FALSE)
  
  ########### haplotyping
  vcfFn <- paste0(param$name, "-all-haplo.vcf")
  if(!file.exists(vcfFn)){
    ## Easier for debug when having local.bam already
    cmd = paste(prepareGATK(), "HaplotypeCaller", "-R", genomeSeq,
                paste("-I", bamFilesClean, collapse=" "),
                "-O", vcfFn,
                "--dont-use-soft-clipped-bases true",
                "--standard-min-confidence-threshold-for-calling 20",
                "--sample-ploidy", 2,
                "--minimum-mapping-quality 20",
                "--output-mode", "EMIT_VARIANTS_ONLY",
                ">", paste0(param$name, "-haplo.stdout"),
                "2>", paste0(param$name, "-haplo.stderr"))
    ezSystem(cmd)
  }
  
  #### Variant filtering
  cmd <- paste(prepareGATK(), "VariantFiltration", "-R", genomeSeq,
               "-V", vcfFn, "-O", basename(vcfOutputFile),
               "--cluster-window-size 35 --cluster-size 3",
               "--filter-name \"my_filter\" --filter-expression \"FS > 30.0 && QD < 2.0\""
               )
  ezSystem(cmd)
  
  ## filter the vcf file
  # ezFilterVcf(vcfFile=vcfFn, 
  #             basename(vcfOutputFile), discardMultiAllelic=FALSE,
  #             bamDataset=bamDataset, param=param)
  gc()
  
  ## create an html report
  setwdNew(reportDir)
  
  chromSizes = ezChromSizesFromVcf(file.path("..", basename(vcfOutputFile)))
  genotype = geno(readVcf(file.path("..", basename(vcfOutputFile)),
                          genome="genomeDummy"))
  gt = genotype$GT
  gt[genotype$DP < param$vcfCall.minReadDepth] = "lowCov"
  nSamples = nrow(bamDataset)
  
  writeIgvSession(genome = getIgvGenome(param), 
                  refBuild=param$ezRef["refBuild"], file="igvSession.xml", 
                  vcfUrls = paste(PROJECT_BASE_URL, vcfOutputFile, sep="/") )
  writeIgvJnlp(jnlpFile="igv.jnlp",
               projectId = sub("\\/.*", "", output$getColumn("Report")),
               sessionUrl = paste(PROJECT_BASE_URL,
                                  output$getColumn("Report"),
                                  "igvSession.xml", sep="/"))
  conds = ezConditionsFromDataset(bamDataset, param=param)
  sampleColors = getSampleColors(conds, colorNames = names(conds))
  idxMat = ezMatrix(match(gt, c("0/0", "0/1", "1/1")) -2,
                    rows=rownames(gt), cols=colnames(gt))
  d = dist(t(idxMat))
  hc = hclust(d, method="ward.D2")
  hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
  hcd = colorClusterLabels(hcd, sampleColors)
  chrom = sub(":.*", "", rownames(gt))
  pos = as.integer(sub("_.*", "", sub(".*:", "", rownames(gt))))
  isRealChrom = !grepl("[\\._]", names(chromSizes)) ## TODO select chromosomes by name
  idxList = split(1:nrow(gt), chrom)
  snpColors = c("0/0"="blue", "0/1"="darkgrey", "1/1"="red")

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "gatkRnaHaplotyper.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="gatkRnaHaplotyper.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  setwd("..")
  
  return("Success")
}
