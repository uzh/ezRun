context("Counting and bam apps with example data")

## this tests do take long therefore we only run them if the environment variable RUN_LONG_TEST is set to TRUE
# Sys.setenv(RUN_LONG_TEST=TRUE)

cwd = getwd()

skipLong = function() {
  if (Sys.getenv("RUN_LONG_TEST") == "TRUE") {
    return()
  } else {
    skip("not running lengthy tests")
  }
}

yeastCommonCountParam = function() {
  param = list()
  param[['cores']] = '1'
  param[['ram']] = '10'
  param[['scratch']] = '10'
  param[['node']] = ''
  param[['process_mode']] = 'SAMPLE'
  param[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
  param[['paired']] = 'true'
  param[['strandMode']] = 'sense'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['featureLevel']] = 'gene'
  param[['specialOptions']] = ''
  param[['mail']] = ''
  param[['dataRoot']] = system.file(package = "ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Count_Result'
  return(param)
}

test_that("Count_RSEM", {
  skipLong()
  ezSystem("rm -fr /scratch/test_rsem/*")
  setwdNew("/scratch/test_rsem")
  param = yeastCommonCountParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = list()
  output[['Name']] = 'wt_1'
  output[['Count [File]']] = 'p1001/Count_RSEM/wt_1.txt'
  output[['Species']] = 'S. cerevisiae'
  output[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[['featureLevel']] = 'isoform'
  output[['refFeatureFile']] = 'genes.gtf'
  output[['strandMode']] = 'sense'
  output[['paired']] = 'true'
  output[['Read Count']] = '9794'
  output[['Genotype [Factor]']] = 'wt'
  param[['trimAdapter']] = 'false'
  param[['trimLeft']] = '0'
  param[['trimRight']] = '0'
  param[['minTailQuality']] = '0'
  param[['bowtie-e']] = '200'
  param[['cmdOptions']] = ' --calc-ci '
  param[['keepBam']] = 'false'
  param[['minAvgQuality']] = '10'
  param[['transcriptFasta']] = ''
  myApp = EzAppRSEM$new()
  myApp$run(
    input = input$copy()$subset(1),
    output = EzDataset$new(metaNew = output, dataRoot = param$dataRoot),
    param = param
  )
  setwd(cwd)
})

test_that("Count_FeatureCounts", {
  skipLong()
  ezSystem("rm -fr /scratch/test_featureCounts/*")
  setwdNew("/scratch/test_featureCounts")
  param = yeastCommonCountParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR_counts/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  param[['gtfFeatureType']] = 'exon'
  param[['allowMultiOverlap']] = 'true'
  param[['countPrimaryAlignmentsOnly']] = 'true'
  param[['minFeatureOverlap']] = '10'
  param[['minMapQuality']] = '10'
  param[['keepMultiHits']] = 'true'
  myApp = EzAppFeatureCounts$new()
  myApp$run(
    input = input$copy()$subset(1),
    output = output$copy()$subset(1),
    param = param
  )
  setwd(cwd)
})

test_that("getAnnotatedSpliceSites() derives the introns from the exons of a transcript", {
  param = yeastCommonCountParam()
  param = ezParam(param)
  gff = ezLoadFeatures(
    param,
    featureFile = system.file(
      "extdata/genes.gtf",
      package = "ezRun",
      mustWork = TRUE
    )
  )
  spliceSites = getAnnotatedSpliceSites(gff)
  expect_named(spliceSites, c("starts", "ends"))
  expect_gt(length(spliceSites$starts), 0)
  ## one intron per gap between the exons of a transcript, so fewer sites than transcripts
  expect_lt(length(spliceSites$starts), length(unique(gff$transcript_id)))
  ## the positions are the first and the last base of the intron, i.e. inside the transcript
  exons = gff[gff$type == "exon" & !is.na(gff$transcript_id), ]
  multiExon = names(which(table(exons$transcript_id) > 1))
  expect_gt(length(multiExon), 0)
  first = exons[exons$transcript_id == multiExon[1], ]
  first = first[order(first$start), ]
  expect_true(paste(first$seqid[1], first$end[1] + 1) %in% spliceSites$starts)
  expect_true(paste(first$seqid[1], first$start[2] - 1) %in% spliceSites$ends)
})

test_that("getJunctionStatsFromSJ() reproduces the RSeQC junction classes", {
  ## RSeQC judges the two splice sites independently: both known is annotated, neither known
  ## is complete_novel and exactly one known is partial_novel
  spliceSites = list(
    starts = c("chr1 100", "chr1 500"),
    ends = c("chr1 200", "chr1 600")
  )
  sjFile = tempfile(fileext = ".SJ.out.tab")
  on.exit(file.remove(sjFile))
  ## annotated, novel combination of known sites (also annotated), partial and complete novel,
  ## plus one junction below the minimum intron size that has to be dropped
  write.table(
    data.frame(
      seqid = "chr1",
      start = c(100L, 100L, 100L, 300L, 700L),
      end = c(200L, 600L, 400L, 800L, 740L),
      strand = 1L,
      motif = 1L,
      annotatedInSjdb = 0L,
      nUnique = c(50L, 10L, 5L, 1L, 99L),
      nMulti = c(0L, 0L, 5L, 0L, 0L),
      maxOverhang = 40L
    ),
    file = sjFile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  result = getJunctionStatsFromSJ(sjFile, spliceSites)
  expect_named(
    result,
    c("splice_junction", "splice_events", "junctionSaturation")
  )
  ## the 41 bp intron is below the default minIntronSize of 50 and must not be counted
  expect_equal(sum(result$splice_junction), 100)
  expect_equal(sum(result$splice_events), 100)
  expect_equal(as.vector(result$splice_junction[["annotated"]]), 50)
  expect_equal(as.vector(result$splice_junction[["partial_novel"]]), 25)
  expect_equal(as.vector(result$splice_junction[["complete_novel"]]), 25)

  ## the saturation curve counts junctions, so it ends at the number of junctions kept and
  ## rises monotonically; the known and novel curves must add up to the total
  saturation = result$junctionSaturation
  expect_named(
    saturation,
    c("all junctions", "known junctions", "novel junctions")
  )
  expect_equal(as.vector(saturation[["all junctions"]][["1"]]), 4)
  expect_equal(
    saturation[["known junctions"]] + saturation[["novel junctions"]],
    saturation[["all junctions"]]
  )
  expect_false(is.unsorted(saturation[["all junctions"]]))
  expect_identical(
    names(saturation[["all junctions"]]),
    as.character(seq(0.05, 1, by = 0.05))
  )

  ## a file holding nothing usable must not break the report
  emptyFile = tempfile(fileext = ".SJ.out.tab")
  writeLines("chr1\t700\t740\t1\t1\t0\t1\t0\t40", emptyFile)
  on.exit(file.remove(emptyFile), add = TRUE)
  expect_identical(getJunctionStatsFromSJ(emptyFile, spliceSites), list())
})

test_that("getUpstreamQcFiles() skips absent, missing and empty upstream files", {
  param = yeastCommonCountParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  ## a column the upstream app did not write at all, e.g. an older version of it
  expect_identical(getUpstreamQcFiles(input, "DupRate"), character(0))

  ## build a dataset that does declare the column, with dataRoot="" so that the
  ## paths are taken as they are
  dupFile = tempfile(fileext = "_dupRate.txt")
  meta = ezRead.table(
    system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    row.names = NULL
  )[1, ]
  meta[["DupRate [File]"]] = dupFile
  datasetFile = tempfile(fileext = ".tsv")
  ezWrite.table(meta, file = datasetFile, row.names = FALSE)
  on.exit(file.remove(datasetFile))
  input2 = EzDataset$new(file = datasetFile, dataRoot = "")

  ## the column is there but the file is not: must not throw, the caller recomputes
  expect_false(file.exists(dupFile))
  expect_length(getUpstreamQcFiles(input2, "DupRate"), 0)

  ## an existing but empty file also counts as missing
  file.create(dupFile)
  on.exit(file.remove(dupFile), add = TRUE)
  expect_length(getUpstreamQcFiles(input2, "DupRate"), 0)

  ## only a file with content is handed on, keeping the sample name because the callers
  ## look the results up by it
  writeLines("ID\tRPK\tdupRate", dupFile)
  usable = getUpstreamQcFiles(input2, "DupRate")
  expect_identical(unname(usable), dupFile)
  expect_identical(names(usable), input2$getNames())
})

test_that("RNA_Bamstats", {
  skipLong()
  ezSystem("rm -fr /scratch/test_RNA_Bamstats/*")
  setwdNew("/scratch/test_RNA_Bamstats")
  param = yeastCommonCountParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = list()
  output[['Name']] = 'RNA_BAM_Statistics'
  output[[
    'Report [File]'
  ]] = 'p1001/QC_RNABamStats_8617_2015-11-17--10-15-39/RNA_BAM_Statistics'
  output[[
    'Html [Link]'
  ]] = 'p1001/QC_RNABamStats_8617_2015-11-17--10-15-39/RNA_BAM_Statistics/00index.html'
  output[['Species']] = ''
  output[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[['refFeatureFile']] = 'genes.gtf'
  param[['process_mode']] = 'DATASET'
  param[['name']] = 'RNA_BAM_Statistics'
  myApp = EzAppRnaBamStats$new()
  myApp$run(input = input, output = output, param = param)
  setwd(cwd)
})

test_that("TEQC", {
  skipLong()
  ezSystem("rm -fr /scratch/test_TEQC/*")
  setwdNew("/scratch/test_TEQC")
  param = yeastCommonCountParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = list()
  output[['Name']] = 'TEQC_Result'
  output[[
    'Report [File]'
  ]] = 'p1001/QC_Teqc_5579_2015-05-04--13-41-58/TEQC_Result'
  output[[
    'Html [Link]'
  ]] = 'p1001/QC_Teqc_5579_2015-05-04--13-41-58/TEQC_Result/00index.html'
  param[['cores']] = 1
  param[['paired']] = "false"
  param[['process_mode']] = 'DATASET'
  param[['name']] = 'TEQC_Result'
  param[['designFile']] = system.file(
    "extdata/genes.bed",
    package = "ezRun",
    mustWork = TRUE
  )
  param[['covUniformityPlot']] = 'true'
  param[['covTargetLengthPlot']] = 'true'
  param[['duplicatesPlot']] = 'true'
  param[['cmdOptions']] = ''
  myApp = EzAppTeqc$new()
  myApp$run(input = input, output = output, param = param)
  setwd(cwd)
})
