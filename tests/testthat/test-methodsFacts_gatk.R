## methods_facts() of EzAppGatkRnaHaplotyper, EzAppJoinGenoTypesRNASeq and EzAppRnaBamStats:
## unconditional facts for list(), gated facts follow the parameters / input, and every literal
## a fact states is still in the code of the function it names.

## deparse() re-wraps long calls when the package was loaded without srcrefs: squash whitespace
codeOf <- function(f) gsub("( ", "(", gsub("\\s+", " ", paste(deparse(f), collapse = " ")), fixed = TRUE)
has <- function(facts, pattern) any(grepl(pattern, facts, fixed = TRUE))
withInput <- function(param, df) { attr(param, "input") <- df; param }

test_that("GatkRnaHaplotyper facts follow addReadGroup and dbsnpFile", {
  app <- EzAppGatkRnaHaplotyper$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_true(has(f0, "exome.bed"))
  expect_true(has(f0, "no exome capture or target file"))
  expect_true(has(f0, "-dont-use-soft-clipped-bases"))
  expect_true(has(f0, "only when dbsnpFile names an existing file"))
  expect_false(has(f0, "AddOrReplaceReadGroups"))
  expect_false(has(f0, "BaseRecalibrator"))

  off <- app$methods_facts(list(addReadGroup = FALSE, markDuplicates = FALSE, dbsnpFile = ""))
  expect_false(has(off, "AddOrReplaceReadGroups"))
  expect_false(has(off, "BaseRecalibrator"))
  expect_true(has(off, "only renamed recal.bam"))
  expect_false(has(off, "only when dbsnpFile names"))

  on <- app$methods_facts(list(addReadGroup = TRUE, dbsnpFile = "/ref/dbsnp.vcf.gz"))
  rg <- grep("AddOrReplaceReadGroups", on, value = TRUE)
  expect_length(rg, 1)
  expect_match(rg, "no read-group ID (-ID) was passed", fixed = TRUE)
  expect_true(has(on, "BaseRecalibrator (--known-sites /ref/dbsnp.vcf.gz"))
  expect_false(has(on, "only renamed recal.bam"))

  code <- codeOf(ezMethodGatkRnaHaplotyper)
  for (lit in c('gtf$gene_biotype %in% c("lncRNA", "protein_coding")', 'gtf$type == "exon"',
                'grepl("chr", seqnames(gtf))', '"exome.bed"', '"-LB"', '"-PU"', '"-SM"', '"ILLUMINA"',
                "stopifnot(param$markDuplicates == FALSE)", "SplitNCigarReads", "--use-original-qualities",
                "ApplyBQSR", "mv splitNtrim.bam recal.bam", "-ERC GVCF", "-dont-use-soft-clipped-bases",
                'paste("--dbsnp", dbsnpFile)'))
    expect_true(grepl(lit, code, fixed = TRUE), label = lit)
  expect_false(grepl('"-ID"', code, fixed = TRUE))
  expect_false(grepl("standard-min-confidence", codeOf(removeSource(ezMethodGatkRnaHaplotyper)), fixed = TRUE))
})

test_that("JoinGenoTypesRNASeq facts carry grouping, minReadDepth and the sample count", {
  app <- EzAppJoinGenoTypesRNASeq$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_true(has(f0, "exactly 500,000"))
  expect_true(has(f0, "ward.D2"))
  expect_true(has(f0, "below minReadDepth were"))

  p <- list(grouping = "Condition", minReadDepth = 20)
  f <- app$methods_facts(p)
  expect_true(has(f, "group of the Condition column"))
  expect_true(has(f, "below 20 (minReadDepth)"))
  big <- app$methods_facts(withInput(p, data.frame(Name = paste0("s", 1:5))))
  small <- app$methods_facts(withInput(p, data.frame(Name = paste0("s", 1:3))))
  expect_true(has(big, "ward.D2"))
  expect_false(has(small, "ward.D2"))
  expect_true(has(small, "No sample clustering"))
  expect_false(has(big, "No sample clustering"))

  code <- codeOf(ezMethodJoinGenoTypesRNASeq)
  for (lit in c("bcftools view -m2 -M2 -v snps", "target_n <- 500000L", "seed <- 1L",
                "stop(\"target_n > n_total\")", "genotype$DP < param$minReadDepth", "allSamples",
                "Mpileup.Rmd"))
    expect_true(grepl(lit, code, fixed = TRUE), label = lit)
  pipe <- codeOf(runGatkPipelineRNASeq)
  expect_true(grepl("CombineGVCFs", pipe, fixed = TRUE))
  expect_true(grepl("nrow(datasetCase) > 1", pipe, fixed = TRUE))
  expect_false(grepl("dbsnp|-L |VariantFiltration|VariantRecalibrator", pipe))
  rmd <- paste(readLines(system.file("templates/Mpileup.Rmd", package = "ezRun")), collapse = "\n")
  for (lit in c('method="ward.D2"', 'c("0/0", "0/1", "1/1")) -2', "eval=length(samples) > 3",
                "if (all(!is.na(d)))", 'grepl("[\\\\._]", names(chromSizes)) & chromSizes > 20000'))
    expect_true(grepl(lit, rmd, fixed = TRUE), label = lit)
})

test_that("RnaBamStats facts follow strandMode, paired, the optional steps and the input", {
  app <- EzAppRnaBamStats$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_true(has(f0, "at least 10 bases"))
  expect_true(has(f0, "101 evenly spaced"))
  expect_true(has(f0, "No IGV session link"))
  expect_false(has(f0, "Fragment sizes"))
  expect_false(has(f0, "Reads were matched to features on"))

  base <- list(strandMode = "antisense", paired = TRUE, fragSizeMax = 500, keepProperPairsOnly = TRUE,
               posErrorRates = TRUE, dupRadar = TRUE)
  f <- app$methods_facts(base)
  expect_true(has(f, "opposite strand (strandMode antisense)"))
  expect_true(has(f, "first reads of proper pairs, up to 500 bp"))
  expect_true(has(f, "at most 100,000"))
  expect_true(has(f, "MarkDuplicates"))
  expect_true(has(app$methods_facts(modifyList(base, list(strandMode = "both"))), "either strand (strandMode both)"))
  se <- app$methods_facts(modifyList(base, list(paired = FALSE)))
  expect_false(has(se, "Fragment sizes"))
  expect_false(has(se, "first and second reads"))
  offs <- app$methods_facts(modifyList(base, list(posErrorRates = FALSE, dupRadar = FALSE)))
  expect_false(has(offs, "mismatch rates"))
  expect_false(has(offs, "dupRadar"))

  ## input: aligner files reused only when the columns exist
  star <- withInput(base, data.frame(`Junctions [File]` = "a_SJ.out.tab", `DupRate [File]` = "a_dupRate.txt",
                                     `StrandFile [Link,File]` = "a_strand.txt", check.names = FALSE))
  bare <- withInput(base, data.frame(`BAM [File]` = "a.bam", check.names = FALSE))
  fs <- app$methods_facts(star); fb <- app$methods_facts(bare)
  expect_true(has(fs, "SJ.out.tab"))
  expect_false(has(fb, "SJ.out.tab"))
  expect_true(has(fb, "junction_annotation.py (--mapq=1)"))
  expect_true(has(fs, "DupRate file"))
  expect_false(has(fb, "DupRate file"))
  expect_true(has(fs, "infer_experiment"))
  expect_false(has(fb, "infer_experiment"))

  code <- codeOf(ezMethodRnaBamStats)
  for (lit in c('getUpstreamQcFiles(input, "Junctions")', 'getUpstreamQcFiles(input, "DupRate")',
                "hasSameUpstreamSettings(input, param)", "ezReadRSeQCStrandness", "Repeats/RNA_repeats.gff"))
    expect_true(grepl(lit, code, fixed = TRUE), label = lit)
  for (lit in c('c("refFeatureFile", "strandMode", "paired")'))
    expect_true(grepl(lit, codeOf(hasSameUpstreamSettings), fixed = TRUE), label = lit)
  par <- codeOf(getStatsFromBamParallel)
  for (lit in c('grepl("^NC_", seqNames)', "nchar(seqNames) <= 6", "head(100)"))
    expect_true(grepl(lit, par, fixed = TRUE), label = lit)
  tt <- codeOf(getTargetTypeCounts)
  for (lit in c("minoverlap = 10", "flank(msgRanges, 2000)", 'ensemblTypes == "protein_coding"', "unannotated"))
    expect_true(grepl(lit, tt, fixed = TRUE), label = lit)
  expect_true(grepl('antisense = flipStrand(strandValues)', codeOf(fixStrand), fixed = TRUE))
  gs <- codeOf(getStatsFromBam)
  for (lit in c('c("protein_coding", "mRNA")', 'gff$tsl %in% "5"', "tCount <= 2", "length.out = 101",
                "breaks = c(0.5, 10, 90, 99.5)", "breaks = c(599, 1200, 2400)", "sum(isInLc) > 40",
                "quantile(avgCov, c(0.25, 0.75))"))
    expect_true(grepl(lit, gs, fixed = TRUE), label = lit)
  pe <- codeOf(getPosErrorFromBam)
  expect_true(grepl("which.max(seqLengths)", pe, fixed = TRUE))
  pr <- codeOf(ezPosSpecErrorRate)
  expect_true(grepl("nMaxReads = 1e+05", pr, fixed = TRUE))
  expect_true(grepl('grepl("N|I|D", bam$cigar)', pr, fixed = TRUE))
  sj <- codeOf(getJunctionStatsFromSJ)
  expect_true(grepl("minIntronSize = 50", sj, fixed = TRUE))
  expect_true(grepl("seq(0.05, 1, by = 0.05)", sj, fixed = TRUE))
  jb <- codeOf(getJunctionPlotsFromBam)
  expect_true(grepl('"--mapq=1"', jb, fixed = TRUE))
  expect_true(grepl("nSim = 10", jb, fixed = TRUE))
  du <- codeOf(dupBam)
  expect_true(grepl("dupDistance = 100", du, fixed = TRUE))
  expect_true(grepl("MarkDuplicates", du, fixed = TRUE))
  expect_true(grepl("MarkDuplicates|markdup", codeOf(bamHasMarkedDuplicates), fixed = TRUE))
  expect_true(grepl("analyzeDuprates", codeOf(getDupRateFromBam), fixed = TRUE))
  expect_true(grepl('tag = "NH"', codeOf(getBamMultiMatchingFromNH), fixed = TRUE))
  expect_true(grepl("nReads - sum(result)", codeOf(getBamMultiMatching), fixed = TRUE))
  rmd <- readLines(system.file("templates/RNABamStats.Rmd", package = "ezRun"))
  expect_true(any(grepl("```{r igv link, include=FALSE, eval=FALSE}", rmd, fixed = TRUE)))
})
