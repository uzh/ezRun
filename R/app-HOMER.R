###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppHomerDiffPeaks <-
  setRefClass(
    "EzAppHomerDiffPeaks",
    contains = "EzApp",
    methods = list(
      ## HOMER behaviour quoted here was checked against the HOMER 5.1 usage text
      ## and getDifferentialPeaksReplicates.pl source; ChIPpeakAnno against 3.46.0.
      methods_facts = function() {
        c(
          ## ezMethodHomerDiffPeaks (app-HOMER.R:74-96); makeTagDirectory (app-HOMER.R:388-415)
          "A HOMER tag directory was made from each BAM file with makeTagDirectory (SAM input, refBuildHOMER genome, no GC check); samples of the sample group were the target and samples of the reference group the background, and fold changes are target over background.",
          ## app-HOMER.R:98 (||), 101-119
          "When the sample group or the reference group had at least two samples, differential peaks were computed with getDifferentialPeaksReplicates.pl using DESeq2 (-DESeq2, fixed regardless of degMethod); otherwise the replicate-free getDifferentialPeaks was used.",
          ## getDifferentialPeaksReplicates.pl Step1/Step2 (no -i, -p or -use given)
          "With replicates, putative peaks were called with findPeaks (-style style) on one pooled tag directory of all target samples, without an input control, and reads per sample were counted in these peaks with annotatePeaks.pl and annotated against the HOMER genome.",
          ## app-HOMER.R:106; getDifferentialPeaksReplicates.pl $norm2total
          "With replicates, when balanced is true, signal was normalized across the peaks by DESeq2 (-balanced) rather than to the total mapped reads per sample, which is used otherwise.",
          ## app-HOMER.R:105 (-all), 132-135
          "With replicates, HOMER reported all peaks (-all) and ezRun kept those with absolute log2 fold change at least log2(repFoldChange) and DESeq2 adjusted p-value at most repFDR, in both directions.",
          ## app-HOMER.R:219-277; mergePeaks -d, getDifferentialPeaks -F/-P usage
          "Without replicates and with peakMode true (set only through specialOptions), peaks were called with findPeaks in each sample separately and merged with mergePeaks at a maximum centre distance of 100 bp before comparison with getDifferentialPeaks (-F 0 -P 1), keeping peaks higher in the target with fold change at least repFoldChange and uncorrected Poisson p-value at most repFDR; this path reads peaks.txt, which findPeaks writes only for style factor or dnase, so with other styles (including the default histone) no differential peaks were produced.",
          ## app-HOMER.R:279-326; annotatePeaks.pl tss default size 4000 bp
          "Without replicates and with peakMode false, transcription start site regions of 4000 bp (plus or minus 2 kb, HOMER annotatePeaks.pl tss default) were compared with getDifferentialPeaks in both directions (-F 0 -P 1), keeping regions with fold change at least repFoldChange and uncorrected Poisson p-value at most repFDR.",
          ## app-HOMER.R:190-218, 328-345; ChIPpeakAnno nearestStart -> nearestLocation, PeakLocForDistance default "start"
          "Without replicates, peaks were assigned to the gene with the nearest transcription start site (ChIPpeakAnno annotatePeakInBatch, one gene per peak, distance from peak start) using GENCODE gene annotation release 37 for hg38, M23 for mm10 and M37 for mm39."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodHomerDiffPeaks
        name <<- "EzAppHomerDiffPeaks"
        appDefaults <<- rbind(
          refBuildHOMER = ezFrame(
            Type = "character",
            DefaultValue = "hg38",
            Description = "Genome version to use from HOMER: hg38, mm10, danRer10, etc."
          ),
          repFoldChange = ezFrame(
            Type = "numeric",
            DefaultValue = 2,
            Description = "Replicate fold change cutoff for peak identification (calculated by DESeq2)"
          ),
          repFDR = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "Replicate FDR cutoff for peak identification (calculated by DESeq2)"
          ),
          balanced = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Do not force the use of normalization factors to match total mapped reads.  This can be useful when analyzing differential peaks between similar data (for example H3K27ac) where we expect similar levels in all experiments. Applying this allows the data to essentially be quantile normalized during the differential calculation."
          ),
          style = ezFrame(
            Type = "character",
            DefaultValue = "histone",
            Description = "Style of peaks found by findPeaks during features selection (factor, histone, super, groseq, tss, dnase, mC)"
          ),
          peakMode = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Call Peaks for differential analysis without replicated. Otherwise check tss regions only."
          ),
          cmdOptions = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "to define batches in the analysis to perform paired test, e.g. -batch 1 2 1 2"
          )
        )
      }
    )
  )

ezMethodHomerDiffPeaks = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  require(parallel)
  require(rtracklayer)
  require(ChIPpeakAnno)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))

  stopifnot(param$sampleGroup != param$refGroup)

  bamFiles <- input$getFullPaths("BAM")
  targetSamples <- input$getNames()[
    input$getColumn(param$grouping) %in%
      param$sampleGroup
  ]
  baseSamples <- input$getNames()[
    input$getColumn(param$grouping) %in%
      param$refGroup
  ]
  targetTagDirs <- mcmapply(
    makeTagDirectory,
    inBam = bamFiles[targetSamples],
    outputDir = paste0(targetSamples, "_tag"),
    MoreArgs = list(genome = param$refBuildHOMER),
    mc.cores = min(param$cores, 4)
  )
  baseTagDirs <- mcmapply(
    makeTagDirectory,
    inBam = bamFiles[baseSamples],
    outputDir = paste0(baseSamples, "_tag"),
    MoreArgs = list(genome = param$refBuildHOMER),
    mc.cores = min(param$cores, 4)
  )

  if (length(targetSamples) >= 2L || length(baseSamples) >= 2L) {
    ## The experiments with replicates
    outputFile <- basename(output$getColumn("DiffPeak"))
    cmd <- paste(
      "getDifferentialPeaksReplicates.pl -DESeq2",
      "-genome",
      param$refBuildHOMER,
      "-all",
      ifelse(param$balanced, "-balanced", ""),
      "-style",
      param$style
    )
    cmd <- paste(
      cmd,
      "-t",
      paste(targetTagDirs, collapse = " "),
      "-b",
      paste(baseTagDirs, collapse = " "),
      param$cmdOptions,
      "> fullResult.tsv"
    )
    ezSystem(cmd)

    homerResult <- ezRead.table("fullResult.tsv", row.names = NULL)
    homerResult[[1]] <- NULL
    if (nrow(homerResult) > 0) {
      colnames(homerResult) <- gsub('Tag.*', '[Signal]', colnames(homerResult))
      colnames(homerResult) <- gsub('bg vs. target ', '', colnames(homerResult))
      colnames(homerResult) <- gsub(
        'adj. p-value',
        'fdr',
        colnames(homerResult)
      )

      doKeep <- abs(homerResult[['Log2 Fold Change']]) >=
        log2(param$repFoldChange) &
        homerResult[['fdr']] <= param$repFDR
      homerResult <- homerResult[doKeep, ]
      if (nrow(homerResult) > 0) {
        resultFile <- paste0(
          param$sampleGroup,
          '_over_',
          param$refGroup,
          '_',
          sub('txt$', 'xlsx', outputFile)
        )
        writexl::write_xlsx(homerResult, resultFile)
        ezWrite.table(homerResult, outputFile, row.names = FALSE)

        peakBedFile <- 'HomerPeaks.bed'
        refFasta <- file.path(
          Sys.getenv("HOMER"),
          'data/genomes',
          param$refBuildHOMER,
          'genome.fa'
        )
        fai <- ezRead.table(paste0(refFasta, ".fai"), header = FALSE)
        colnames(fai) <- c("LENGTH", "OFFSET", "LINEBASES", "LINEWDITH")
        homerResult <- homerResult[homerResult$Chr %in% rownames(fai), ]
        bed <- data.frame(
          chr = homerResult$Chr,
          start = homerResult$Start,
          end = homerResult$End,
          name = homerResult[['Gene Name']],
          score = homerResult[['Peak Score']],
          strand = homerResult$Strand
        )
        bed <- bed[bed$start > 0, ]
        ezWrite.table(bed, peakBedFile, row.names = FALSE, col.names = FALSE)

        peakSeqFile <- 'HomerPeaks.fa'

        cmd <- paste(
          "bedtools",
          " getfasta -fi",
          refFasta,
          "-bed ",
          peakBedFile,
          "-name -fo ",
          peakSeqFile
        )
        ezSystem(cmd)
      }
    } else {
      ezWrite.table(homerResult, outputFile, row.names = FALSE)
    }
  } else {
    ## The experiments without replicates;
    ## focus on tss regions
    #cmd <- paste("annotatePeaks.pl tss", param$ezRef["refFastaFile"],
    #             "-gtf", param$ezRef["refFeatureFile"], "> tss.txt")

    param$refBuild <- switch(
      param$refBuildHOMER,
      'hg38' = 'Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_37-2021-05-04',
      'mm10' = 'Mus_musculus/GENCODE/GRCm38.p6/Annotation/Release_M23-2019-11-05',
      'mm39' = 'Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03',
      stop("unsupported refBuildHOMER")
    )
    param$ezRef <- NULL
    param <- ezParam(param)
    localAnnotation <- ezFeatureAnnotation(param, dataFeatureType = "gene")
    localAnnotation <- unique(localAnnotation[, grep(
      '^gene_id$|^description$|name$|symbol$|^type$',
      colnames(localAnnotation),
      ignore.case = TRUE
    )])

    gtfFile = param$ezRef["refFeatureFile"]
    gtf = rtracklayer::import(gtfFile)
    idx = gtf$type == 'gene'
    if (!any(idx)) {
      idx = gtf$type == 'start_codon'
    }
    gtf = gtf[idx]
    if (grepl('gtf$', gtfFile)) {
      names_gtf = make.unique(gtf$'gene_id')
    } else {
      names_gtf = make.unique(gtf$'ID')
    }
    names(gtf) = names_gtf
    if (param$peakMode) {
      cmd <- paste("findPeaks", targetTagDirs, "-style", param$style, "-o auto")
      system(cmd)
      cmd <- paste("findPeaks", baseTagDirs, "-style", param$style, "-o auto")
      system(cmd)
      cmd <- paste(
        "mergePeaks -d 100",
        file.path(targetTagDirs, "peaks.txt"),
        file.path(baseTagDirs, "peaks.txt"),
        "> mergedPeakFile.txt"
      )
      system(cmd)
      cmd <- paste(
        "getDifferentialPeaks mergedPeakFile.txt",
        targetTagDirs,
        baseTagDirs,
        param$cmdOptions,
        '-F 0',
        '-P 1',
        "> fullResult.tsv"
      )
      system(cmd)
      cmd <- paste(
        "annotatePeaks.pl fullResult.tsv",
        param$refBuildHOMER,
        "-annStats annoStats.txt",
        "> annotatedPeaks.txt"
      )
      system(cmd)
      peakAnnot <- ezRead.table('annotatedPeaks.txt', row.names = NULL)
      colnames(peakAnnot)[1] <- '#PeakID'
      toRemove <- c(
        'Chr',
        'Start',
        'End',
        'Strand',
        'Peak Score',
        'Focus Ratio/Region Size'
      )
      toKeep <- colnames(peakAnnot)[!(colnames(peakAnnot) %in% toRemove)]
      peakAnnot <- peakAnnot[, toKeep]
      skip <- grep("^#PeakID", readLines('fullResult.tsv', n = 200))
      peakStats <- ezRead.table(
        'fullResult.tsv',
        row.names = NULL,
        skip = skip - 1
      )
      homerResult <- merge(
        peakAnnot,
        peakStats,
        by.x = '#PeakID',
        by.y = '#PeakID'
      )
      homerResult <- homerResult[order(homerResult[['p-value']]), ]
      homerResult <- unique(homerResult)
      homerResult <- homerResult[
        homerResult[['Fold Change vs. Background']] >= param$repFoldChange,
      ]
      homerResult <- homerResult[homerResult[['p-value']] <= param$repFDR, ]
      file.remove("mergedPeakFile.txt")
    } else {
      cmd <- paste("annotatePeaks.pl tss", param$refBuildHOMER, "> tss.txt")
      ezSystem(cmd)
      cmd <- paste(
        "getDifferentialPeaks",
        "tss.txt",
        targetTagDirs,
        baseTagDirs,
        param$cmdOptions,
        '-F 0',
        '-P 1',
        "> fullResult.tsv"
      )
      ezSystem(cmd)

      cmd <- paste(
        "getDifferentialPeaks",
        "tss.txt",
        baseTagDirs,
        targetTagDirs,
        param$cmdOptions,
        '-F 0',
        '-P 1',
        "> fullResult_reverse.tsv"
      )
      ezSystem(cmd)

      skip = grep("^#PeakID", readLines('fullResult.tsv', n = 200))
      homerResult <- ezRead.table('fullResult.tsv', skip = skip - 1)
      homerResult <- unique(homerResult)

      homerResult <- homerResult[
        homerResult[['Fold Change vs. Background']] >= param$repFoldChange,
      ]
      homerResult <- homerResult[homerResult[['p-value']] <= param$repFDR, ]

      homerResultRev <- ezRead.table('fullResult_reverse.tsv', skip = skip - 1)
      homerResultRev <- unique(homerResultRev)
      homerResultRev <- homerResultRev[
        homerResultRev[['Fold Change vs. Background']] >= param$repFoldChange,
      ]
      homerResultRev <- homerResultRev[
        homerResultRev[['p-value']] <= param$repFDR,
      ]
      homerResultRev[['Fold Change vs. Background']] <- 1 /
        homerResultRev[['Fold Change vs. Background']]
      homerResult <- rbind(homerResult, homerResultRev)
      file.remove("tss.txt")
    }
    if (nrow(homerResult) > 1) {
      peaksRD = makeGRangesFromDataFrame(homerResult, keep.extra.columns = TRUE)
      names(peaksRD) = mcols(peaksRD)$name
      annotatedPeaks <- annotatePeakInBatch(
        peaksRD,
        AnnotationData = gtf,
        output = 'nearestStart',
        multiple = FALSE,
        FeatureLocForDistance = 'TSS'
      )
      annotatedPeaks = as.data.frame(annotatedPeaks)
      homerResult = merge(
        annotatedPeaks,
        localAnnotation,
        by.x = 'feature',
        by.y = 'gene_id',
        all.x = T
      )
      comparison <- paste0(param$sampleGroup, '_over_', param$refGroup)
      resultFile <- paste0(comparison, '_diffPeaks.xlsx')
      homerResult <- homerResult[order(homerResult[['p.value']]), ]
      writexl::write_xlsx(homerResult, resultFile)
      ezWrite.table(
        homerResult,
        basename(output$getColumn("DiffPeak")),
        row.names = FALSE
      )

      peakBedFile <- paste0(comparison, '_diffPeaks.bed')
      bed <- data.frame(
        chr = homerResult$seqnames,
        start = homerResult$start,
        end = homerResult$end,
        name = rownames(homerResult),
        score = homerResult$score,
        strand = homerResult$strand
      )
      bed <- bed[bed$start > 0, ]
      ezWrite.table(bed, peakBedFile, row.names = FALSE, col.names = FALSE)

      peakSeqFile <- paste0(comparison, '_diffPeaks.fa')
      cmd <- paste(
        "bedtools",
        " getfasta -fi",
        param$ezRef['refFastaFile'],
        "-bed ",
        peakBedFile,
        "-name -fo ",
        peakSeqFile
      )
      ezSystem(cmd)
    } else {
      cmd <- paste('touch', basename(output$getColumn("DiffPeak")))
      ezSystem(cmd)
    }
  }

  return("Success")
}

makeTagDirectory <- function(
  inBam,
  outputDir,
  genome = NULL,
  checkGC = FALSE,
  isAntisense = FALSE,
  strandedPaired = FALSE
) {
  localSamFile <- sub(".bam$", ".sam", basename(inBam))
  stopifnot(!file.exists(localSamFile))
  ezSystem(paste('samtools view -h', inBam, '>', localSamFile))
  cmd <- paste("makeTagDirectory", outputDir, localSamFile, "-format sam")
  if (!is.null(genome)) {
    cmd <- paste(cmd, "-genome", genome)
  }
  if (isTRUE(checkGC)) {
    cmd <- paste(cmd, "-checkGC")
  }
  if (isTRUE(isAntisense)) {
    cmd <- paste(cmd, "-flip")
  }
  if (isTRUE(strandedPaired)) {
    cmd <- paste(cmd, "-sspe")
  }
  ezSystem(cmd)
  file.remove(localSamFile)
  return(outputDir)
}
