###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMacs3 = function(input = NA, output = NA, param = NA) {
  opt = paste(param$cmdOptions, '-q', param$qValue)
  if (param$paired) {
    opt = paste(opt, '-f BAMPE')
  }
  ## With BAMPE file, --shift cannot be set.
  dataset = input$meta

  ## -g option: mappable genome size
  ## TODO: the MACS3 defaults should be used or the user should be asked
  if (param$genomeSize == 0) {
    gsize <- sum(as.numeric(fasta.seqlengths(param$ezRef["refFastaFile"])))
    gsize <- round(gsize * 0.8)
    ezLog("Use calculated gsize: ", gsize)
    opt <- paste(opt, "-g", as.character(gsize))
  }
  ## --keep-dup: behavior towards duplicate tags at the exact same location
  if (!grepl("--keep-dup", opt)) {
    opt <- paste(opt, "--keep-dup all")
  }

  if (param$mode == "ChIP-seq") {
    ## --extsize: extend reads in 5'->3' direction to fix-sized fragments when model building is deactivated.
    ## This size is taken from sushi app.
    if (grepl("--nomodel", opt) && !grepl("--extsize", opt)) {
      opt <- paste(opt, "--extsize 147")
    }
    bamFile <- input$getFullPaths("BAM")
    outBam <- basename(output$getColumn("BAM"))
    if (param$removeDuplicates) {
      removeDuplicatesFromBam(
        inBam = bamFile,
        outBam = outBam,
        ram = param$ram,
        cores = param$cores
      )
    } else {
      file.copy(from = bamFile, to = outBam, overwrite = TRUE)
      Rsamtools::indexBam(outBam)
    }

    if (isTRUE(param$useControl)) {
      if (!any(grepl("Control", input$colNames, ignore.case = TRUE))) {
        stop(
          "The parameter 'useControl' is 'true' but no column named 'Control [File]' is available."
        )
      }

      cmd = paste(
        "macs3",
        "callpeak -t",
        outBam,
        "-c",
        input$getFullPaths("Control"),
        "-B",
        opt,
        "-n",
        output$getNames()
      )
      ezSystem(cmd)
      bedgraphFileTreat = paste0(output$getNames(), '_treat_pileup.bdg')
      bedgraphFileControl = paste0(output$getNames(), '_control_lambda.bdg')
      cmd = paste(
        "macs3",
        " bdgcmp -t",
        bedgraphFileTreat,
        "-c",
        bedgraphFileControl,
        "-o",
        paste0(output$getNames(), "_FE.bdg"),
        "-m FE"
      )
      ezSystem(cmd)
      bdgSorted = "sorted.bdg"
      ezSystem(paste(
        "bedSort",
        paste0(output$getNames(), "_FE.bdg"),
        bdgSorted
      ))
      ezSystem(paste("samtools idxstats", outBam, "> idxStats.txt"))
      idxStats <- ezRead.table("idxStats.txt", header = FALSE)
      idxStats <- idxStats[rownames(idxStats) != "*", ] ## remove the unmapeed
      ezWrite.table(
        idxStats %>% dplyr::select(1),
        file = "chromSizes.txt",
        col.names = FALSE
      )
      cmd = paste(
        "bedGraphToBigWig",
        bdgSorted,
        "chromSizes.txt",
        paste0(output$getNames(), "_processed.bw")
      )
      ezSystem(cmd)
      ezSystem("rm *.bdg")
    } else {
      cmd = paste("macs3", "callpeak -t", outBam, opt, "-n", output$getNames())
      ezSystem(cmd)
      bam2bw(
        file = outBam,
        destination = basename(output$getColumn("BigWigFile")),
        paired = param$paired,
        method = "deepTools",
        cores = param$cores
      )
    }
  } else if (param$mode == "ATAC-seq") {
    if (!param$paired) {
      stop("For ATAC-seq, we only support paired-end data.")
    }

    ## --extsize: extend reads in 5'->3' direction to fix-sized fragments.
    if (!grepl("--extsize", opt)) {
      ## https://github.com/taoliu/MACS/issues/145
      opt <- paste(opt, "--extsize 200")
    } else {
      opt <- gsub("--extsize 147", "--extsize 200", opt)
    }

    ## Preprocess ATAC-seq bam file
    atacBamProcess(input = input, output = output, param = param)

    cmd = paste(
      "macs3",
      "callpeak -t",
      basename(output$getColumn("BAM")),
      opt,
      "-n",
      output$getNames()
    )
    ezSystem(cmd)
    bam2bw(
      file = basename(output$getColumn("BAM")),
      destination = basename(output$getColumn("BigWigFile")),
      paired = param$paired,
      method = "deepTools",
      cores = param$cores
    )
  } else {
    stop("MACS3 only supports ChIP-seq or ATAC-seq data.")
  }

  peakBedFile = basename(output$getColumn("BED"))
  if (grepl('broad', opt)) {
    file.rename(
      from = paste0(output$getNames(), "_peaks.broadPeak"),
      to = peakBedFile
    )
  } else {
    file.rename(
      from = paste0(output$getNames(), "_peaks.narrowPeak"),
      to = peakBedFile
    )
  }
  peakSeqFile = basename(output$getColumn("PeakSequences"))
  cmd = paste(
    "bedtools",
    " getfasta -fi",
    param$ezRef["refFastaFile"],
    " -bed ",
    peakBedFile,
    " -name -fo ",
    peakSeqFile
  )
  ezSystem(cmd)
  peakXlsFile <- basename(output$getColumn("CalledPeaks"))
  annotatePeaks(peakXlsFile, peakSeqFile, param)

  ## Build a per-sample IGV (igv.js) HTML page: coverage bigWig + called peaks + gene annotation.
  ## The bigWig/peaks track URLs point at the gstore project path and become reachable once the
  ## results are delivered; the genome/gene tracks are served from REF_HOST and always resolve.
  localPeakBed <- basename(output$getColumn("BED"))
  igvJsonFile <- sub("\\.html$", ".json", basename(output$getColumn("IGV")))
  igvJson <- writeMacs3IgvSession(
    param,
    output,
    jsonFileName = igvJsonFile,
    bigwigFile = output$getColumn("BigWigFile"),
    peakBedFile = if (file.exists(localPeakBed)) output$getColumn("BED") else NULL,
    baseUrl = PROJECT_BASE_URL
  )
  writeNfCoreIgvHtml(
    param,
    igvJson,
    title = paste(output$getNames(), "coverage"),
    htmlTemplate = "templates/igvNfCoreTemplate.html",
    htmlFileName = basename(output$getColumn("IGV"))
  )
  return("Success")
}

##' @title Writes a per-sample IGV session for MACS3
##' @description Builds an igv.js session (JSON) for a single MACS3 sample with a genome
##'   sequence track, the coverage bigWig, the called peaks (BED, optional), and gene
##'   annotation (transcripts GTF + exons BED). Modeled on \code{writeIgvSessionFile} in
##'   the PeakCombiner app.
##' @param param the parameter list, providing \code{param$ezRef}.
##' @param output the output \code{EzDataset}, providing the sample name.
##' @param jsonFileName file name to write the JSON session to.
##' @param bigwigFile the gstore-relative path of the coverage bigWig track.
##' @param peakBedFile the gstore-relative path of the called-peaks BED track, or
##'   \code{NULL} to omit the peaks track (e.g. when no peaks were called).
##' @param baseUrl the project base URL the bigWig/peaks paths are relative to.
##' @return the JSON session content as a string (also written to \code{jsonFileName}).
##' @template roxygen-template
writeMacs3IgvSession <- function(
  param,
  output,
  jsonFileName,
  bigwigFile,
  peakBedFile = NULL,
  baseUrl
) {
  refBuildName <- param$ezRef@refBuildName
  refUrlBase <- file.path(REF_HOST, param$ezRef@refBuild)
  fastaUrl <- sub(
    "Annotation.*",
    "Sequence/WholeGenomeFasta/genome.fa",
    refUrlBase
  )
  faiUrl <- paste0(fastaUrl, ".fai")

  sampleName <- output$getNames()
  tracks <- list()
  tracks[[length(tracks) + 1]] <- list(type = "sequence")
  tracks[[length(tracks) + 1]] <- list(
    id = sampleName,
    url = file.path(baseUrl, bigwigFile),
    format = "bigWig",
    name = paste(sampleName, "coverage")
  )
  if (!is.null(peakBedFile)) {
    tracks[[length(tracks) + 1]] <- list(
      id = "peaks",
      url = file.path(baseUrl, peakBedFile),
      format = "bed",
      type = "annotation",
      name = "peaks"
    )
  }
  tracks[[length(tracks) + 1]] <- list(
    id = "genes",
    url = file.path(
      REF_HOST,
      param$ezRef@refBuild,
      "Genes/transcripts.only.gtf"
    ),
    format = "gtf",
    type = "annotation",
    name = "genes"
  )
  tracks[[length(tracks) + 1]] <- list(
    id = "exons",
    url = file.path(REF_HOST, param$ezRef@refBuild, "Genes/genes.bed"),
    format = "bed",
    type = "annotation",
    name = "exons"
  )
  jsonLines <- list(
    version = "3.5.3",
    showSampleNames = FALSE,
    reference = list(id = refBuildName, fastaUrl = fastaUrl, indexURL = faiUrl),
    tracks = tracks
  )
  jsonFile <- rjson::toJSON(jsonLines, indent = 5, method = "C")
  write(jsonFile, jsonFileName)
  return(jsonFile)
}

##' @template app-template
##' @templateVar method ezMethodMacs3(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppMacs3 <-
  setRefClass(
    "EzAppMacs3",
    contains = "EzApp",
    methods = list(
      ## MACS3, samtools, Rsamtools and bedtools unconditional; Picard MarkDuplicates unless
      ## removeDuplicates is false (default TRUE; skipped when the BAM is already marked, data-
      ## dependent); deepTools for ATAC-seq, or ChIP-seq without useControl; UCSC bigWig tools for
      ## ChIP-seq with useControl (R default TRUE); ChIPseeker/ChIPpeakAnno/GenomicFeatures/
      ## rtracklayer unless annotatePeaks is false (default TRUE).
      citation = function(param = list()) {
        atac <- identical(param$mode, "ATAC-seq")
        control <- !isFALSE(as.logical(param$useControl))
        annotate <- !isFALSE(as.logical(param$annotatePeaks))
        c(
          "Zhang, Y. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biology 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137",
          "Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25(16), 2078-2079 (2009). https://doi.org/10.1093/bioinformatics/btp352",
          "Morgan, M. & Pagès, H. Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. R package. https://doi.org/10.18129/B9.bioc.Rsamtools",
          if (!isFALSE(as.logical(param$removeDuplicates))) "Picard Toolkit. Broad Institute. https://broadinstitute.github.io/picard/",
          if (atac || !control) "Ramírez, F. et al. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research 44, W160-W165 (2016). https://doi.org/10.1093/nar/gkw257",
          if (!atac && control) "Kent, W.J. et al. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics 26, 2204-2207 (2010). https://doi.org/10.1093/bioinformatics/btq351",
          "Quinlan, A.R. & Hall, I.M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842 (2010). https://doi.org/10.1093/bioinformatics/btq033",
          if (annotate) "Yu, G. et al. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics 31, 2382-2383 (2015). https://doi.org/10.1093/bioinformatics/btv145",
          if (annotate) "Zhu, L.J. et al. ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 11, 237 (2010). https://doi.org/10.1186/1471-2105-11-237",
          if (annotate) "Lawrence, M. et al. Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology 9(8), e1003118 (2013). https://doi.org/10.1371/journal.pcbi.1003118",
          if (annotate) "Lawrence, M., Gentleman, R. & Carey, V. rtracklayer: an R package for interfacing with genome browsers. Bioinformatics 25(14), 1841-1842 (2009). https://doi.org/10.1093/bioinformatics/btp328"
        )
      },
      methods_facts = function(param = list()) {
        known <- length(param) > 0
        opt <- if (known && ezIsSpecified(param$cmdOptions)) param$cmdOptions else ""
        mode <- if (known) param$mode
        atac <- identical(mode, "ATAC-seq")
        chip <- identical(mode, "ChIP-seq")
        paired <- isTRUE(as.logical(param$paired))
        control <- chip && isTRUE(as.logical(param$useControl))
        dedup <- known && !isFALSE(as.logical(param$removeDuplicates))
        gsize0 <- known && isTRUE(as.numeric(param$genomeSize) == 0)
        c(
          ## ezMethodMacs3: opt = cmdOptions + -q qValue (+ -f BAMPE when paired)
          if (!known) "MACS3 callpeak was run with the q-value cutoff -q qValue, plus -f BAMPE for paired-end data.",
          if (known) paste0("MACS3 callpeak was run with the q-value cutoff -q ", param$qValue,
                            if (paired) " and the paired-end BAM input format -f BAMPE", "."),
          ## ezMethodMacs3: --keep-dup all unless cmdOptions has --keep-dup
          if (!grepl("--keep-dup", opt)) paste0("ezRun added --keep-dup all to MACS3, so MACS3 itself removed no duplicate reads",
                                                if (known) (if (dedup) "; duplicates were removed before peak calling as described below" else "; duplicates were not removed before either (removeDuplicates false)"), "."),
          ## ezMethodMacs3: genomeSize == 0 -> -g round(0.8 * sum(fasta.seqlengths))
          if (gsize0) "The MACS3 genome size (-g) was set by ezRun to 80% of the summed sequence lengths of the reference genome FASTA (genomeSize 0).",
          if (known && !gsize0 && !grepl("(^|\\s)-g\\s", opt)) paste0("The genomeSize ", param$genomeSize, " was not passed to MACS3: ezRun adds -g only when genomeSize is 0, and cmdOptions had no -g."),
          ## ezMethodMacs3 / atacBamProcess -> removeDuplicatesFromBam (bamHasMarkedDuplicates -> samtools -F 1024, else dupBam)
          if (dedup && (atac || chip)) "Before peak calling, duplicate reads were removed: when the BAM header showed duplicates already marked (a MarkDuplicates or markdup program record), the flagged reads were dropped with samtools view -F 1024; otherwise Picard MarkDuplicates was run with REMOVE_DUPLICATES=true.",
          if (known && !dedup && (atac || chip)) "The input BAM was used without duplicate removal (removeDuplicates false), so duplicate reads were not removed.",
          ## atacBamProcess -> filteroutBam(chrs = M, MT, chrM, chrMT; mapQ = 10)
          if (atac) "In ATAC-seq mode (paired-end only), the reads on the mitochondrial chromosomes (M, MT, chrM, chrMT) were removed and only alignments with a MAPQ of at least 10 were kept (samtools view -q 10) before peak calling.",
          ## atacBamProcess: shiftATAC -> alignmentSieve --ATACshift, then sortBam + indexBam
          if (atac && isTRUE(as.logical(param$shiftATAC))) "The filtered ATAC-seq reads were then shifted with deepTools alignmentSieve --ATACshift and sorted again.",
          ## ezMethodMacs3 ChIP-seq branch: no filteroutBam
          if (chip) "In ChIP-seq mode there was no mapping-quality or mitochondrial filter: peaks were called on the input BAM after the duplicate step only.",
          ## ezMethodMacs3 ATAC-seq: --extsize 200 added, or --extsize 147 replaced by 200
          if (atac && !grepl("--extsize", opt)) "In ATAC-seq mode --extsize 200 was added to the MACS3 options.",
          if (atac && grepl("--extsize 147", opt, fixed = TRUE)) "In ATAC-seq mode the --extsize 147 given in cmdOptions was replaced by --extsize 200.",
          ## ezMethodMacs3 ChIP-seq: --nomodel without --extsize -> --extsize 147
          if (chip && grepl("--nomodel", opt) && !grepl("--extsize", opt)) "Because cmdOptions had --nomodel without --extsize, --extsize 147 was added to the MACS3 options.",
          ## ezMethodMacs3 ChIP-seq with useControl: callpeak -c -B, bdgcmp -m FE, bedSort, bedGraphToBigWig
          if (control) "Peaks were called against the Control BAM (-c) with pileup tracks (-B); the coverage bigWig is the fold enrichment of treatment over control (macs3 bdgcmp -m FE), converted with bedSort and bedGraphToBigWig.",
          ## bam2bw(method = "deepTools"): bamCoverage --binSize 10 --normalizeUsing CPM
          if (atac || (chip && !control)) "The coverage bigWig was made from the peak-calling BAM with deepTools bamCoverage in 10 bp bins normalized to counts per million (--binSize 10, --normalizeUsing CPM).",
          ## ezMethodMacs3: BED = broadPeak when the options contain 'broad', else narrowPeak
          if (known && grepl("broad", opt)) "The BED output holds MACS3's broadPeak file.",
          ## ezMethodMacs3: bedtools getfasta -name on the peaks BED
          "The DNA sequence of every peak was extracted with bedtools getfasta -name from the reference genome FASTA.",
          ## annotatePeaks(): ChIPseeker annotatePeak(tssRegion = c(-1000, 1000)), ChIPpeakAnno annotatePeakInBatch(nearestStart, TSS), dustyScore
          if (known && !isFALSE(as.logical(param$annotatePeaks))) "Peaks were annotated with ChIPseeker annotatePeak (TxDb built from the GTF, tssRegion -1000 to 1000) and with ChIPpeakAnno annotatePeakInBatch (output nearestStart, distance to the TSS, one gene per peak), plus a low-complexity dustyScore of each peak sequence (ShortRead); the peak table is sorted by -log10 p-value."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMacs3
        name <<- "EzAppMacs3"
        appDefaults <<- rbind(
          useControl = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "should control samples be used"
          ),
          shiftATAC = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "should all reads aligning to + strand were offset by +4bp, all reads aligning to the - strand are offset -5 bp"
          ),
          annotatePeaks = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "use gtf to annotate peaks"
          ),
          genomeSize = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "genome size"
          ),
          qValue = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "The q-value (minimum FDR) cutoff to call significant regions."
          )
        )
      }
    )
  )

##' @title Annotates peaks
##' @description Annotates peaks and writes them into a separate table.
##' @template input-template
##' @template output-template
##' @param param a list of parameters to extract the \code{ezRef@@refFeatureFile} and the \code{ezRef@@refAnnotationFile} from.
##' @template roxygen-template
annotatePeaks = function(peakFile, peakSeqFile, param) {
  require(rtracklayer)
  require(ChIPpeakAnno)
  require(ShortRead)
  require(ChIPseeker)
  require(GenomicFeatures)

  data <- c()
  tryCatch(
    expr = {
      data = ezRead.table(
        sub('.xlsx$', '.xls', peakFile),
        comment.char = "#",
        row.names = NULL
      )
    },
    error = function(e) {
      ezLog(paste("No peaks detected. Skip peak annotation"))
    }
  )
  if (is.null(data) | nrow(data) == 0) {
    return('No peaks detected. Skip peak annotation')
  }
  data = data[order(data$chr, data$start), ]

  if (!param$annotatePeaks) {
    writexl::write_xlsx(data, peakFile)
    return('done')
  }

  gtfFile = param$ezRef@refFeatureFile
  myTxDB <- txdbmaker::makeTxDbFromGFF(file = gtfFile, format = 'gtf')

  gtf <- rtracklayer::import(gtfFile)
  if ('gene' %in% unique(gtf$type)) {
    idx = gtf$type == 'gene'
  } else if ('transcript' %in% unique(gtf$type)) {
    idx = gtf$type == 'transcript'
  } else if ('start_codon' %in% unique(gtf$type)) {
    idx = gtf$type == 'start_codon'
  } else {
    ezLog('gtf is incompatabible. Peak annotation skipped!')
    return(NULL)
  }
  gtf = gtf[idx]
  if (grepl('gtf$', gtfFile)) {
    names_gtf = make.unique(gtf$'gene_id')
  } else {
    names_gtf = make.unique(gtf$'ID')
  }
  names(gtf) = names_gtf
  peaksRD = makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
  names(peaksRD) = mcols(peaksRD)$name
  annot_ChIPseeker <- annotatePeak(
    peaksRD,
    TxDb = myTxDB,
    tssRegion = c(-1000, 1000),
    verbose = FALSE
  )
  tryCatch(
    {
      annotatedPeaks <- annotatePeakInBatch(
        peaksRD,
        AnnotationData = gtf,
        output = 'nearestStart',
        multiple = FALSE,
        FeatureLocForDistance = 'TSS'
      )
    },
    error = function(e) {
      writexl::write_xlsx(data, peakFile)
      return('no valid peaks for annotation')
    }
  )
  if (!exists('annotatedPeaks')) {
    writexl::write_xlsx(data, peakFile)
    return('no valid peaks for annotation')
  }
  annotatedPeaks = as.data.frame(annotatedPeaks)
  annotatedPeaks = annotatedPeaks[, c(
    "peak",
    "feature",
    "feature_strand",
    "start_position",
    "end_position",
    "insideFeature",
    "distancetoFeature"
  )]
  colnames(annotatedPeaks) = c(
    "peak",
    "feature",
    "feature_strand",
    "feature_start",
    "feature_end",
    "insideFeature",
    "distancetoFeature"
  )

  annotatedPeaks = merge(
    data,
    annotatedPeaks,
    by.x = 'name',
    by.y = 'peak',
    all.x = T
  )
  localAnnotation <- ezFeatureAnnotation(param, dataFeatureType = "gene")
  localAnnotation = unique(localAnnotation[, grep(
    '^gene_id$|^description$|name$|symbol$|^type$',
    colnames(localAnnotation),
    ignore.case = TRUE
  )])
  if (!is.null(ncol(localAnnotation))) {
    annotatedPeaks = merge(
      annotatedPeaks,
      localAnnotation,
      by.x = 'feature',
      by.y = 'gene_id',
      all.x = T
    )
  }
  colnames(annotatedPeaks) = gsub('-log10', '_-log10', colnames(annotatedPeaks))
  seqs = readDNAStringSet(peakSeqFile)
  dustyScores = data.frame(
    ID = names(seqs),
    dustyScore_peakSequence = dustyScore(seqs),
    stringsAsFactors = F
  )
  annotatedPeaks = merge(
    annotatedPeaks,
    dustyScores,
    by.x = 'name',
    by.y = 'ID',
    all.x = TRUE
  )

  annot_ChIPseeker <- data.frame(annot_ChIPseeker@anno)
  keepCol_ChIPSeeker <- c(
    "name",
    "annotation",
    "geneId",
    "transcriptId",
    "distanceToTSS",
    "geneChr",
    "geneStart",
    "geneEnd",
    "geneLength",
    "geneStrand"
  )
  annot_ChIPseeker <- annot_ChIPseeker[, keepCol_ChIPSeeker]
  colnames(annot_ChIPseeker)[2:ncol(annot_ChIPseeker)] <- paste0(
    'ChIPSeeker_',
    colnames(annot_ChIPseeker)[2:ncol(annot_ChIPseeker)]
  )
  annotatedPeaks <- merge(
    annotatedPeaks,
    annot_ChIPseeker,
    by.x = 'name',
    by.y = 'name',
    all.x = TRUE
  )
  annotatedPeaks <- annotatedPeaks[!duplicated(annotatedPeaks$name), ]
  annotatedPeaks <- annotatedPeaks[
    order(annotatedPeaks[['_-log10(pvalue)']], decreasing = TRUE),
  ]
  writexl::write_xlsx(annotatedPeaks, peakFile)
  return('done')
}

### import Macs3's BED6+4 file: narrowPeak or broadPeak
import.Macs3Peaks <- function(con) {
  bed <- ezRead.table(con, header = FALSE, row.names = NULL)
  colnames(bed) <- c(
    "chr",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "fold-change",
    "p-value",
    "q-value",
    "summit"
  )
  bed <- transform(bed, start = start + 1)
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  return(bed)
}
