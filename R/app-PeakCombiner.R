###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCombinePeaks <- function(input = NA,
                                 output = NA,
                                 param = NA) {
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(Rsubread)
    param$paired <- as.logical(input$meta$paired[1])
    outFolder = output$getColumn("PeakCountResult") |> basename()
    setwdNew(outFolder)
    samples <- input$getNames()
    bamFiles   <- input$getFullPathsList("BAM")
    peakFiles   <- input$getFullPathsList("BED")
    peaksList   <- lapply(peakFiles, import_narrowPeak)
    peaksGrl          <- GRangesList(peaksList)
    all_peaks_combined <- unlist(peaksGrl, use.names = FALSE)
    
    support <- rowSums(sapply(peaksList, function(gr)
        countOverlaps(all_peaks_combined, gr) > 0L))
    mcols(all_peaks_combined)$sampleCount <- support
    
    filtered <- all_peaks_combined[support >= param$minSamples]
    filtered <- sortSeqlevels(filtered)
    filtered <- sort(filtered)
    
    # extract sample name from original peak name
    filtered$sample <- sub("_peak_.*$", "", filtered$name)
    
    # collapse overlapping peaks and keep revmap
    collapsed <- GenomicRanges::reduce(filtered, with.revmap = TRUE)
    idx <- mcols(collapsed)$revmap
    
    # annotate collapsed peaks
    collapsed$collapsed_samples <- sapply(idx, function(i)
        paste(sort(unique(filtered$sample[i])), collapse = ","))
    
    collapsed$n_samples <- sapply(idx, function(i)
        length(unique(filtered$sample[i])))
    
    collapsed$max_fold_enrichment <- sapply(idx, function(i)
        max(filtered$fold_enrichment[i], na.rm = TRUE))
    
    collapsed$mean_fold_enrichment <- sapply(idx, function(i)
        mean(filtered$fold_enrichment[i], na.rm = TRUE))
    
    # optionally remove extra chromosomes
    if (param$skipExtraChr) {
        collapsed <- collapsed[nchar(as.character(seqnames(collapsed))) <= 5]
    }
    
    #ADD PEAK NAMES (shown in IGV)
    collapsed$name <- paste0("peak_", seq_along(collapsed))
    
    # export to BED
    outBed <- "combined_filtered_peaks.bed"
    export(collapsed, outBed, format = "BED")
    
    collapsed_df <- data.frame(collapsed)
    collapsed_df$revmap <- NULL
    ezWrite.table(collapsed_df,
                  file = sub("\\.bed$", ".txt", outBed),
                  row.names = FALSE)
    message("Wrote ",
            length(collapsed),
            " peaks (in ≥",
            param$minSamples,
            " samples) to ",
            outBed)
    
    #Count peaks per sample via Rsubread
    regions_SAF <- data.frame(
        GeneID = collapsed_df$name,
        Chr = collapsed_df$seqnames,
        Start = collapsed_df$start,
        End = collapsed_df$end,
        Strand = collapsed_df$strand,
        ID = collapsed_df$name,
        stringsAsFactors = F
    )
    
    countStats = c()
    for (i in 1:length(bamFiles)) {
        bamFile = bamFiles[[i]]
        result = featureCounts(
            bamFile,
            annot.inbuilt = NULL,
            annot.ext = regions_SAF,
            isGTFAnnotationFile = FALSE,
            useMetaFeatures = FALSE,
            allowMultiOverlap = TRUE,
            isPairedEnd = param$paired,
            requireBothEndsMapped = FALSE,
            checkFragLength = FALSE,
            nthreads = param$cores,
            strandSpecific = 0,
            minMQS = 10,
            readExtension5 = 0,
            readExtension3 = 0,
            read2pos = NULL,
            minOverlap = 5,
            nonOverlap = 5,
            ignoreDup = FALSE,
            splitOnly = FALSE,
            countMultiMappingReads = TRUE,
            fraction = FALSE,
            primaryOnly = TRUE,
            countChimericFragments = FALSE,
            chrAliases = NULL,
            reportReads = NULL,
            reportReadsPath = NULL
        )
        countStats = rbind(countStats, c(result$stat[1, 2], sum(result$stat[, 2])))
        
        peakCounts = data.frame(
            PeakID = rownames(result$counts),
            matchCounts = result$counts[, 1],
            stringsAsFactors = F
        )
        peakCounts_annot <- merge(regions_SAF,
                                  peakCounts,
                                  by.x = 'GeneID',
                                  by.y = 'PeakID')
        peakCounts_annot$ID <- NULL
        colnames(peakCounts_annot)[1] <- 'Geneid'
        peakCounts_annot$Length <- peakCounts_annot$End - peakCounts_annot$Start + 1
        peakCounts_annot <- peakCounts_annot[, c('Geneid',
                                                 'Chr',
                                                 'Start',
                                                 'End',
                                                 'Strand',
                                                 'Length',
                                                 'matchCounts')]
        ezWrite.table(
            peakCounts_annot,
            file = sub(
                '_processed.bam$',
                '_peak_counts.txt',
                basename(bamFile)
            ),
            row.names = FALSE
        )
    }
    countStats <- data.frame(
        Name = samples,
        AssignedReads = countStats[, 1],
        TotalReads = countStats[, 2],
        peakFraction = countStats[, 1] / countStats[, 2]
    )
    ezWrite.table(countStats, file = 'peakCountStats.txt', row.names = FALSE)
    
    jsonFileName <- paste0(outFolder, "/igv_session.json")
    htmlFileName <- paste0(outFolder, "/igv_session.html")
    bwFiles   <- input$getColumn("BigWigFile")
    ##TODO: bed peakFile
    jsonFile <- writeIgvSessionFile(param, outFolder, jsonFileName = jsonFileName, bigwigFiles = bwFiles,
                                    baseUrl = file.path(PROJECT_BASE_URL))
    writeNfCoreIgvHtml(param, jsonFile, title = "MultiSample Coverage Tracks", htmlTemplate = "templates/igvNfCoreTemplate.html", htmlFileName = basename(htmlFileName))
    
    return('success')
}


##' @description write IGV session in json format
writeIgvSessionFile <- function(param, outFolder, jsonFileName, bigwigFiles, baseUrl){
    refBuildName = param$ezRef@refBuildName
    refUrlBase = file.path(REF_HOST, param$ezRef@refBuild)
    fastaUrl = sub("Annotation.*", "Sequence/WholeGenomeFasta/genome.fa", refUrlBase)
    faiUrl = paste0(fastaUrl, ".fai")
    
    trackNames <- basename(bigwigFiles) |> str_replace("\\..*", "")
    tracks <- list()
    tracks[[1]] <- list(type=	"sequence")
    for (i in 1:length(bigwigFiles)){
        tracks[[i+1]] <- list(id = trackNames[[i]],
                              url = paste(baseUrl,bigwigFiles[[i]],sep = "/"),
                              format =	"bigWig",
                              name	= trackNames[[i]])
        
    }
    tracks[[i+2]] <- list(
        id = "peaks",
        url = file.path(baseUrl, outFolder,'combined_filtered_peaks.bed'),
        format =	"bed",
        type = "annotation",
        name = "consensusPeaks")
    tracks[[i+3]] <- list(
        id = "genes",
        url = file.path(REF_HOST, param$ezRef@refBuild,'Genes/transcripts.only.gtf'),
        format =	"gtf",
        type = "annotation",
        name = "genes")
    
    tracks[[i+4]] <- list(
        id = "exons",
        url = file.path(REF_HOST, param$ezRef@refBuild,'Genes/genes.bed'),
        format =	"bed",
        type = "annotation",
        name = "exons")
    jsonLines <- list( version =	"3.5.3",
                       showSampleNames = FALSE,
                       reference = list(id = refBuildName , fastaUrl = fastaUrl, indexURL = faiUrl),
                       tracks = tracks
    )
    jsonFile <- rjson::toJSON(jsonLines, indent=5, method="C")
    write(jsonFile, basename(jsonFileName))
    return(jsonFile)
}




#' @template app-template
EzAppPeakCombiner <-
    setRefClass(
        "EzAppPeakCombiner",
        contains = "EzApp",
        methods = list(
            initialize = function() {
                "Initializes the application using its specific defaults."
                runMethod <<- ezMethodCombinePeaks
                name <<- "EzAppPeakCombiner"
                appDefaults <<- rbind(
                    minSamples = ezFrame(
                        Type = "integer",
                        DefaultValue = 2L,
                        Description = "minimum number of samples a peak must be present in to be retained"
                    ),
                    skipExtraChr = ezFrame(
                        Type = "logical",
                        DefaultValue = TRUE,
                        Description = "whether to skip chromosomes with names longer than 5 characters"
                    )
                )
            }
        )
    )

import_narrowPeak <- function(f) {
    import(
        f,
        format    = "BED",
        extraCols = c(
            fold_enrichment = "numeric",
            negLog10_pvalue = "numeric",
            negLog10_qvalue = "numeric",
            peakScore       = "integer"
        )
    )
}