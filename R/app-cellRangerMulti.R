###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCellRangerMulti <- function(input = NA, output = NA, param = NA) {
  sampleName <- input$getNames()

  dirList <- prepareFastqData(input, param)
  sampleDirs <- dirList$sampleDirs

  #3. Create the multi configuration file and command
  conf <- buildMultiConfigFile(input, param, dirList)
  configFileName <- conf$configFileName
  refDir <- conf$refDir

  cellRangerFolder <- str_sub(sampleName, 1, 45) %>% str_c("-cellRanger")
  #3. Build command
  cmd <- paste(
    "cellranger multi",
    paste0("--id=", cellRangerFolder),
    paste0("--localmem=", param$ram),
    paste0("--localcores=", param$cores),
    paste0("--csv=", configFileName)
  )

  #4. Add additional cellranger options if specified
  if (ezIsSpecified(param$cmdOptions)) {
    cmd <- paste(cmd, param$cmdOptions)
  }

  #5. Execute the command
  ezSystem(cmd)

  #7. Delete temp files and rename the final cellranger output folder
  unlink(basename(sampleDirs), recursive = TRUE)
  if (exists("featureDirs")) {
    unlink(basename(featureDirs))
  }
  file.rename(file.path(cellRangerFolder, "outs"), sampleName)
  unlink(cellRangerFolder, recursive = TRUE)
  if (ezIsSpecified(param$controlSeqs) || ezIsSpecified(param$secondRef)) {
    futile.logger::flog.info(sprintf("Removing %s", refDir))
    unlink(refDir, recursive = TRUE)
  }

  #8. Optional removal of the bam files
  if (!param$keepBam) {
    futile.logger::flog.info(ezSystem(
      'find . -name "*_alignments.bam*" -type f'
    ))
    ezSystem('find . -name "*_alignments.bam*" -type f -delete')
  } else if (ezIsSpecified(param$secondRef) && param$secondRef != '') {
    bamFiles <- ezSystem(
      'find . -name "*_alignments.bam" -type f',
      intern = TRUE
    )
    refDir <- param$ezRef["refFastaFile"]
    for (bamFile in bamFiles) {
      out <- tryCatch(
        ezSystem(paste(
          'samtools view',
          '-T',
          refDir,
          '-@',
          param$cores,
          '-o',
          sub('.bam$', '.cram', bamFile),
          '-C',
          bamFile
        )),
        error = function(e) NULL
      )
      ezSystem(paste('rm', bamFile))
    }
  }

  #9. Generate expanded dataset.tsv:
  libraryTypes <- as.vector(str_split(param$TenXLibrary, ",", simplify = TRUE))
  ds <- output$meta
  samplePath <- file.path(sampleName, 'per_sample_outs')
  subSampleDirs <- list.dirs(samplePath, recursive = FALSE)
  subSamples <- basename(subSampleDirs)

  # Detect CellRanger version from param$CellRangerVersion
  # Format: "Aligner/CellRanger/X.Y.Z" - extract major version
  # v10+ changed output structure: removed count/ subdirectory and multi/count/ path
  crVersion <- as.numeric(
    str_extract(param$CellRangerVersion, "(?<=/)[0-9]+(?=\\.)")
  )
  isV9OrBelow <- crVersion <= 9

  expandedDS <- data.frame(
    Name = subSamples,
    Species = ds$Species,
    refBuild = ds$refBuild,
    refFeatureFile = ds$refFeatureFile,
    featureLevel = ds$featureLevel,
    transcriptTypes = ds$transcriptTypes,
    SCDataOrigin = '10X'
  )
  expandedDS[['ResultDir [File]']] <- file.path(
    param[['resultDir']],
    samplePath,
    subSamples
  )
  expandedDS[['Report [Link]']] <- file.path(
    expandedDS[['ResultDir [File]']],
    'web_summary.html'
  )
  if (isV9OrBelow) {
    # v9 and below: count/ subdirectory exists
    expandedDS[['CountMatrix [Link]']] <- file.path(
      expandedDS[['ResultDir [File]']],
      'count',
      'sample_filtered_feature_bc_matrix'
    )
  } else {
    # v10+: matrices directly in sample dir
    expandedDS[['CountMatrix [Link]']] <- file.path(
      expandedDS[['ResultDir [File]']],
      'sample_filtered_feature_bc_matrix'
    )
  }
  if ("Multiplexing" %in% libraryTypes) {
    if (isV9OrBelow) {
      # v9 and below
      expandedDS[['UnfilteredCountMatrix [Link]']] <- file.path(
        expandedDS[['ResultDir [File]']],
        'count',
        'sample_raw_feature_bc_matrix'
      )
    } else {
      # v10+
      expandedDS[['UnfilteredCountMatrix [Link]']] <- file.path(
        expandedDS[['ResultDir [File]']],
        'sample_raw_feature_bc_matrix'
      )
    }
  } else {
    if (isV9OrBelow) {
      # v9 and below: multi/count/ exists
      expandedDS[['UnfilteredCountMatrix [Link]']] <- file.path(
        param[['resultDir']],
        sampleName,
        'multi/count',
        'raw_feature_bc_matrix'
      )
    } else {
      # v10+: aggregated at top level
      expandedDS[['UnfilteredCountMatrix [Link]']] <- file.path(
        param[['resultDir']],
        sampleName,
        'raw_feature_bc_matrix'
      )
    }
  }
  expandedDS[['Condition [Factor]']] = c('')
  expandedDS[['Order Id [B-Fabric]']] = ds[['Order Id [B-Fabric]']]
  commonPath <- file.path(
    '/srv/GT/analysis/CM_datasets',
    basename(param[['resultDir']])
  )
  if (!dir.exists(commonPath)) {
    dir.create(commonPath)
  }
  dsPath <- file.path(commonPath, 'expanded_dataset.tsv')
  if (!file.exists(dsPath)) {
    ezWrite.table(expandedDS, dsPath, row.names = FALSE)
  } else {
    ezWrite.table(
      expandedDS,
      dsPath,
      row.names = FALSE,
      append = TRUE,
      col.names = FALSE
    )
  }
  ezSystem(paste(
    '/usr/local/ngseq/bin/g-req copynow -f',
    dsPath,
    file.path(param$dataRoot, param[['resultDir']])
  ))
  return("Success")
}

prepareFastqData <- function(input, param) {
  sampleName <- input$getNames()

  #1. Prepare GEX data
  sampleDirs <- getFastqDirs(input, "RawDataDir", sampleName)

  #1.1. decompress tar files if they are in tar format
  if (all(grepl("\\.tar$", sampleDirs))) {
    sampleDirs <- tarExtract(sampleDirs)
  }

  #1.2. Subsample if chosen
  if (ezIsSpecified(param$nReads) && param$nReads > 0) {
    sampleDirs <- sapply(sampleDirs, subsample, param)
  }

  sampleDirs <- normalizePath(sampleDirs)

  #1.3. Fix FileNames if sampleName in dataset was changed
  fileLevelDirs <- list.files(sampleDirs)
  if (length(fileLevelDirs) == 1L) {
    if (fileLevelDirs != sampleName) {
      setwd(sampleDirs)
      ezSystem(paste('mv', fileLevelDirs, sampleName))
      cmd <- paste(
        'rename',
        paste0('s/', fileLevelDirs, '/', sampleName, '/'),
        paste0(sampleName, '/*.gz')
      )
      ezSystem(cmd)
      setwd('..')
    }
  } else if (length(fileLevelDirs) > 1L) {
    if (all(fileLevelDirs %in% sampleName)) {
      #...
    } else {
      stop('multiple runs and renaming samples is an unsupported case')
    }
  }
  dirList <- list(sampleName = sampleName, sampleDirs = sampleDirs)

  #2. Check the dataset for the other modalities and get the fastq files
  libraryTypes <- as.vector(str_split(param$TenXLibrary, ",", simplify = TRUE))
  otherModColNames <- c(
    "MultiDataDir",
    "FeatureDataDir",
    "VdjTDataDir",
    "VdjBDataDir"
  )
  #2.1 VDJ-T
  if ("VDJ-T" %in% libraryTypes) {
    dataInfo <- getCellRangerMultiData(input, "VdjTDataDir", sampleName)
    dirList <- c(
      dirList,
      list(
        vdjtName = dataInfo[["multiName"]],
        vdjtDirs = dataInfo[["multiDirs"]]
      )
    )
  }
  #2.2 VDJ-B
  if ("VDJ-B" %in% libraryTypes) {
    dataInfo <- getCellRangerMultiData(input, "VdjBDataDir", sampleName)
    dirList <- c(
      dirList,
      list(
        vdjbName = dataInfo[["multiName"]],
        vdjbDirs = dataInfo[["multiDirs"]]
      )
    )
  }
  #2.3 Feature Barcoding
  if (
    "FeatureBarcoding" %in%
      libraryTypes ||
      (param$MultiplexingType == "antibody")
  ) {
    dataInfo <- getCellRangerMultiData(input, "FeatureDataDir", sampleName)
    dirList <- c(
      dirList,
      list(
        featureName = dataInfo[["multiName"]],
        featureDirs = dataInfo[["multiDirs"]]
      )
    )
  }
  #2.4 Multiplexing (CMO only - OCM multiplexing info is embedded in GEX reads)
  if (
    "Multiplexing" %in%
      libraryTypes &&
      !("fixedRNA" %in% libraryTypes) &&
      param$MultiplexingType != "antibody" &&
      param$MultiplexingType != "ocm"
  ) {
    dataInfo <- getCellRangerMultiData(input, "MultiDataDir", sampleName)
    dirList <- c(
      dirList,
      list(
        multiplexName = dataInfo[["multiName"]],
        multiplexDirs = dataInfo[["multiDirs"]]
      )
    )
  }

  return(dirList)
}

getSampleMultiplexFiles <- function(input) {
  sampleName <- rownames(input$meta)
  projectId <- strsplit(dirname(input$meta[['RawDataDir']]), '/')[[1]][1]
  sampleMultiplexFolder <- file.path(
    input$dataRoot,
    projectId,
    paste0('o', input$meta[['Order Id']], '_metaData')
  )
  sampleMultiplexFiles <- list.files(sampleMultiplexFolder, full.names = TRUE)
  if (length(sampleMultiplexFiles) == 0) {
    stop(paste0(
      'No multiplexing files found in ',
      sampleMultiplexFolder,
      '. Please add this information to gstore!'
    ))
  }
  return(sampleMultiplexFiles)
}

getCellRangerMultiData <- function(input, multiColName, sampleName) {
  #2.1. Locate the multiplex sample
  multiDirs <- getFastqDirs(input, multiColName, sampleName)
  multiName <- gsub(".tar", "", basename(multiDirs))

  #2.2. Decompress the sample that contains the antibodies reads if they are in tar format
  if (all(grepl("\\.tar$", multiDirs))) {
    multiDirs <- tarExtract(multiDirs)
  }

  multiDirs <- normalizePath(multiDirs)
  return(list(multiName = multiName, multiDirs = multiDirs))
}

buildMultiConfigFile <- function(input, param, dirList) {
  configFileName <- tempfile(
    pattern = "multi_config",
    tmpdir = ".",
    fileext = ".csv"
  )
  fileConn <- file(configFileName)
  fileContents <- c()
  refDir <- NULL

  libraryTypes <- as.vector(str_split(param$TenXLibrary, ",", simplify = TRUE))

  # Add 'GEX' to libraryTypes if 'fixedRNA' is in libraryTypes and 'GEX' is not
  libraryTypes <- if (
    "fixedRNA" %in% libraryTypes && !"GEX" %in% libraryTypes
  ) {
    union(libraryTypes, "GEX")
  } else {
    libraryTypes
  }

  hasFb <- "FeatureBarcoding" %in% libraryTypes
  hasMult <- "Multiplexing" %in% libraryTypes
  isFixed <- "fixedRNA" %in% libraryTypes

  # Write [Gene Expression] section
  if ("GEX" %in% libraryTypes) {
    refDir <- getCellRangerGEXReference(param)
    fileContents <- append(fileContents, "[gene-expression]")
    fileContents <- append(fileContents, sprintf("reference,%s", refDir))
    fileContents <- append(fileContents, paste("create-bam,true"))
    if (ezIsSpecified(param$expectedCells)) {
      fileContents <- append(
        fileContents,
        sprintf("expect-cells,%s", param$expectedCells)
      )
    }
    if (!(isFixed)) {
      # include introns option only available with plain GEX, not fixedRNA
      includeIntronsLine <-
        ifelse(
          ezIsSpecified(param$includeIntrons) && param$includeIntrons,
          "include-introns,true",
          "include-introns,false"
        )
      fileContents <- append(fileContents, includeIntronsLine)
    }
    fileContents <- append(fileContents, c(""))
  }
  if (isFixed) {
    refDir <- getCellRangerGEXReference(param)
    myProbesetFile <- file.path(
      '/srv/GT/databases/10x_Probesets/Chromium',
      param$probesetFile
    )
    outputFile <- sub('.csv', '_filtered.csv', basename(myProbesetFile))
    maxHeaderLine <- max(grep('#', readLines(myProbesetFile)))
    headerSection <- readLines(myProbesetFile, n = maxHeaderLine)
    headerSection[grep('reference_genome', headerSection)] = paste0(
      '#reference_genome=',
      basename(refDir)
    )
    probeInfo <- ezRead.table(
      myProbesetFile,
      sep = ',',
      row.names = NULL,
      skip = maxHeaderLine
    )
    if (ezIsSpecified(param$customProbesFile) && param$customProbesFile != '') {
      customProbes <- ezRead.table(
        file.path(param$dataRoot, param$customProbesFile),
        sep = ',',
        row.names = NULL
      ) %>%
        mutate(
          gene_id = ifelse(
            startsWith(gene_id, "Gene_"),
            gene_id,
            paste0("Gene_", gene_id)
          ),
          probe_id = ifelse(
            startsWith(probe_id, "Gene_"),
            probe_id,
            paste0("Gene_", probe_id)
          )
        )
      probeInfo <- bind_rows(list(probeInfo, customProbes))
    }
    annotation <- ezRead.table(
      file.path(refDir, 'star', 'geneInfo.tab'),
      row.names = NULL,
      skip = 1,
      header = FALSE
    )
    intersectionGeneIDs <- intersect(annotation$V1, probeInfo$gene_id)
    probeInfo <- probeInfo[probeInfo$gene_id %in% intersectionGeneIDs, ]
    if (length(probeInfo$gene_name) == 0) {
      probeInfo_gene_name <- limma::strsplit2(probeInfo$probe_id, '\\|')[, 2]
    } else {
      probeInfo_gene_name <- probeInfo$gene_name
    }
    intersectionGeneNames <- intersect(annotation$V2, probeInfo_gene_name)
    probeInfo <- probeInfo[probeInfo_gene_name %in% intersectionGeneNames, ]
    probeInfo_gene_name <- probeInfo_gene_name[
      probeInfo_gene_name %in% intersectionGeneNames
    ]
    myID_probes <- paste0(probeInfo_gene_name, '--', probeInfo$gene_id)
    annotationIDs <- paste0(annotation$V2, '--', annotation$V1)

    probeInfo <- probeInfo[which(myID_probes %in% annotationIDs), ]

    writeLines(headerSection, outputFile)
    ezWrite.table(
      probeInfo,
      outputFile,
      sep = ',',
      row.names = FALSE,
      append = TRUE
    )
    fileContents <- append(
      fileContents,
      sprintf("probe-set,%s", file.path(getwd(), outputFile))
    )

    # Determine Chemistry for Flex v1/v2
    # Priority 1: User specified chemistry in parameters (not "auto")
    # Priority 2: Auto-detect based on probe set version (v2.x.x vs v1.x.x)
    # Cell Ranger 10.0+ requires explicit chemistry for reliable Flex v2 detection
    if (ezIsSpecified(param$chemistry) && param$chemistry != "auto") {
      chemistry <- param$chemistry
    } else if (grepl("v2\\.", param$probesetFile)) {
      # Flex v2: Matches v2.0.0, v2.1.0, etc.
      chemistry <- ifelse(hasMult, "Flex-v2-RNA-R2", "Flex-v2-singleplex")
    } else {
      # Flex v1: Default legacy behavior (MFRP/SFRP)
      chemistry <- ifelse(hasMult, "MFRP", "SFRP")
    }
    fileContents <- append(fileContents, sprintf("chemistry,%s", chemistry))
    fileContents <- append(fileContents, c(""))
  }
  if (any(c("VDJ-T", "VDJ-B") %in% libraryTypes)) {
    vdjRefDir <- getCellRangerVDJReference(param)
    fileContents <- append(fileContents, "[vdj]")
    fileContents <- append(fileContents, sprintf("reference,%s", vdjRefDir))
    fileContents <- append(fileContents, c(""))
  }
  if (hasFb) {
    featureRefFile <- file.path(param$dataRoot, param$FeatureBarcodeFile)
    fileContents <- append(fileContents, "[feature]")
    fileContents <- append(
      fileContents,
      sprintf("reference,%s", featureRefFile)
    )
    fileContents <- append(fileContents, c(""))
  }
  if (hasMult) {
    multiplexBarcodeFile <- NULL # Initialize as NULL; only created for CMO/HTO
    if (param$MultiplexingType == "antibody") {
      # HTO/Hashtag multiplexing - needs [feature] section
      multiplexBarcodeFile <- tempfile(
        pattern = "multi_barcode_set",
        tmpdir = ".",
        fileext = ".csv"
      )
      multiplexBarcodeFile <- file.path(getwd(), multiplexBarcodeFile)
      fileContents <- append(fileContents, "[feature]")
      fileContents <- append(
        fileContents,
        sprintf("reference,%s", multiplexBarcodeFile)
      )
      fileContents <- append(fileContents, c(""))
    } else if (!isFixed && param$MultiplexingType != "ocm") {
      # CMO/CellPlex multiplexing - needs cmo-set reference
      multiplexBarcodeFile <- tempfile(
        pattern = "multi_barcode_set",
        tmpdir = ".",
        fileext = ".csv"
      )
      multiplexBarcodeFile <- file.path(getwd(), multiplexBarcodeFile)
      fileContents <- append(
        fileContents,
        sprintf("cmo-set,%s", multiplexBarcodeFile)
      )
      fileContents <- append(fileContents, c(""))
    }
    # OCM: No barcode reference file needed - chemistry handles it internally
  }

  # Fastq Files
  fileContents <- append(
    fileContents,
    c("[libraries]", "fastq_id,fastqs,feature_types")
  )
  if (any(c("GEX", "fixedRNA") %in% libraryTypes)) {
    fileContents <- append(
      fileContents,
      sprintf(
        "%s,%s,%s",
        dirList$sampleName,
        dirList$sampleDirs,
        "Gene Expression"
      )
    )
  }
  if ("VDJ-T" %in% libraryTypes) {
    fileContents <- append(
      fileContents,
      sprintf("%s,%s,%s", dirList$vdjtName, dirList$vdjtDirs, "VDJ-T")
    )
  }
  if ("VDJ-B" %in% libraryTypes) {
    fileContents <- append(
      fileContents,
      sprintf("%s,%s,%s", dirList$vdjbName, dirList$vdjbDirs, "VDJ-B")
    )
  }
  if (hasFb) {
    fileContents <- append(
      fileContents,
      sprintf(
        "%s,%s,%s",
        dirList$featureName,
        dirList$featureDirs,
        "Antibody Capture"
      )
    )
  }
  # Multiplexing library entry - skip for OCM (multiplexing info embedded in GEX reads)
  if (hasMult && !isFixed && param$MultiplexingType != "ocm") {
    if (param$MultiplexingType == "antibody") {
      fileContents <- append(
        fileContents,
        sprintf(
          "%s,%s,%s",
          dirList$featureName,
          dirList$featureDirs,
          "Antibody Capture"
        )
      )
    } else {
      fileContents <- append(
        fileContents,
        sprintf(
          "%s,%s,%s",
          dirList$multiplexName,
          dirList$multiplexDirs,
          "Multiplexing Capture"
        )
      )
    }
  }
  fileContents <- append(fileContents, "")

  # sample mapping
  if (hasMult) {
    sampleName <- rownames(input$meta)
    sampleMultiplexFiles <- getSampleMultiplexFiles(input)
    # we match according to just the beginning ^ and the sample names as a
    # prefix, since sometimes parts of the library information are
    # as postfixes. This may potentially cause collisions but we risk it
    names(sampleMultiplexFiles) <- paste0(
      '^',
      sub('_Sample2Barcode.csv', '', basename(sampleMultiplexFiles))
    )
    sampleMultiplexFile <- sampleMultiplexFiles[sapply(
      names(sampleMultiplexFiles),
      grepl,
      pattern = sampleName,
      ignore.case = TRUE
    )]
    sampleMultiplexMapping <- read_csv(
      sampleMultiplexFile,
      show_col_types = FALSE
    )
    concatCols <- function(y) {
      return(paste(as.character(y), collapse = ","))
    }
    if (!("fixedRNA" %in% libraryTypes)) {
      if (param$MultiplexingType == "ocm") {
        # OCM: No barcode reference file needed - use ocm_barcode_ids column directly
        fileContents <- append(
          fileContents,
          c("[samples]", "sample_id,ocm_barcode_ids")
        )
        fileContents <- append(
          fileContents,
          apply(sampleMultiplexMapping, 1, concatCols)
        )
      } else {
        # CMO or HTO: Load and filter barcode reference file
        multiplexBarcodeSet <- read_csv(
          file.path(
            "/srv/GT/databases/10x/CMO_files",
            param$MultiplexBarcodeSet
          ),
          show_col_types = FALSE
        )
        # Handle pipe-separated IDs for double-hashing (e.g., "B0301|B0304")
        all_barcode_ids <- unique(unlist(strsplit(
          sampleMultiplexMapping$cmo_ids,
          "\\|"
        )))
        multiplexBarcodeSet <- multiplexBarcodeSet %>%
          filter(id %in% all_barcode_ids)
        data.table::fwrite(
          multiplexBarcodeSet,
          file = multiplexBarcodeFile,
          sep = ","
        )

        if (param$MultiplexingType == "antibody") {
          fileContents <- append(
            fileContents,
            c("[samples]", "sample_id,hashtag_ids")
          )
        } else {
          fileContents <- append(
            fileContents,
            c("[samples]", "sample_id,cmo_ids")
          )
        }
        fileContents <- append(
          fileContents,
          apply(sampleMultiplexMapping, 1, concatCols)
        )
      }
    } else if ("fixedRNA" %in% libraryTypes) {
      sampleMultiplexMapping$description = sampleMultiplexMapping$sample_id
      fileContents <- append(
        fileContents,
        c("[samples]", "sample_id,probe_barcode_ids,description")
      )
      fileContents <- append(
        fileContents,
        apply(sampleMultiplexMapping, 1, concatCols)
      )
    }
  }

  # write outputs
  writeLines(fileContents, fileConn)
  close(fileConn)

  return(list(configFileName = configFileName, refDir = refDir))
}

##' @author Noé, Falko
##' @template app-template
##' @templateVar method ezMethodCellRanger(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppCellRangerMulti <-
  setRefClass(
    "EzAppCellRangerMulti",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCellRangerMulti
        name <<- "EzAppCellRangerMulti"
        appDefaults <<- rbind(
          TenXLibrary = ezFrame(
            Type = "charVector",
            DefaultValue = "GEX",
            Description = "Which 10X library?"
          ),
          chemistry = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "Assay configuration."
          ),
          expectedCells = ezFrame(
            Type = "numeric",
            DefaultValue = 10000,
            Description = "Expected number of cells."
          ),
          includeIntrons = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Count reads on introns."
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          ),
          keepBam = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "keep bam file produced by CellRanger"
          )
        )
      }
    )
  )
