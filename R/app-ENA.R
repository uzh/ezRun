## Download an ENA XML record (run/sample/experiment) with retries and validation.
## ENA occasionally returns a truncated or empty body that still parses as XML but
## lacks the expected content; retrying on a too-small/unparseable response avoids
## poisoning the downstream extraction. Returns a parsed XMLInternalDocument.
.enaFetchXml <- function(accession, tries = 4, minBytes = 400) {
  destfile <- paste0(accession, ".xml")
  for (attempt in seq_len(tries)) {
    cmd <- paste0(
      "curl -s -f -o ",
      destfile,
      " -X GET \'https://www.ebi.ac.uk/ena/browser/api/xml/",
      accession,
      "?download=true\'"
    )
    try(ezSystem(cmd), silent = TRUE)
    if (file.exists(destfile) && file.info(destfile)$size >= minBytes) {
      doc <- tryCatch(xmlParse(destfile), error = function(e) NULL)
      if (!is.null(doc)) {
        return(doc)
      }
    }
    Sys.sleep(5 * attempt)
  }
  stop(
    "Could not retrieve a valid ENA XML for ",
    accession,
    " after ",
    tries,
    " attempts"
  )
}

## Download one ENA fastq with retries. wget treats an FTP "No such directory" reply as
## fatal and does not retry it, although ENA returns it transiently for existing files.
## A file already present with the expected md5 is reused, so a crashed job can be
## resumed in its old scratch directory.
.enaDownloadFastq <- function(ftpPath, md5 = NA, tries = 5, retryWait = 60) {
  destfile <- basename(ftpPath)
  isComplete <- function() {
    file.exists(destfile) &&
      (is.na(md5) ||
        md5 == "" ||
        unname(tools::md5sum(destfile)) == md5)
  }
  if (isComplete()) {
    ezLog("Reusing existing download: ", destfile)
    return(destfile)
  }
  for (attempt in seq_len(tries)) {
    unlink(destfile)
    try(
      ezSystem(paste0("wget --no-verbose --tries=3 --wait=10 ftp://", ftpPath)),
      silent = TRUE
    )
    if (isComplete()) {
      return(destfile)
    }
    Sys.sleep(retryWait * attempt)
  }
  stop("Could not download ", ftpPath, " after ", tries, " attempts")
}

## Python env with bfabricPy; the job user's ~/.bfabricpy.yml and
## BFABRICPY_CONFIG_ENV select the B-Fabric instance and account.
ENA_BFABRIC_PYTHON <- "/usr/local/ngseq/miniforge3/envs/gi_py3.12.8/bin/python"

## Register the downloaded files in B-Fabric as short-term storage (application 372),
## so they are deleted automatically after a while. Runs before the files are copied
## to gStore: sizes and md5 sums come from the copies in the working directory.
## A failure is logged but does not fail the app; the data are still delivered.
.enaRegisterShortTermStorage <- function(datasetFile, containerId, name) {
  script <- system.file("python/register_bfabric_storage.py", package = "ezRun")
  cmd <- paste(
    ENA_BFABRIC_PYTHON, "-I", shQuote(script),
    "--dataset", shQuote(datasetFile),
    "--container-id", containerId,
    "--name", shQuote(name),
    "--application 372 --local-dir .",
    "--threads", min(8, ezThreads())
  )
  tryCatch(
    {
      if (!nzchar(script) || !file.exists(ENA_BFABRIC_PYTHON)) {
        stop("missing ", if (nzchar(script)) ENA_BFABRIC_PYTHON else "register_bfabric_storage.py")
      }
      ezSystem(cmd)
      TRUE
    },
    error = function(e) {
      ezLog(
        "B-Fabric short-term storage registration FAILED: ", conditionMessage(e),
        level = "error"
      )
      FALSE
    }
  )
}

## Extract the ENA-SPOT-COUNT (read count) from a run's RUN_ATTRIBUTES list without
## crashing on a missing/empty TAG. Returns NA when the tag is absent, in which case
## the caller falls back to countReadsInFastq().
.enaSpotCount <- function(runAttr) {
  if (is.null(runAttr) || !is.list(runAttr)) {
    return(NA)
  }
  for (k in seq_along(runAttr)) {
    tag <- runAttr[[k]]$TAG
    if (
      !is.null(tag) &&
        length(tag) == 1 &&
        !is.na(tag) &&
        tag == "ENA-SPOT-COUNT"
    ) {
      return(runAttr[[k]]$VALUE)
    }
  }
  NA
}

## ENA lists paired runs with 3 files when some reads lost their mate:
## <run>.fastq.gz (unpaired) plus <run>_1 / <run>_2. Keep the two mates (and
## their md5 sums) so that a run has 1 (single-end) or 2 (paired) files.
## Returns list(ftp, md5) as ';'-separated strings.
.enaNormaliseFastqFiles <- function(ftp, md5) {
  files <- strsplit(ftp, ';', fixed = TRUE)[[1]]
  md5s <- strsplit(ifelse(is.na(md5), '', md5), ';', fixed = TRUE)[[1]]
  if (length(files) <= 2) {
    return(list(ftp = ftp, md5 = md5))
  }
  idx1 <- grep('_1\\.f(ast)?q(\\.gz)?$', files)
  idx2 <- grep('_2\\.f(ast)?q(\\.gz)?$', files)
  if (length(idx1) != 1 || length(idx2) != 1) {
    stop(
      "Cannot identify the read 1 / read 2 files among ",
      length(files), " files: ", ftp
    )
  }
  keep <- c(idx1, idx2)
  ezLog(
    "Ignoring unpaired file(s) of a paired run: ",
    paste(basename(files[-keep]), collapse = ", ")
  )
  list(
    ftp = paste(files[keep], collapse = ';'),
    md5 = if (length(md5s) == length(files)) {
      paste(md5s[keep], collapse = ';')
    } else {
      md5
    }
  )
}

## Shell command that registers the delivered dataset in production SUSHI as a
## child of the ENA order dataset. BatchMode makes ssh fail instead of waiting
## for a password when the key is not accepted.
.enaSushiRegisterCmd <- function(datasetPath, containerId, datasetName, datasetId) {
  sprintf(
    paste0(
      "ssh -o BatchMode=yes trxcopy@fgcz-h-082 ",
      "\"bash -lc 'cd /srv/sushi/production/master && ",
      "bundle exec sushi_fabric --dataset %s --run --input_dataset_application 372 ",
      "--project %s --dataset_name %s --dataset_id %s'\""
    ),
    datasetPath,
    containerId,
    datasetName,
    datasetId
  )
}

ezMethodGetEnaData <- function(input = NA, output = NA, param = NA) {
  require(XML)
  for (i in seq_along(5)) {
    containerId <- sub(
      'p',
      '',
      strsplit(output$getColumn('ENA Result'), '/')[[1]][1]
    )
  }
  cmd = paste0(
    "curl -o fastqLinks.txt -X GET ",
    "\'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
    param[['projectID']],
    "&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,fastq_md5&format=tsv&download=true\'"
  )
  ezSystem(cmd)
  fastqInfo = ezRead.table('fastqLinks.txt', row.names = NULL)
  fastqInfo[['Name']] = ''
  fastqInfo[['ReadCount']] = 0
  colnames(fastqInfo)[grep('scientific_name', colnames(fastqInfo))] = 'Species'

  ##Filter for included or excluded samples:
  if (param$excludedSamples != '') {
    excludedSamples <- as.vector(gsub(
      ' ',
      '',
      limma::strsplit2(param$excludedSamples, split = ',')
    ))
    toExclude <- which(fastqInfo[['sample_accession']] %in% excludedSamples)
    if (length(toExclude) > 0) {
      ezLog(paste(
        'Samples excluded:',
        paste(fastqInfo[toExclude, 'sample_accession'], collapse = ", ")
      ))
      fastqInfo <- fastqInfo[-c(toExclude), ]
    }
  }

  if (param$includedSamples != '') {
    includedSamples <- as.vector(gsub(
      ' ',
      '',
      limma::strsplit2(param$includedSamples, split = ',')
    ))
    includedSamples <- includedSamples[includedSamples != '']
    toInclude <- which(fastqInfo[['sample_accession']] %in% includedSamples)
    if (length(toInclude) == 0) {
      stop(
        "includedSamples matches no sample_accession of project ",
        param[['projectID']], ": ", paste(includedSamples, collapse = ", ")
      )
    }
    notFound <- setdiff(includedSamples, fastqInfo[['sample_accession']])
    if (length(notFound) > 0) {
      ezLog(
        "includedSamples not found in the project and ignored: ",
        paste(notFound, collapse = ", "),
        level = "warn"
      )
    }
    ezLog(paste(
      'Samples included:',
      paste(unique(fastqInfo[toInclude, 'sample_accession']), collapse = ", ")
    ))
    fastqInfo <- fastqInfo[toInclude, ]
  }

  ###Filter fastqInfo for empty fastq_ftp entries
  toRemove <- which(fastqInfo[['fastq_ftp']] == '')

  if (length(toRemove) > 0) {
    ezLog(paste(
      'Samples removed because of missing raw data:',
      paste(capture.output(print(fastqInfo[toRemove, ])), collapse = "\n")
    ))
    fastqInfo <- fastqInfo[-c(toRemove), ]
  }

  ## Paired runs can list an extra file with the unpaired reads; keep the mates only
  for (i in seq_len(nrow(fastqInfo))) {
    files <- .enaNormaliseFastqFiles(
      fastqInfo$fastq_ftp[i],
      fastqInfo$fastq_md5[i]
    )
    fastqInfo$fastq_ftp[i] <- files$ftp
    fastqInfo$fastq_md5[i] <- files$md5
  }

  ##Check sequencing mode and keep only supported mode in mixed cases
  pairedDataPos <- grep(';', fastqInfo$fastq_ftp)
  if (length(pairedDataPos) > 0 & length(pairedDataPos) != nrow(fastqInfo)) {
    if (param$supportedMode == 'paired') {
      fastqInfo <- fastqInfo[pairedDataPos, ]
    } else if (param$supportedMode == 'single') {
      fastqInfo <- fastqInfo[-c(pairedDataPos), ]
    }
  }

  setwdNew(rownames(input$meta))

  fastqInfo$Name <- fastqInfo$sample_accession ## use the accession as a fallback
  for (i in 1:nrow(fastqInfo)) {
    #download ERR xml File (retried/validated against truncated ENA responses)
    runInfo <- .enaFetchXml(fastqInfo$run_accession[i])
    #Extract read number from xml, extract sampleID from xml
    runAttr = xmlToList(runInfo)$RUN$RUN_ATTRIBUTES
    readCount <- .enaSpotCount(runAttr)
    if (!is.na(readCount)) {
      fastqInfo[['ReadCount']][i] <- readCount
    }
    sampleID <- xmlToList(runInfo)$RUN$RUN_LINKS[[2]]$XREF_LINK$ID

    if (is.list(sampleID) | any(grepl('^ERP', sampleID))) {
      sampleID <- fastqInfo$sample_accession[i]
    }
    xml <- .enaFetchXml(sampleID)
    sampleInfo <- xmlToList(xml)

    if (!is.null(sampleInfo$SAMPLE$TITLE)) {
      cleanName <- gsub("[\\# \\(\\):,;]", "_", sampleInfo$SAMPLE$TITLE)
      fastqInfo[['Name']][i] <- paste(cleanName, sampleID, sep = '_')
      sampleAttributes <- xmlToDataFrame(xmlRoot(xml)[[1]][[5]])
      sampleAttributes <- sampleAttributes[
        grep('ENA', sampleAttributes$TAG, invert = TRUE),
      ]

      if (!is.null(nrow(sampleAttributes)) && nrow(sampleAttributes) > 0) {
        for (j in 1:nrow(sampleAttributes)) {
          attrName <- as.character(sampleAttributes[j, 1])
          if (is.null(fastqInfo[[attrName]])) {
            fastqInfo[[attrName]] = '-'
          }
          fastqInfo[i, attrName] = as.character(sampleAttributes[j, 2])
        }
      }
    } else {
      out <- tryCatch(
        {
          experimentID <- xmlToList(runInfo)[[1]]$EXPERIMENT_REF$.attrs
          if (length(experimentID) > 1) {
            experimentID <- experimentID['accession']
          }
          xml <- .enaFetchXml(experimentID)
          experimentInfo <- xmlToList(xml)
          cleanName <- gsub(
            "[\\# \\(\\):,;]",
            "_",
            experimentInfo$EXPERIMENT$TITLE
          )
          fastqInfo[['Name']][i] <- paste(cleanName, sampleID, sep = '_')
        },
        error = function(e) {
          ezLog(
            'sample title not available in experiment accession file',
            level = "error"
          )
          return("")
        }
      )
    }

    if (grepl(';', fastqInfo$fastq_ftp[i])) {
      paired = TRUE
      fastqFiles <- strsplit(fastqInfo$fastq_ftp[i], ';')[[1]]
      fastqMd5s <- strsplit(fastqInfo$fastq_md5[i], ';')[[1]]
      stopifnot(length(fastqFiles) == 2)
      for (k in seq_along(fastqFiles)) {
        .enaDownloadFastq(fastqFiles[k], fastqMd5s[k])
      }
    } else {
      paired = FALSE
      .enaDownloadFastq(fastqInfo$fastq_ftp[i], fastqInfo$fastq_md5[i])
    }
  }
  myPath = output$meta[['ENA Result [File]']]
  dataset <- createDataset(fastqInfo, myPath, paired = paired)

  if (sum(dataset[['Read Count']]) == 0) {
    dataset[['Read Count']] <- countReadsInFastq(basename(fastqInfo$fastq_ftp))
  }

  if (length(dataset$Name) == length(unique(dataset$Name))) {
    ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
  } else {
    ezLog('Data pooling needed')
    readColumns <- intersect(c('Read1 [File]', 'Read2 [File]'), colnames(dataset))
    for (sampleName in unique(dataset$Name)) {
      ## select by the original name; the sanitised stem is only for file names
      isSample <- dataset[['Name']] == sampleName
      fileStem <- gsub('[/ ;]', '_', sampleName)
      for (readColumn in readColumns) {
        readTag <- if (readColumn == 'Read1 [File]') '_R1' else '_R2'
        pooledFile <- paste0(fileStem, readTag, '.fastq.gz')
        sampleFiles <- basename(dataset[[readColumn]][isSample])
        ezSystem(paste(
          'cat', paste(shQuote(sampleFiles), collapse = ' '),
          '>', shQuote(pooledFile)
        ))
        ezSystem(paste('rm', paste(shQuote(sampleFiles), collapse = ' ')))
        dataset[[readColumn]][isSample] <- file.path(
          dirname(dataset[[readColumn]][isSample]),
          pooledFile
        )
      }
      dataset[['Read Count']][isSample] <- sum(as.numeric(
        dataset[['Read Count']][isSample]
      ))
    }
    dataset = dataset[, -grep('md5sum', colnames(dataset))]
    dataset = unique(dataset)
    ##TODO: Recompute md5sums after pooling & handle PairedEnd data
    ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
  }

  if (param$tarOutput) {
    readColumns <- intersect(c('Read1 [File]', 'Read2 [File]'), colnames(dataset))
    ## folder and tar names must not contain '/' (sample titles may)
    tarStems <- gsub('[/ ;]', '_', dataset$Name)
    for (i in seq_len(nrow(dataset))) {
      mySample <- tarStems[i]
      dir.create(mySample)
      ## move each read file into the sample folder with an Illumina-style name
      for (readColumn in readColumns) {
        readTag <- if (readColumn == 'Read1 [File]') 'R1' else 'R2'
        ezSystem(paste(
          'mv',
          shQuote(basename(dataset[[readColumn]][i])),
          shQuote(file.path(
            mySample,
            paste0(mySample, '_S1_L001_', readTag, '_001.fastq.gz')
          ))
        ))
      }
      ezSystem(paste('tar cf', shQuote(paste0(mySample, '.tar')), shQuote(mySample)))
      ezSystem(paste('rm -Rf', shQuote(mySample)))
    }

    ##Update dataset
    dataset[['RawDataDir [File]']] <- file.path(
      dirname(dataset[['Read1 [File]']])[1],
      paste0(tarStems, '.tar')
    )
    dataset[['Read1 [File]']] <- NULL
    dataset[['Read2 [File]']] <- NULL
    ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
  }

  ##Register dataset/resources in BF as short term storage
  datasetName = paste0('ENA_App_', output$getColumn('projectID'))

  .enaRegisterShortTermStorage('dataset.tsv', containerId, datasetName)
  datasetPath <- paste0(
    '/srv/GT/analysis/datasets/ENA_App/dataset_',
    containerId,
    '.tsv'
  )
  ezWrite.table(dataset, datasetPath, row.names = FALSE)

  ## Register the dataset in SUSHI. A failure fails the job: the data are on
  ## gStore already, but without the SUSHI dataset the order looks delivered
  ## while it is not.
  registerCmd <- .enaSushiRegisterCmd(
    datasetPath,
    containerId,
    datasetName,
    param$datasetId
  )
  scriptPath <- tempfile(fileext = ".sh")
  writeLines(c("#!/bin/bash", "set -e", registerCmd), con = scriptPath)
  res <- system(paste('bash', shQuote(scriptPath)))
  if (res != 0) {
    stop(
      "SUSHI registration failed (exit code ", res, "). The FASTQ files are ",
      "delivered; register the dataset by running:\n", registerCmd
    )
  }
  ezLog("Registered ", datasetName, " in SUSHI (project ", containerId, ")")
  return("Success")
}

createDataset <- function(fastqInfo, myPath, paired = FALSE) {
  if (ncol(fastqInfo) > 10) {
    xtraCols <- fastqInfo[, 11:ncol(fastqInfo)]
    colnames(fastqInfo[, 11:ncol(fastqInfo)])
  } else {
    xtraCols <- data.frame(sample_accession = fastqInfo$sample_accession)
  }
  if (!paired) {
    dataset = data.frame(
      Name = gsub(' ', '_', fastqInfo$Name),
      Read1 = file.path(myPath, basename(fastqInfo$fastq_ftp)),
      md5sum = fastqInfo$fastq_md5,
      Species = fastqInfo$Species,
      ReadCount = fastqInfo$ReadCount,
      xtraCols,
      stringsAsFactors = FALSE
    )
    colnames(dataset) = c(
      'Name',
      'Read1 [File]',
      'md5sum',
      'Species',
      'Read Count',
      colnames(xtraCols)
    )
  } else {
    dataset = data.frame(
      Name = gsub(' ', '_', fastqInfo$Name),
      Read1 = file.path(
        myPath,
        sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[1, ]
      ),
      Read2 = file.path(
        myPath,
        sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[2, ]
      ),
      md5sum = fastqInfo$fastq_md5,
      Species = fastqInfo$Species,
      ReadCount = fastqInfo$ReadCount,
      xtraCols,
      stringsAsFactors = FALSE
    )
    colnames(dataset) = c(
      'Name',
      'Read1 [File]',
      'Read2 [File]',
      'md5sum',
      'Species',
      'Read Count',
      colnames(xtraCols)
    )
  }
  dataset[['Read Count']] <- as.numeric(dataset[['Read Count']])
  return(dataset)
}

##' @template app-template
##' @templateVar method ezMethodGetEnaData(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppENA <-
  setRefClass(
    "EzAppENA",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodGetEnaData
        name <<- "EzAppENA"
        appDefaults <<- rbind(
          excludedSamples = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = 'sample_accession SAMN-ids for samples to exclude'
          ),
          includedSamples = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = 'sample_accession SAMN-ids for samples to include'
          )
        )
      }
    )
  )
