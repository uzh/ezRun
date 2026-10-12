###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## BRB-seq (Alithea MERCURIUS V5): pooled 3' RNA-seq, 96 wells per plate.
## R1 = 14 nt well barcode + 14 nt UMI + polyT (151 nt), R2 = cDNA.
## One input row = one sequencing pool (already split by Illumina index).
## The well barcode is resolved here by STARsolo, then every library of the
## sample sheet is written as a FeatureCounts-format count file so that
## CountQC / DESeq2 / EdgeR run on the output unchanged.

BRB_KIT_FILE <- "/srv/GT/databases/BRB/whiteList_SampleBarcodes/BRB_V5_barcodes.tsv"
BRB_STAT_ROWS <- c(
  "Assigned", "Unassigned_Unmapped", "Unassigned_Read_Type", "Unassigned_Singleton",
  "Unassigned_MappingQuality", "Unassigned_Chimera", "Unassigned_FragmentLength",
  "Unassigned_Duplicate", "Unassigned_MultiMapping", "Unassigned_Secondary",
  "Unassigned_NonSplit", "Unassigned_NoFeatures", "Unassigned_Overlapping_Length",
  "Unassigned_Ambiguity"
)

ezMethodBRBseq <- function(input = NA, output = NA, param = NA) {
  require(Matrix)
  poolNames <- input$getNames()
  resultDir <- param$name
  dir.create(resultDir, showWarnings = FALSE)

  ## 1. Sample sheet + barcode kit
  sheet <- readBRBSampleSheet(getBRBSampleSheetPath(input, param))
  kit <- ezRead.table(param$barcodeFile, row.names = NULL)
  plateVersion <- checkBRBSampleSheet(sheet, kit, param$plateVersion, poolNames)
  whitelistFile <- file.path(getwd(), "brb_whitelist.txt")
  writeLines(kit$Barcode[kit$Plate == plateVersion], whitelistFile)

  ## 2. STARsolo per pool (or adopt finished runs from param$starsoloDir)
  feature <- if (isTRUE(as.logical(param$countIntrons))) "GeneFull" else "Gene"
  poolDirs <- if (ezIsSpecified(param$starsoloDir)) {
    adoptBRBStarsolo(param$starsoloDir, poolNames, feature)
  } else {
    runBRBStarsolo(input, param, getSTARSoloReference(param), whitelistFile, poolNames, feature)
  }

  ## 3. Libraries -> pools, per-library counts
  pools <- lapply(poolDirs, readBRBPool, feature = feature)
  asg <- assignBRBLibraries(sheet, pools)
  genes <- brbGenesToWrite(rownames(pools[[1]]$m), param)
  libs <- writeBRBCounts(asg$libs, pools, genes, resultDir)
  writeBRBDataset(libs, input, output, param, resultDir)

  ## 4. Keep the per-pool STARsolo summaries + raw matrices, then the report
  poolSummary <- keepBRBPoolOutputs(poolDirs, resultDir, feature, param)
  ezWrite.table(sheet, file.path(resultDir, "sampleSheet_used.tsv"), row.names = FALSE)
  ezWrite.table(libs, file.path(resultDir, "library_assignment.tsv"), row.names = FALSE)
  ezWrite.table(asg$umiByPool, file.path(resultDir, "barcode_umis_by_pool.tsv"), row.names = FALSE)
  makeBRBseqReport(resultDir, param,
    libs = libs, umiByPool = asg$umiByPool, poolSummary = poolSummary,
    plateVersion = plateVersion,
    ySignal = brbYGeneSignal(pools, libs, param)
  )
  return("Success")
}

##' Sample sheet location: an uploaded sheet (param$sampleSheet) wins; otherwise
##' the B-Fabric-derived sheet at <project>/o<order>_metaData/o<order>_BRB_SampleSheet.tsv
##' written by inst/python/brbseq_samplesheet.py.
getBRBSampleSheetPath <- function(input, param) {
  if (ezIsSpecified(param$sampleSheet)) {
    path <- param$sampleSheet
    if (!file.exists(path)) path <- file.path(param$dataRoot, path)
    if (!file.exists(path)) stop("sample sheet not found: ", param$sampleSheet)
    return(path)
  }
  read1 <- strsplit(as.character(input$getColumn("Read1"))[1], ",")[[1]][1]
  projectId <- strsplit(read1, "/")[[1]][1]
  if (!input$hasColumn("Order Id")) stop("no 'Order Id' column: upload a sample sheet")
  orderId <- as.character(input$getColumn("Order Id"))[1]
  path <- file.path(input$dataRoot, projectId, paste0("o", orderId, "_metaData"),
                    paste0("o", orderId, "_BRB_SampleSheet.tsv"))
  if (!file.exists(path)) {
    stop("No BRB sample sheet: upload one or create ", path,
         " with inst/python/brbseq_samplesheet.py --order ", orderId, " --deliver")
  }
  path
}

readBRBSampleSheet <- function(path) {
  sheet <- read.delim(path, check.names = FALSE, colClasses = "character",
                      na.strings = character(0), strip.white = TRUE)
  need <- c("Name", "Barcode", "Well")
  if (!all(need %in% colnames(sheet))) {
    stop("sample sheet ", path, " needs the columns ", paste(need, collapse = ", "))
  }
  if (!"Role" %in% colnames(sheet)) sheet$Role <- "sample"
  if (!"Pool" %in% colnames(sheet)) sheet$Pool <- ""
  sheet$Role[sheet$Role == ""] <- "sample"
  sheet$Barcode <- toupper(sheet$Barcode)
  sheet$Well <- sprintf("%s%02d", toupper(substr(sheet$Well, 1, 1)),
                        as.integer(substring(sheet$Well, 2)))
  sheet[sheet$Name != "", , drop = FALSE]
}

##' Validate the sheet and resolve the plate version ("auto" = best match of the
##' sheet's well/barcode pairs against the kit; V5A and V5D share sequences).
checkBRBSampleSheet <- function(sheet, kit, plateVersion, poolNames) {
  dupName <- unique(sheet$Name[duplicated(sheet$Name)])
  if (length(dupName) > 0) stop("duplicate sample sheet names: ", paste(dupName, collapse = ", "))
  bad <- sheet$Name[!grepl("^[A-Za-z0-9_-]+$", sheet$Name)]
  if (length(bad) > 0) stop("sample sheet names may only use letters, digits, _ and -: ", paste(bad, collapse = ", "))
  poolBc <- paste(sheet$Pool, sheet$Barcode)
  if (any(duplicated(poolBc))) {
    stop("barcode used twice in the same pool: ", paste(unique(sheet$Barcode[duplicated(poolBc)]), collapse = ", "),
         " (several plates of one version must be split with the Pool column)")
  }
  unknownPool <- setdiff(unique(sheet$Pool[sheet$Pool != ""]), poolNames)
  if (length(unknownPool) > 0) stop("sample sheet Pool not in the input dataset: ", paste(unknownPool, collapse = ", "))

  if (is.null(plateVersion) || plateVersion %in% c("", "auto")) {
    key <- paste(sheet$Well, sheet$Barcode)
    score <- sapply(split(paste(kit$Well, kit$Barcode), kit$Plate), function(k) sum(key %in% k))
    plateVersion <- names(which.max(score))
    ezLog(sprintf("plate version auto-detected: %s (%d/%d wells match)", plateVersion, max(score), nrow(sheet)))
  }
  wl <- kit$Barcode[kit$Plate == plateVersion]
  if (length(wl) != 96) stop("barcode kit has no plate ", plateVersion)
  offList <- sheet$Name[!sheet$Barcode %in% wl]
  if (length(offList) > 0) stop("barcodes not on plate ", plateVersion, ": ", paste(offList, collapse = ", "))
  wellBc <- setNames(kit$Barcode[kit$Plate == plateVersion], kit$Well[kit$Plate == plateVersion])
  mism <- sheet$Name[wellBc[sheet$Well] != sheet$Barcode]
  if (length(mism) > 0) {
    ezLog(paste("WARNING: well label does not match the barcode on plate", plateVersion, "for:",
                paste(mism, collapse = ", ")))
  }
  plateVersion
}

##' STARsolo for each pool into pools/<pool>/. Settings follow the Alithea
##' BRB-seq layout and are not exposed: changing them breaks the demultiplexing.
runBRBStarsolo <- function(input, param, refDir, whitelistFile, poolNames, feature) {
  read1 <- input$getFullPathsList("Read1")
  read2 <- input$getFullPathsList("Read2")
  keepBam <- isTRUE(as.logical(param$keepBam))
  poolDirs <- setNames(file.path(getwd(), "pools", poolNames), poolNames)
  for (nm in poolNames) {
    if (length(read1[[nm]]) != length(read2[[nm]])) {
      stop(sprintf("unequal number of Read1/Read2 files for pool %s", nm))
    }
    dir.create(poolDirs[[nm]], recursive = TRUE, showWarnings = FALSE)
    cmd <- paste(
      "STAR --runMode alignReads",
      "--runThreadN", param$cores,
      "--genomeDir", refDir,
      ## cDNA (R2) first, barcode read (R1) last
      "--readFilesIn", paste(read2[[nm]], collapse = ","), paste(read1[[nm]], collapse = ","),
      "--readFilesCommand zcat",
      "--soloType CB_UMI_Simple",
      "--soloCBwhitelist", whitelistFile,
      "--soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 14",
      "--soloBarcodeReadLength 0", ## R1 is 151 nt: barcode + UMI + polyT
      "--soloCBmatchWLtype 1MM", ## min Hamming distance between V5 barcodes is 6
      "--soloUMIdedup 1MM_All",
      "--soloCellFilter None", ## keep every well; samples come from the sheet
      "--soloStrand Forward",
      "--soloFeatures", paste(unique(c("Gene", feature)), collapse = " "),
      "--soloCellReadStats Standard",
      if (keepBam) {
        paste("--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM CB UB",
              "--limitBAMsortRAM", format(param$ram * 0.4 * 1e9, scientific = FALSE))
      } else {
        "--outSAMtype None"
      },
      "--outFileNamePrefix", paste0(poolDirs[[nm]], "/"),
      if (ezIsSpecified(param$cmdOptions)) param$cmdOptions else ""
    )
    ezSystem(cmd)
    if (!file.exists(file.path(poolDirs[[nm]], "Solo.out", feature, "raw", "matrix.mtx"))) {
      stop("STARsolo produced no ", feature, " matrix for pool ", nm)
    }
    if (keepBam) {
      ezSystem(paste("samtools index -@", param$cores,
                     file.path(poolDirs[[nm]], "Aligned.sortedByCoord.out.bam")))
    }
  }
  poolDirs
}

##' Use finished STARsolo runs (one <dir>/<pool>/Solo.out per input row), e.g. to
##' redo the sample assignment after a sheet correction without remapping.
adoptBRBStarsolo <- function(starsoloDir, poolNames, feature) {
  poolDirs <- setNames(file.path(starsoloDir, poolNames), poolNames)
  miss <- poolNames[!file.exists(file.path(poolDirs, "Solo.out", feature, "raw"))]
  if (length(miss) > 0) stop("no Solo.out/", feature, "/raw under ", starsoloDir, " for pool(s): ", paste(miss, collapse = ", "))
  poolDirs
}

readBRBPool <- function(poolDir, feature) {
  rawDir <- file.path(poolDir, "Solo.out", feature, "raw")
  rd <- function(f) {
    fn <- file.path(rawDir, f)
    if (!file.exists(fn)) fn <- paste0(fn, ".gz")
    fn
  }
  m <- Matrix::readMM(rd("matrix.mtx"))
  feat <- read.delim(rd("features.tsv"), header = FALSE, colClasses = "character")
  rownames(m) <- feat$V1
  colnames(m) <- readLines(rd("barcodes.tsv"))
  cr <- read.delim(file.path(poolDir, "Solo.out", feature, "CellReads.stats"), check.names = FALSE)
  cr <- cr[cr$CB != "CBnotInPasslist", , drop = FALSE]
  rownames(cr) <- cr$CB
  list(m = as(m, "CsparseMatrix"), cr = cr, geneNames = setNames(feat$V2, feat$V1))
}

##' Each sheet library is taken from its Pool, or, when Pool is empty, from the
##' pool where its barcode has the most UMIs. UMIs of a sheet barcode in any
##' other pool are reported as contamination.
assignBRBLibraries <- function(sheet, pools) {
  umi <- sapply(pools, function(p) {
    v <- setNames(numeric(nrow(sheet)), sheet$Barcode)
    hit <- intersect(sheet$Barcode, colnames(p$m))
    v[hit] <- Matrix::colSums(p$m[, hit, drop = FALSE])
    v
  })
  umi <- matrix(umi, nrow = nrow(sheet), dimnames = list(sheet$Name, names(pools)))
  libs <- sheet
  libs$AssignedFrom <- ifelse(libs$Pool == "", "data", "sheet")
  libs$Pool[libs$Pool == ""] <- colnames(umi)[apply(umi[libs$Pool == "", , drop = FALSE], 1, which.max)]
  libs$UMIs <- umi[cbind(libs$Name, libs$Pool)]
  libs$UMIsOtherPools <- rowSums(umi) - libs$UMIs
  ## a data-assigned barcode with UMIs in two pools is ambiguous: flag, keep the larger
  libs$PoolAmbiguous <- libs$AssignedFrom == "data" & libs$UMIsOtherPools > 0.1 * libs$UMIs
  if (any(libs$PoolAmbiguous)) {
    ezLog(paste("WARNING: barcode found in more than one pool:", paste(libs$Name[libs$PoolAmbiguous], collapse = ", ")))
  }
  umiByPool <- data.frame(Name = sheet$Name, Well = sheet$Well, Barcode = sheet$Barcode,
                          Role = sheet$Role, AssignedPool = libs$Pool, umi, check.names = FALSE)
  list(libs = libs, umiByPool = umiByPool)
}

##' Genes of the count files: the STARsolo features, restricted to transcriptTypes
##' (a no-op for an index built with them; matters for adopted runs), in the
##' order of the reference GTF, as FeatureCounts writes them.
brbGenesToWrite <- function(features, param) {
  genes <- features
  if (ezIsSpecified(param$transcriptTypes)) {
    anno <- ezFeatureAnnotation(param$ezRef@refAnnotationFile, dataFeatureType = "gene")
    genes <- intersect(genes, anno$gene_id[anno$type %in% param$transcriptTypes])
  }
  gtfOrder <- unique(data.table::fread(
    cmd = paste("grep -o 'gene_id \"[^\"]*\"'", shQuote(param$ezRef@refFeatureFile)),
    header = FALSE, sep = "\t", quote = "", col.names = "x"
  )$x)
  gtfOrder <- sub('^gene_id "(.*)"$', "\\1", gtfOrder)
  c(intersect(gtfOrder, genes), setdiff(genes, gtfOrder))
}

##' <Name>.txt (Identifier / matchCounts, deduplicated UMIs) and <Name>-stats.txt
##' (Rsubread status rows, filled with READ counts from STARsolo CellReads.stats).
writeBRBCounts <- function(libs, pools, genes, resultDir) {
  libs$ReadCount <- 0
  libs$GenesDetected <- 0
  for (i in seq_len(nrow(libs))) {
    p <- pools[[libs$Pool[i]]]
    bc <- libs$Barcode[i]
    cnt <- setNames(integer(length(genes)), genes)
    if (bc %in% colnames(p$m)) {
      v <- p$m[, bc]
      cnt[] <- as.integer(v[genes])
    }
    counts <- matrix(cnt, ncol = 1, dimnames = list(genes, "matchCounts"))
    ezWrite.table(counts, file = file.path(resultDir, paste0(libs$Name[i], ".txt")))

    cr <- if (bc %in% rownames(p$cr)) p$cr[bc, ] else NULL
    st <- setNames(numeric(length(BRB_STAT_ROWS)), BRB_STAT_ROWS)
    if (!is.null(cr)) {
      st["Assigned"] <- cr$featureU
      st["Unassigned_MultiMapping"] <- cr$genomeM
      st["Unassigned_NoFeatures"] <- cr$genomeU - cr$featureU
      st["Unassigned_Unmapped"] <- cr$cbMatch - cr$genomeU - cr$genomeM
      libs$ReadCount[i] <- cr$cbMatch
    }
    ezWrite.table(data.frame(Status = names(st), Count = st),
                  file = file.path(resultDir, paste0(libs$Name[i], "-stats.txt")), row.names = FALSE)
    libs$GenesDetected[i] <- sum(cnt > 0)
  }
  libs
}

##' Per-library dataset.tsv in the FeatureCounts column layout, including the
##' Read Count that is only known after the run (fallback for a manual import).
writeBRBDataset <- function(libs, input, output, param, resultDir) {
  rel <- if (output$hasColumn("ResultDir")) as.character(output$getColumn("ResultDir"))[1] else resultDir
  ds <- data.frame(
    Name = libs$Name,
    `Count [File,Link]` = file.path(rel, paste0(libs$Name, ".txt")),
    `Stats [File,Link]` = file.path(rel, paste0(libs$Name, "-stats.txt")),
    Species = as.character(input$getColumn("Species"))[1],
    refBuild = param$refBuild,
    featureLevel = "gene",
    refFeatureFile = param$refFeatureFile,
    strandMode = "sense",
    paired = "false",
    `Read Count` = libs$ReadCount,
    transcriptTypes = paste(param$transcriptTypes, collapse = ","),
    check.names = FALSE
  )
  facCols <- grep("\\[Factor\\]$", colnames(libs), value = TRUE)
  for (cc in facCols) ds[[cc]] <- libs[[cc]]
  ds[["Role [Factor]"]] <- libs$Role
  ds[["Pool [Factor]"]] <- libs$Pool
  for (cc in c("Sample Id [B-Fabric]")) if (cc %in% colnames(libs)) ds[[cc]] <- libs[[cc]]
  if (input$hasColumn("Order Id")) ds[["Order Id [B-Fabric]"]] <- as.character(input$getColumn("Order Id"))[1]
  ezWrite.table(ds, file.path(resultDir, "dataset.tsv"), row.names = FALSE)
}

##' Copy the per-pool STARsolo summaries and gzipped raw matrices (and the BAM if
##' kept) into the result; returns one row of Summary.csv metrics per pool.
keepBRBPoolOutputs <- function(poolDirs, resultDir, feature, param) {
  rows <- list()
  for (nm in names(poolDirs)) {
    src <- file.path(poolDirs[[nm]], "Solo.out", feature)
    dst <- file.path(resultDir, "pools", nm)
    dir.create(file.path(dst, "raw"), recursive = TRUE, showWarnings = FALSE)
    for (f in c("Summary.csv", "Features.stats", "CellReads.stats")) {
      file.copy(file.path(src, f), dst, overwrite = TRUE)
    }
    for (f in c("matrix.mtx", "features.tsv", "barcodes.tsv")) {
      fn <- file.path(src, "raw", f)
      if (file.exists(fn)) {
        ezSystem(paste("pigz -p 4 -c", fn, ">", file.path(dst, "raw", paste0(f, ".gz"))))
      } else {
        file.copy(paste0(fn, ".gz"), file.path(dst, "raw"), overwrite = TRUE)
      }
    }
    logFinal <- file.path(poolDirs[[nm]], "Log.final.out")
    if (file.exists(logFinal)) file.copy(logFinal, dst, overwrite = TRUE)
    bam <- file.path(poolDirs[[nm]], "Aligned.sortedByCoord.out.bam")
    if (isTRUE(as.logical(param$keepBam)) && file.exists(bam)) {
      file.rename(bam, file.path(dst, paste0(nm, ".bam")))
      file.rename(paste0(bam, ".bai"), file.path(dst, paste0(nm, ".bam.bai")))
    }
    s <- read.csv(file.path(src, "Summary.csv"), header = FALSE, col.names = c("metric", "value"))
    rows[[nm]] <- data.frame(Pool = nm, t(setNames(s$value, s$metric)), check.names = FALSE)
  }
  do.call(rbind, rows)
}

##' Mean log2 CPM of Y-chromosome marker genes per library (sex check), for
##' species with a known marker set; NULL otherwise.
brbYGeneSignal <- function(pools, libs, param) {
  yGenes <- c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "RPS4Y1")
  out <- lapply(seq_len(nrow(libs)), function(i) {
    p <- pools[[libs$Pool[i]]]
    bc <- libs$Barcode[i]
    if (!bc %in% colnames(p$m)) return(NA_real_)
    v <- p$m[, bc]
    idx <- names(p$geneNames)[p$geneNames %in% yGenes]
    if (length(idx) == 0 || sum(v) == 0) return(NA_real_)
    log2(mean(v[idx] / sum(v) * 1e6) + 1)
  })
  y <- unlist(out)
  if (all(is.na(y))) return(NULL)
  data.frame(Name = libs$Name, Role = libs$Role, yLog2CPM = y)
}

makeBRBseqReport <- function(resultDir, param, ...) {
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  setwd(resultDir)
  makeQuartoReport(
    param = param,
    ...,
    qmdFile = "BRBseq.qmd",
    reportTitle = paste("BRB-seq -", param$name),
    buttons = TRUE,
    number = TRUE,
    colour = TRUE
  )
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodBRBseq(input=NA, output=NA, param=NA)
##' @description Use this reference class to demultiplex and count BRB-seq pools
##' with STARsolo into per-library FeatureCounts-format count files.
EzAppBRBseq <-
  setRefClass(
    "EzAppBRBseq",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBRBseq
        name <<- "EzAppBRBseq"
        appDefaults <<- rbind(
          plateVersion = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "BRB V5 plate version (V5A..V5D); 'auto' detects it from the sample sheet's well/barcode pairs."
          ),
          sampleSheet = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Uploaded sample sheet (TSV: Name, Barcode, Well, Role, Pool, ...). Empty: use <project>/o<order>_metaData/o<order>_BRB_SampleSheet.tsv."
          ),
          barcodeFile = ezFrame(
            Type = "character",
            DefaultValue = BRB_KIT_FILE,
            Description = "Barcode kit table (Plate, Well, Barcode)."
          ),
          countIntrons = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Count reads on full gene bodies (STARsolo GeneFull) instead of exons only."
          ),
          keepBam = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Keep one coordinate-sorted BAM per pool (CB/UB tagged)."
          ),
          starsoloDir = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Optional directory with one finished STARsolo run per pool (<dir>/<pool>/Solo.out); skips mapping."
          ),
          transcriptTypes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "Transcript types in the STAR index and the count files."
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          )
        )
      }
    )
  )
