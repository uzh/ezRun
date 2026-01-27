###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodVisiumQC <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  ezLoadPackage('SpotClean')
  ezLoadPackage('S4Vectors')
  ezLoadPackage('Seurat')

  setwdNew('VisiumQC')
  dataset <- input$meta
  if (all(grepl('^H', input$getColumn('Slide')))) {
    param$visiumType <- 'HD'
  } else {
    param$visiumType <- 'SD'
  }
  param$sizeFactors <- as.numeric(unlist(strsplit(param$sizeFactors, ',')))
  ###CollectStats
  for (j in 1:nrow(dataset)) {
    sampleName <- rownames(dataset)[j]
    samplePath <- file.path(param$dataRoot, dirname(dataset[['Count [Link]']][j]))
    umiCounts <- sum(
      ezRead.table(file.path(
        samplePath,
        paste0(sampleName, '-counts.txt')
      ))$matchCount
    )
    sampleStats <- data.frame(
      sampleName = sampleName,
      umiCounts = umiCounts,
      ezRead.table(file.path(samplePath, 'metrics_summary.csv'), sep = ','),
      check.names = FALSE
    )
    if (j == 1) {
      stats <- sampleStats
    } else {
      commonCols <- intersect(colnames(stats),colnames(sampleStats))
      sampleStats <- sampleStats[,commonCols]
      stats <- stats[,commonCols]
      stats <- rbind(stats, sampleStats)
    }
  }

  for (i in 2:ncol(stats)) {
    stats[, i] <- round(stats[, i], 3)
  }

  rownames(stats) <- sub('-spaceRanger', '', rownames(stats))
  colnames(stats) <- gsub(' ', '_', colnames(stats))
  colnames(stats) <- gsub('-', '__', colnames(stats))
  ezWrite.table(stats, file = 'metrics_summary.tsv', row.names = FALSE)

  if (param$visiumType == 'SD') {
    myPlots <- list()
    for (j in 1:nrow(dataset)) {
      sampleName <- rownames(dataset)[j]
      samplePath <- file.path(param$dataRoot, dirname(dataset[['Count [Link]']][j]))
      data_raw <- read10xRaw(file.path(samplePath, "raw_feature_bc_matrix"))
      if (
        file.exists(file.path(samplePath, "spatial", "tissue_positions.csv"))
      ) {
        tissueFile <- file.path(samplePath, "spatial", "tissue_positions.csv")
      } else {
        tissueFile <- file.path(
          samplePath,
          "spatial",
          "tissue_positions_list.csv"
        )
      }
      if (
        file.exists(file.path(samplePath, "spatial", "tissue_hires_image.png"))
      ) {
        imageFile <- file.path(samplePath, "spatial", "tissue_hires_image.png")
      } else {
        imageFile <- file.path(samplePath, "spatial", "tissue_lowres_image.png")
      }
      scaleFile <- file.path(samplePath, "spatial", "scalefactors_json.json")
      data_slide_info <- read10xSlide(tissueFile, imageFile, scaleFile)
      slide_obj <- createSlide(data_raw, data_slide_info)
      metadata(slide_obj)$slide$total_counts <- Matrix::colSums(data_raw)
      metadata(slide_obj)$slide$tissue <- as.numeric(as.character(
        metadata(slide_obj)$slide$tissue
      ))

      myPlots[[sampleName]] <- list()
      myPlots[[sampleName]][['Overview']] <- visualizeSlide(
        slide_obj,
        title = paste(sampleName, "overview")
      )
      myPlots[[sampleName]][['Tissue']] <- visualizeHeatmap(
        slide_obj,
        "tissue",
        title = paste(sampleName, "tissue")
      )
      myPlots[[sampleName]][['Signal']] <- visualizeHeatmap(
        slide_obj,
        "total_counts",
        title = paste(sampleName, "total_counts")
      )
    }
    write_rds(myPlots, 'myPlots.rds')

    sampleName <- rownames(dataset)[1]
    samplePath <- samplePath <- file.path(param$dataRoot, dirname(dataset[['Count [Link]']][1]))
    img = Read10X_Image(
      file.path(samplePath, "spatial"),
      image.name = "tissue_hires_image.png"
    )
    img@scale.factors$lowres <- img@scale.factors$hires
    scData <- Load10X_Spatial(samplePath, image = img)
    write_rds(scData, 'scData.rds')
  }
  write_rds(param, 'param.rds')

  reportTitle <- 'VisiumQC - MultipleSample QC Metrics'
  makeRmdReport(rmdFile = "VisiumQC.Rmd", reportTitle = reportTitle)
  #CleanUp
  if (param$visiumType == 'SD') {
    ezSystem('rm scData.rds myPlots.rds')
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodVisiumQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
EzAppVisiumQC <-
  setRefClass(
    "EzAppVisiumQC",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodVisiumQC
        name <<- "EzAppVisiumQC"
        appDefaults <<- rbind(
          sizeFactors = ezFrame(
            Type = "character",
            DefaultValue = '1, 3, 5, 10, 20, 30',
            Description = "size factors to try for spatialFeaturePlots"
          )
        )
      }
    )
  )
