###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodXeniumQC <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
    ezLoadPackage('jsonlite')
    ezLoadPackage('dplyr')

    setwdNew('XeniumQC')
    dataset <- input$meta
    xeniumPaths <- input$getColumn("XeniumPath")

    ###Collect Stats from each Xenium output directory
    stats <- data.frame()

    for (j in 1:nrow(dataset)){
        sampleName <- rownames(dataset)[j]
        xeniumPath <- file.path(param$dataRoot, xeniumPaths[j])

        # Read metrics_summary.csv
        metricsFile <- file.path(xeniumPath, 'metrics_summary.csv')
        if (!file.exists(metricsFile)) {
            warning(paste("Metrics file not found for", sampleName, ":", metricsFile))
            next
        }

        sampleMetrics <- ezRead.table(metricsFile, sep = ',', header = TRUE)

        # Extract panel information from gene_panel.json
        panelFile <- file.path(xeniumPath, 'gene_panel.json')
        panelName <- NA
        panelDesignId <- NA
        if (file.exists(panelFile)) {
            panelData <- tryCatch({
                jsonlite::fromJSON(panelFile)
            }, error = function(e) {
                list()
            })
            # Correct path in the JSON structure
            if (!is.null(panelData$payload$panel$identity$name)) {
                panelName <- panelData$payload$panel$identity$name
            }
            if (!is.null(panelData$payload$panel$identity$design_id)) {
                panelDesignId <- panelData$payload$panel$identity$design_id
            }
        }

        # Remove panel columns from sampleMetrics if they exist to avoid duplicates
        sampleMetrics$panel_name <- NULL
        sampleMetrics$panel_design_id <- NULL

        # Add sample name and panel info
        sampleStats <- data.frame(
            sampleName = sampleName,
            panel_name = ifelse(is.na(panelName), "NA", panelName),
            panel_design_id = ifelse(is.na(panelDesignId), "NA", panelDesignId),
            sampleMetrics,
            check.names = FALSE
        )

        if (nrow(stats) == 0){
            stats <- sampleStats
        } else {
            stats <- rbind(stats, sampleStats)
        }
    }

    # Clean up column names and round numeric values
    for (i in 2:ncol(stats)){
        if (is.numeric(stats[,i])){
            stats[,i] <- round(stats[,i], 3)
        }
    }

    # Make unique rownames by combining name with row index if there are duplicates
    if (anyDuplicated(stats$sampleName)) {
        rownames(stats) <- make.unique(as.character(stats$sampleName), sep = "_")
    } else {
        rownames(stats) <- stats$sampleName
    }

    colnames(stats) <- gsub(' ', '_', colnames(stats))
    colnames(stats) <- gsub('-', '__', colnames(stats))

    # Reorder columns to put interesting metrics first
    priority_cols <- c("sampleName", "panel_name", "panel_design_id",
                       "num_cells_detected", "median_genes_per_cell", "median_transcripts_per_cell",
                       "total_high_quality_decoded_transcripts", "fraction_transcripts_decoded_q20",
                       "fraction_transcripts_assigned", "fraction_empty_cells",
                       "segmented_cell_stain_frac", "decoded_transcripts_per_100um2")

    # Get remaining columns
    remaining_cols <- setdiff(colnames(stats), priority_cols)

    # Reorder: priority columns first, then the rest
    available_priority <- intersect(priority_cols, colnames(stats))
    stats <- stats[, c(available_priority, remaining_cols)]

    # Write summary table
    ezWrite.table(stats, file = 'metrics_summary.tsv', row.names = FALSE)

    # Save param for Rmd
    write_rds(param, 'param.rds')

    # Generate report
    reportTitle <- 'XeniumQC - Multiple Sample QC Metrics'
    makeRmdReport(rmdFile = "XeniumQC.Rmd", reportTitle = reportTitle)

    return("Success")
}

##' @template app-template
##' @templateVar method ezMethodXeniumQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run XeniumQC
##' @section Functions:
EzAppXeniumQC <-
  setRefClass("EzAppXeniumQC",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodXeniumQC
        name <<- "EzAppXeniumQC"
        appDefaults <<- data.frame()
      }
    )
  )
