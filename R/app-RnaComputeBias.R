###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


EzAppRnaComputeBias <-
    setRefClass("EzAppRnaComputeBias",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodRnaComputeBias
                        name <<- "EzAppRnaComputeBias"
                        appDefaults <<- rbind(minReadsPerSample=ezFrame(Type="numeric",
                                                                  DefaultValue=30000,
                                                                  Description=""),
                                              maxReadsPerSample=ezFrame(Type="numeric",
                                                                        DefaultValue=5e6,
                                                                        Description=""),
                                              minReadsPerGene=ezFrame(Type="numeric",
                                                                        DefaultValue=3,
                                                                        Description=""),
                                              minPresentFraction=ezFrame(Type="numeric",
                                                                         DefaultValue=0.2,
                                                                         Description="")
                                              )
                    }
                )
    )


ezMethodRnaComputeBias <- function(input, output, param){
    ezWrite.table(input$meta, 'input_dataset.tsv')
    dsFile <- 'input_dataset.tsv'
    qcSummaryDir = getwd()
    dsName = 'report-RNAseq'
    dir.create(dsName)
    plateId <- unique(input$getColumn("PlateName"))
    orderId <- unique(input$getColumn("Order Id"))
    tag <- paste0("Plate_",plateId, '_o', orderId)
    ezComputeBias(dsFile = dsFile, dsName = dsName, param = param,  qcSummaryDir = qcSummaryDir, 
                  minReadsPerSample=param$minReadsPerSample, maxReadsPerSample=param$maxReadsPerSample,
                  minReadsPerGene=param$minReadsPerGene, minPresentFraction=param$minPresentFraction, tag = tag)
    return('success')
}
