###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @template app-template
##' @templateVar method ezMethodVeloCyto(input=NA, output=NA, param=NA)
##' @description Use this reference class to run velocyto on CellRanger outputs
##' @author Lennart Opitz
EzAppVeloCyto <-
    setRefClass("EzAppVeloCyto",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodVeloCyto
                        name <<- "EzAppVeloCyto"
                        appDefaults <<- rbind(
                            outputDir = ezFrame(
                                Type = "character",
                                DefaultValue = ".",
                                Description = "Output directory"
                            )
                        )
                    }
                )
    )


ezMethodVeloCyto <- function(input=NA, output=NA, param=NA){
    gtfFile <- param$ezRef["refFeatureFile"]
    sampleName <- input$getNames()
    
    ##Copy data to scratch
    cellRangerPath <- file.path(input$dataRoot,input$getColumn("ResultDir"))
    cmd <- paste('cp -R',  cellRangerPath, '.')
    ezSystem(cmd)

    sampleDir <- basename(cellRangerPath)
    cmd <- paste('rsync -av --remove-source-files', paste0(sampleDir,'/*'), paste0(sampleDir,'/outs'))
    ezSystem(cmd)

    cwd <- getwd()
    sampleBam <- list.files('.', pattern = 'sample_alignments.bam$', recursive=TRUE)
    if(length(sampleBam) == 1L) { #CellRanger Multi Output
        setwd(dirname(sampleBam))
        system('mv sample_alignments.bam possorted_genome_bam.bam')
        system('samtools index possorted_genome_bam.bam')
        system('mv sample_filtered_feature_bc_matrix filtered_feature_bc_matrix')
        system(sprintf('mv * %s', file.path(cwd, sampleDir, "outs")))
        setwd(cwd)
    }
    out <- tryCatch(local_CondaEnv("gi_velocyto", pathToMiniConda = "/usr/local/ngseq/miniforge3"), error = function(e) NULL)
    cmd <- paste('velocyto run10x', sampleDir, gtfFile, '-@', param$cores)
    ezSystem(cmd)
    file.copy(file.path(sampleName, 'velocyto', paste0(sampleName,'.loom')), '.')
    ezSystem(paste('rm -Rf ', sampleName))
    return('success')
}
