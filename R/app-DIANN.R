###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDIANN = function(input=NA, output=NA, param=NA,
                           htmlFile="00index.html"){
  library(yaml)
  setwdNew(basename(output$getColumn("Result")))
  ## todo write inputs / params / dataset
  dir.create("work")
  write_yaml(list(param=param), "work/params.yml")
  ## get the input
  rawDir <- "work/input/raw"
  dir.create(rawDir, recursive = TRUE)
  rawFiles <- input$getColumn("RAW")
  rawUrl <- file.path("fgcz-ms.uzh.ch:/srv/www/htdocs", rawFiles)
  for (x in rawUrl){
    ezSystem(paste("scp", x, rawDir))
  }
  
  
  cmd <- paste0("scp ", "fgcz-r-033:", param$`03b_additional_fasta_database_path`,
  " work/input/", basename(param$`03b_additional_fasta_database_path`))
  ezSystem(cmd)
  cmd <- paste0("scp ", "fgcz-r-036:", param$`order_fasta`,
                " work/input/", "order.fasta")
  ezSystem(cmd)
  
  ds <- input$meta |> 
    rownames_to_column("Resource") |> 
    rename("Relative Path"=RAW) |>
    mutate(File=basename(`Relative Path`))
  write_csv(ds, "work/dataset.csv")
  
  
  snakeFile <- system.file("templates/Snakefile.DIANN3step.smk", package = "ezRun", mustWork = TRUE)
  file.copy(snakeFile, basename(snakeFile))
  ## example cmd: snakemake uses the on-the-fly generated environment
  ## snakemake -s /scratch/cache/bfabric/bfabric_app_runner/ephemeral/env_dbox3jko/lib/python3.13/site-packages/diann_runner/Snakefile.DIANN3step.smk --cores 64 -p all -d /scratch/5858/A386_DIANN_v23/WU339908/work
  cmd <- paste("snakemake",
               "-s", basename(snakeFile),
               "--cores 64 -p all",
               "-d ./work")
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppDIANN <-
  setRefClass("EzAppDIANN",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDIANN
                  name <<- "EzAppDIANN"
                  appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=TRUE, Description="whether to run the GO analysis"),
                                        nSampleClusters=ezFrame(Type="numeric", DefaultValue=6, Description="Number of SampleClusters, default value 6"),
                                        selectByFtest=ezFrame(Type="logical", DefaultValue=FALSE, Description="select topGenes by Test instead of SD"),
                                        topGeneSize=ezFrame(Type="numeric", DefaultValue=100, Description="number of genes to consider in gene clustering, mds etc"))
                }
              )
  )
