###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

# DIA-NN SUSHI app. Wraps wolski/diann-runner's snakemake pipeline.
# Container runtime: apptainer (the runner is patched to dispatch to apptainer
# when docker is unavailable, which is the case for trxcopy@fgcz-c-050).
# Params input: bundled DIA / DDA YAML templates under inst/templates/, with
# an optional `customParamsYml` path override per-run.

# Resolve the params.yml file the snakemake workflow will read.
# Falls back to the bundled DIA template when paramsTemplate/customParamsYml
# are not set, so existing job submissions stay working.
ezDiannResolveParamsYml <- function(param) {
  if (!is.null(param$customParamsYml) &&
      nzchar(param$customParamsYml) &&
      toupper(param$customParamsYml) != "NONE") {
    if (!file.exists(param$customParamsYml)) {
      stop("customParamsYml does not exist: ", param$customParamsYml)
    }
    return(param$customParamsYml)
  }
  template <- if (!is.null(param$paramsTemplate) && nzchar(param$paramsTemplate)) {
    param$paramsTemplate
  } else {
    "default-DIA"
  }
  path <- system.file(
    file.path("templates", paste0("DIANN_params_", template, ".yml")),
    package = "ezRun",
    mustWork = FALSE
  )
  if (!nzchar(path)) {
    stop("Unknown paramsTemplate '", template,
         "'. Ship DIANN_params_", template, ".yml in ezRun/inst/templates/.")
  }
  path
}

# Stage a single file from a local NFS path or remote SSH source into destDir.
# Treats anything under /srv/, /misc/ or /scratch/ as already-local; everything
# else is SCPed from `fgcz-ms.uzh.ch:/srv/www/htdocs`, matching the existing
# bfabric-app-runner inputs.yml SSH convention.
ezDiannStageOne <- function(relativeOrAbsolute, destDir) {
  dir.create(destDir, recursive = TRUE, showWarnings = FALSE)
  if (startsWith(relativeOrAbsolute, "/srv/") ||
      startsWith(relativeOrAbsolute, "/misc/") ||
      startsWith(relativeOrAbsolute, "/scratch/")) {
    ezSystem(paste("cp", shQuote(relativeOrAbsolute), shQuote(destDir)))
  } else {
    src <- file.path("fgcz-ms.uzh.ch:/srv/www/htdocs", relativeOrAbsolute)
    ezSystem(paste("scp", src, shQuote(destDir)))
  }
}

ezMethodDIANN <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)

  dir.create("work", showWarnings = FALSE)
  rawDir <- "work/input/raw"
  inputDir <- "work/input"
  dir.create(rawDir, recursive = TRUE, showWarnings = FALSE)

  ## 1. Resolve params.yml: bundled template OR user override ---------------
  paramsYmlSrc <- ezDiannResolveParamsYml(param)
  yamlIn <- yaml::read_yaml(paramsYmlSrc)
  if (is.null(yamlIn$params)) {
    stop("params.yml at ", paramsYmlSrc, " has no top-level `params:` block.")
  }
  yamlParams <- yamlIn$params

  ## 2. Inject registration block from SUSHI / B-Fabric job context ---------
  registration <- list(
    workunit_id = as.character(param$workunit_id %||% "0"),
    container_id = as.character(param$container_id %||% "0"),
    application_id = "386",
    application_name = "DIANN_v23",
    container_type = "order",
    storage_id = "2",
    storage_output_folder = as.character(param$resultDir %||% "")
  )
  yaml::write_yaml(
    list(params = yamlParams, registration = registration),
    "work/params.yml"
  )

  ## 3. Stage raw files (one SCP per row of the SUSHI dataset) --------------
  rawFiles <- input$getColumn("RAW")
  for (rf in rawFiles) {
    ezDiannStageOne(rf, rawDir)
  }

  ## 4. Stage FASTA(s) referenced in params.yml -----------------------------
  ## The runner expects FASTAs under work/input/ keyed by basename.
  fastaPath <- yamlParams[["03_fasta_database_path"]]
  if (!is.null(fastaPath) && nzchar(fastaPath) && toupper(fastaPath) != "NONE") {
    ezDiannStageOne(fastaPath, inputDir)
  }
  altFastaPath <- yamlParams[["03b_additional_fasta_database_path"]]
  if (!is.null(altFastaPath) && nzchar(altFastaPath) &&
      toupper(altFastaPath) != "NONE") {
    ezDiannStageOne(altFastaPath, inputDir)
  }
  if (!is.null(param$order_fasta) && nzchar(param$order_fasta)) {
    ## SUSHI's per-order FASTA lands in input/order.fasta regardless of source.
    if (startsWith(param$order_fasta, "/")) {
      ezSystem(paste("cp", shQuote(param$order_fasta),
                     file.path(inputDir, "order.fasta")))
    } else {
      ezSystem(paste("scp",
                     file.path("fgcz-r-036:", param$order_fasta),
                     file.path(inputDir, "order.fasta")))
    }
  }

  ## 5. Build work/dataset.csv from the SUSHI dataset (3-col format) --------
  ds <- input$meta |>
    tibble::rownames_to_column("Name") |>
    dplyr::rename(`Relative Path` = "RAW")
  outCols <- intersect(c("Relative Path", "Name", "Grouping Var"), colnames(ds))
  readr::write_csv(ds[, outCols, drop = FALSE], "work/dataset.csv")

  ## 6. Run snakemake via the patched, pre-installed diann-runner -----------
  runnerSrc <- Sys.getenv(
    "DIANN_RUNNER_PATH",
    unset = "/usr/local/ngseq/opt/diann-runner"
  )
  if (!dir.exists(runnerSrc)) {
    stop("diann-runner not found at ", runnerSrc,
         ". Install via fgcz-r-029 or set DIANN_RUNNER_PATH.")
  }
  file.symlink(runnerSrc, "diann-runner")
  Sys.setenv(DIANN_RUNNER_SIF_CACHE_DIR = "/misc/ngseq12/opt/sif")
  snakeFile <- "diann-runner/src/diann_runner/Snakefile.DIANN3step.smk"
  stopifnot(file.exists(snakeFile))
  ezSystem(paste(
    "snakemake -s", snakeFile,
    "--cores 64 -p all -d ./work"
  ))

  ## 7. Stage outputs back to SUSHI's expected dataset columns --------------
  ezSystem(paste(
    "mv", "work/out-DIANN_quantC",
    shQuote(output$getColumn("DIANN Quant"))
  ))
  ezSystem(paste(
    "mv", "work/qc_result",
    shQuote(output$getColumn("qc_result"))
  ))

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDIANN(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description SUSHI app for DIA-NN proteomics quantification + prolfqua QC.
EzAppDIANN <-
  setRefClass(
    "EzAppDIANN",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDIANN
        name <<- "EzAppDIANN"
        appDefaults <<- rbind(
          paramsTemplate = ezFrame(
            Type = "character",
            DefaultValue = "default-DIA",
            Description = paste(
              "Name of a bundled DIANN params template",
              "('default-DIA' or 'default-DDA'); ignored when",
              "customParamsYml is set."
            )
          ),
          customParamsYml = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = paste(
              "Absolute path to a user-supplied params.yml (under /srv/gstore",
              "or /misc); overrides paramsTemplate when set."
            )
          ),
          order_fasta = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Optional per-order custom FASTA, staged as input/order.fasta."
          )
        )
      }
    )
  )
