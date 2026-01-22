ezSetupProjectInGitRepo <- function(
  project = NULL,
  gitRepo = paste0("~/git/", system("whoami", intern = TRUE), "-scripts")
) {
  stopifnot(!is.null(project))
  stopifnot(file.exists(gitRepo))
  setwdNew(file.path(gitRepo, project))
  projectMasterScript <- paste0("00-", project, "-custom-analyses.R")
  if (!file.exists(projectMasterScript)) {
    file.copy(
      system.file("templates/00-manage-project-template.R", package = "ezRun"),
      projectMasterScript
    )
  }
  rstudioapi::navigateToFile(projectMasterScript)
}

ezSetupAnalysis <- function(
  project = NULL,
  analysisName = "custom-analysis",
  subDir = paste0(analysisDate, "_", analysisName),
  parentAnalysis = NULL,
  user = system("whoami", intern = TRUE),
  analysisDate = format(Sys.time(), "%Y-%m-%d"),
  gitRepo = paste0("~/git/", system("whoami", intern = TRUE), "-scripts"),
  templateFile = system.file("templates/report-template.Rmd", package = "ezRun")
) {
  stopifnot(!is.null(project))
  stopifnot(file.exists(gitRepo))
  if (!file.exists(file.path(gitRepo, project))) {
    dir.create(file.path(gitRepo, project))
  }
  rmdFile <- paste0(
    gitRepo,
    "/",
    project,
    "/",
    analysisDate,
    "_",
    analysisName,
    ".Rmd"
  )
  stopifnot(!file.exists(rmdFile))
  file.copy(templateFile, rmdFile)
  analysisDir <- paste0("/srv/GT/analysis/", user, "/", project)
  setwdNew(analysisDir)
  if (!is.null(subDir)) {
    setwdNew(subDir)
  }
  stopifnot(!file.exists(basename(rmdFile)))
  file.symlink(rmdFile, ".")
  #rstudioapi::navigateToFile(basename(rmdFile)) ## that does not work, because rstudio will follow the link and open the original file
  rstudioapi::filesPaneNavigate(normalizePath("."))
}


ezPublishAnalysis <- function() {
  ## copy analysis to gstore
  ## register in B-Fabric
  ## write a B-Fabric comment
  ## commit/push the analysis-script to gitlab
  ##
}

## Nice to have:
## - storage duration for analyses
## - deletion of analyses
