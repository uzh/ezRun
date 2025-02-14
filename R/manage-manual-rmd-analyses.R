
ezSetupProjectInGitRepo <- function(project=NULL,
                                    gitRepo=paste0("~/git/", system("whoami", intern = TRUE), "-scripts")){
  stopifnot(!is.null(project))
  stopifnot(file.exists(gitRepo))
  setwdNew(file.path(gitRepo, project))
  file.copy(system.file("templates/00-manage-project-template.R", package="ezRun"), paste("00-", project, "-custom-analyses.R"))
  rstudioapi::navigateToFile(paste("00-", project, "-custom-analyses.R"))
}

ezSetupAnalysis <- function(project=NULL, analysisName="custom-analysis",
                            parentAnalysis=NULL,
                            user=system("whoami", intern = TRUE),
                            analysisDate=format(Sys.time(), "%Y-%m-%d"),
                            gitRepo=paste0("~/git/", system("whoami", intern = TRUE), "-scripts"),
                            templateFile=system.file("templates/report-template.Rmd", package="ezRun")
){
  stopifnot(!is.null(project))
  stopifnot(file.exists(gitRepo))
  if (!file.exists(file.path(gitRepo, project))){
    dir.create(file.path(gitRepo, project))
  }
  rmdFile <- paste0(gitRepo, "/", project, "/", analysisDate, "_", analysisName, ".Rmd")
  file.copy(templateFile, rmdFile)
  analysisDir <- paste0("/srv/GT/analysis/", user, "/", project, "/", analysisDate, "_",
                        analysisName)
  setwdNew(analysisDir)
  file.symlink(rmdFile, ".")
  rstudioapi::navigateToFile(file.path(getwd(), basename(rmdFile)))
}


ezPublishAnalysis <- function(){
  ## copy analysis to gstore
  ## register in B-Fabric
  ## write a B-Fabric comment
  ## commit/push the analysis-script to gitlab
  ## 
}

## TBD:
## - storage duration for analyses
## - deletion of analyses





