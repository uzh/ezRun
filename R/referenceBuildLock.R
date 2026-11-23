###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Build a shared reference exactly once
##' @description Many SUSHI jobs can start at the same moment and all need the
##'   same reference (e.g. an index built on demand under /srv/GT/reference).
##'   This function makes sure that exactly one of them builds it, while the
##'   others wait and then use the result.
##'
##'   The lock is a directory: \code{dir.create()} creates it, or fails if it
##'   already exists, in one atomic step, so only one job can hold it, also
##'   across nodes. A lock file written after a \code{file.exists()} check is
##'   not safe: jobs that start together all pass the check and all build.
##'
##'   Only the job that created the lock removes it, also when the build fails.
##'   A failing build therefore never leaves a stale lock, and a failing job
##'   never removes the lock of another job.
##'
##'   \code{build} should write into a private directory and move the finished
##'   result into place as its last step (\code{file.rename()} within the same
##'   filesystem), so that \code{isComplete} never sees a half-written
##'   reference.
##' @param lockDir path of the lock directory, usually
##'   \code{paste0(refDir, ".lock")}. A lock file of the same name written by an
##'   older ezRun version is respected as well.
##' @param isComplete a function without arguments that returns \code{TRUE}
##'   when the reference is complete and can be used.
##' @param build a function without arguments that builds the reference. It is
##'   only called by the job holding the lock.
##' @param what a short label of the reference, used in log and error messages.
##' @param timeout minutes to wait for another job before giving up.
##' @param pollSeconds seconds between two checks while waiting.
##' @return \code{TRUE} (invisibly) if this job built the reference,
##'   \code{FALSE} if it was complete already or built by another job.
##' @keywords internal
ez_build_reference_once <- function(
  lockDir,
  isComplete,
  build,
  what = "reference",
  timeout = INDEX_BUILD_TIMEOUT,
  pollSeconds = 60
) {
  dir.create(dirname(lockDir), recursive = TRUE, showWarnings = FALSE)
  deadline <- Sys.time() + timeout * 60
  hasLock <- FALSE
  waitLogged <- FALSE
  while (!isComplete()) {
    if (dir.create(lockDir, showWarnings = FALSE)) {
      hasLock <- TRUE
      break
    }
    ## somebody else builds and we wait
    lockOwner <- ez_build_lock_owner(lockDir)
    if (Sys.time() > deadline) {
      stop(
        what, " building still in progress after ", timeout, " min; lock: ",
        lockDir, " (", lockOwner, "). If that job no longer runs, ",
        "remove the lock and resubmit."
      )
    }
    if (!waitLogged) {
      ezLog("waiting for ", what, " built by another job: ", lockOwner)
      waitLogged <- TRUE
    }
    Sys.sleep(pollSeconds)
  }
  if (!hasLock) {
    return(invisible(FALSE))
  }
  on.exit(unlink(lockDir, recursive = TRUE), add = TRUE)
  writeLines(
    paste0(
      "host=", Sys.info()[["nodename"]],
      " pid=", Sys.getpid(),
      " slurmJob=", Sys.getenv("SLURM_JOB_ID", "none"),
      " started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    file.path(lockDir, "owner.txt")
  )
  if (isComplete()) {
    ## another job finished between our last check and taking the lock
    return(invisible(FALSE))
  }
  build()
  if (!isComplete()) {
    stop(what, " build finished, but the result is incomplete")
  }
  invisible(TRUE)
}

##' @title Create a private build directory next to a reference
##' @description Returns a new, empty directory next to \code{target} (same
##'   filesystem, so its content can later be moved into place with an atomic
##'   \code{file.rename()}). The name contains host and process id, so it is
##'   unique per job. The caller removes it, e.g. with \code{on.exit()}.
##' @param target the final reference directory or file.
##' @return the path of the build directory.
##' @keywords internal
ez_new_build_dir <- function(target) {
  buildDir <- paste0(
    sub("/+$", "", target), ".build_", Sys.info()[["nodename"]], "_",
    Sys.getpid()
  )
  unlink(buildDir, recursive = TRUE)
  dir.create(buildDir, recursive = TRUE)
  buildDir
}

##' @title Move the content of a build directory into place
##' @description Moves every file and directory of \code{fromDir} into
##'   \code{toDir} with \code{file.rename()}, replacing entries of the same
##'   name. The entries named in \code{last} are moved at the very end, in the
##'   given order: use the file whose presence marks a complete reference, so
##'   that a complete-looking reference is never half-written.
##' @param fromDir the build directory.
##' @param toDir the final directory; created if missing.
##' @param last names of entries to move last.
##' @return \code{toDir}, invisibly.
##' @keywords internal
ez_move_into_place <- function(fromDir, toDir, last = character(0)) {
  dir.create(toDir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(fromDir, all.files = TRUE, no.. = TRUE)
  entries <- c(setdiff(entries, last), last[last %in% entries])
  for (entry in entries) {
    target <- file.path(toDir, entry)
    ## file.rename() replaces files and links, but not non-empty directories
    if (dir.exists(target)) {
      unlink(target, recursive = TRUE)
    }
    if (!file.rename(file.path(fromDir, entry), target)) {
      stop("could not move ", file.path(fromDir, entry), " to ", target)
    }
  }
  invisible(toDir)
}

##' @title Describe who holds a reference build lock
##' @param lockDir path of the lock directory
##' @return a single character string
##' @keywords internal
ez_build_lock_owner <- function(lockDir) {
  ownerFile <- file.path(lockDir, "owner.txt")
  if (dir.exists(lockDir) && file.exists(ownerFile)) {
    return(paste(readLines(ownerFile, warn = FALSE), collapse = " "))
  }
  if (file.exists(lockDir) && !dir.exists(lockDir)) {
    ## lock file written by an older ezRun version
    return(paste(
      "lock file of an older ezRun version, modified",
      format(file.mtime(lockDir), "%Y-%m-%d %H:%M:%S")
    ))
  }
  "unknown"
}
