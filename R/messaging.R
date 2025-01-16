###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Modified version of Sys.time()
##' @description Displays the system time in a different format.
##' @return Returns the system time.
##' @template roxygen-template
##' @examples 
##' ezTime()
ezTime = function(){
  format(Sys.time(), "%Y-%m-%d--%H-%M-%S")
}

##' @title Writes start time of job
##' @description Writes the starting date and time of \code{jobName} into a separate file.
##' @param jobName a character representing the job name.
##' @return Returns a list containing \code{jobName} and the running time of R.
##' @template roxygen-template
##' @template connection-template
##' @examples
##' con = file("example.txt", "w")
##' ezJobStart("a job", con=con)
##' close(con)
ezJobStart = function(jobName, con=stderr()){
  ezWrite(ezTime(), " ", jobName, " start", con=con)
  job = list(name=jobName, start=proc.time())
  return(job)
}

##' @title Writes elapsed time of job
##' @description Writes the current date and time, the status and how much time elapsed since the jobs' start.
##' @param job a list with fields \code{$name} and \code{$start}.
##' @param status a character describing the status of the job.
##' @param x corresponds to \code{job$start}.
##' @template roxygen-template
##' @template connection-template
##' @examples
##' con = file("example.txt", "w")
##' job = ezJobStart("a job",con=con)
##' ezWriteElapsed(job,con=con)
##' close(con)
ezWriteElapsed = function(job, status="done", con=stderr()){
  ezWrite(ezTime(), " ", job$name, " ", status, ": ", getElapsed(job$start), con=con)
  return(invisible(list(name=job$name, start=proc.time())))
}

##' @describeIn ezWriteElapsed Gets the elapsed time in minutes since a reference time. 
getElapsed = function(x){
  paste(signif((proc.time() - x)[ "elapsed"]/60, digits=4), "min")
}


ezSessionInfo <- function(){
  ezRunDetails = library(help = ezRun)
  RemoteSha <- sub('.*\\s+','',ezRunDetails$info[[1]][grep('RemoteSha', ezRunDetails$info[[1]])])
  githubUrl <- file.path('https://github.com/uzh/ezRun/tree', RemoteSha)
  cat('ezRun tag:', RemoteSha, '\n')
  cat('ezRun github link:', githubUrl, '\n \n')
  
  print(sessionInfo())
}



##' @title Mails information, the nodename and working directory
##' @description Mails the nodename and working directory to a specified email-address. A text message and the subject can also be defined.
##' @param text a message to send.
##' @param subject the subject of the email.
##' @param to the recipient of the email.
##' @template roxygen-template
##' @examples 
##' \dontrun{
##' ezMail(to="example@@example.com")
##' }
ezMail = function(text="done", subject="r-mail", to=""){
  stopifnot(to != "")
  text = paste(c(text, 
                 paste("Host:", Sys.info()["nodename"]),
                 paste("Dir:", getwd())
  ), collapse="\n")
  res = system(paste("ssh fgcz-c-044 \" mail -s '", subject, "'", to, "\""), input=text)
  stopifnot(res == 0)
}

##' @title Is the email-address valid?
##' @description Checks, whether \code{addressString} is a valid email-address.
##' @param addressString the character to check.
##' @return Returns TRUE if \code{addressString} contains an @@.
##' @template roxygen-template
##' @examples 
##' ezValidMail("example@@example.com")
ezValidMail = function(addressString){
  return(!is.null(addressString) && grepl("@", addressString))
}

.waitUntilFileExists = function(myFile, maxWaitSeconds=300, interval=30){
  waitSeconds = 0
  while(!file.exists(myFile) && waitSeconds < maxWaitSeconds){
    Sys.sleep( interval)
    waitSeconds = waitSeconds + interval
  }
  if (!file.exists(myFile)){
    stop(paste("waited", maxWaitSeconds, "but file is still missing:", myFile))
  }
  return()
}

##' @title Logs time and base directory together with a parameter
##' @description Logs a time with an additional context provided in \code{message} and the base directory. Also displays a defined method name and the parameters name and Id.
##' @param methodName the name of the method to display.
##' @param param the parameter to display.
##' @param message a character to start the expression and provide context information.
##' @return Returns a message about \code{param} with a context.
##' @template roxygen-template
##' @examples 
##' param = ezParam()
##' logMessage("a method", param, "Starting")
##' logMessage("a method", param, "Finished")
logMessage = function(methodName, param, message){
  msg = paste(message, methodName, param$name, param$projectId, 
              basename(getwd()), format(Sys.time(), "%Y-%m-%d %H:%M:%S\n"), sep="\t")
  cat(msg)
}
