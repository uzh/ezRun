##' @title The R5 class representing a runnable app using \code{<%= method %>()}
##' @field runMethod the function that will be executed in the \code{run} method.
##' @field name the name of the application.
##' @field appDefaults the defaults to run the application with.
##' @section Functions:
##' \itemize{
##'   \item{\code{<%= method %> }}:
##'   {The function to run this application.}
##' }
##' @param input a list, file path or an object of the class EzDataset containing the input.
##' @param output a list, file path or an object of the class EzDataset containing the output information.
##' @param param a list of parameters to customize the application run.
##' @author Rehrauer, Hubert
##' @author Schmid, Peter
##' @seealso \code{\link{EzApp}}
##' @seealso \code{\link{EzDataset}}
