#' Infix operator to run background jobs
#'
#' This infix operator can be used to create a background job in RStudio/Posit and, once completed, the value of rhs is
#' assigned to lhs.
#'
#' @param lhs the object that the rhs value is assigned to
#' @param rhs the value you want to assign to lhs
#'
#' @return prints the ID of the background job in the console and, once completed, the value of lhs is assigned to rhs
#'
#' @examples
#' # Can only be executed in Rstudio
#' \dontrun{x %<=% rnorm(1e7)}
`%<=%` <- function(lhs, rhs) {
  Pkgs = names(sessionInfo()$otherPkgs)
  lhs  = as.character(enquote(substitute(lhs))[2])
  rhs  = as.character(enquote(substitute(rhs))[2])
  Job  = c(if(is.null(Pkgs)) NULL else paste0("LibraryM(", paste0(Pkgs, collapse = ", "), ")"),
           paste0(lhs, " <- ", rhs))
  tmpR = tempfile(fileext = ".R")
  writeLines(Job, tmpR)
  if(!file.exists(tmpR))
    stop("Temporary R script not created")
  jobRunScript(tmpR, exportEnv = "R_GlobalEnv", importEnv = TRUE)
}

#' Infix operator to run background jobs
#'
#' This infix operator can be used to create a background job for a block of code in RStudio/Posit and, once completed,
#' all objects created in the block of code are imported into the global environment.
#'
#' @param lhs not used, see details and examples
#' @param rhs the block of code that you want to run
#'
#' @return prints the ID of the background job in the console and, once completed, the objects created in the block of code
#' are imported into the global environment
#'
#' @details
#' You can use this infix operator in two different ways. Either you set the left-hand side to \code{NULL} or you use the syntax
#' `` `%{}%` ``\code{({BlockOfCode})}
#'
#' @examples
#' # Can only be executed in Rstudio
#' \dontrun{
#' NULL %{}% {
#'  x = rnorm(1e7)
#'  y = rnorm(1e7)
#' }
#' `%{}%` ({
#'  x = rnorm(1e7)
#'  y = rnorm(1e7)
#' })
#' }
`%{}%` <- function(lhs, rhs) {
  Pkgs = names(sessionInfo()$otherPkgs)
  Job  = c(if(is.null(Pkgs)) NULL else paste0("LibraryM(", paste0(Pkgs, collapse = ", "), ")"),
           as.character(enquote(substitute(rhs))[2]))
  tmpR = tempfile(fileext = ".R")
  writeLines(Job, tmpR)
  if(!file.exists(tmpR))
    stop("Temporary R script not created")
  jobRunScript(tmpR, exportEnv = "R_GlobalEnv", importEnv = TRUE)
}


#' Function to load multiple packages at once
#'
#'
#' @param ... the packages that you want to load
#'
#' @return invisible \code{NULL}
#'
#' @examples
#' LibraryM(CalibrationCurves)
LibraryM <- function(...) {
  libs = as.list(substitute(list(...)))[-1L]
  if(!is.character(libs)) {
    if(length(libs[[1]]) > 1)
      libs = sapply(substitute(libs)[-1], as.character)
    else
      libs = as.character(substitute(libs))
  }
  lapply(libs, library, character.only = TRUE)
  invisible(NULL)
}
