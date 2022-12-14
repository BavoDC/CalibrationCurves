#' Print function for a CalibrationCurve object
#'
#' Prints the call, confidence level and values for the performance measures.
#'
#' @param x an object of type CalibrationCurve, resulting from \code{\link{val.prob.ci.2}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{val.prob.ci.2}}
#' @return The original \code{CalibrationCurve} object is returned.
print.CalibrationCurve <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(
    paste(
      "A ",
      x$cl.level * 100,
      "% confidence interval is given for the calibration intercept, calibration slope and c-statistic. \n\n",
      sep = ""
    )
  )
  print(x$stats, ...)
  invisible(x)
}

#' Print function for a ggplotCalibrationCurve object
#'
#' Prints the ggplot, call, confidence level and values for the performance measures.
#'
#' @param x an object of type ggplotCalibrationCurve, resulting from \code{\link{valProbggplot}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{valProbggplot}}
#' @return The original \code{ggplotCalibrationCurve} object is returned.
print.ggplotCalibrationCurve <- function(x, ...) {
  print(x$ggPlot)
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(
    paste(
      "A ",
      x$cl.level * 100,
      "% confidence interval is given for the calibration intercept, calibration slope and c-statistic. \n\n",
      sep = ""
    )
  )
  print(x$stats, ...)
  invisible(x)
}
