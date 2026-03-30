#' Print function for a CalibrationCurve object
#'
#' Prints the call, confidence level and values for the performance measures.
#'
#' @param x an object of type CalibrationCurve, resulting from \code{\link{val.prob.ci.2}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{val.prob.ci.2}}
#' @return The original \code{CalibrationCurve} object is returned.
#' @export
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
  if(!is.null(x$warningMessages))
    for(w in x$warningMessages)
      warning(paste0(w, "\n"), immediate. = TRUE)
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
#' @export
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
  if(!is.null(x$warningMessages))
    for(w in x$warningMessages)
      warning(paste0(w, "\n"), immediate. = TRUE)
  invisible(x)
}


#' Print function for a GeneralizedCalibrationCurve object
#'
#' Prints the call, confidence level and values for the performance measures.
#'
#' @param x an object of type GeneralizedCalibrationCurve, resulting from \code{\link{genCalCurve}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{genCalCurve}}
#' @return The original \code{GeneralizedCalibrationCurve} object is returned.
#' @export
print.GeneralizedCalibrationCurve <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(
    paste(
      "A ",
      x$cl.level * 100,
      "% confidence interval is given for the calibration intercept and calibration slope. \n\n",
      sep = ""
    )
  )
  print(x$stats, ...)
  if(!is.null(x$warningMessages))
    for(w in x$warningMessages)
      warning(paste0(w, "\n"), immediate. = TRUE)
  invisible(x)
}


#' Print function for a SurvivalCalibrationCurve object
#'
#' @param x an object of type SurvivalCalibrationCurve, resulting from \code{\link{valProbSurvival}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{valProbSurvival}}
#' @return The original \code{SurvivalCalibrationCurve} object is returned.
#' @export
print.SurvivalCalibrationCurve <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(
    paste(
      "A ",
      (1 - x$alpha) * 100,
      "% confidence interval is given for the statistics. \n\n",
      sep = ""
    )
  )
  cat("Calibration performance:\n")
  cat("------------------------\n\n")
  cat("In the large\n\n")
  print(x$stats$Calibration$InTheLarge, ...)
  cat("\nSlope\n\n")
  print(x$stats$Calibration$Slope, ...)
  cat("\nAdditional statistics\n\n")
  print(x$stats$Calibration$Statistics, ...)
  print(x$stats$Calibration$BrierScore, ...)

  cat("\n\nDiscrimination performance:\n")
  cat("-------------------------------\n\n")
  cat("Concordance statistic\n\n")
  print(x$stats$Concordance, ...)
  cat("\n\nTime-dependent AUC\n\n")
  print(x$stats$TimeDependentAUC)
  invisible(x)
}


#' Print function for a ClusteredCalibrationCurve object
#'
#' Prints the ggplot, call, confidence level and values for the performance measures.
#'
#' @param x an object of type ggplotCalibrationCurve, resulting from \code{\link{valProbggplot}}.
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso \code{\link{valProbggplot}}
#' @return The original \code{ggplotCalibrationCurve} object is returned.
#' @export
print.ClusteredCalibrationCurve <- function(x, ...) {
  print(x$ggPlot)
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(
    paste(
      "A ",
      x$cl.level * 100,
      "% confidence interval is used. \n\n",
      sep = ""
    )
  )
  # print(x$stats, ...)
  # if(!is.null(x$warningMessages))
  #   for(w in x$warningMessages)
  #     warning(paste0(w, "\n"), immediate. = TRUE)
  # invisible(x)
}

#' Print function for a CompRisksCalibrationCurve object
#'
#' @param x an object of class \code{CompRisksCalibrationCurve}, resulting
#'   from \code{\link{valProbCompRisks}}.
#' @param ... arguments passed to \code{\link{print}}
#' @seealso \code{\link{valProbCompRisks}}
#' @return The original object, invisibly.
#' @export
print.CompRisksCalibrationCurve <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(sprintf("Competing-risks calibration (cause %s, time horizon %s)\n\n",
              x$cause, x$timeHorizon))
  if (!is.null(x$stats$Calibration$InTheLarge)) {
    cat("Calibration in the large (O/E ratio):\n")
    print(x$stats$Calibration$InTheLarge, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Slope)) {
    cat("Calibration slope:\n")
    print(x$stats$Calibration$Slope, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Statistics)) {
    cat("Calibration statistics:\n")
    print(x$stats$Calibration$Statistics, ...)
    cat("\n")
  }
  if (!is.null(x$stats$AUC)) {
    cat("AUC:\n")
    print(x$stats$AUC, ...)
    cat("\n")
  }
  if (!is.null(x$stats$BrierScore)) {
    cat("Brier score:\n")
    print(x$stats$BrierScore, ...)
    cat("\n")
  }
  invisible(x)
}

#' Print function for a MulticlassCalibrationCurve object
#'
#' @param x an object of class \code{MulticlassCalibrationCurve}, resulting
#'   from \code{\link{valProbMulticat}}.
#' @param ... arguments passed to \code{\link{print}}
#' @seealso \code{\link{valProbMulticat}}
#' @return The original object, invisibly.
#' @export
print.MulticlassCalibrationCurve <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(sprintf("Multiclass calibration (%s outcome)\n\n", x$type))
  cat("Proper scoring rules:\n")
  cat(sprintf("  Brier score (overall): %.4f\n", x$stats$BrierScore["Overall"]))
  cat(sprintf("  Log-loss:              %.4f\n", x$stats$LogLoss))
  cat("\nCalibration slopes (one-vs-rest):\n")
  print(x$Calibration$Slopes, ...)
  cat("\nCalibration intercepts (one-vs-rest):\n")
  print(x$Calibration$Intercepts, ...)
  if (!is.null(x$stats$MSEC)) {
    cat("\nCalibration statistics per category (ICI, E50, E90, Emax):\n")
    for (nm in names(x$stats$MSEC)) {
      cat(sprintf("  %s: ", nm))
      print(x$stats$MSEC[[nm]], ...)
    }
  }
  invisible(x)
}

#' Print function for a RecalibratedPredictions object
#'
#' @param x an object of class \code{RecalibratedPredictions}, resulting
#'   from \code{\link{recalibrate}}.
#' @param ... arguments passed to \code{\link{print}}
#' @seealso \code{\link{recalibrate}}
#' @return The original object, invisibly.
#' @export
print.RecalibratedPredictions <- function(x, ...) {
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(sprintf("Recalibration method: %s\n\n", x$method))
  tab <- rbind(
    "Brier score"  = c(x$before$Brier,     x$after$Brier),
    "Log-loss"     = c(x$before$LogLoss,    x$after$LogLoss),
    "Cal. slope"   = c(x$before$Slope,      x$after$Slope),
    "Cal. intercept" = c(x$before$Intercept, x$after$Intercept)
  )
  colnames(tab) <- c("Before", "After")
  print(round(tab, 4), ...)
  if (x$method == "platt") {
    cat("\nPlatt scaling transform: intercept =",
        round(x$transform["intercept"], 4),
        ", slope =", round(x$transform["slope"], 4), "\n")
  }
  invisible(x)
}

