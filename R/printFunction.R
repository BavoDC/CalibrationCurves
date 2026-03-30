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
  cat(sprintf("A %s%% confidence interval is given for the statistics.\n\n",
              (1 - x$alpha) * 100))

  cat("Calibration performance:\n")
  cat("------------------------\n\n")
  if (!is.null(x$stats$Calibration$InTheLarge)) {
    cat("  Calibration in the large (O/E ratio):\n")
    print(x$stats$Calibration$InTheLarge, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Slope)) {
    cat("  Calibration slope:\n")
    print(x$stats$Calibration$Slope, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Statistics)) {
    cat("  Calibration statistics (ICI, E50, E90, Emax):\n")
    print(x$stats$Calibration$Statistics, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$BrierScore)) {
    bs        <- x$stats$Calibration$BrierScore
    model_row <- bs[!tolower(as.character(bs$model)) %in% c("null model", "null", "reference"), , drop = FALSE]
    if (nrow(model_row) > 0 && "Brier" %in% names(model_row)) {
      cat(sprintf("  Brier score: %.4f", model_row$Brier[1]))
      if ("IPA" %in% names(model_row) && !is.na(model_row$IPA[1]))
        cat(sprintf("    IPA (scaled Brier): %.4f", model_row$IPA[1]))
      cat("\n\n")
    }
  }

  cat("Discrimination performance:\n")
  cat("---------------------------\n\n")
  if (!is.null(x$stats$Concordance)) {
    cat("  Concordance statistic:\n")
    print(x$stats$Concordance, ...)
    cat("\n")
  }
  if (!is.null(x$stats$TimeDependentAUC)) {
    cat("  Time-dependent AUC:\n")
    print(x$stats$TimeDependentAUC, ...)
    cat("\n")
  }
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
  cat(sprintf("Competing-risks calibration (cause %s, time horizon %s)\n",
              x$cause, x$timeHorizon))
  cat(sprintf("A %s%% confidence interval is given for the statistics.\n\n",
              (1 - x$alpha) * 100))

  cat("Calibration performance:\n")
  cat("------------------------\n\n")
  if (!is.null(x$stats$Calibration$InTheLarge)) {
    cat("  Calibration in the large (O/E ratio):\n")
    print(x$stats$Calibration$InTheLarge, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Slope)) {
    cat("  Calibration slope:\n")
    print(x$stats$Calibration$Slope, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Intercept)) {
    cat("  Calibration intercept:\n")
    print(x$stats$Calibration$Intercept, ...)
    cat("\n")
  }
  if (!is.null(x$stats$Calibration$Statistics)) {
    cat("  Calibration statistics (ICI, E50, E90, Emax):\n")
    print(x$stats$Calibration$Statistics, ...)
    cat("\n")
  }
  if (!is.null(x$stats$BrierScore)) {
    bs        <- x$stats$BrierScore
    model_row <- bs[!tolower(as.character(bs$model)) %in% c("null model", "null", "reference"), , drop = FALSE]
    if (nrow(model_row) > 0 && "Brier" %in% names(model_row)) {
      cat(sprintf("  Brier score: %.4f", model_row$Brier[1]))
      if ("IPA" %in% names(model_row) && !is.na(model_row$IPA[1]))
        cat(sprintf("    IPA (scaled Brier): %.4f", model_row$IPA[1]))
      cat("\n\n")
    }
  }

  cat("Discrimination performance:\n")
  cat("---------------------------\n\n")
  if (!is.null(x$stats$AUC)) {
    cat("  AUC:\n")
    print(x$stats$AUC, ...)
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
  cat(sprintf("Multiclass calibration (%s outcome)\n", x$type))
  cat(sprintf("A %s%% confidence interval is given for the statistics.\n\n",
              x$cl.level * 100))

  cat("Calibration performance:\n")
  cat("------------------------\n\n")
  if (!is.null(x$Calibration$Slopes)) {
    cat("  Calibration slopes (one-vs-rest):\n")
    print(x$Calibration$Slopes, ...)
    cat("\n")
  }
  if (!is.null(x$Calibration$Intercepts)) {
    cat("  Calibration intercepts (one-vs-rest):\n")
    print(x$Calibration$Intercepts, ...)
    cat("\n")
  }
  if (!is.null(x$stats$MSEC)) {
    cat("  Calibration statistics per category (ICI, E50, E90, Emax):\n")
    for (nm in names(x$stats$MSEC)) {
      cat(sprintf("    %s: ", nm))
      print(x$stats$MSEC[[nm]], ...)
    }
    cat("\n")
  }
  if (!is.null(x$stats$BrierScore))
    cat(sprintf("  Brier score (overall): %.4f\n", x$stats$BrierScore["Overall"]))
  if (!is.null(x$stats$LogLoss))
    cat(sprintf("  Log-loss:              %.4f\n", x$stats$LogLoss))
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
  cat("Calibration performance:\n")
  cat("------------------------\n\n")
  tab <- rbind(
    "Brier score"      = c(x$before$Brier,      x$after$Brier),
    "Log-loss"         = c(x$before$LogLoss,     x$after$LogLoss),
    "Cal. slope"       = c(x$before$Slope,       x$after$Slope),
    "Cal. intercept"   = c(x$before$Intercept,   x$after$Intercept)
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

