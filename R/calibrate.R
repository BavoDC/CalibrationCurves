#' Automatic calibration assessment for multiple model types
#'
#' A convenience wrapper that detects the class of a fitted prediction model and
#' routes to the appropriate calibration function: \code{\link{val.prob.ci.2}} /
#' \code{\link{valProbggplot}} for binary (logistic) models, \code{\link{genCalCurve}}
#' for generalized models, \code{\link{valProbSurvival}} for Cox models, and
#' \code{\link{valProbCompRisks}} for competing-risks models (CSC or FGR).
#' For ordinal/multinomial models, \code{\link{valProbMulticat}} is called.
#'
#' @param fit a fitted prediction model.  Supported classes:
#'   \itemize{
#'     \item \code{glm} with \code{family = binomial} — binary calibration;
#'     \item \code{glm} with other families — generalized calibration via
#'       \code{\link{genCalCurve}};
#'     \item \code{coxph} — time-to-event calibration via
#'       \code{\link{valProbSurvival}};
#'     \item \code{CSC} or \code{FGR} (from \pkg{riskRegression}) — competing-risks
#'       calibration via \code{\link{valProbCompRisks}};
#'     \item \code{polr} (from \pkg{MASS}) or \code{multinom} (from \pkg{nnet}) —
#'       multiclass calibration via \code{\link{valProbMulticat}}.
#'   }
#' @param newdata data frame containing the validation data.
#' @param y character string or vector.  For binary/generalized models, the
#'   observed outcome vector or column name in \code{newdata}.  For survival
#'   and competing-risks models this is ignored (extracted from the formula).
#' @param timeHorizon numeric, time horizon for survival / competing-risks
#'   calibration. Default \code{5}.
#' @param cause integer, cause of interest for competing-risks models. Default
#'   \code{1}.
#' @param plotCal character, passed down to the called function.  Default
#'   \code{"ggplot"}.
#' @param ... additional arguments passed to the underlying calibration function.
#'
#' @return The return value of the underlying calibration function.
#'
#' @seealso \code{\link{val.prob.ci.2}}, \code{\link{valProbggplot}},
#'   \code{\link{genCalCurve}}, \code{\link{valProbSurvival}},
#'   \code{\link{valProbCompRisks}}, \code{\link{valProbMulticat}}
#'
#' @examples
#' library(CalibrationCurves)
#' data("traindata")
#' data("testdata")
#' glmFit <- glm(y ~ ., data = traindata, family = binomial)
#' # The following detects that 'glmFit' is a binary logistic model and
#' # calls valProbggplot automatically.
#' \dontrun{
#' calibrationCurve(glmFit, newdata = testdata, y = "y")
#' }
#'
#' \dontrun{
#' library(survival)
#' data(trainDataSurvival)
#' data(testDataSurvival)
#' sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
#'               data = trainDataSurvival, x = TRUE, y = TRUE)
#' # Automatically calls valProbSurvival
#' calibrationCurve(sFit, newdata = testDataSurvival, timeHorizon = 5)
#' }
#' @export
calibrationCurve <- function(fit,
                      newdata     = NULL,
                      y           = NULL,
                      timeHorizon = 5,
                      cause       = 1,
                      plotCal     = "ggplot",
                      ...) {

  # ---- Binary / GLM --------------------------------------------------------
  if (inherits(fit, "glm")) {
    fam <- family(fit)$family
    if (fam == "binomial") {
      pHat <- predict(fit, newdata = newdata, type = "response")
      yObs <- .extract_outcome(y, newdata, fit)
      return(valProbggplot(pHat, yObs, ...))
    } else {
      yObs <- .extract_outcome(y, newdata, fit)
      yHat <- predict(fit, newdata = newdata, type = "response")
      return(genCalCurve(yObs, yHat, family = family(fit), ...))
    }
  }

  # ---- Cox PH --------------------------------------------------------------
  if (inherits(fit, "coxph")) {
    if (is.null(newdata))
      stop("'newdata' must be supplied for survival models.")
    return(valProbSurvival(fit, valdata = newdata,
                           timeHorizon = timeHorizon,
                           plotCal = plotCal, ...))
  }

  # ---- Competing risks: CSC / FGR ------------------------------------------
  compRiskClasses <- c("CSC", "FGR", "selectCox", "riskRegression")
  if (any(sapply(compRiskClasses, function(cl) inherits(fit, cl)))) {
    if (is.null(newdata))
      stop("'newdata' must be supplied for competing-risks models.")
    return(valProbCompRisks(fit, valdata = newdata,
                            cause = cause,
                            timeHorizon = timeHorizon,
                            plotCal = plotCal, ...))
  }

  # ---- Ordinal / multinomial -----------------------------------------------
  if (inherits(fit, "polr") || inherits(fit, "multinom")) {
    if (is.null(newdata))
      stop("'newdata' must be supplied for ordinal/multinomial models.")
    return(valProbMulticat(fit, valdata = newdata,
                           plotCal = plotCal, ...))
  }

  # ---- Fallback: attempt to extract predicted probabilities ----------------
  pHat <- tryCatch(
    predict(fit, newdata = newdata, type = "response"),
    error = function(e) NULL
  )
  if (is.null(pHat))
    pHat <- tryCatch(
      predict(fit, newdata = newdata),
      error = function(e) NULL
    )

  if (is.null(pHat))
    stop(paste("Cannot determine calibration function for model class:",
               paste(class(fit), collapse = ", ")))

  yObs <- .extract_outcome(y, newdata, fit)
  if (is.null(yObs))
    stop("'y' or 'newdata' with the outcome column must be provided.")

  # Treat as binary if values in [0,1]
  if (all(pHat >= 0 & pHat <= 1)) {
    message("Model class not recognised; assuming binary outcome. ",
            "Call val.prob.ci.2/valProbggplot directly for full control.")
    return(valProbggplot(pHat, yObs, ...))
  }

  stop(paste("Unsupported model class:", paste(class(fit), collapse = ", ")))
}


# Internal helper: extract the outcome variable from newdata or the user
.extract_outcome <- function(y, newdata, fit) {
  if (is.null(y)) {
    # Try to infer response name from the model formula
    respName <- tryCatch(
      as.character(formula(fit)[[2]]),
      error = function(e) NULL
    )
    if (!is.null(respName) && !is.null(newdata) && respName %in% names(newdata))
      return(newdata[[respName]])
    return(NULL)
  }
  if (is.character(y) && length(y) == 1) {
    if (!is.null(newdata) && y %in% names(newdata))
      return(newdata[[y]])
    stop(paste0("Column '", y, "' not found in 'newdata'."))
  }
  return(y)
}
