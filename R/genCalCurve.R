#' Calibration performance using the generalized calibration framework
#'
#' Function to assess the calibration performance of a prediction model where the outcome's distribution is a member of the exponential family (De Cock Campo, 2023).
#' The function plots the generalized calibration curve and computes the generalized calibration slope and intercept.
#'
#'
#' @param y a vector with the values for the response variable
#' @param yHat a vector with the predicted values
#' @param family a description of the type of distribution and link function in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' (See family for details of family functions.)
#' @param plot logical, indicating if a plot should be made or not.
#' @param Smooth logical, indicating if the flexible calibration curve should be estimated.
#' @param GLMCal logical, indicating if the GLM calibration curve has to be estimated.
#' @param lwdIdeal the line width of the ideal line.
#' @param colIdeal the color of the ideal line.
#' @param ltyIdeal the line type of the ideal line.
#' @param lwdSmooth the line width of the flexible calibration curve.
#' @param colSmooth the color of the flexible calibration curve.
#' @param ltySmooth the line type of the flexible calibration curve.
#' @param argzSmooth arguments passed to \code{\link{loess}}.
#' @param lwdGLMCal the line width of the GLM calibration curve.
#' @param colGLMCal the color of the GLM calibration curve.
#' @param ltyGLMCal the line type of the GLM calibration curve.
#' @param AddStats logical, indicating whether to add the values of the generalized calibration slope and intercept to the plot.
#' @param Digits the number of digits of the generalized calibration slope and intercept.
#' @param cexStats the font size of the statistics shown on the plot.
#' @param lwdLeg the line width in the legend.
#' @param Legend logical, indicating whether the legend has to be added.
#' @param legendPos the position of the legend on the plot.
#' @param xLim,yLim numeric vectors of length 2, giving the x and y coordinates ranges (see \code{\link{plot.window}})
#' @param posStats numeric vector of length 2, specifying the x and y coordinates of the statistics (generalized calibration curve and intercept) printed on the plot. Default is \code{NULL}
#' which places the statistics in the top left corner of the plot.
#' @param confLimitsSmooth character vector to indicate if and how the confidence limits for the flexible calibration curve have to be computed. \code{"none"} omits the confidence limits,
#' \code{"bootstrap"} uses 2000 bootstrap samples to calculate the 95\% confidence limits and \code{"pointwise"} uses the pointwise confidence limits.
#' @param confLevel the confidence level for the calculation of the pointwise confidence limits of the flexible calibration curve.
#' @param Title the title of the plot
#' @param xlab x-axis label, default is \code{"Predicted value"}.
#' @param ylab y-axis label, default is \code{"Empirical average"}.
#' @param EmpiricalDistribution logical, indicating if the empirical distribution of the predicted values has to be added to the bottom of the plot.
#' @param length.seg controls the length of the histogram lines. Default is \code{1}.
#' @param ... arguments to be passed to \code{\link{plot}}, see \code{\link{par}}
#'
#' @return An object of type \code{GeneralizedCalibrationCurve} with the following slots:
#' @return \item{call}{the matched call.}
#' @return \item{ggPlot}{the ggplot object.}
#' @return \item{stats}{a vector containing performance measures of calibration.}
#' @return \item{cl.level}{the confidence level used.}
#' @return \item{Calibration}{contains the calibration intercept and slope, together with their confidence intervals.}
#' @return \item{Cindex}{the value of the c-statistic, together with its confidence interval.}
#' @return \item{warningMessages}{if any, the warning messages that were printed while running the function.}
#' @return \item{CalibrationCurves}{The coordinates for plotting the calibration curves. }
#' @export
#'
#' @references De Cock Campo, B. (2023). Towards reliable predictive analytics: a generalized calibration framework. arXiv:2309.08559, available at \url{https://arxiv.org/abs/2309.08559}.
#'
#' @examples
#' library(CalibrationCurves)
#' library(mgcv)
#' data("poissontraindata")
#' data("poissontestdata")
#'
#' glmFit = glm(Y ~ ., data = poissontraindata, family = poisson)
#'
#' # Example of a well calibrated poisson prediction model
#' yOOS = poissontestdata$Y
#' yHat = predict(glmFit, newdata = poissontestdata, type = "response")
#' genCalCurve(yOOS, yHat, family = "poisson", plot = TRUE)
#'
#' # Example of an overfit poisson prediction model
#' gamFit = gam(Y ~ x1 + x3 + x1:x3 + s(x5), data = poissontraindata, family = poisson)
#' yHat = as.vector(predict(gamFit, newdata = poissontestdata, type = "response"))
#' genCalCurve(yOOS, yHat, family = "poisson", plot = TRUE)
#'
#' # Example of an underfit poisson prediction model
#' glmFit = glm(Y ~ x2, data = poissontraindata, family = poisson)
#' yOOS = poissontestdata$Y
#' yHat = predict(glmFit, newdata = poissontestdata, type = "response")
#' genCalCurve(yOOS, yHat, family = "poisson", plot = TRUE)
genCalCurve <- function(y, yHat, family, plot = TRUE, Smooth = FALSE, GLMCal = TRUE, lwdIdeal = 2, colIdeal = "gray", ltyIdeal = 1,
                        lwdSmooth = 1, colSmooth = "blue", ltySmooth = 1, argzSmooth = alist(degree = 2),
                        lwdGLMCal = 1, colGLMCal = "red", ltyGLMCal = 1,
                        AddStats = T, Digits = 3, cexStats = 1, lwdLeg = 1.5, Legend = TRUE, legendPos = "bottomright",
                        xLim = NULL, yLim = NULL, posStats = NULL,
                        confLimitsSmooth = c("none", "bootstrap", "pointwise"), confLevel = 0.95,
                        Title = "Calibration plot",
                        xlab = "Predicted value", ylab = "Empirical average",
                        EmpiricalDistribution = TRUE, length.seg = 1, ...) {
  bootSamples <- BT.samples
  call   = match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if(!is.null(posStats)) {
    if(!is.vector(posStats))
      stop("Has to be of type vector.")
    if(length(posStats) != 2)
      stop("Length of the vector has to be equal to 2.")
  }
  confLimitsSmooth = match.arg(confLimitsSmooth)
  a     = 1 - confLevel
  wmess = NULL

  Eta    = family$linkfun(yHat)
  ClInt = tryCatch(
    glm(
      y ~ offset(Eta),
      family = family,
      control = glm.control(maxit = 1e2)
    ),
    error = function(e)
      T,
    warning = function(w)
      T
  )
  if(is.logical(ClInt)) {
    # https://stackoverflow.com/questions/8212063/glm-starting-values-not-accepted-log-link
    ClInt = glm(I(y + .Machine$double.eps) ~ offset(Eta), family = family, control = glm.control(maxit = 1e2))
  }
  ClSl  =
    tryCatch(
      glm(y ~ Eta, family = family, control = glm.control(maxit = 1e2)),
      error = function(e)
        T,
      warning = function(w)
        T
    )
  if(is.logical(ClSl)) {
    lmFit = lm(y ~ Eta)
    ClSl  = glm(y ~ Eta, family = family, control = glm.control(maxit = 1e2), start = coef(lmFit))
  }
  ClSl2 = tryCatch(
    glm(y ~ Eta - 1, family = family),
    error = function(e)
      T,
    warning = function(w)
      T
  )
  if(is.logical(ClSl2)) {
    lmFit = lm(y ~ Eta - 1)
    ClSl2  = glm(y ~ Eta - 1, family = family, control = glm.control(maxit = 1e2), start = coef(lmFit))
  }
  CalibrStats = c("Calibration intercept" = unname(coef(ClInt)), "Calibration slope" = unname(coef(ClSl)[2]))

  ClIntCL = confint(ClInt, level = confLevel)
  ClSlCL  = confint(ClSl, level = confLevel)[2, ]

  y    = y[order(yHat)]
  Eta  = Eta[order(yHat)]
  yHat = sort(yHat)
  calCurves = list()

  if(GLMCal) {
    glmFit = glm(y ~ Eta, family = family)
    rangeY = range(glmFit$fitted)
    glmCal = data.frame(x = yHat, y = fitted(glmFit))
    calCurves$GLMCalibration = glmCal
  }
  if(Smooth) {
    argzSmooth$formula = y ~ yHat
    SmFit <- Sm <- do.call("loess", argzSmooth)
    Sm     = data.frame(Sm$x, Sm$fitted)
    rangeY = if(GLMCal) c(min(rangeY, SmFit$fitted), max(rangeY, SmFit$fitted)) else range(SmFit$fitted)
    calCurves$FlexibleCalibration = Sm
  }
  xLim = if(is.null(xLim)) range(yHat) else xLim
  yLim = if(is.null(yLim)) c(min(c(xLim, rangeY)), max(c(xLim, rangeY))) else yLim
  yLim[1] =
    if(yLim[1] <= 0.5) {
      0 - 0.1 * diff(range(yLim))
    } else {
      yLim[1] * 0.9
    }
  if(plot) {
    plot(mean(xLim), mean(yLim), col = "white", pch = 1, xlab = xlab, ylab = ylab,
         xlim = xLim, ylim = yLim, main = Title, ...)
    clip(min(c(xLim, yHat)), max(c(xLim, yHat)), min(c(yLim, rangeY)), max(c(yLim, rangeY)))
    abline(0, 1, col = colIdeal, lwd = lwdIdeal, lty = ltyIdeal)
  }

  labLeg = "Ideal"
  colLeg = colIdeal
  ltyLeg = ltyIdeal
  lwdLeg = lwdIdeal

  if(Smooth) {
    if(plot)
      lines(Sm, lty = ltySmooth, lwd = lwdSmooth, col = colSmooth)
    if(confLimitsSmooth != "none") {
      if(confLimitsSmooth == "bootstrap") {
        yHatGrid = seq(min(yHat), max(yHat), length = 200)
        resBoot  = replicate(2000, bootSamples(y, yHat, yHatGrid))
        clBoot   = apply(resBoot, 1, quantile, c(0.025, 0.975))
        dfCL     = data.frame(x = yHatGrid, ymin = clBoot[1, ], ymax = clBoot[2, ])
        rownames(dfCL) = NULL
      } else {
        cl.loess = predict(SmFit, type = "fitted", se = TRUE)
        dfCL     = data.frame(x = yHat, ymin = with(cl.loess, fit - qnorm(1 - a / 2) * se.fit),
                              ymax = with(cl.loess, fit + qnorm(1 - a / 2) * se.fit))
      }
      if(plot)
        with(dfCL,
             polygon(
               x = c(x, rev(x)),
               y = c(ymax,
                     rev(ymin)),
               col = rgb(177, 177, 177, 177, maxColorValue = 255),
               border = NA
             )
      )
    }
    labLeg = c(labLeg, "Flexible calibration")
    colLeg = c(colLeg, colSmooth)
    ltyLeg = c(ltyLeg, ltySmooth)
    lwdLeg = c(lwdLeg, lwdSmooth)
  }
  if(GLMCal) {
    if(plot)
      lines(yHat, fitted(glmFit), lty = ltyGLMCal, lwd = lwdGLMCal, col = colGLMCal)
    labLeg = c(labLeg, "GLM calibration")
    colLeg = c(colLeg, colGLMCal)
    ltyLeg = c(ltyLeg, ltyGLMCal)
    lwdLeg = c(lwdLeg, lwdGLMCal)
  }
  if(plot)
    do.call("clip", as.list(par()$usr))
  if(EmpiricalDistribution) {
    x     <- yHat
    bins  <- seq(min(x), max(x), length = 101)
    f0	  <- table(cut(x, bins))
    bins  <- (bins[-101])
    maxf  <- max(f0)
    f0	  <- (0.1 * f0) / maxf

    if(plot) {
      segments(bins, yLim[1], bins, yLim[1] + length.seg * f0)
      lines(c(min(bins) - 0.01, max(bins) + 0.01), c(yLim[1], yLim[1]))
    }
  }
  if(AddStats & plot) {
    StatsPlot = paste0('Calibration\n',
                       '...intercept: ',
                       sprintf(paste0("%.", Digits, "f"), CalibrStats[1]), '\n',
                       '...slope: ',
                       sprintf(paste0("%.", Digits, "f"), CalibrStats[2]), '\n')
    if(is.null(posStats))
      text(xLim[1], xLim[2] * 0.85, StatsPlot, pos = 4, cex = cexStats)
    else
      text(posStats[1], posStats[2], StatsPlot, pos = 4, cex = cexStats)
  }
  if(plot) {
    if(Legend)
      if(is.character(legendPos))
        legend(legendPos, legend = labLeg, col = colLeg, lty = ltyLeg, bty = "n", lwd = lwdLeg)
    else
      legend(legendPos[1], legendPos[2], legend = labLeg, col = colLeg, lty = ltyLeg, bty = "n", lwd = lwdLeg)
  }
  Results =
    structure(
      list(
        call   = call,
        stats  = CalibrStats,
        cl.level = confLevel,
        Calibration = list(
          Intercept = c("Point estimate" = CalibrStats[1],
                        "Lower confidence limit" = ClIntCL[1],
                        "Upper confidence limit" = ClIntCL[2]),
          Slope = c("Point estimate" = CalibrStats[2],
                    "Lower confidence limit" = ClSlCL[1],
                    "Upper confidence limit" = ClSlCL[2])
        ),
        warningMessages = wmess,
        CalibrationCurves = calCurves
      ), class = "GeneralizedCalibrationCurve"
    )
  return(Results)
  return(CalibrStats)
}
