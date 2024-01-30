#' Calibration performance: ggplot version
#'
#' The function \code{valProbggplot} is an adaptation of \code{\link{val.prob}} from Frank Harrell's rms package,
#' \url{https://cran.r-project.org/package=rms}. Hence, the description of some of the functions of \code{valProbggplot}
#' come from the the original \code{\link{val.prob}}.
#' \cr \cr The key feature of \code{valProbggplot} is the generation of logistic and flexible calibration curves and related statistics.
#' When using this code, please cite: Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M.J., Steyerberg,
#' E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology},
#' \bold{74}, pp. 167-176
#'
#' @inheritParams rms::val.prob
#' @param smooth \code{"loess"} generates a flexible calibration curve based on \code{\link{loess}},
#'  \code{"rcs"} generates a calibration curves based on restricted cubic splines (see \code{\link{rcs}} and
#'  \code{\link[Hmisc]{rcspline.plot}}), \code{"none"} suppresses the flexible curve. We recommend to use loess unless N is large,
#'   for example N>5000. Default is \code{"loess"}.
#' @param CL.smooth \code{"fill"} shows pointwise 95\% confidence limits for the flexible calibration curve with a gray
#' area between the lower and upper limits, \code{TRUE} shows pointwise 95\% confidence limits for the flexible calibration curve
#'  with dashed lines, \code{FALSE} suppresses the confidence limits. Default is \code{"fill"}.
#' @param CL.BT \code{TRUE} uses confidence limits based on 2000 bootstrap samples, \code{FALSE} uses closed form confidence limits.
#' Default is \code{FALSE}.
#' @param nr.knots specifies the number of knots for rcs-based calibration curve. The default as well as the highest allowed value is 5.
#' In case the specified number of knots leads to estimation problems, then the number of knots is automatically reduced to the closest
#'  value without estimation problems.
#' @param dostats specifies whether and which performance measures are shown in the figure.
#' \code{TRUE} shows the \code{"abc"} of model performance (Steyerberg et al., 2011): calibration intercept, calibration slope,
#'  and c-statistic. \code{TRUE} is default.
#'  \code{FALSE} suppresses the presentation of statistics in the figure. A \code{c()} list of specific stats shows the specified
#'  stats. The key stats which are also mentioned in this paper are \code{"C (ROC)"} for the c statistic, \code{"Intercept"} for the
#'  calibration intercept, \code{"Slope"} for the calibration slope, and \code{"ECI"} for the estimated calibration index
#'  (Van Hoorde et al, 2015). The full list of possible statistics is taken from \code{\link{val.prob}}
#'  and augmented with the estimated calibration index: \code{"Dxy", "C (ROC)", "R2", "D", "D:Chi-sq", "D:p", "U", "U:Chi-sq",
#'   "U:p", "Q", "Brier", "Intercept", "Slope", "Emax", "Brier scaled", "Eavg", "ECI"}. These statistics are always returned by the function.
#' @param xlim,ylim numeric vectors of length 2, giving the x and y coordinates ranges (see \code{\link{xlim}} and \code{\link{ylim}}).
#' @param size,size.leg controls the font size of the statistics (\code{size}) or plot legend (\code{size.leg}). Default is 3 and 5, respectively.
#' @param roundstats specifies the number of decimals to which the statistics are rounded when shown in the plot. Default is 2.
#' @param d0lab,d1lab controls the labels for events and non-events (i.e. outcome y) for the histograms.
#' Defaults are \code{d1lab="1"} for events and \code{d0lab="0"} for non-events.
#' @param size.d01 controls the size of the labels for events and non-events. Default is 5.
#' @param dist.label controls the horizontal position of the labels for events and non-events. Default is 0.01.
#' @param dist.label2 controls the vertical distance between the labels for events and non-events. Default is 0.03.
#' @param line.bins controls the horizontal (y-axis) position of the histograms. Default is -0.05.
#' @param cutoff puts an arrow at the specified risk cut-off(s). Default is none.
#' @param length.seg controls the length of the histogram lines. Default is \code{0.85}.
#' @param col.ideal controls the color of the ideal line on the plot. Default is \code{"red"}.
#' @param lwd.ideal controls the line width of the ideal line on the plot. Default is \code{1}.
#' @param lty.ideal linetype of the ideal line. Default is \code{1}.
#' @param logistic.cal \code{TRUE} plots the logistic calibration curve, \code{FALSE} suppresses this curve. Default is \code{FALSE}.
#' @param xlab x-axis label, default is \code{"Predicted Probability"}.
#' @param ylab y-axis label, default is \code{"Observed proportion"}.
#' @param statloc the "abc" of model performance (Steyerberg et al., 2011)-calibration intercept, calibration slope,
#' and c statistic-will be added to the plot, using statloc as the upper left corner of a box (default is c(0,.85).
#' You can specify a list or a vector. Use locator(1) for the mouse, \code{FALSE} to suppress statistics. This is plotted after
#' the curve legends.
#' @param  pl \code{TRUE} to plot the calibration curve(s). If \code{FALSE} no calibration curves will be plotted,
#'  but statistics will still be computed and outputted.
#' @param connect.smooth Defaults to \code{TRUE} to draw smoothed estimates using a line. Set to \code{FALSE} to instead use dots at individual estimates
#' @param legendloc if \code{pl=TRUE}, list with components \code{x,y} or vector \code{c(x,y)} for bottom right corner of legend for
#' curves and points. Default is \code{c(.50, .27)} scaled to lim. Use \code{locator(1)} to use the mouse, \code{FALSE} to suppress legend.
#' @param col.log if \code{logistic.cal=TRUE}, the color of the logistic calibration curve. Default is \code{"black"}.
#' @param lty.log if \code{logistic.cal=TRUE}, the linetype of the logistic calibration curve. Default is \code{1}.
#' @param lwd.log if \code{logistic.cal=TRUE}, the line width of the logistic calibration curve. Default is \code{1}.
#' @param col.smooth the color of the flexible calibration curve. Default is \code{"black"}.
#' @param lty.smooth the linetype of the flexible calibration curve. Default is \code{1}.
#' @param lwd.smooth the line width of the flexible calibration curve. Default is \code{1}.
#' @param cl.level if \code{dostats=TRUE}, the confidence level for the calculation of the confidence intervals of the calibration intercept,
#'  calibration slope and c-statistic. Default is \code{0.95}.
#' @param method.ci method to calculate the confidence interval of the c-statistic. The argument is passed to \code{\link{auc.nonpara.mw}} from
#' the auRoc-package and possible methods to compute the confidence interval are \code{"newcombe"}, \code{"pepe"}, \code{"delong"} or
#' \code{"jackknife"}. Bootstrap-based methods are not available. The default method is \code{"pepe"} and here, the confidence interval is
#' the logit-transformation-based confidence interval as documented in Qin and Hotilovac (2008). See \code{\link{auc.nonpara.mw}} for
#' more information on the other methods.
#' @param allowPerfectPredictions Logical, indicates whether perfect predictions (i.e. values of either 0 or 1) are allowed. Default is \code{FALSE}, since we transform
#' the predictions using the logit transformation to calculate the calibration measures. In case of 0 and 1, this results in minus infinity and infinity, respectively. if
#' \code{allowPerfectPredictions = TRUE}, 0 and 1 are replaced by 1e-8 and 1 - 1e-8, respectively.
#' @param argzLoess a list with arguments passed to the \code{\link{loess}} function
#'
#' @return An object of type \code{ggplotCalibrationCurve} with the following slots:
#' @return \item{call}{the matched call.}
#' @return \item{ggPlot}{the ggplot object.}
#' @return \item{stats}{a vector containing performance measures of calibration.}
#' @return \item{cl.level}{the confidence level used.}
#' @return \item{Calibration}{contains the calibration intercept and slope, together with their confidence intervals.}
#' @return \item{Cindex}{the value of the c-statistic, together with its confidence interval.}
#' @return \item{warningMessages}{if any, the warning messages that were printed while running the function.}
#' @return \item{CalibrationCurves}{The coordinates for plotting the calibration curves. }
#'
#' @note In order to make use (of the functions) of the package auRoc, the user needs to install JAGS. However, since our package only uses the
#' \code{auc.nonpara.mw} function which does not depend on the use of JAGS, we therefore copied the code and slightly adjusted it when
#' \code{method="pepe"}.
#'
#' @details When using the predicted probabilities of an uninformative model (i.e. equal probabilities for all observations), the model has no predictive value.
#'  Consequently, where applicable, the value of the performance measure corresponds to the worst possible theoretical value. For the ECI, for example, this equals 1 (Edlinger et al., 2022).
#'
#' @references Edlinger, M, van Smeden, M, Alber, HF, Wanitschek, M, Van Calster, B. (2022). Risk prediction models for discrete ordinal outcomes: Calibration and the impact of the proportional odds assumption. \emph{Statistics in Medicine}, \bold{41( 8)}, pp. 1334– 1360
#' @references Qin, G., & Hotilovac, L. (2008). Comparison of non-parametric confidence intervals for the area under the ROC curve of a continuous-scale diagnostic test. \emph{Statistical Methods in Medical Research}, \bold{17(2)}, pp. 207-21
#' @references Steyerberg, E.W., Van Calster, B., Pencina, M.J. (2011). Performance measures for prediction models and markers : evaluation of predictions and classifications. \emph{Revista Espanola de Cardiologia}, \bold{64(9)}, pp. 788-794
#' @references Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, \bold{74}, pp. 167-176
#' @references Van Hoorde, K., Van Huffel, S., Timmerman, D., Bourne, T., Van Calster, B. (2015). A spline-based tool to assess and visualize the calibration of multiclass risk predictions. \emph{Journal of Biomedical Informatics}, \bold{54}, pp. 283-93
#'
#' @importFrom Hmisc cut2
#' @import ggplot2
#'
#' @examples
#'
#' # Load package
#' library(CalibrationCurves)
#' set.seed(1783)
#'
#' # Simulate training data
#' X      = replicate(4, rnorm(5e2))
#' p0true = binomial()$linkinv(cbind(1, X) %*% c(0.1, 0.5, 1.2, -0.75, 0.8))
#' y      = rbinom(5e2, 1, p0true)
#' Df     = data.frame(y, X)
#'
#' # Fit logistic model
#' FitLog = lrm(y ~ ., Df)
#'
#' # Simulate validation data
#' Xval   = replicate(4, rnorm(5e2))
#' p0true = binomial()$linkinv(cbind(1, Xval) %*% c(0.1, 0.5, 1.2, -0.75, 0.8))
#' yval   = rbinom(5e2, 1, p0true)
#' Pred   = binomial()$linkinv(cbind(1, Xval) %*% coef(FitLog))
#'
#' # Default calibration plot
#' valProbggplot(Pred, yval)
#'
#' # Adding logistic calibration curves and other additional features
#' valProbggplot(Pred, yval, CL.smooth = TRUE, logistic.cal = TRUE, lty.log = 2,
#'  col.log = "red", lwd.log = 1.5)
#'
#' valProbggplot(Pred, yval, CL.smooth = TRUE, logistic.cal = TRUE, lty.log = 9,
#' col.log = "red", lwd.log = 1.5, col.ideal = colors()[10], lwd.ideal = 0.5)

valProbggplot <- function(p, y, logit, group,
                          weights = rep(1, length(y)), normwt = FALSE, pl = TRUE,
                          smooth = c("loess", "rcs", "none"), CL.smooth = "fill",
                          CL.BT = FALSE, lty.smooth = 1, col.smooth = "black", lwd.smooth = 1,
                          nr.knots = 5, logistic.cal = FALSE, lty.log = 1,
                          col.log = "black", lwd.log = 1, xlab = "Predicted probability", ylab = "Observed proportion",
                          xlim = c(-0.02, 1), ylim = c(-0.15, 1), m, g, cuts, emax.lim = c(0, 1),
                          legendloc =  c(0.50 , 0.27), statloc = c(0, .85), dostats = TRUE, cl.level = 0.95, method.ci = "pepe",
                          roundstats = 2, riskdist = "predicted", size = 3, size.leg = 5, connect.group = FALSE, connect.smooth = TRUE,
                          g.group = 4, evaluate = 100, nmin = 0, d0lab = "0", d1lab = "1", size.d01 = 5,
                          dist.label = 0.01, line.bins = -.05, dist.label2 = .04, cutoff, length.seg = 0.85,
                          lty.ideal = 1, col.ideal = "red", lwd.ideal = 1, allowPerfectPredictions = FALSE, argzLoess = alist(degree = 2))
{
  call   = match.call()
  smooth = match.arg(smooth)
  if (smooth == "none")
    smooth <- "F"
  if(!is.logical(connect.smooth))
    stop("Argument connect.smooth has to be of type logical.")
  if (!is.numeric(nr.knots))
    stop("Nr.knots must be numeric.")
  if(nr.knots > 5 | nr.knots < 3)
    stop(paste(
      "Number of knots = ",
      nr.knots,
      sep = "",
      ", only 5 >= nk >=3 is allowed."
    ))
  if (!missing(p))
    if(allowPerfectPredictions & any(!(p > 0 | p < 1)))
      stop("Probabilities can not be > 1 or < 0.")
  else if (any(!(p >= 0 | p <= 1)))
    stop("Probabilities can not be >= 1 or <= 0.")
  if(allowPerfectPredictions) {
    if(all(p %in% 0:1))
      stop("All predicted values are equal to 0 or 1, implying that the underlying process is deterministic. Please check your model or the input.")
    if(any(p %in% c(0, 1))) {
      p = sapply(p, function(x) {
        if(is.na(x) | is.nan(x))
          x
        else if(x == 0)
          1e-8
        else if(x == 1)
          1 - 1e-8
        else
          x
      })
      wmess = paste0("There are predictions with value 0 or 1! These are replaced by values 1e-8 and 1 - 1e-8, respectively. ",
                     "Take this into account when interpreting the performance measures, as these are not calculated with the original values.",
                     "\n\nPlease check your model, as this may be an indication of overfitting. Predictions of 0 or 1 imply that these predicted values are deterministic.\n\n",
                     "We observe this in the following cases:\n - logistic regression: with quasi-complete separation, the coefficients tend to infinity;\n",
                     " - tree-based methods: one of the leaf nodes contains only observations with either 0 or 1;\n",
                     " - neural networks: the weights tend to infinity and this is known as weight/gradient explosion.")
      warning(wmess, immediate. = TRUE)
    } else {
      wmess = NULL
    }
  } else {
    wmess = NULL
  }

  a = 1 - cl.level

  if (missing(p))
    p <- 1 / (1 + exp(-logit))
  else
    logit <- log(p / (1 - p))
  if (!all(y %in% 0:1)) {
    stop("The vector with the binary outcome can only contain the values 0 and 1.")
  }
  if (length(p) != length(y))
    stop("lengths of p or logit and y do not agree")
  names(p) <- names(y) <- names(logit) <- NULL
  if (!missing(group)) {
    if (length(group) == 1 && is.logical(group) && group)
      group <- rep("", length(y))
    if (!is.factor(group))
      group <-
        if (is.logical(group) || is.character(group))
          as.factor(group)
    else
      cut2(group, g = g.group)
    names(group) <- NULL
    nma <- !(is.na(p + y + weights) | is.na(group))
    ng <- length(levels(group))
  } else {
    nma <- !is.na(p + y + weights)
    ng <- 0
  }
  if(any(nma == FALSE)) {
    tmpmess = "There are observations with missing values. These are removed."
    warning(tmpmess, immediate. = TRUE)
    wmess = c(wmess, tmpmess)
  }
  logit <- logit[nma]
  y <- y[nma]
  p <- p[nma]
  if(ng > 0) {
    group <- group[nma]
    weights <- weights[nma]
    return(val.probg(p, y, group, evaluate, weights, normwt, nmin))
  }

  # Fix 'No visible global binding for global variable' note
  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  ymin <- ymax <- xbeta <- lower <- upper <- xend <- yend <- NULL

  # Sort vector with probabilities
  y     <- y[order(p)]
  logit <- logit[order(p)]
  p     <- p[order(p)]


  if (length(p) > 5000 & smooth == "loess") {
    warning("Number of observations > 5000, RCS is recommended.",
            immediate. = TRUE)
  }
  if (length(p) > 1000 & CL.BT == TRUE) {
    warning("Number of observations is > 1000, this could take a while...",
            immediate. = TRUE)
  }

  if(length(unique(p)) == 1) {
    # Adjusted 2022-09-26
    P       <- mean(y)
    Intc    <- log(P/(1 - P))
    n       <- length(y)
    D       <- -1/n
    L01     <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = TRUE)
    L.cal   <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = TRUE)
    U.chisq <- L01 - L.cal
    U.p     <- 1 - pchisq(U.chisq, 1)
    U       <- (U.chisq - 1)/n
    Q       <- D - U
    cl.auc  <- ci.auc(y, p, cl.level, method.ci)


    stats   <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y - p[1])^2), Intc, 0, rep(abs(p[1] - P), 2), 1)
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                      "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier",
                      "Intercept", "Slope", "Emax", "Eavg", "ECI")
    Results =
      structure(
        list(
          call = call,
          stats = stats,
          cl.level = cl.level,
          Calibration = list(
            Intercept = c("Point estimate" = unname(stats["Intercept"]),
                          "Lower confidence limit" = NA,
                          "Upper confidence limit" = NA),
            Slope = c("Point estimate" = unname(stats["Slope"]),
                      "Lower confidence limit" = NA,
                      "Upper confidence limit" = NA)
          ),
          Cindex = c("Point estimate" = unname(stats["C (ROC)"]),
                     "Lower confidence limit" = cl.auc[2],
                     "Upper confidence limit" = cl.auc[3])
        ), class = "CalibrationCurve"
      )
    return(Results)
  }
  i <- !is.infinite(logit)
  nm <- sum(!i)
  if(nm > 0)
    warning(paste(nm, "observations deleted from logistic calibration due to probs. of 0 or 1"))
  i.2  <- i
  f.or <- glm(y[i] ~ logit[i], family = binomial) # lrm(y[i] ~ logit[i])
  f    <- lrm.fit(logit[i], y[i])
  # glm(y ~ offset(Eta), family = binomial, control = glm.control(maxit = 1e2))
  # glm(y ~ Eta, family = binomial)

  cl.slope <- confint(f, level = cl.level)[2, ]
  f2   <-	lrm.fit(offset = logit[i], y = y[i])
  if(f2$fail){
    warning("The lrm function did not converge when computing the calibration intercept!",immediate.=TRUE)
    f2 <- list()
    f2$coef <- NA
    cl.interc <- rep(NA,2)
  } else{
    cl.interc <- confint(f2, level = cl.level)
  }
  stats <- f$stats
  cl.auc <- ci.auc(y, p, cl.level, method.ci)

  n <- stats["Obs"]
  predprob <- seq(emax.lim[1], emax.lim[2], by = 0.0005)
  lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
  calp <- 1/(1 + exp( - lt))
  emax <- max(abs(predprob - calp))
  if (pl) {
    gg = ggplot(data.frame()) +
      xlim(xlim) + ylim(ylim) +
      geom_line(data = data.frame(x = 0:1, y = 0:1), aes(x = x, y = y, colour = "Ideal"), linewidth = lwd.ideal, show.legend = TRUE) +
      labs(x = xlab, y = ylab)

    legCol    = c("Ideal" = col.ideal)
    lt        = lty.ideal
    lw.d      = lwd.ideal
    marks     = NA
    calCurves = list()


    if (logistic.cal) {
      if (min(p) > plogis(-7) | max(p) < plogis(7)) {
        lrm.fit.1 = lrm(y[i.2] ~ qlogis(p[i.2]))
        logCal    = data.frame(x = p[i.2], y = plogis(lrm.fit.1$linear.predictors))
        gg = gg + geom_line(data = logCal, show.legend = TRUE,
                            aes(x = x, y = y, color = "Logistic calibration"), linewidth = lwd.log, linetype = lty.log)
      } else {
        logit  = seq(-7, 7, length = 200)
        prob   = 1 / (1 + exp(-logit))
        pHat   = binomial()$linkinv(cbind(1, logit) %*% coef(f))
        logCal = data.frame(x = prob, y = pHat)
        gg = gg + geom_line(data = logCal, aes(x = x, y = y, color = "Logistic calibration"), linewidth = lwd.log, linetype = lty.log)
      }
      calCurves$LogisticCalibration = logCal
      legCol = c(legCol, "Logistic calibration" = col.log)
      lt      <- c(lt, lty.log)
      lw.d    <- c(lw.d, lwd.log)
      marks   <- c(marks, NA)
    }
    if (smooth == "loess") {
      argzLoess$formula = y ~ p
      SmFit = do.call("loess", argzLoess)
      Sm    = data.frame(x = unname(SmFit$x), y = SmFit$fitted)
      Sm.01 = Sm

      if(any(Sm$y < 0)) {
        sel = which(Sm$y < 0)
        sel = c(sel[length(sel)], sel[length(sel)] + 1)
        tmp = Sm[sel, ]
        Sm    = Sm[Sm$y >= 0 & Sm$y <= 1, ]
        Sm  = rbind.data.frame(
          data.frame(x = predict(lm(x ~ y, data = tmp), data.frame(y = 0)), y = 0),
          Sm
        )
      }

      if (connect.smooth) {
        gg = gg + geom_line(data = Sm, aes(x = x, y = y, color = "Flexible calibration (Loess)"), linetype = lty.smooth, linewidth = lwd.smooth)

        legCol = c(legCol, "Flexible calibration (Loess)" = col.smooth)
        lt    <- c(lt, lty.smooth)
        lw.d  <- c(lw.d, lwd.smooth)
        marks <- c(marks, NA)
      } else {
        gg = gg + geom_point(data = Sm, aes(x = x, y = y, color = "Flexible calibration (Loess)"))
        legCol = c(legCol, "Flexible calibration (Loess)")
        lt    <- c(lt, 0)
        lw.d  <- c(lw.d, 1)
        marks <- c(marks, 1)
      }

      if(CL.smooth != FALSE) {
        if(CL.BT) {
          to.pred = seq(min(p), max(p), length = 200)
          res.BT  = replicate(2000, BT.samples(y, p, to.pred))
          CL.BT   = apply(res.BT, 1, quantile, c(0.025, 0.975))
          colnames(CL.BT) = to.pred
          dfCL    = data.frame(x = to.pred, y = apply(res.BT, 1, quantile, 0.5), ymin = CL.BT[1, ], ymax = CL.BT[2, ])
          rownames(dfCL) = NULL
        } else {
          cl.loess = predict(SmFit, type = "fitted", se = TRUE)
          dfCL     = data.frame(x = p, ymin = with(cl.loess, fit - qnorm(1 - a / 2) * se.fit), ymax = with(cl.loess, fit + qnorm(1 - a / 2) * se.fit))
        }
        dfCL[dfCL$ymax < 0, "ymax"] <- dfCL[dfCL$ymin < 0, "ymin"] <- 0
        dfCL[dfCL$ymax > 1, "ymax"] <- dfCL[dfCL$ymin > 1, "ymin"] <- 1
        if (CL.smooth == "fill") {
          gg = gg + geom_ribbon(data = dfCL, aes(x = x, ymin = ymin, ymax = ymax),
                                fill = rgb(177, 177, 177, 177, maxColorValue = 255))
        } else{
          gg =
            gg +
            geom_line(data = dfCL[dfCL$ymin > 0, ], aes(x = x, y = ymin, color = "CL flexible"), linetype = 2, linewidth = 1) +
            geom_line(data = dfCL[dfCL$ymax < 1, ], aes(x = x, y = ymax), linetype = 2, linewidth = 1, col = col.smooth)

          legCol = c(legCol, "CL flexible" = col.smooth)
          lt      <- c(lt, 2)
          lw.d    <- c(lw.d, 1)
          marks   <- c(marks, NA)
        }
      }
      if(any(Sm$y < 0)) {
        sel = which(Sm$y < 0)
        sel = c(sel[length(sel)], sel[length(sel)] + 1)
        tmp = Sm[sel, ]
        Sm    = Sm[Sm$y >= 0 & Sm$y <= 1, ]
        Sm  = rbind.data.frame(
          data.frame(x = predict(lm(x ~ y, data = tmp), data.frame(y = 0)), y = 0),
          Sm
        )
      }
      colnames(Sm) = c("x", "y")
      if(exists("dfCL", envir = environment())) {
        flexCal = if("CL.BT" %in% names(call) && call$CL.BT) list(loessFit = Sm, BootstrapConfidenceLimits = dfCL) else merge(Sm, dfCL, by = "x")
      } else {
        flexCal = Sm
      }
      calCurves$FlexibleCalibration = flexCal

      cal.smooth <- approx(Sm.01, xout = p, ties = "ordered")$y
      eavg       <- mean(abs(p - cal.smooth))
      ECI        <- mean((p - cal.smooth) ^ 2) * 100
    } else if (smooth == "rcs") {
      argzRCS = alist(x = p,
                      y = y,
                      model = "logistic",
                      nk = nr.knots,
                      show = "prob",
                      statloc = "none",
                      plot = FALSE,
                      showknots = FALSE,
                      xrange = c(min(na.omit(p)), max(na.omit(p))),
                      lty = lty.smooth)
      nkDecrease <- function(Argz) {
        tryCatch(
          do.call(".rcspline.plot", Argz),
          error = function(e) {
            nk = eval(Argz$nk)
            warning(paste0("The number of knots led to estimation problems, nk will be set to ", nk - 1), immediate. = TRUE)
            if(nk < 3)
              stop("Nk = 3 led to estimation problems.")
            Argz$nk = nk - 1
            nkDecrease(Argz)
          }
        )
      }
      rcsFit        = nkDecrease(argzRCS)
      rcsDf         = as.data.frame(rcsFit)
      calCurves$RCS = rcsDf
      gg = gg +
        geom_line(data = rcsDf, aes(x = x, y = xbeta, color = "Flexible calibration (RCS)"), linetype = lty.smooth, linewidth = lwd.smooth) +
        geom_line(data = rcsDf, aes(x = x, y = lower, color = "CL flexible"), linetype = 2, linewidth = 1) +
        geom_line(data = rcsDf, aes(x = x, y = upper), linetype = 2, linewidth = 1, col = col.smooth)

      legCol = c(legCol, "Flexible calibration (RCS)" = col.smooth, "CL flexible" = col.smooth)
      lt      <- c(lt, lty.smooth, 2)
      lw.d    <- c(lw.d, rep(lwd.smooth, 2))
      marks   <- c(marks, NA, NA)
    }
    if (!missing(m) | !missing(g) | !missing(cuts)) {
      if (!missing(m))
        q <- cut2(p,
                  m = m,
                  levels.mean = TRUE,
                  digits = 7)
      else if (!missing(g))
        q <- cut2(p,
                  g = g,
                  levels.mean = TRUE,
                  digits = 7)
      else if (!missing(cuts))
        q <- cut2(p,
                  cuts = cuts,
                  levels.mean = TRUE,
                  digits = 7)
      means <- as.single(levels(q))
      prop <- tapply(y, q, function(x) mean(x, na.rm = TRUE))
      gg = gg + geom_point(data = data.frame(x = means, y = prop), aes(x = x, y = y, color = "Grouped observations"), shape = 2, size = 3)
      #18.11.02: CI triangles
      ng	<- tapply(y, q, length)
      og	<- tapply(y, q, sum)
      ob	<- og / ng
      se.ob	<- sqrt(ob * (1 - ob) / ng)
      g		<- length(as.single(levels(q)))

      for (i in 1:g)
        gg = gg + geom_line(data = data.frame(x = c(means[i], means[i]), y = c(prop[i], min(1, prop[i] + 1.96 * se.ob[i]))), aes(x = x, y = y))
      for (i in 1:g)
        gg = gg + geom_line(data = data.frame(x = c(means[i], means[i]), y = c(prop[i], max(0, prop[i] - 1.96 * se.ob[i]))), aes(x = x, y = y))

      if (connect.group) {
        gg = gg + geom_line(data = data.frame(x = means, y = prop), aes(x = x, y = y))
        lt   <- c(lt, 1)
        lw.d <- c(lw.d, 1)
      } else {
        lt   <- c(lt, 0)
        lw.d <- c(lw.d, 0)
      }
      legCol = c(legCol, "Grouped observations" = "black")
      marks   <- c(marks, 2)
    }
  }
  lr      = stats["Model L.R."]
  p.lr    = stats["P"]
  D       = (lr - 1) / n
  L01     = -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
  U.chisq = L01 - f$deviance[2]
  p.U     = 1 - pchisq(U.chisq, 2)
  U       = (U.chisq - 2) / n
  Q       = D - U
  Dxy     = stats["Dxy"]
  C       = stats["C"]
  R2      = stats["R2"]
  B       = sum((p - y) ^ 2) / n
  Bmax    = mean(y) * (1 - mean(y)) ^ 2 + (1 - mean(y)) * mean(y) ^ 2
  Bscaled = 1 - B / Bmax

  stats = c(Dxy,
             C,
             R2,
             D,
             lr,
             p.lr,
             U,
             U.chisq,
             p.U,
             Q,
             B,
             f2$coef[1],
             f$coef[2],
             emax,
             Bscaled)
  names(stats) = c(
    "Dxy",
    "C (ROC)",
    "R2",
    "D",
    "D:Chi-sq",
    "D:p",
    "U",
    "U:Chi-sq",
    "U:p",
    "Q",
    "Brier",
    "Intercept",
    "Slope",
    "Emax",
    "Brier scaled"
  )
  if (smooth == "loess")
    stats <- c(stats, c(Eavg = eavg), c(ECI = ECI))

  # Cut off definition
  if(!missing(cutoff)) {
    gg = gg + geom_segment(aes(x = cutoff, y = .1, xend = cutoff, yend = -0.025), arrow = arrow(length = unit(.15, "npc")))
  }
  if(pl) {
    if(!is.logical(statloc)) {
      if(dostats[1] == TRUE){
        stats.2 <- paste('Calibration\n',
                         '...intercept: '
                         , sprintf(paste("%.", roundstats, "f", sep = ""), stats["Intercept"]), " (",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.interc[1]), " to ",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.interc[2]), ")", '\n',
                         '...slope: '
                         , sprintf(paste("%.", roundstats, "f", sep = ""), stats["Slope"]), " (",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.slope[1]), " to ",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.slope[2]), ")", '\n',
                         'Discrimination\n',
                         '...c-statistic: '
                         , sprintf(paste("%.", roundstats, "f", sep = ""), stats["C (ROC)"]), " (",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.auc[2]), " to ",
                         sprintf(paste("%.", roundstats, "f", sep = ""), cl.auc[3]), ")"
                         , sep = '')
        gg = gg + annotate("text", x = statloc[1], y = statloc[2], label = stats.2, hjust = 0, size = size)
      } else {
        dostats <- dostats
        leg <- format(names(stats)[dostats])	#constant length
        leg <- paste0(leg, ":", format(stats[dostats], digits=roundstats))
        if(!is.list(statloc))
          statloc <- list(x = statloc[1], y = statloc[2])
        gg =
          gg +
          annotate("text", x = statloc$x, y = statloc$y, label = paste(format(names(stats[dostats])), collapse = "\n"),
                        hjust = 0, size = size) +
          annotate("text", x = statloc$x + (xlim[2] - xlim[1])/3, y = statloc$y, label = paste(format(round(stats[dostats], digits = roundstats)), collapse = "\n"),
                   hjust = 0, size = size)
      }
    }
    if(is.character(riskdist)) {
      if (riskdist == "calibrated") {
        x <- f$coef[1] + f$coef[2] * log(p / (1 - p))
        x <- 1 / (1 + exp(-x))
        x[p == 0] <- 0
        x[p == 1] <- 1
      } else {
        x <- p
      }
      bins <- seq(0, min(1, max(xlim)), length = 101)
      x <- x[x >= 0 & x <= 1]
      #08.04.01,yvon: distribution of predicted prob according to outcome
      f0	<- table(cut(x[y == 0], bins))
      f1	<- table(cut(x[y == 1], bins))
      j0	<- f0 > 0
      j1	<- f1 > 0
      bins0 <- (bins[-101])[j0]
      bins1 <- (bins[-101])[j1]
      f0	<- f0[j0]
      f1	<- f1[j1]
      maxf <- max(f0, f1)
      f0	<- (0.1 * f0) / maxf
      f1	<- (0.1 * f1) / maxf

      gg =
        gg +
        geom_segment(data = data.frame(x = bins1, xend = bins1, y = rep(line.bins, length(bins1)), yend = c(length.seg * f1 + line.bins)),
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_segment(data = data.frame(x = bins0, xend = bins0, y = rep(line.bins, length(bins0)), yend = c(length.seg * -f0 + line.bins)),
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_line(data = data.frame(x = c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 0.01), y = c(line.bins, line.bins)),
                  aes(x = x, y = y)) +
        annotate(geom = "text", x = max(bins0, bins1) + dist.label, y = line.bins + dist.label2, label = d1lab, size = size.d01) +
        annotate(geom = "text", x = max(bins0, bins1) + dist.label, y = line.bins - dist.label2, label = d0lab, size = size.d01)

    }
  }
  gg =
    gg +
    scale_color_manual("", values = legCol, breaks = names(legCol)) +
    guides(colour = guide_legend(override.aes = list(linetype = lt, shape = marks, linewidth = lw.d * 0.5, size = 7))) +
    theme_bw() +
    theme(plot.background=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.margin = margin(11, 11, 5.5, 5.5, "points"), legend.position = "bottom")
  Results =
    structure(
      list(
        call   = call,
        ggPlot = gg,
        stats  = stats,
        cl.level = cl.level,
        Calibration = list(
          Intercept = c("Point estimate" = unname(stats["Intercept"]),
                        "Lower confidence limit" = cl.interc[1],
                        "Upper confidence limit" = cl.interc[2]),
          Slope = c("Point estimate" = unname(stats["Slope"]),
                    "Lower confidence limit" = cl.slope[1],
                    "Upper confidence limit" = cl.slope[2])
        ),
        Cindex = c("Point estimate" = unname(stats["C (ROC)"]),
                   "Lower confidence limit" = cl.auc[2],
                   "Upper confidence limit" = cl.auc[3]),
        warningMessages = wmess,
        CalibrationCurves = calCurves
      ), class = "ggplotCalibrationCurve"
    )
  return(Results)
}
