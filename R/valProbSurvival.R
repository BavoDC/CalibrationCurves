#' Plot a calibration curve for a Cox Proportional Hazards model
#'
#' @inheritParams val.prob.ci.2
#' @param fit the model fit, has to be of type \code{\link[survival]{coxph}}
#' @param valdata the validation data set
#' @param alpha the significance level
#' @param timeHorizon the time point at which the predictions have to be evaluated
#' @param nk the number of knots, for the restricted cubic splines fit
#' @param plotCal indicates if and how the calibration curve has to be plotted.
#' \code{plotCal = "none"} plots no calibration curve, \code{plotCal = "base"} plots
#' the calibration curve using base R (see \code{\link{plot}}) and \code{plotCal = "ggplot"}
#' creates a plot using \code{\link[ggplot2]{ggplot}}
#' @param addCox logical, indicates if the Cox's estimated calibration curve has to be added
#' to the plot
#' @param addRCS logical, indicates if the restricted cubic splines' (RCS) estimated calibration curve has to be added
#' to the plot
#' @param CL.cox \code{"fill"} shows pointwise 95\% confidence limits for the Cox calibration curve with a gray
#' area between the lower and upper limits and \code{"line"} shows the confidence limits with a dotted line
#' @param CL.rcs \code{"fill"} shows pointwise 95\% confidence limits for the RCS calibration curve with a gray
#' area between the lower and upper limits and \code{"line"} shows the confidence limits with a dotted line
#' @param lty.cox if \code{addCox = TRUE}, the linetype of the Cox calibration curve
#' @param col.cox if \code{addCox = TRUE}, the color of the Cox calibration curve
#' @param lwd.cox if \code{addCox = TRUE}, the linewidth of the Cox calibration curve
#' @param fill.cox if \code{addCox = TRUE} and \code{CL.cox = "fill"}, the fill of the Cox calibration curve
#' @param lty.rcs if \code{addRCS = TRUE}, the linetype of the RCS calibration curve
#' @param col.rcs if \code{addRCS = TRUE}, the color of the RCS calibration curve
#' @param lwd.rcs if \code{addRCS = TRUE}, the linewidth of the RCS calibration curve
#' @param fill.rcs if \code{addRCS = TRUE} and \code{CL.rcs = "fill"}, the fill of the RCS calibration curve
#' @param size.d01 controls the size of the labels for events and non-events. Default is 5 and
#' this value is multiplied by 0.25 when \code{plotCal = "base"}.
#'
#' @return An object of type \code{SurvivalCalibrationCurves} with the following slots:
#' @return \item{call}{the matched call.}
#' @return \item{stats}{a list containing performance measures of calibration.}
#' @return \item{alpha}{the significance level used.}
#' @return \item{Calibration}{contains the estimated calibration slope, together with their confidence intervals.}
#' @return \item{CalibrationCurves}{The coordinates for plotting the calibration curves. }
#'
#' @references van Geloven N, Giardiello D, Bonneville E F, Teece L, Ramspek C L, van Smeden M et al. (2022). Validation of prediction models in the presence of competing risks: a guide through modern methods. \emph{BMJ}, \bold{377:e069249}, doi:10.1136/bmj-2021-069249
#'
#' @examples
#' \dontrun{
#' library(CalibrationCurves)
#' data(trainDataSurvival)
#' data(testDataSurvival)
#' sFit = coxph(Surv(ryear, rfs) ~ csize + cnode + grade3, data = trainDataSurvival,
#'  x = TRUE, y = TRUE)
#' calPerf = valProbSurvival(sFit, gbsg5, plotCal = "base", nk = 5)
#' }

valProbSurvival <- function(fit, valdata, alpha = 0.05, timeHorizon = 5, nk = 3,
                            plotCal = c("none", "base", "ggplot"),
                            addCox = FALSE, addRCS = TRUE,
                            CL.cox = c("fill", "line"),
                            CL.rcs = c("fill", "line"),
                            xlab = "Predicted probability", ylab = "Observed proportion",
                            xlim = c(-0.02, 1), ylim = c(-0.15, 1),
                            lty.ideal = 1, col.ideal = "red", lwd.ideal = 1,
                            lty.cox = 1, col.cox = "grey", lwd.cox = 1, fill.cox = "lightgrey",
                            lty.rcs = 1, col.rcs = "black", lwd.rcs = 1, fill.rcs = rgb(177, 177, 177, 177, maxColorValue = 255),
                            riskdist = "predicted", d0lab = "0", d1lab = "1", size.d01 = 5, dist.label = 0.01,
                            line.bins = -0.05, dist.label2 = 0.04, length.seg = 0.85,
                            legendloc =  c(0.50 , 0.27)) {
  callFn = match.call()

  # Fix 'No visible global binding for global variable' note
  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  pred <- obs <- lower <- upper <- xend <- yend <- NULL

  ### Validation of the original model -----------------------
  # Discrimination ---------------------------------------
  if(!inherits(fit, "coxph"))
    stop("Only model fits of class coxph are allowed")
  plotCal   = match.arg(plotCal)
  CL.cox    = match.arg(CL.cox)
  CL.rcs    = match.arg(CL.rcs)
  stats     = list()

  valdata$LP = predict(fit, newdata = valdata, type = "lp")

  argzConc =
    alist(
      data = valdata,
      reverse = TRUE
    )
  argzConc$obj = update(fit$formula, "~ - . + LP")

  HarrellC = do.call("concordance", argzConc)
  argzConc$timewt = "n/G2"
  UnoC     = do.call("concordance", argzConc)

  res_C <- matrix(
    c(
      with(HarrellC, {
        c(concordance,
          concordance - qnorm(1 - alpha/2) * sqrt(var),
          concordance + qnorm(1 - alpha/2) * sqrt(var))
      }),
      with(UnoC, {
        c(concordance,
          concordance - qnorm(1 - alpha/2) * sqrt(var),
          concordance + qnorm(1 - alpha/2) * sqrt(var))
      })
    )
    ,
    nrow = 2,
    ncol = 3,
    byrow = T,
    dimnames = list(c("Harrell C", "Uno C"),
                    c("Estimate", "2.5 %", "97.5 %"))
  )

  stats$Concordance = res_C

  # Uno's time dependent AUC
  UnoTDAUC =
    with(
      timeROC(
        T      = valdata[[all.vars(fit$formula)[1]]],
        delta  = valdata[[all.vars(fit$formula)[2]]],
        marker = valdata$LP,
        cause = 1,
        weighting = "marginal",
        times = max(as.numeric(fit$y)) - 0.01,
        iid = TRUE
      ),
     c(
        "Uno AUC" = AUC[[2]],
        "2.5 %"   = AUC[[2]] - qnorm(1 - alpha / 2) * inference$vect_sd_1[[2]],
        "97. 5 %" = AUC[[2]] + qnorm(1 - alpha / 2) * inference$vect_sd_1[[2]]
      )
     )

  stats$TimeDependentAUC = UnoTDAUC

  # Calibration -----------------------------------------

  # Observed / Expected ratio

  # Observed
  adjFormula = update(fit$formula, ~ - . + 1)
  obj        = summary(survfit(adjFormula, data = valdata), times = timeHorizon)
  obs_t      = 1 - obj$surv

  # Predicted risk
  valdata$pred <- predictRisk(fit, newdata = valdata, times = timeHorizon)

  # Expected
  exp_t = mean(valdata$pred)
  OE_t  = obs_t / exp_t

  OE_summary <- c(
    "OE" = OE_t,
    "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
    "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  )

  stats$Calibration$InTheLarge = OE_summary

  # Calibration plot ----------------------------------
  #valdata$pred.cll <- predictCox(fit, times = timeHorizon, newdata = valdata, type = "lp")


  # Estimate actual risk
  calCox = cph(Surv(ryear, rfs) ~ LP,
             x    = TRUE,
             y    = TRUE,
             surv = TRUE,
             data = valdata
    )

  vcal <- eval(substitute(
    cph(Surv(ryear, rfs) ~ rcs(LP, nk),
                   x    = TRUE,
                   y    = TRUE,
                   surv = TRUE,
                   data = valdata
  ), list(nk = nk)))

  datCox <- cbind.data.frame(
    "obs" = 1 - survest(calCox, times = timeHorizon, newdata = valdata)$surv,

    "lower" = 1 - survest(calCox, times = timeHorizon, newdata = valdata)$upper,

    "upper" = 1 - survest(calCox, times = timeHorizon, newdata = valdata)$lower,
    "pred" = valdata$pred
  )

  datCox <- datCox[order(datCox$pred), ]

  dat_cal <- cbind.data.frame(
    "obs" = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$surv,

    "lower" = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$upper,

    "upper" = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$lower,
    "pred" = valdata$pred
  )

  dat_cal <- dat_cal[order(dat_cal$pred), ]

  # Numerical measures
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  stats$Calibration$Statistics <- c(
    "ICI" = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
    "Emax" = max(absdiff_cph)
  )

  # calibration slope (fixed time point)-------------------------------------
  gval <- coxph(Surv(ryear, rfs) ~ LP, data = valdata)

  stats$Calibration$Slope <- c(
    "calibration slope" = unname(gval$coef),
    "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
    "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
  )

  # Overall performance ---------------------------------------
  stats$Calibration$BrierScore <-
    Score(list("cox" = fit),
                formula  = adjFormula,
                data     = valdata,
                conf.int = TRUE,
                times    = timeHorizon - 0.01,
                cens.model = "km",
                metrics  = "brier",
                summary  = "ipa"
    )$Brier$score

  # Save calibration curves
  calCurves = list(
    CoxCalibration = datCox,
    RCS = dat_cal
    )

  ## Double check!!!
  if(is.character(riskdist)) {
    if (riskdist == "calibrated") {
      x = datCox$obs
      x[datCox$pred == 0] <- 0
      x[datCox$pred == 1] <- 1
    } else {
      x = datCox$pred
    }
    y     = calCox$y[, 2]
    bins  <- seq(0, min(1, max(xlim)), length = 101)
    x     <- x[x >= 0 & x <= 1]
    f0	  <- table(cut(x[y == 0], bins))
    f1	  <- table(cut(x[y == 1], bins))
    j0	  <- f0 > 0
    j1	  <- f1 > 0
    bins0 <- (bins[-101])[j0]
    bins1 <- (bins[-101])[j1]
    f0	  <- f0[j0]
    f1	  <- f1[j1]
    maxf  <- max(f0, f1)
    f0	  <- (0.1 * f0) / maxf
    f1	  <- (0.1 * f1) / maxf

  ## Plotting
  if(plotCal == "base") {
    par(xaxs = "i", yaxs = "i", las = 1)
    plot(
      dat_cal$pred,
      dat_cal$obs,
      type = "l",
      col = "white",
      lty = 1,
      xlim = xlim,
      ylim = ylim,
      lwd = 2,
      xlab = xlab,
      ylab = ylab
    )
    abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
    legCol    = c("Ideal" = col.ideal)
    lt        = lty.ideal
    lw.d      = lwd.ideal
    marks     = NA
    if(addCox) {
      legCol = c(legCol, "Cox calibration" = col.cox)
      lt     = c(lt, lty.cox)
      lw.d   = c(lw.d, lwd.cox)
      marks  = c(marks, NA)
      if(CL.cox == "line") {
        legCol = c(legCol, "CL Cox" = col.cox)
        lt     = c(lt, 2)
        lw.d   = c(lw.d, 1)
        marks  = c(marks, NA)
        lines(datCox$pred,
              datCox$lower,
              type = "l",
              lty = 2,
              lwd = 2)
        lines(datCox$pred,
              datCox$upper,
              type = "l",
              lty = 2,
              lwd = 2)
      } else {
        polygon(
          x = with(datCox, c(pred, rev(pred))),
          y = with(datCox, c(upper, rev(lower))),
          col = fill.cox,
          border = NA
        )
      }
      lines(datCox$pred, datCox$obs, type = "l",
            lwd = lwd.cox, lty = lty.cox)
    }

    if(addRCS) {
      legCol = c(legCol, "Flexible calibration (rcs)" = col.rcs)
      lt     = c(lt, lty.rcs)
      lw.d   = c(lw.d, lwd.rcs)
      marks  = c(marks, NA)
      if(CL.rcs == "line") {
        legCol = c(legCol, "CL flexible (rcs)" = col.rcs)
        lt     = c(lt, 2)
        lw.d   = c(lw.d, 1)
        marks  = c(marks, NA)
        lines(dat_cal$pred,
              dat_cal$lower,
              type = "l",
              lty = 2,
              lwd = 2)
        lines(dat_cal$pred,
              dat_cal$upper,
              type = "l",
              lty = 2,
              lwd = 2)
      } else {
        polygon(
          x = with(dat_cal, c(pred, rev(pred))),
          y = with(dat_cal, c(upper, rev(lower))),
          col = fill.cox,
          border = NA
        )
      }
      lines(dat_cal$pred, dat_cal$obs, type = "l", lwd = lwd.rcs, lty = lty.rcs)
    }
    ll <- legendloc
    if (!is.logical(ll))
      if (!is.list(ll))
        ll <- list(x = ll[1], y = ll[2])
    legend(ll,
           c("Ideal calibration",
             "Cox calibration",
             "95% confidence interval"),
           col = c(2, 1, 1),
           lty = c(2, 1, 2),
           lwd = c(2, 2, 2),
           bty = "n",
           cex = 0.85)
    segments(bins1, line.bins, bins1, length.seg * f1 + line.bins)
    segments(bins0, line.bins, bins0, length.seg * -f0 + line.bins)
    lines(c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 0.01), c(line.bins, line.bins))
    text(max(bins0, bins1) + dist.label,
         line.bins + dist.label2,
         d1lab,
         cex = size.d01 * 0.25)
    text(max(bins0, bins1) + dist.label,
         line.bins - dist.label2,
         d0lab,
         cex = size.d01 * 0.25)
  } else if(plotCal == "ggplot") {
    gg =
      ggplot(data.frame()) +
      geom_line(data = data.frame(x = 0:1, y = 0:1),
                aes(x = x, y = y, colour = "Ideal"),
                linewidth = lwd.ideal,
                show.legend = TRUE) +
      labs(x = xlab, y = ylab)
    legCol    = c("Ideal" = col.ideal)
    lt        = lty.ideal
    lw.d      = lwd.ideal
    marks     = NA
    if(addCox) {
      gg =
        gg +
        geom_line(data = datCox, aes(x = pred, y = obs, color = "Cox calibration"), linetype = lty.cox, linewidth = lwd.cox)
      legCol = c(legCol, "Cox calibration" = col.cox)
      lt     = c(lt, lty.cox)
      lw.d   = c(lw.d, lwd.cox)
      marks  = c(marks, NA)
      gg =
        if(CL.cox == "line") {
          legCol = c(legCol, "CL Cox" = col.cox)
          lt     = c(lt, 2)
          lw.d   = c(lw.d, 1)
          marks  = c(marks, NA)
          gg +
            geom_line(data = datCox, aes(x = pred, y = lower, color = "CL Cox"), linetype = 2, linewidth = 1) +
            geom_line(data = datCox, aes(x = pred, y = upper, color = "CL Cox"), linetype = 2, linewidth = 1)
        } else {
          gg +
            geom_ribbon(data = datCox, aes(x = pred, ymin = lower, ymax = upper),
                        fill = fill.cox)
        }
    }

    if(addRCS) {
      gg =
        gg +
        geom_line(data = dat_cal, aes(x = pred, y = obs, color = "Flexible calibration (rcs)"), linetype = lty.cox, linewidth = lwd.cox)
      legCol = c(legCol, "Flexible calibration (rcs)" = col.rcs)
      lt     = c(lt, lty.rcs)
      lw.d   = c(lw.d, lwd.rcs)
      marks  = c(marks, NA)
      gg =
        if(CL.rcs == "line") {
          legCol = c(legCol, "CL flexible (rcs)" = col.rcs)
          lt     = c(lt, 2)
          lw.d   = c(lw.d, 1)
          marks  = c(marks, NA)
          gg +
            geom_line(data = dat_cal, aes(x = pred, y = lower, color = "CL flexible (rcs)"), linetype = 2, linewidth = 1) +
            geom_line(data = dat_cal, aes(x = pred, y = upper, color = "CL flexible (rcs)"), linetype = 2, linewidth = 1)
        } else {
          gg +
            geom_ribbon(data = dat_cal, aes(x = pred, ymin = lower, ymax = upper),
                        fill = fill.rcs)
        }
    }

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
      gg =
        gg +
        scale_color_manual("", values = legCol, breaks = names(legCol)) +
        guides(colour = guide_legend(override.aes = list(linetype = lt, shape = marks, linewidth = lw.d * 0.5, size = 7))) +
        theme_bw() +
        theme(plot.background  = element_blank(),
              panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text        = element_text(size = 12),
              axis.title       = element_text(size = 14),
              plot.margin      = margin(11, 11, 5.5, 5.5, "points"), legend.position = "bottom")
      gg =
        gg + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  }
  Results = structure(
    list(call = callFn,
         stats = stats,
         Calibration = list(
           Slope = c("Point estimate" = unname(stats$Calibration$Slope[1]),
                     "Lower confidence limit" = unname(stats$Calibration$Slope[2]),
                     "Upper confidence limit" = unname(stats$Calibration$Slope[3]))
         ),
         alpha = alpha,
         CalibrationCurves = calCurves),
    class = "SurvivalCalibrationCurve"
  )
  if(plotCal == "ggplot")
    Results$ggPlot = gg
  return(Results)
}

