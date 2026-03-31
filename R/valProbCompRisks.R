#' Plot a calibration curve for a competing risks model
#'
#' Evaluates the calibration performance of a competing risks prediction model
#' (either Cause-Specific Cox or Fine-Gray subdistribution hazards) at a specified
#' time horizon.  The calibration curve is estimated via a flexible regression
#' approach (restricted cubic splines on the complementary log-log scale) and,
#' optionally, using pseudo-observations smoothed with loess.
#'
#' @inheritParams valProbSurvival
#' @param fit the model fit.  Must be of class \code{\link[riskRegression]{CSC}}
#'   (cause-specific hazards) or \code{\link[riskRegression]{FGR}} (Fine-Gray
#'   subdistribution hazards).
#' @param cause integer or character, the cause of interest (default \code{1}).
#' @param nk integer, number of **internal** knots for the natural cubic spline
#'   used in the FGR calibration model on the complementary log-log scale.
#'   \code{nk} internal knots produce \code{nk + 1} natural-spline basis functions
#'   (as in \code{\link[splines]{ns}(x, df = nk + 1)}).
#'   Austin et al. advise using between 3 (more smoothing) and 5 (less smoothing)
#'   internal knots; the default is \code{nk = 3} which is robust for typical
#'   validation set sizes.  If the requested number of knots leads to numerical
#'   problems, \code{nk} is automatically reduced until the model converges.
#' @param pseudo logical.  If \code{TRUE}, an additional calibration curve based
#'   on pseudo-observations (loess-smoothed) is computed and plotted.  Requires the
#'   \pkg{geepack} package for intercept/slope estimation. Default is \code{FALSE}.
#' @param bandwidth numeric, the loess bandwidth for smoothing pseudo-observations
#'   when \code{pseudo = TRUE}. Default is \code{0.10}.
#' @param addRCS logical, indicates whether the RCS-based calibration curve should
#'   be added to the plot. Default is \code{TRUE}.
#' @param addPseudo logical, if \code{pseudo = TRUE}, indicates whether to add the
#'   pseudo-observation loess curve to the plot. Default is \code{TRUE}.
#' @param col.pseudo color of the pseudo-observation loess curve. Default \code{"blue"}.
#' @param lty.pseudo linetype of the pseudo-observation loess curve. Default \code{2}.
#' @param lwd.pseudo line width of the pseudo-observation loess curve. Default \code{1}.
#' @param fill.pseudo fill color for the confidence band around the pseudo
#'   curve (only for \code{plotCal = "ggplot"}). Default uses a blue-tinted transparent fill.
#'
#' @details
#' Two calibration curve approaches are provided:
#' \itemize{
#'   \item \strong{RCS (Flexible regression)}: The complementary log-log transform of
#'     the predicted risks is used as covariate in a Fine-Gray subdistribution hazards
#'     model with restricted cubic splines (\code{\link[riskRegression]{FGR}}).
#'     This approach closely follows van Geloven et al. (2022).
#'   \item \strong{Pseudo-observations} (optional, requires \pkg{geepack}): Pseudo-
#'     observations are computed via \code{\link[riskRegression]{Score}} and smoothed
#'     with loess.  Calibration intercept and slope are estimated via a GEE model
#'     (\code{geepack::geese}) on the pseudo-observations.
#' }
#'
#' @return An object of class \code{CompRisksCalibrationCurve} with the following slots:
#' \item{call}{the matched call.}
#' \item{stats}{a list containing performance measures of calibration.}
#' \item{alpha}{the significance level used.}
#' \item{Calibration}{calibration slope together with confidence intervals.}
#' \item{CalibrationCurves}{The coordinates for plotting the calibration curves.}
#' \item{ggPlot}{the ggplot object when \code{plotCal = "ggplot"}.}
#'
#' @references van Geloven N, Giardiello D, Bonneville E F, Teece L, Ramspek C L,
#' van Smeden M et al. (2022). Validation of prediction models in the presence of
#' competing risks: a guide through modern methods. \emph{BMJ}, \bold{377:e069249},
#' doi:10.1136/bmj-2021-069249
#'
#' @examples
#' \dontrun{
#' library(CalibrationCurves)
#' library(riskRegression)
#'
#' # Rotterdam breast cancer data: rdata = development, vdata = validation
#' data(rdata)
#' data(vdata)
#'
#' # Fit a cause-specific hazards model on the development cohort
#' # cause 1 = breast-cancer recurrence/death; cause 2 = competing death
#' cscFit <- CSC(Hist(time, status_num) ~ age + size + ncat + hr_status,
#'               data = rdata)
#'
#' # Assess calibration on the external validation cohort at 5 years
#' calComp <- valProbCompRisks(cscFit, vdata, cause = 1,
#'                             timeHorizon = 5, plotCal = "ggplot")
#' calComp
#' }
#' @importFrom stats gaussian
#' @importFrom geepack geese
#' @importFrom prodlim Hist
#' @importFrom survival survfit Surv
#' @export
valProbCompRisks <- function(fit,
                             valdata,
                             cause       = 1,
                             weights     = NULL,
                             alpha       = 0.05,
                             timeHorizon = 5,
                             nk          = 3,
                             pseudo      = FALSE,
                             bandwidth   = 0.10,
                             plotCal     = c("none", "base", "ggplot"),
                             addRCS      = TRUE,
                             addPseudo   = TRUE,
                             CL.rcs      = c("fill", "line"),
                             xlab        = "Predicted probability",
                             ylab        = "Observed proportion",
                             xlim        = c(-0.02, 1),
                             ylim        = c(-0.15, 1),
                             lty.ideal   = 1, col.ideal = "red", lwd.ideal = 1,
                             lty.rcs     = 1, col.rcs   = "black", lwd.rcs = 1,
                             fill.rcs    = rgb(177, 177, 177, 177, maxColorValue = 255),
                             col.pseudo  = "blue", lty.pseudo = 2, lwd.pseudo = 1,
                             fill.pseudo = rgb(173, 216, 230, 100, maxColorValue = 255),
                             riskdist    = "predicted",
                             d0lab = "0", d1lab = "1", size.d01 = 5,
                             dist.label  = 0.01, line.bins = -0.05,
                             dist.label2 = 0.04, length.seg = 0.85,
                             legendloc   = c(0.50, 0.27)) {

  callFn  <- match.call()
  plotCal <- match.arg(plotCal)
  CL.rcs  <- match.arg(CL.rcs)

  # Suppress global-binding notes
  pred <- obs <- lower <- upper <- xend <- yend <- NULL

  # --------------------------------------------------------------------------
  # Input checks
  # --------------------------------------------------------------------------
  allowedClasses <- c("CauseSpecificCox", "CSC", "FGR", "selectCox", "riskRegression")
  if (!inherits(fit, allowedClasses))
    stop("'fit' must be a CSC or FGR model from the riskRegression package.")

  stats <- list()

  # --------------------------------------------------------------------------
  # Predicted risks at timeHorizon
  # --------------------------------------------------------------------------
  valdata$pred <- predictRisk(fit, newdata = valdata,
                              times = timeHorizon, cause = cause)

  # --------------------------------------------------------------------------
  # Discrimination: AUC and C-index via Score()
  # --------------------------------------------------------------------------
  # Extract formula robustly: newer riskRegression versions may not store
  # fit$formula directly; fall back to formula() S3 dispatch or the call.
  fitFormula <- tryCatch(formula(fit), error = function(e) NULL)
  if (is.null(fitFormula)) fitFormula <- fit$formula
  if (is.null(fitFormula) && !is.null(fit$call))
    fitFormula <- tryCatch(fit$call$formula, error = function(e) NULL)
  vars      <- all.vars(fitFormula)
  timeVar   <- vars[1L]
  statusVar <- vars[2L]

  adjFormula <- stats::as.formula(
    paste0("Hist(", timeVar, ", ", statusVar, ") ~ 1"),
    env = asNamespace("prodlim")
  )

  score_res <- tryCatch(
    Score(
      list("model" = fit),
      formula     = adjFormula,
      cens.model  = "km",
      cens.method = c("ipcw", "pseudo"),
      data        = valdata,
      conf.int    = TRUE,
      times       = timeHorizon,
      metrics     = c("auc", "brier"),
      summary     = "ipa",
      cause       = cause,
      plots       = "calibration"
    ),
    error = function(e) NULL
  )

  if (!is.null(score_res)) {
    stats$AUC        <- score_res$AUC$score
    stats$BrierScore <- score_res$Brier$score
  }

  # --------------------------------------------------------------------------
  # Calibration in the large: O/E via Aalen-Johansen (survfit multi-state)
  # Matches the reference approach: survfit with factor status, pstate column
  # --------------------------------------------------------------------------
  z <- qnorm(1 - alpha / 2)

  aj_sfit <- tryCatch(
    summary(
      survival::survfit(
        survival::Surv(valdata[[timeVar]], factor(valdata[[statusVar]])) ~ 1
      ),
      times = timeHorizon
    ),
    error = function(e) NULL
  )

  if (!is.null(aj_sfit) &&
      !is.null(aj_sfit$pstate) &&
      ncol(aj_sfit$pstate) >= cause + 1L) {
    aj_obs <- aj_sfit$pstate[1L, cause + 1L]
    aj_se  <- aj_sfit$std.err[1L, cause + 1L]
    exp_t  <- mean(valdata$pred, na.rm = TRUE)
    OE_t   <- aj_obs / exp_t
    stats$Calibration$InTheLarge <- c(
      "OE"     = OE_t,
      "2.5 %"  = exp(log(OE_t) - z * aj_se / aj_obs),
      "97.5 %" = exp(log(OE_t) + z * aj_se / aj_obs)
    )
  }

  # --------------------------------------------------------------------------
  # Calibration curve: RCS on cloglog scale via FGR
  # --------------------------------------------------------------------------
  valdata$cll_pred <- log(-log(1 - pmax(pmin(valdata$pred, 1 - 1e-8), 1e-8)))

  fitFGRwithKnots <- function(current_nk) {
    rcs_mat           <- splines::ns(valdata$cll_pred, df = current_nk + 1L)
    colnames(rcs_mat) <- paste0("cll_basis_", seq_len(ncol(rcs_mat)))
    valdata_rcs       <- cbind.data.frame(valdata, rcs_mat)
    rcs_formula <- stats::reformulate(
      termlabels = colnames(rcs_mat),
      response   = paste0("Hist(", timeVar, ", ", statusVar, ")"),
      env        = asNamespace("prodlim")
    )
    result <- tryCatch(
      list(fit = FGR(formula = rcs_formula, cause = cause, data = valdata_rcs),
           data = valdata_rcs),
      error = function(e) {
        if (current_nk > 2) {
          warning(paste0("RCS calibration via FGR failed with nk = ", current_nk,
                         ", retrying with nk = ", current_nk - 1L), immediate. = TRUE)
          fitFGRwithKnots(current_nk - 1L)
        } else {
          warning("RCS calibration via FGR failed: ", conditionMessage(e))
          NULL
        }
      }
    )
    result
  }

  fgr_result  <- fitFGRwithKnots(nk)
  calib_fgr   <- if (!is.null(fgr_result)) fgr_result$fit  else NULL
  valdata_rcs <- if (!is.null(fgr_result)) fgr_result$data else NULL

  dat_rcs <- NULL
  if (!is.null(calib_fgr)) {
    obs_rcs  <- predict(calib_fgr, times = timeHorizon, newdata = valdata_rcs)
    dat_rcs  <- data.frame(pred = valdata$pred, obs = obs_rcs)
    dat_rcs  <- dat_rcs[order(dat_rcs$pred), ]

    absdiff  <- abs(dat_rcs$pred - dat_rcs$obs)
    stats$Calibration$Statistics <- c(
      "ICI"  = mean(absdiff),
      setNames(quantile(absdiff, c(0.5, 0.9)), c("E50", "E90")),
      "Emax" = max(absdiff)
    )
  }

  # --------------------------------------------------------------------------
  # Calibration intercept and slope via GEE on pseudo-observations
  # Always computed from Score()$Calibration$plotframe (Giardiello reference
  # approach, section 2.1.4). The `pseudo` flag controls only the LOESS curve.
  # --------------------------------------------------------------------------
  dat_pseudo <- NULL

  pseudos <- tryCatch({
    if (!is.null(score_res) && !is.null(score_res$Calibration$plotframe))
      data.frame(score_res$Calibration$plotframe)
    else
      NULL
  }, error = function(e) NULL)

  if (!is.null(pseudos) && nrow(pseudos) > 0) {
    pseudos        <- pseudos[order(pseudos$risk), ]
    pseudos$cll_pred <- log(-log(1 - pmax(pmin(pseudos$risk, 1 - 1e-8), 1e-8)))

    fit_cal_int_main <- tryCatch(
      geepack::geese(
        pseudovalue ~ offset(cll_pred),
        data      = pseudos,
        id        = pseudos$riskRegression_ID,
        scale.fix = TRUE,
        family    = gaussian(),
        mean.link = "cloglog",
        corstr    = "independence",
        jack      = TRUE
      ),
      error = function(e) NULL
    )
    fit_cal_slope_main <- tryCatch(
      geepack::geese(
        pseudovalue ~ offset(cll_pred) + cll_pred,
        data      = pseudos,
        id        = pseudos$riskRegression_ID,
        scale.fix = TRUE,
        family    = gaussian(),
        mean.link = "cloglog",
        corstr    = "independence",
        jack      = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_cal_slope_main)) {
      sm      <- summary(fit_cal_slope_main)$mean
      est_raw <- sm["cll_pred", "estimate"]
      san_se  <- sm["cll_pred", "san.se"]
      stats$Calibration$Slope <- c(
        "calibration slope" = 1 + est_raw,
        "2.5 %"  = 1 + (est_raw - z * san_se),
        "97.5 %" = 1 + (est_raw + z * san_se)
      )
    }
    if (!is.null(fit_cal_int_main)) {
      sm_int  <- summary(fit_cal_int_main)$mean
      est_int <- sm_int[1L, "estimate"]
      san_int <- sm_int[1L, "san.se"]
      stats$Calibration$Intercept <- c(
        "calibration intercept" = est_int,
        "2.5 %"  = est_int - z * san_int,
        "97.5 %" = est_int + z * san_int
      )
    }
  }

  if (pseudo) {
    pseudo_df <- if (!is.null(pseudos)) {
      data.frame(
        pseudovalue       = pseudos$pseudovalue,
        risk              = pseudos$risk,
        riskRegression_ID = pseudos$riskRegression_ID,
        cll_pred          = pseudos$cll_pred
      )
    } else {
      NULL
    }

    if (!is.null(pseudo_df)) {

      # Smooth pseudo-obs with loess for the curve
      sm_loess <- predict(
        stats::loess(pseudovalue ~ risk, data = pseudo_df,
                     degree = 1, span = bandwidth),
        se = TRUE
      )
      dat_pseudo <- data.frame(
        pred  = pseudo_df$risk,
        obs   = sm_loess$fit,
        lower = pmax(sm_loess$fit - qt(1 - alpha / 2, sm_loess$df) * sm_loess$se, 0),
        upper = sm_loess$fit + qt(1 - alpha / 2, sm_loess$df) * sm_loess$se
      )
      dat_pseudo <- dat_pseudo[order(dat_pseudo$pred), ]
    }
  }

  # --------------------------------------------------------------------------
  # Spike histogram data (Giardiello reference approach: single rescaled
  # spike set for all predicted risks, placed below y = 0)
  # --------------------------------------------------------------------------
  spike_bounds  <- c(ylim[1] + 0.02, 0)   # e.g. c(-0.13, 0) for default ylim
  bin_breaks    <- seq(0, max(xlim), length.out = 101)
  spike_freqs   <- table(cut(valdata$pred,
                             breaks = bin_breaks, include.lowest = TRUE))
  spike_bins    <- bin_breaks[-1L]         # right edge of each bin
  spike_fv      <- spike_freqs[spike_freqs > 0]
  spike_pos     <- unname(spike_bins[spike_freqs > 0])
  if (length(spike_fv) > 1L) {
    spike_rescaled <- unname(
      spike_bounds[1L] +
        (spike_bounds[2L] - spike_bounds[1L]) *
        (spike_fv - min(spike_fv)) / (max(spike_fv) - min(spike_fv))
    )
  } else {
    spike_rescaled <- rep(spike_bounds[2L], length(spike_fv))
  }

  spike_df <- if (length(spike_pos) > 0L) {
    data.frame(
      x    = spike_pos,
      xend = spike_pos,
      y    = spike_bounds[1L],
      yend = spike_rescaled
    )
  } else NULL

  # --------------------------------------------------------------------------
  # Plotting
  # --------------------------------------------------------------------------
  if (plotCal == "base") {
    par(xaxs = "i", yaxs = "i", las = 1)
    plot(c(0, 1), c(0, 1), type = "n",
         xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)

    legCol <- c("Ideal" = col.ideal)
    lt     <- lty.ideal; lw.d <- lwd.ideal; marks <- NA

    if (addRCS && !is.null(dat_rcs)) {
      legCol <- c(legCol, "Flexible calibration (rcs)" = col.rcs)
      lt     <- c(lt, lty.rcs); lw.d <- c(lw.d, lwd.rcs); marks <- c(marks, NA)
      lines(dat_rcs$pred, dat_rcs$obs, lwd = lwd.rcs, lty = lty.rcs, col = col.rcs)
    }

    if (pseudo && addPseudo && !is.null(dat_pseudo)) {
      legCol <- c(legCol, "Pseudo-obs (loess)" = col.pseudo)
      lt     <- c(lt, lty.pseudo); lw.d <- c(lw.d, lwd.pseudo); marks <- c(marks, NA)
      lines(dat_pseudo$pred, dat_pseudo$obs, lwd = lwd.pseudo,
            lty = lty.pseudo, col = col.pseudo)
    }
    ll <- legendloc
    if (!is.logical(ll) && !is.list(ll))
      ll <- list(x = ll[1], y = ll[2])
    legend(ll, names(legCol), col = legCol,
           lty = lt, lwd = lw.d, bty = "n", cex = 0.85)
    if (!is.null(spike_df)) {
      segments(spike_df$x, spike_df$y, spike_df$xend, spike_df$yend)
      lines(c(min(spike_df$x) - 0.01, max(spike_df$x) + 0.01),
            c(spike_bounds[1L], spike_bounds[1L]))
    }

  } else if (plotCal == "ggplot") {
    gg <- ggplot(data.frame()) +
      geom_line(data = data.frame(x = 0:1, y = 0:1),
                aes(x = x, y = y, colour = "Ideal"),
                linewidth = lwd.ideal) +
      labs(x = xlab, y = ylab)

    legCol <- c("Ideal" = col.ideal)
    lt     <- lty.ideal; lw.d <- lwd.ideal; marks <- NA

    if (addRCS && !is.null(dat_rcs)) {
      gg     <- gg + geom_line(data = dat_rcs,
                               aes(x = pred, y = obs,
                                   color = "Flexible calibration (rcs)"),
                               linetype = lty.rcs, linewidth = lwd.rcs)
      legCol <- c(legCol, "Flexible calibration (rcs)" = col.rcs)
      lt     <- c(lt, lty.rcs); lw.d <- c(lw.d, lwd.rcs); marks <- c(marks, NA)
    }

    if (pseudo && addPseudo && !is.null(dat_pseudo)) {
      gg <- gg +
        geom_ribbon(data = dat_pseudo,
                    aes(x = pred, ymin = lower, ymax = upper),
                    fill = fill.pseudo, alpha = 0.4) +
        geom_line(data = dat_pseudo,
                  aes(x = pred, y = obs, color = "Pseudo-obs (loess)"),
                  linetype = lty.pseudo, linewidth = lwd.pseudo)
      legCol <- c(legCol, "Pseudo-obs (loess)" = col.pseudo)
      lt     <- c(lt, lty.pseudo); lw.d <- c(lw.d, lwd.pseudo); marks <- c(marks, NA)
    }

    if (!is.null(spike_df))
      gg <- gg + geom_segment(
        data = spike_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE)

    gg <- gg +
      scale_color_manual("", values = legCol, breaks = names(legCol)) +
      guides(colour = guide_legend(
        override.aes = list(linetype = lt, shape = marks,
                            linewidth = lw.d * 0.5, size = 7))) +
      theme_bw() +
      theme(plot.background  = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 12),
            axis.title       = element_text(size = 14),
            plot.margin      = margin(11, 11, 5.5, 5.5, "points"),
            legend.position  = "bottom") +
      coord_cartesian(xlim = xlim, ylim = ylim)
    print(gg)
  }

  # --------------------------------------------------------------------------
  # Return object
  # --------------------------------------------------------------------------
  calCurves <- list(RCS = dat_rcs)
  if (!is.null(dat_pseudo))
    calCurves$Pseudo <- dat_pseudo

  Results <- structure(
    list(
      call             = callFn,
      stats            = stats,
      cause            = cause,
      timeHorizon      = timeHorizon,
      alpha            = alpha,
      Calibration      = list(
        Slope = if (!is.null(stats$Calibration$Slope)) {
          c("Point estimate"          = unname(stats$Calibration$Slope[1]),
            "Lower confidence limit"  = unname(stats$Calibration$Slope[2]),
            "Upper confidence limit"  = unname(stats$Calibration$Slope[3]))
        } else NULL
      ),
      CalibrationCurves = calCurves
    ),
    class = "CompRisksCalibrationCurve"
  )

  if (plotCal == "ggplot")
    Results$ggPlot <- gg

  return(Results)
}
