#' Recalibrate a prediction model using Platt scaling or isotonic regression
#'
#' Applies post-hoc recalibration to a set of predicted probabilities using
#' either Platt scaling (logistic regression on the logit scale) or isotonic
#' regression. Optionally plots (1) a before/after calibration curve comparison
#' and (2) a transformation plot showing the recalibration mapping.
#'
#' @param p numeric vector of predicted probabilities (before recalibration).
#' @param y binary outcome vector (0/1).
#' @param method character, recalibration method: \code{"platt"} (default) or
#'   \code{"isotonic"}.
#' @param plot logical, whether to produce plots.  Default \code{TRUE}.
#' @param plottype character, which plots to produce: \code{"both"} (default),
#'   \code{"calibration"} (before/after calibration comparison only), or
#'   \code{"transform"} (recalibration mapping only).
#' @param alpha numeric, significance level for confidence bands. Default \code{0.05}.
#' @param smooth character, smoother for the calibration curve: \code{"loess"}
#'   (default) or \code{"rcs"}.
#' @param nr.knots number of knots when \code{smooth = "rcs"}. Default \code{3}.
#' @param xlab x-axis label for the calibration plot. Default
#'   \code{"Predicted probability"}.
#' @param ylab y-axis label for the calibration plot. Default
#'   \code{"Observed proportion"}.
#' @param col.before color of the calibration curve before recalibration.
#'   Default \code{"steelblue"}.
#' @param col.after color of the calibration curve after recalibration.
#'   Default \code{"darkgreen"}.
#' @param verbose logical, if \code{TRUE} (default) prints a summary of the
#'   recalibration transform.
#'
#' @details
#' \strong{Platt scaling} fits the logistic regression model
#' \deqn{\text{logit}(P(y=1)) = \alpha + \beta \cdot \text{logit}(\hat{p})}
#' on the validation data and uses the fitted probabilities as recalibrated
#' predictions.  The calibration slope \eqn{\beta} and intercept \eqn{\alpha}
#' are reported.
#'
#' \strong{Isotonic regression} fits an isotone (monotone non-decreasing)
#' mapping from predicted to observed probabilities using the pool adjacent
#' violators algorithm via \code{\link[stats]{isoreg}}.  The resulting
#' step-function maps each original prediction to a calibrated value.
#'
#' For both methods, the following proper scoring rules are computed and
#' compared before and after recalibration:
#' \itemize{
#'   \item Brier score;
#'   \item Log-loss (negative mean log-likelihood);
#'   \item Calibration intercept and slope (from logistic regression on
#'         recalibrated predictions).
#' }
#'
#' @return An object of class \code{RecalibratedPredictions} with:
#' \item{call}{the matched call.}
#' \item{method}{the recalibration method used.}
#' \item{p_original}{the original predicted probabilities.}
#' \item{p_recal}{the recalibrated predicted probabilities.}
#' \item{before}{list of calibration statistics computed on the original
#'   predictions.}
#' \item{after}{list of calibration statistics computed on the recalibrated
#'   predictions.}
#' \item{transform}{for \code{"platt"}, a named vector
#'   \code{c(intercept, slope)}; for \code{"isotonic"}, the fitted
#'   \code{isoreg} object.}
#' \item{ggPlot}{list with ggplot objects (\code{calibration} and
#'   \code{transform}) when \code{plot = TRUE}.}
#'
#' @references
#' Platt J C (1999). Probabilistic outputs for support vector machines and
#' comparisons to regularized likelihood methods. In \emph{Advances in Large
#' Margin Classifiers}, MIT Press.
#'
#' Fawcett T and Niculescu-Mizil A (2007). PAV and the ROC convex hull.
#' \emph{Machine Learning}, \bold{68(1)}, pp. 97--106.
#'
#' Van Calster B, McLernon D J, van Smeden M, et al. (2019). Calibration:
#' the Achilles heel of predictive analytics. \emph{BMC Medicine},
#' \bold{17}, 230.
#'
#' @examples
#' library(CalibrationCurves)
#' data("traindata")
#' data("testdata")
#' glmFit <- glm(y ~ ., data = traindata, family = binomial)
#' p_test <- predict(glmFit, newdata = testdata, type = "response")
#' y_test <- testdata$y
#' # Recalibrate using Platt scaling
#' recal <- recalibrate(p_test, y_test, method = "platt")
#' recal
#' # Recalibrate using isotonic regression
#' recal_iso <- recalibrate(p_test, y_test, method = "isotonic")
#' @export
recalibrate <- function(p,
                        y,
                        method     = c("platt", "isotonic"),
                        plot       = TRUE,
                        plottype   = c("both", "calibration", "transform"),
                        alpha      = 0.05,
                        smooth     = c("loess", "rcs"),
                        nr.knots   = 3,
                        xlab       = "Predicted probability",
                        ylab       = "Observed proportion",
                        col.before = "steelblue",
                        col.after  = "darkgreen",
                        verbose    = TRUE) {

  callFn   <- match.call()
  method   <- match.arg(method)
  plottype <- match.arg(plottype)
  smooth   <- match.arg(smooth)

  # Input validation
  if (!is.numeric(p) || any(p < 0 | p > 1, na.rm = TRUE))
    stop("'p' must be a numeric vector of probabilities in [0, 1].")
  if (!all(y %in% c(0, 1)))
    stop("'y' must be a binary vector with values 0 and 1.")
  if (length(p) != length(y))
    stop("'p' and 'y' must have the same length.")

  n   <- length(p)
  lp  <- qlogis(pmax(pmin(p, 1 - 1e-8), 1e-8))

  # --------------------------------------------------------------------------
  # Helper: compute calibration summary stats
  # --------------------------------------------------------------------------
  .cal_stats <- function(phat, ytrue) {
    eps      <- 1e-15
    brier    <- mean((phat - ytrue)^2)
    logloss  <- -mean(ytrue * log(pmax(phat, eps)) +
                      (1 - ytrue) * log(pmax(1 - phat, eps)))
    lp_hat   <- qlogis(pmax(pmin(phat, 1 - 1e-8), 1e-8))
    fit_cal  <- tryCatch(glm(ytrue ~ lp_hat, family = binomial),
                         error = function(e) NULL)
    slope <- intercept <- NA
    if (!is.null(fit_cal)) {
      slope     <- unname(coef(fit_cal)[2])
      intercept <- unname(coef(fit_cal)[1])
    }
    list(Brier     = brier,
         LogLoss   = logloss,
         Intercept = intercept,
         Slope     = slope)
  }

  # --------------------------------------------------------------------------
  # Helper: flexible calibration curve data
  # --------------------------------------------------------------------------
  .cal_curve <- function(phat, ytrue) {
    if (smooth == "loess") {
      tryCatch({
        lo  <- loess(ytrue ~ phat, degree = 2, span = 0.75)
        od  <- order(phat)
        data.frame(pred = phat[od], obs = fitted(lo)[od])
      }, error = function(e) NULL)
    } else {
      tryCatch({
        lp_h  <- qlogis(pmax(pmin(phat, 1 - 1e-8), 1e-8))
        kn    <- min(nr.knots, 5)
        lp_rcs <- Hmisc::rcspline.eval(lp_h, knots = kn, inclx = TRUE)
        fit_r  <- glm(ytrue ~ lp_rcs, family = binomial)
        od     <- order(phat)
        data.frame(pred = phat[od], obs = fitted(fit_r)[od])
      }, error = function(e) NULL)
    }
  }

  # --------------------------------------------------------------------------
  # Before stats
  # --------------------------------------------------------------------------
  before <- .cal_stats(p, y)
  curve_before <- .cal_curve(p, y)

  # --------------------------------------------------------------------------
  # Recalibration
  # --------------------------------------------------------------------------
  transform_obj <- NULL

  if (method == "platt") {
    fit_platt  <- glm(y ~ lp, family = binomial)
    p_recal    <- fitted(fit_platt)
    transform_obj <- c(intercept = unname(coef(fit_platt)[1]),
                       slope     = unname(coef(fit_platt)[2]))
    if (verbose) {
      message("Platt scaling: intercept = ",
              round(transform_obj["intercept"], 4),
              ", slope = ", round(transform_obj["slope"], 4))
    }
  } else {
    # Isotonic regression via pool-adjacent-violators
    # isoreg() requires x to be ordered; sort then build an approxfun
    ord        <- order(p)
    p_s        <- p[ord]
    y_s        <- as.numeric(y)[ord]
    iso_fit    <- isoreg(p_s, y_s)
    # approxfun with method="constant" handles ties in p robustly
    p_recal_fn <- approxfun(
      x      = p_s,
      y      = iso_fit$yf,
      method = "constant",
      yleft  = iso_fit$yf[1L],
      yright = iso_fit$yf[length(iso_fit$yf)],
      ties   = "ordered"
    )
    p_recal    <- pmin(pmax(p_recal_fn(p), 0), 1)
    transform_obj <- iso_fit
    if (verbose)
      message("Isotonic regression applied (pool adjacent violators).")
  }

  # --------------------------------------------------------------------------
  # After stats
  # --------------------------------------------------------------------------
  after <- .cal_stats(p_recal, y)
  curve_after <- .cal_curve(p_recal, y)

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  gg_cal <- gg_trans <- NULL

  if (plot && (plottype %in% c("both", "calibration"))) {
    plot_df <- rbind(
      if (!is.null(curve_before))
        cbind(curve_before, Group = "Before recalibration") else NULL,
      if (!is.null(curve_after))
        cbind(curve_after,  Group = "After recalibration")  else NULL
    )
    if (!is.null(plot_df)) {
      gg_cal <- ggplot(plot_df, aes(x = pred, y = obs, color = Group)) +
        geom_abline(slope = 1, intercept = 0, color = "red",
                    linetype = 1, linewidth = 0.8) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = c(
          "Before recalibration" = col.before,
          "After recalibration"  = col.after
        )) +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = xlab, y = ylab,
             title = "Calibration curves before and after recalibration",
             color = NULL) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position  = "bottom",
              axis.text        = element_text(size = 11),
              axis.title       = element_text(size = 13))
      print(gg_cal)
    }
  }

  if (plot && (plottype %in% c("both", "transform"))) {
    trans_df <- data.frame(p_orig = sort(p),
                           p_recal_sorted = p_recal[order(p)])

    gg_trans <- ggplot(trans_df, aes(x = p_orig, y = p_recal_sorted)) +
      geom_abline(slope = 1, intercept = 0, color = "grey60",
                  linetype = 2, linewidth = 0.7) +
      (if (method == "platt") {
        geom_line(color = col.after, linewidth = 1)
      } else {
        geom_step(color = col.after, linewidth = 1)
      }) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Original predicted probability",
           y = "Recalibrated predicted probability",
           title = paste(
             if (method == "platt") "Platt scaling" else "Isotonic regression",
             "transformation")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 11),
            axis.title       = element_text(size = 13))
    if (plottype == "both") {
      # Already printed cal plot; print transform now
      print(gg_trans)
    } else {
      print(gg_trans)
    }
  }

  # --------------------------------------------------------------------------
  # Return
  # --------------------------------------------------------------------------
  Results <- structure(
    list(
      call       = callFn,
      method     = method,
      p_original = p,
      p_recal    = p_recal,
      before     = before,
      after      = after,
      transform  = transform_obj,
      ggPlot     = list(calibration = gg_cal, transform = gg_trans)
    ),
    class = "RecalibratedPredictions"
  )
  return(Results)
}
