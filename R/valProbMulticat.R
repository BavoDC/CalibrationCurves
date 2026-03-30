#' Plot calibration curves for ordinal and multinomial prediction models
#'
#' Evaluates calibration performance of prediction models for ordered or
#' nominal categorical outcomes, following Van Hoorde et al. (2012, 2015)
#' and Van Calster et al. (2019). For each outcome category a calibration
#' curve is produced.
#'
#' @param fit the fitted model.  Supported classes:
#'   \itemize{
#'     \item \code{polr} from \pkg{MASS} (proportional-odds / ordinal outcome);
#'     \item \code{multinom} from \pkg{nnet} (nominal / multinomial outcome);
#'     \item a named list of predicted probability matrices (see Details).
#'   }
#' @param valdata the validation data frame containing the response and the
#'   predictors used to fit the model.
#' @param y character string naming the response column in \code{valdata}, or a
#'   factor/integer vector.  Required when \code{fit} is a probability matrix.
#' @param alpha numeric, significance level for confidence intervals. Default
#'   \code{0.05}.
#' @param plot logical, indicates whether calibration curves should be plotted.
#'   Default \code{TRUE}.
#' @param plotCal character, \code{"base"} for base-R graphics, \code{"ggplot"}
#'   for \pkg{ggplot2}. Default \code{"ggplot"}.
#' @param smooth character, \code{"loess"} or \code{"rcs"} for the smoother
#'   used to estimate individual category curves. Default \code{"loess"}.
#' @param nr.knots number of knots when \code{smooth = "rcs"}. Default \code{3}.
#' @param xlab,ylab axis labels.
#' @param xlim,ylim axis limits.  Default \code{c(0,1)} for both.
#' @param nrow integer, number of rows in the faceted ggplot. Default
#'   \code{NULL} (determined automatically).
#' @param cl.level confidence level for calibration slope/intercept CIs. Default
#'   \code{0.95}.
#'
#' @details
#' \strong{Ordinal outcome (polr)}
#'
#' For an ordered response with \eqn{K} categories, the model produces
#' \eqn{K-1} cumulative probabilities \eqn{\hat{F}_k = P(Y \le k)}.  Per
#' category, calibration is assessed by plotting observed cumulative frequencies
#' against predicted cumulative probabilities (Van Hoorde et al., 2015).  The
#' overall calibration slope and intercept are obtained from a proportional-odds
#' regression of the actual outcome on the logit of the predicted probabilities.
#'
#' \strong{Multinomial outcome (multinom)}
#'
#' For a nominal response, one-versus-rest binary calibration curves are produced
#' for each category (Van Hoorde et al., 2012).  The calibration slope is
#' estimated via logistic regression of the binary indicator for each category on
#' the logit of the corresponding predicted probability. A multivariate Brier score
#' is returned as the primary proper scoring rule.
#'
#' \strong{Proper scoring rules returned:}
#' \itemize{
#'   \item Brier score (both overall multiclass and per-category);
#'   \item Log-loss (negative log-likelihood);
#'   \item Mean Squared Error of calibration (MSEC) per category.
#' }
#'
#' @return An object of class \code{MulticlassCalibrationCurve} with:
#' \item{call}{the matched call.}
#' \item{type}{\code{"ordinal"} or \code{"multinomial"}.}
#' \item{stats}{list of scoring rules and calibration statistics.}
#' \item{Calibration}{per-category calibration slope and intercept.}
#' \item{CalibrationCurves}{data frames of curve coordinates, one per category.}
#' \item{ggPlot}{ggplot object when \code{plotCal = "ggplot"}.}
#'
#' @references
#' Van Hoorde K, Vergouwe Y, Timmerman D, Van Huffel S, Steyerberg E W,
#' Van Calster B (2014). Simple polytomous logistic regression calibration for
#' external validation of a multinomial clinical prediction model.
#' \emph{Statistics in Medicine}, \bold{33(2)}, pp. 203--216.
#'
#' Van Hoorde K, Van Huffel S, Timmerman D, Bourne T, Van Calster B (2015).
#' A continuous calibration belt for logistic calibration assessment.
#' \emph{Statistics in Medicine}, \bold{34(10)}, pp. 1598--1617.
#'
#' Van Calster B, McLernon D J, van Smeden M, et al. (2019). Calibration: the
#' Achilles heel of predictive analytics. \emph{BMC Medicine}, \bold{17}, 230.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' library(CalibrationCurves)
#' # Ordinal example using Wine dataset
#' data(Pima.tr, package = "MASS")
#' data(Pima.te, package = "MASS")
#' oPima <- polr(type ~ ., data = Pima.tr, Hess = TRUE)
#' calOrd <- valProbMulticat(oPima, Pima.te)
#' calOrd
#' }
#' @export
valProbMulticat <- function(fit,
                            valdata  = NULL,
                            y        = NULL,
                            alpha    = 0.05,
                            plot     = TRUE,
                            plotCal  = c("ggplot", "base", "none"),
                            smooth   = c("loess", "rcs"),
                            nr.knots = 3,
                            xlab     = "Predicted probability",
                            ylab     = "Observed proportion",
                            xlim     = c(0, 1),
                            ylim     = c(0, 1),
                            nrow     = NULL,
                            cl.level = 0.95) {

  callFn  <- match.call()
  plotCal <- match.arg(plotCal)
  smooth  <- match.arg(smooth)
  z       <- qnorm(1 - alpha / 2)

  # --------------------------------------------------------------------------
  # Determine type and extract predictions + outcome
  # --------------------------------------------------------------------------
  if (inherits(fit, "polr")) {
    type     <- "ordinal"
    pred_mat <- predict(fit, newdata = valdata, type = "probs")
    if (is.null(valdata))
      stop("'valdata' must be provided.")
    yRespName <- as.character(fit$call$formula[[2]])
    yObs      <- valdata[[yRespName]]
  } else if (inherits(fit, "multinom")) {
    type      <- "multinomial"
    pred_mat  <- predict(fit, newdata = valdata, type = "probs")
    if (is.null(valdata))
      stop("'valdata' must be provided.")
    yRespName <- as.character(fit$call$formula[[2]])
    yObs      <- valdata[[yRespName]]
  } else if (is.matrix(fit) || is.data.frame(fit)) {
    # User supplies a matrix of predicted probabilities directly
    pred_mat  <- as.matrix(fit)
    if (is.null(y))
      stop("'y' must be provided when 'fit' is a probability matrix.")
    yObs      <- y
    type      <- "multinomial"
  } else {
    stop("'fit' must be of class 'polr', 'multinom', or a probability matrix.")
  }

  if (is.null(dim(pred_mat)) || ncol(pred_mat) < 2)
    stop("Prediction matrix must have at least 2 columns (categories).")

  if (!is.factor(yObs))
    yObs <- factor(yObs)

  cats     <- levels(yObs)
  K        <- length(cats)
  n        <- nrow(pred_mat)
  colnames(pred_mat) <- cats

  # --------------------------------------------------------------------------
  # Proper scoring rules
  # --------------------------------------------------------------------------
  # Indicator matrix
  Y_ind <- model.matrix(~ yObs - 1)
  colnames(Y_ind) <- cats

  # Multiclass Brier score  (avg over obs and categories)
  brier_overall <- mean(rowSums((pred_mat - Y_ind)^2))

  # Per-category Brier score
  brier_cat <- colMeans((pred_mat - Y_ind)^2)

  # Log-loss (negative mean log-likelihood)
  eps       <- 1e-15
  logloss   <- -mean(rowSums(Y_ind * log(pmax(pred_mat, eps))))

  stats <- list(
    BrierScore       = c("Overall" = brier_overall, brier_cat),
    LogLoss          = logloss
  )

  # --------------------------------------------------------------------------
  # Per-category calibration curves and slope/intercept
  # --------------------------------------------------------------------------
  calCurves   <- vector("list", K)
  names(calCurves) <- cats
  calSlopes   <- matrix(NA, nrow = K, ncol = 3,
                        dimnames = list(cats,
                                        c("Slope", paste0(c((alpha/2)*100, (1-alpha/2)*100), " %"))))
  calIntercepts <- calSlopes

  for (k in seq_len(K)) {
    pk    <- pred_mat[, k]
    yk    <- as.integer(yObs == cats[k])

    # Logit of predicted probability
    lp_k  <- qlogis(pmax(pmin(pk, 1 - 1e-8), 1e-8))

    # Calibration slope and intercept (logistic regression one-vs-rest)
    fitCal <- tryCatch(
      glm(yk ~ lp_k, family = binomial),
      error = function(e) NULL
    )
    if (!is.null(fitCal)) {
      ci_slope <- confint.default(fitCal, level = cl.level)
      calSlopes[k, ]     <- c(coef(fitCal)[2], ci_slope[2, ])
      calIntercepts[k, ] <- c(coef(fitCal)[1], confint.default(fitCal, level = cl.level)[1, ])
    }

    # ICI / E50 / E90 per category
    obs_fit <- if (smooth == "loess") {
      tryCatch({
        lo   <- loess(yk ~ pk, degree = 2)
        data.frame(pred = pk, obs = fitted(lo))
      }, error = function(e) NULL)
    } else {
      tryCatch({
        dd  <- data.frame(yk = yk, p = pk)
        lp  <- Hmisc::rcspline.eval(lp_k, knots = nr.knots, inclx = TRUE)
        glm_rcs <- glm(yk ~ lp, family = binomial, data = dd)
        data.frame(pred = pk, obs = fitted(glm_rcs))
      }, error = function(e) NULL)
    }

    if (!is.null(obs_fit)) {
      obs_fit         <- obs_fit[order(obs_fit$pred), ]
      absdiff_k       <- abs(obs_fit$pred - obs_fit$obs)
      stats$MSEC[[cats[k]]] <- c(
        "ICI"  = mean(absdiff_k),
        setNames(quantile(absdiff_k, c(0.5, 0.9)), c("E50", "E90")),
        "Emax" = max(absdiff_k)
      )
      calCurves[[k]] <- obs_fit
    }
  }

  stats$CalibrationSlopes     <- calSlopes
  stats$CalibrationIntercepts <- calIntercepts

  # --------------------------------------------------------------------------
  # Plotting
  # --------------------------------------------------------------------------
  gg <- NULL
  if (plotCal != "none" && plot) {

    # Build a long data frame for ggplot
    plot_data <- do.call(rbind, lapply(cats, function(cat) {
      cd <- calCurves[[cat]]
      if (is.null(cd)) return(NULL)
      cbind(cd, Category = cat)
    }))

    if (plotCal == "ggplot" && !is.null(plot_data)) {
      gg <- ggplot(plot_data, aes(x = pred, y = obs, color = Category)) +
        geom_abline(slope = 1, intercept = 0,
                    linetype = 1, color = "red", linewidth = 0.8) +
        geom_line(linewidth = 1) +
        facet_wrap(~ Category, nrow = nrow) +
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(limits = ylim) +
        labs(x = xlab, y = ylab,
             title = paste0(
               if (type == "ordinal") "Ordinal" else "Multinomial",
               " calibration curves")) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position  = "none",
          strip.background = element_rect(fill = "grey95"),
          axis.text        = element_text(size = 11),
          axis.title       = element_text(size = 13)
        )
      print(gg)

    } else if (plotCal == "base" && !is.null(plot_data)) {
      old_par <- par(mfrow = c(ceiling(K / 2), 2), mar = c(4, 4, 2, 1))
      on.exit(par(old_par), add = TRUE)
      for (cat in cats) {
        cd <- calCurves[[cat]]
        if (is.null(cd)) next
        plot(cd$pred, cd$obs, type = "l",
             xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab,
             main = paste("Category:", cat),
             lwd = 2)
        abline(0, 1, col = "red", lty = 1, lwd = 1)
      }
    }
  }

  # --------------------------------------------------------------------------
  # Return
  # --------------------------------------------------------------------------
  Results <- structure(
    list(
      call              = callFn,
      type              = type,
      stats             = stats,
      alpha             = alpha,
      cl.level          = cl.level,
      Calibration       = list(
        Slopes     = calSlopes,
        Intercepts = calIntercepts
      ),
      CalibrationCurves = calCurves
    ),
    class = "MulticlassCalibrationCurve"
  )
  if (plotCal == "ggplot")
    Results$ggPlot <- gg
  return(Results)
}
