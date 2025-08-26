#' Meta-Analytical Calibration Curve (MAC2)
#'
#' Computes meta-analytical calibration curves using multiple methods (logistic regression,
#' loess, splines, or kernel density estimation) and performs meta-analysis across clusters
#' to generate aggregated calibration curves with confidence and prediction intervals.
#'
#' @param data Optional data frame containing the columns for `preds`, `y`, and `cluster`.
#'   If provided, `preds`, `y`, and `cluster` should be column names (unquoted).
#'   Default is `NULL`.
#' @param preds A numeric vector of predicted probabilities.
#' @param y A numeric vector of binary outcomes (0 or 1).
#' @param cluster A factor or character vector identifying cluster memberships.
#' @param grid_length Integer; length of the grid for calibration curve evaluation. Default is `100`.
#' @param methods Character vector; methods to use for calibration. Options are:
#'   `"log"` (logistic regression), `"loess"`, `"splines"`, and `"kde"` (kernel density estimation).
#' @param plot Logical; whether to plot the calibration curves. Default is `TRUE`.
#' @param cluster_curves Logical; whether to include cluster-specific curves in the plot. Default is `FALSE`.
#' @param knots Integer; number of knots for splines. Default is `3`.
#' @param transf Character; transformation for predictions: `"logit"` or `"identity"`. Default is `"logit"`.
#' @param method_choice Character; which method to use for meta-analysis. Options are:
#'   `"log"`, `"loess"`, `"splines"`, or `"kde"`. Default is `"splines"`.
#' @param method.tau Character; method for between-study heterogeneity estimation. Default is `"REML"`.
#' @param prediction Logical; whether to compute prediction intervals. Default is `TRUE`.
#' @param random Logical; whether to use random-effects model. Default is `TRUE`.
#' @param sm Character; summary measure for meta-analysis. Default is `"PLOGIT"`.
#' @param hakn Logical; whether to use Hartung-Knapp adjustment. Default is `FALSE`.
#' @param linewidth Numeric; line width for the meta-curve. Default is `1`.
#' @param method.predict Character; method for prediction intervals. Default is `"HTS"`.
#'
#' @details
#' This function calculates calibration curves for multiple methods and aggregates them
#' using meta-analysis. The `method_choice` argument determines which method is used
#' for the meta-analytical aggregation.
#'
#' @return A list containing:
#' \describe{
#'   \item{cluster_data}{Data frame with linear predictors and standard errors for each method per cluster}
#'   \item{plot_data}{Data frame with meta-analysis results including predictions and intervals}
#'   \item{plot}{A `ggplot2` object if `plot = TRUE`, otherwise `NULL`}
#' }
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   preds = runif(100),
#'   y = rbinom(100, 1, 0.5),
#'   cluster = sample(1:5, 100, replace = TRUE)
#' )
#' result <- MAC2(
#'   data = data, preds = preds, y = y, cluster = cluster,
#'   methods = c("log", "splines"), method_choice = "splines"
#' )
#' result$plot
#'
#' @importFrom dplyr filter mutate group_by ungroup
#' @importFrom meta metagen
#' @importFrom ggplot2 ggplot geom_abline geom_ribbon geom_line xlab ylab theme_classic
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_alpha_manual
#' @importFrom ggplot2 scale_fill_manual coord_cartesian theme unit
#' @importFrom rms lrm rcs
#' @importFrom fANCOVA loess.as
#' @importFrom zoo na.approx
#' @importFrom stats plogis qlogis
#'
#' @export
MAC2 <- function(data = NULL,
                 preds,
                 y,
                 cluster,
                 grid_length = 100,
                 methods = c("log", "loess", "splines", "kde"),
                 plot = TRUE,
                 cluster_curves = FALSE,
                 knots = 3,
                 transf = "logit",
                 method_choice = "splines",
                 method.tau = "REML",
                 prediction = TRUE,
                 random = TRUE,
                 sm = "PLOGIT",
                 hakn = FALSE,
                 linewidth = 1,
                 method.predict = "HTS") {
  # --- Extract from data if provided ---
  if (!is.null(data)) {
    preds <- data[[deparse(substitute(preds))]]
    y <- data[[deparse(substitute(y))]]
    cluster <- data[[deparse(substitute(cluster))]]
  }

  # --- Base dataframe ---
  df <- data.frame(
    predictions = as.numeric(preds),
    outcome = as.numeric(y),
    cluster = as.factor(cluster)
  )

  # --- Grid computation ---
  grid <- seq(0.01, 0.99, length.out = grid_length)
  transform_function <- if (transf == "logit") qlogis else identity
  data_all_lp <- data.frame()

  # --- Process each cluster ---
  for (subcluster in unique(df$cluster)) {
    risk_cluster <- df %>% filter(cluster == subcluster)
    observed_grid <- data.frame(
      x = grid,
      cluster = subcluster,
      nsample = nrow(risk_cluster)
    )

    risk_cluster$transf_preds <- transform_function(risk_cluster$predictions)

    # --- Logistic regression method ---
    if ("log" %in% methods) {
      log_model <- rms::lrm(data = risk_cluster, outcome ~ transf_preds)
      log_data <- predict(log_model,
        newdata = data.frame(transf_preds = transform_function(grid)),
        type = "lp",
        se.fit = TRUE
      )
      log_data <- data.frame(
        log = log_data$linear.predictors,
        log_se = log_data$se.fit
      )
      observed_grid <- cbind(observed_grid, log_data)
    }

    # --- Loess method ---
    if ("loess" %in% methods) {
      tryCatch(
        {
          loess_model <- fANCOVA::loess.as(
            x = risk_cluster$transf_preds,
            y = risk_cluster$outcome,
            degree = 2,
            criterion = "aicc",
            plot = FALSE,
            control = loess.control(surface = "direct")
          )

          loess_data <- predict(loess_model,
            newdata = data.frame(x = transform_function(grid)),
            se = TRUE
          )

          if (any(is.na(loess_data$fit)) || any(is.na(loess_data$se.fit))) {
            loess_data$fit <- zoo::na.approx(loess_data$fit, rule = 2)
            loess_data$se.fit <- zoo::na.approx(loess_data$se.fit, rule = 2)
          }

          loess_data <- data.frame(
            loess = transform_function(loess_data$fit),
            loess_se = abs(loess_data$se.fit / (loess_data$fit * (1 - loess_data$fit)))
          )
          observed_grid <- cbind(observed_grid, loess_data)
        },
        error = function(e) {
          message("LOESS was not computed because: ", e$message)
        }
      )
    }

    # --- Splines method ---
    if ("splines" %in% methods) {
      knots_sub <- knots
      splines_model <- suppressWarnings(
        rms::lrm(data = risk_cluster, outcome ~ rcs(transf_preds, knots_sub))
      )

      # --- Model selection logic ---
      while (splines_model$fail && knots_sub != 3) {
        knots_sub <- knots_sub - 1
        splines_model <- suppressWarnings(
          rms::lrm(data = risk_cluster, outcome ~ rcs(transf_preds, knots_sub))
        )
      }

      if (knots_sub > 3) {
        splines_model3 <- suppressWarnings(
          rms::lrm(data = risk_cluster, outcome ~ rcs(transf_preds, 3))
        )
        test <- rms::lrtest(splines_model3, splines_model)

        if (test$stats["P"] > 0.05) {
          if (knots_sub > 4) {
            splines_model4 <- suppressWarnings(
              rms::lrm(data = risk_cluster, outcome ~ rcs(transf_preds, 4))
            )
            if (!splines_model4$fail) {
              test <- rms::lrtest(splines_model4, splines_model3)
            } else {
              test$stats["P"] <- 0
            }

            if (test$stats["P"] > 0.05) {
              splines_model <- splines_model4
              knots_sub <- 4
            } else {
              splines_model <- splines_model3
              knots_sub <- 3
            }
          } else {
            splines_model <- splines_model3
            knots_sub <- 3
          }
        }
      } else {
        splines_model <- suppressWarnings(
          rms::lrm(data = risk_cluster, outcome ~ transf_preds)
        )
        knots_sub <- 1 # effectively no spline
      }

      splines_data <- predict(splines_model,
        newdata = data.frame(transf_preds = transform_function(grid)),
        type = "lp",
        se.fit = TRUE
      )

      splines_data <- data.frame(
        splines = splines_data$linear.predictors,
        splines_se = splines_data$se.fit,
        knots_used = knots_sub
      )
      observed_grid <- cbind(observed_grid, splines_data)
      message("Spline model for cluster ", subcluster, " fitted with ", knots_sub, " knots.")
    }

    # --- KDE method ---
    # if ("kde" %in% methods) {
    #   tryCatch(
    #     {
    #       kde_data <- kde_curve(risk_cluster$outcome, risk_cluster$predictions)
    #       observed_grid <- cbind(observed_grid, kde = kde_data$y)
    #     },
    #     error = function(e) {
    #       message("KDE was not computed because: ", e$message)
    #     }
    #   )
    # }

    data_all_lp <- rbind(data_all_lp, observed_grid)
  }

  # --- Meta-analysis ---
  adhoc.hakn.ci <- if (hakn) "IQWiG6" else NULL
  adhoc.hakn.pi <- if (hakn) "se" else NULL
  curve <- data.frame()

  for (value in unique(data_all_lp$x)) {
    data_v <- data_all_lp %>% filter(x == value)

    meta_inputs <- switch(method_choice,
      "log" = list(TE = data_v$log, seTE = data_v$log_se),
      "loess" = list(TE = data_v$loess, seTE = data_v$loess_se),
      "splines" = list(TE = data_v$splines, seTE = data_v$splines_se),
      # "kde" = list(TE = data_v$kde, seTE = NA),
      stop("Invalid method choice: ", method_choice)
    )

    meta <- meta::metagen(
      TE = meta_inputs$TE,
      seTE = meta_inputs$seTE,
      studlab = data_v$cluster,
      random = random,
      prediction = prediction,
      method.tau = method.tau,
      sm = sm,
      backtransf = TRUE,
      method.random.ci = hakn,
      method.predict = method.predict
    )

    data_v_b <- data.frame(
      y = plogis(meta$TE.random),
      upper = plogis(meta$upper.random),
      lower = plogis(meta$lower.random),
      up_pre = plogis(meta$upper.predict),
      low_pre = plogis(meta$lower.predict),
      tau = meta$tau2,
      x = value
    )
    curve <- rbind(curve, data_v_b)
  }

  # --- Plotting ---
  plot_obj <- NULL
  if (plot) {
    plot_obj <- ggplot(data = curve) +
      geom_abline(linetype = "dashed", alpha = 0.1) +
      geom_ribbon(aes(
        x = x, ymin = pmax(0, pmin(lower, 1)),
        ymax = pmax(0, pmin(upper, 1)),
        fill = "CI 95%", alpha = "CI 95%"
      )) +
      geom_ribbon(aes(
        x = x, ymin = pmax(0, pmin(low_pre, 1)),
        ymax = pmax(0, pmin(up_pre, 1)),
        fill = "PI 95%", alpha = "PI 95%"
      )) +
      geom_line(aes(x = x, y = y),
        color = "black",
        linewidth = linewidth, linetype = "dashed"
      ) +
      xlab("Estimated probability") +
      ylab("Observed proportion") +
      theme_classic(base_size = 8, base_family = "serif") +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      scale_alpha_manual(values = c(0.4, 0.2), name = "Heterogeneity") +
      scale_fill_manual(values = c("red", "red"), name = "Heterogeneity") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme(
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (cluster_curves) {
      plot_obj <- plot_obj +
        geom_line(
          data = data_all_lp,
          aes(x = x, y = plogis(data_all_lp[, method_choice]), group = cluster),
          linewidth = linewidth / 2,
          show.legend = FALSE,
          linetype = "solid"
        )
    }
  }

  # --- Return ---
  list(
    cluster_data = data_all_lp,
    plot_data = curve,
    plot = plot_obj
  )
}
