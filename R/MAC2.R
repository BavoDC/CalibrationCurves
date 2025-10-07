#' Internal function for the Meta-Analytical Calibration Curve (MAC2)
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
#' @param grid the grid for the calibration curve evaluation
#' @param methods Character vector; methods to use for calibration. Options are:
#'   `"log"` (logistic regression), `"loess"` and `"splines"`.
#' @param plot Logical; whether to plot the calibration curves. Default is `TRUE`.
#' @param cluster_curves Logical; whether to include cluster-specific curves in the plot. Default is `FALSE`.
#' @param knots Integer; number of knots for splines. Default is `3`.
#' @param transf Character; transformation for predictions: `"logit"` or `"identity"`. Default is `"logit"`.
#' @param method_choice Character; which method to use for meta-analysis. Options are:
#'   `"log"`, `"loess"` or `"splines"`. Default is `"splines"`.
#' @param method.tau Character; method for between-study heterogeneity estimation. Default is `"REML"`.
#' @param prediction Logical; whether to compute prediction intervals. Default is `TRUE`.
#' @param random Logical; whether to use random-effects model. Default is `TRUE`.
#' @param sm Character; summary measure for meta-analysis. Default is `"PLOGIT"`.
#' @param hakn Logical; whether to use Hartung-Knapp adjustment. Default is `FALSE`.
#' @param linewidth Numeric; line width for the meta-curve. Default is `1`.
#' @param method.predict Character; method for prediction intervals. Default is `"HTS"`.
#' @param verbose logical, indicates whether progress has to be printed in the console.
#' @param cl.level the confidence level for the calculation of the confidence interval. Default is \code{0.95}.
#' @param alpha.lr the alpha-level used for the likelihood ratio test, selecting the number of knots for the
#' restricted cubic splines
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
#'
MAC2 <- function(data = NULL,
                 preds,
                 y,
                 cluster,
                 grid,
                 methods = c("log", "loess", "splines"),
                 cl.level = 0.95,
                 alpha.lr = 0.05 / 3,
                 plot = TRUE,
                 cluster_curves = FALSE,
                 knots = 3,
                 transf = "logit",
                 method_choice = c("splines", "log", "loess"),
                 method.tau = "REML",
                 prediction = TRUE,
                 random = TRUE,
                 sm = "PLOGIT",
                 hakn = FALSE,
                 linewidth = 1,
                 method.predict = "HTS",
                 verbose = FALSE) {
  # --- Extract from data if provided ---
  callFn        = match.call()
  method_choice = match.arg(method_choice)

  if (!is.null(data)) {
    if(!all(sapply(c("preds", "y", "cluster"), function(a) as.character(callFn[a])) %in% colnames(data)))
      stop(paste("Variables", paste0(
        callFn[c("preds", "y", "cluster")], collapse = ", "
      ), "were not found in the data.frame."))
    preds   = eval(callFn$preds, data)
    logit   = Logit(preds)
    y       = eval(callFn$y, data)
    cluster = eval(callFn$cluster, data)
  }

  # --- Base dataframe ---
  df <- data.frame(
    predictions = as.numeric(preds),
    outcome = as.numeric(y),
    cluster = as.factor(cluster)
  )

  # --- Grid computation ---
  transform_function <- if (transf == "logit") Logit else identity
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
      log_model <- lrm(data = risk_cluster, outcome ~ transf_preds)
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
          loess_model <- loess.as(
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



          # Handle NA values in loess predictions
          # To-do: why do these result in NA?!
          if (any(is.na(loess_data$fit)) || any(is.na(loess_data$se.fit))) {
            loess_data$fit    = na.approx(loess_data$fit, rule = 2)
            loess_data$se.fit = na.approx(loess_data$se.fit, rule = 2)
          }
          loess_data$fit <- ifelse(loess_data$fit >= 1, 1 - 1e-6, loess_data$fit)
          loess_data$fit <- ifelse(loess_data$fit <= 0, 1e-6, loess_data$fit)
          loess_data <- data.frame(
            loess = transform_function(loess_data$fit),
            loess_se = abs(loess_data$se.fit / (loess_data$fit * (1 - loess_data$fit)))
          )
          observed_grid <- cbind(observed_grid, loess_data)
        },
        error = function(e) {
          message("Fitting loess resulted in the following error message: ", e$message)
        }
      )
    }

    # --- Splines method ---
    if ("splines" %in% methods) {
      nkDecrease <- function(Argz) {
        tryCatch(
          do.call("lrm", Argz),
          error = function(e) {
            nk = Argz$formula[[3]][[3]]
            warning(paste0("The number of knots led to estimation problems, nk will be set to ", nk - 1), immediate. = TRUE)
            nk = nk - 1
            cat(paste("fitting with", nk, "knots"))
            Argz = list(
              formula = eval(substitute(outcome ~ rcs(transf_preds, k), list(k = nk))),
              data    = risk_cluster
            )
            nkDecrease(Argz)
          },
          warning = function(w) {
            nk = Argz$formula[[3]][[3]]
            warning(paste0("The number of knots led to estimation problems, nk will be set to ", nk - 1), immediate. = TRUE)
            nk = nk - 1
            cat(paste("fitting with", nk, "knots"))
            Argz = list(
              formula = eval(substitute(outcome ~ rcs(transf_preds, k), list(k = nk))),
              data    = risk_cluster
            )
            nkDecrease(Argz)
          })
      }
      knots_sub   = knots
      argzSplines = list(
        formula = eval(substitute(outcome ~ rcs(transf_preds, k), list(k = knots_sub))),
        data    = risk_cluster
      )
      splines_model = nkDecrease(argzSplines)
      knots_sub     = splines_model$sformula[[3]][[3]]

      if (knots_sub > 3) {
        splines_model3 = lrm(data = risk_cluster, outcome ~ rcs(transf_preds, 3))
        splines_model4 = lrm(data = risk_cluster, outcome ~ rcs(transf_preds, 4))
        splines_model5 = lrm(data = risk_cluster, outcome ~ rcs(transf_preds, 5))
        test3          = lrtest(splines_model3, splines_model)
        test4          = lrtest(splines_model4, splines_model)
        test5          = lrtest(splines_model5, splines_model)

        if(test3$stats["P"] > alpha.lr) {
          splines_model = splines_model3
          knots_sub     = 3
        } else if(test4$stats["P"] > alpha.lr) {
          splines_model = splines_model4
          knots_sub     = 4
        } else if(test5$stats["P"] > alpha.lr) {
          splines_model = splines_model5
          knots_sub     = 5
        }
      }

      splines_data <- predict(splines_model,
        newdata = data.frame(transf_preds = transform_function(grid)),
        type = "lp",
        se.fit = TRUE
      )

      splines_data <- data.frame(
        splines    = splines_data$linear.predictors,
        splines_se = splines_data$se.fit,
        knots_used = knots_sub
      )
      observed_grid <- cbind(observed_grid, splines_data)
      if (verbose) {
        message("Spline model for cluster ", subcluster, " fitted with ", knots_sub, " knots.")
      }
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
      "log"     = list(TE = data_v$log, seTE = data_v$log_se),
      "loess"   = list(TE = data_v$loess, seTE = data_v$loess_se),
      "splines" = list(TE = data_v$splines, seTE = data_v$splines_se),
      # "kde" = list(TE = data_v$kde, seTE = NA),
      stop("Invalid method choice: ", method_choice)
    )

    meta <- metagen(
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
      y = Ilogit(meta$TE.random),
      upper = Ilogit(meta$upper.random),
      lower = Ilogit(meta$lower.random),
      up_pre = Ilogit(meta$upper.predict),
      low_pre = Ilogit(meta$lower.predict),
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
          aes(x = x, y = Ilogit(data_all_lp[, method_choice]), group = cluster),
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
