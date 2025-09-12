#' Internal function for the Mixed-Effects Model Calibration Curve (MIXC)
#'
#' Generates calibration plots for mixed-effects models, assessing the agreement between
#' predicted probabilities and observed binary outcomes while accounting for clustering.
#' Supports random intercept and random slope models with confidence and prediction intervals.
#'
#' @param data Optional data frame containing the columns for `preds`, `y`, and `cluster`.
#'   If provided, `preds`, `y`, and `cluster` should be column names (unquoted).
#'   Default is `NULL`.
#' @param preds A numeric vector of predicted probabilities.
#' @param y A numeric vector of binary outcomes (0 or 1).
#' @param cluster A factor or character vector of cluster identifiers.
#' @param grid_length Integer; number of points for the prediction grid. Default is `100`.
#' @param method Character; type of mixed-effects model: `"intercept"` (random intercept)
#'   or `"slope"` (random slope). Default is `"slope"`.
#' @param plot Logical; whether to generate a calibration plot. Default is `TRUE`.
#' @param cluster_curves Logical; whether to include cluster-specific curves in the plot.
#'   Default is `FALSE`.
#' @param nsims_pi Integer; number of simulations for prediction intervals. Default is `10000`.
#' @param CI Logical; whether to calculate confidence intervals. Default is `TRUE`.
#' @param CI_method Character; method for confidence intervals: `"delta"` or empty string.
#'   Default is `""`.
#'
#' @details
#' This function fits mixed-effects logistic regression models to account for clustering
#' and generates calibration curves with optional confidence and prediction intervals.
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted mixed-effects model object}
#'   \item{cluster_data}{Data frame with calibration data for each cluster}
#'   \item{plot_data}{Data frame with calibration data for the average cluster}
#'   \item{observed_data}{Data frame with calibration data for individual patients}
#'   \item{plot}{A `ggplot2` object if `plot = TRUE`, otherwise `NULL`}
#' }
#'
#'
#' @importFrom dplyr filter group_by summarise
#' @importFrom lme4 glmer
#' @importFrom merTools predictInterval
#' @importFrom ggplot2 ggplot geom_abline geom_ribbon geom_line xlab ylab theme_classic
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous coord_cartesian theme
#' @importFrom ggplot2 scale_fill_manual unit
#' @importFrom Matrix tcrossprod
#' @importFrom rms rcs
#'
#' @export
MIXC <- function(data = NULL,
                 preds,
                 y,
                 cluster,
                 grid_length = 100,
                 method = c("slope", "intercept"),
                 plot = TRUE,
                 cluster_curves = FALSE,
                 nsims_pi = 10000,
                 CI = TRUE,
                 CI_method = "") {
  # --- Extract from data if provided ---
  callFn = match.call()
  method = match.arg(method)
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
    logit_preds = Logit(as.numeric(preds)),
    y           = as.factor(y),
    cluster     = as.factor(cluster)
  )

  # --- Grid setup ---
  grid <- seq(0.01, 0.99, length.out = grid_length)

  # --- Model fitting ---
  if (method == "intercept") {
    fit_model <- glmer(
      y ~ rcs(logit_preds, 3) + (1 | cluster),
      data = df,
      family = "binomial",
      verbose = 0
    )
  } else if (method == "slope") {
    fit_model <- glmer(
      y ~ rcs(logit_preds, 3) + (rcs(logit_preds, 3) | cluster),
      data = df,
      family = "binomial",
      verbose = 0
    )
  }

  # --- Predictions for original data ---
  df$re_preds <- predict(fit_model, df)
  df$obs_prob <- Ilogit(df$re_preds)

  # --- Cluster range calculation ---
  cluster_ranges <- df %>%
    group_by(cluster) %>%
    summarise(
      min_logit = min(logit_preds),
      max_logit = max(logit_preds),
      .groups = "drop"
    )

  # --- Cluster-specific calibration data ---
  cluster_cal_data <- data.frame()

  for (i in seq_len(nrow(cluster_ranges))) {
    current_cluster <- cluster_ranges$cluster[i]
    cluster_logit_preds <- Logit(grid)

    temp_data <- data.frame(
      cluster = current_cluster,
      logit_preds = cluster_logit_preds,
      x = grid
    )

    temp_data$p_obs_logit <- predict(fit_model, newdata = temp_data)
    temp_data$pred_prob <- Ilogit(temp_data$logit_preds)
    temp_data$obs_prob <- Ilogit(temp_data$p_obs_logit)

    cluster_cal_data <- rbind(cluster_cal_data, temp_data)
  }

  # --- Average calibration data ---
  avg_logit_preds <- Logit(grid)
  avg_cal_data <- data.frame(logit_preds = avg_logit_preds)

  avg_cal_data$p_obs_logit <- predict(fit_model, newdata = avg_cal_data, re.form = NA)
  avg_cal_data$pred_prob <- Ilogit(avg_cal_data$logit_preds)
  avg_cal_data$obs_prob <- Ilogit(avg_cal_data$p_obs_logit)
  avg_cal_data$x <- grid

  # --- Confidence intervals ---
  if (CI) {
    deltaCI <- function(eta, se, alpha) {
      mu_hat = Ilogit(eta)
      se_mu  = binomial()$mu.eta(eta) * se
      data.frame(
        y_hat = mu_hat,
        lower_delta = pmax(0, pmin(mu_hat + qnorm(alpha / 2) * se_mu, 1)),
        upper_delta = pmax(0, pmin(mu_hat + qnorm(1 - alpha / 2) * se_mu, 1))
      )
    }

    mm <- model.matrix(~ rcs(logit_preds, 3), avg_cal_data)
    # pvar1 <- diag(mm %*% Matrix::tcrossprod(vcov(fit_model), mm))
    V <- vcov(fit_model) # dense
    pvar1 <- rowSums(as.array((mm %*% V) * mm))

    if (CI_method == "delta") {
      tvar1 <- pvar1 + as.numeric(VarCorr(fit_model)$cluster[1])
      delta_results <- deltaCI(avg_cal_data$p_obs_logit, sqrt(tvar1), 0.05)

      avg_cal_data$p_lower_ci <- delta_results$lower_delta
      avg_cal_data$p_upper_ci <- delta_results$upper_delta
      avg_cal_data$logit_lower_ci <- Logit(avg_cal_data$p_lower_ci)
      avg_cal_data$logit_upper_ci <- Logit(avg_cal_data$p_upper_ci)
    } else {
      avg_cal_data$logit_lower_ci <- avg_cal_data$p_obs_logit - 1.96 * sqrt(pvar1)
      avg_cal_data$logit_upper_ci <- avg_cal_data$p_obs_logit + 1.96 * sqrt(pvar1)
      avg_cal_data$p_lower_ci <- Ilogit(avg_cal_data$logit_lower_ci)
      avg_cal_data$p_upper_ci <- Ilogit(avg_cal_data$logit_upper_ci)
    }

    avg_cal_data$cluster <- "average"

    # --- Prediction intervals ---
    PI <- suppressWarnings(predictInterval(
      merMod = fit_model,
      newdata = avg_cal_data,
      level = 0.95,
      n.sims = nsims_pi,
      stat = "median",
      type = "probability",
      include.resid.var = TRUE
    ))

    avg_cal_data <- cbind(avg_cal_data, PI)
  }

  # --- Plotting ---
  plot_obj <- NULL
  if (plot) {
    plot_obj <- ggplot(avg_cal_data) +
      geom_abline(linetype = "dashed", alpha = 0.1) +
      geom_ribbon(aes(
        x = pred_prob, ymin = p_lower_ci, ymax = p_upper_ci,
        fill = "CI 95%"
      ), alpha = 0.3) +
      geom_ribbon(aes(
        x = pred_prob, ymin = lwr, ymax = upr,
        fill = "PI 95%"
      ), alpha = 0.2) +
      geom_line(aes(x = pred_prob, y = obs_prob),
        linewidth = 1, linetype = "dashed", color = "black"
      ) +
      xlab("Estimated probability") +
      ylab("Observed proportion") +
      theme_classic(base_size = 8, base_family = "serif") +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      scale_fill_manual(name = "Uncertainty", values = c("green4", "green")) +
      theme(
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (cluster_curves) {
      plot_obj <- plot_obj +
        geom_line(
          data = cluster_cal_data,
          aes(x = pred_prob, y = obs_prob, color = cluster),
          linewidth = 0.2, linetype = "dotted", show.legend = FALSE
        )
    }
  }

  # --- Patient-level calibration data ---
  patient_cal_data <- data.frame()

  for (current_cluster in unique(df$cluster)) {
    cluster_subset <- df %>% filter(cluster == current_cluster)

    temp_patient_data <- data.frame(
      cluster = cluster_subset$cluster,
      logit_preds = cluster_subset$logit_preds
    )

    temp_patient_data$p_obs_logit <- predict(fit_model, newdata = temp_patient_data)
    temp_patient_data$pred_prob <- Ilogit(temp_patient_data$logit_preds)
    temp_patient_data$obs_prob <- Ilogit(temp_patient_data$p_obs_logit)

    patient_cal_data <- rbind(patient_cal_data, temp_patient_data)
  }

  # --- Return ---
  list(
    model = fit_model,
    cluster_data = cluster_cal_data,
    plot_data = avg_cal_data,
    observed_data = patient_cal_data,
    plot = plot_obj
  )
}
