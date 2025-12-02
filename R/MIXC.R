#' Internal function for the Mixed-Effects Model Calibration Curve (MIXC)
#'
#' Estimates the calibration curve using a logistic generalized linear mixed model.
#'
#' @param data optional data frame containing the variables \code{p}, \code{y},
#'   and \code{cluster}. If supplied, variable names should be given without
#'   quotation marks.
#' @param p predicted probabilities (numeric vector) or name of the column in
#'   \code{data}.
#' @param y binary outcome variable or the name of the column in \code{data}.
#' @param cluster Cluster identifier (factor, character, or integer) or name of
#'   the column in \code{data}.
#' @param grid the grid for the calibration curve evaluation
#' @param method character, type of mixed-effects model: \code{"intercept"} (random intercept)
#'   or \code{"slope"} (random slope). Default is \code{"slope"}.
#' @param plot logical, indicating whether to generate a calibration plot. Default is \code{TRUE}.
#' @param cluster_curves logical, whether to include cluster-specific curves in the plot.
#'   Default is \code{FALSE}.
#' @param nsims_pi integer, number of simulations for prediction intervals. Default is \code{10000}.
#' @param CI logical, whether to calculate confidence intervals. Default is \code{TRUE}.
#' @param CI_method character, method for computing the confidence intervals of the observed proportions.
#'  If \code{"delta"}, the delta method is applied. Conversely, when \code{CI_method == "naive"}, no correction is applied.
#'  Default is \code{"naive"}.
#' @param cl.level the confidence level for the calculation of the confidence interval. Default is \code{0.95}.
#'
#' @details
#' This function estimates the calibration curves using a logistic generalized linear mixed model.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{model}}{The fitted mixed-effects model object}
#'   \item{\code{cluster_data}}{Data frame with calibration data for each cluster}
#'   \item{\code{plot_data}}{Data frame with calibration data for the average cluster}
#'   \item{\code{observed_data}}{Data frame with calibration data for individual observations}
#'   \item{\code{plot}}{A \code{ggplot2} object if \code{plot = TRUE}, otherwise \code{NULL}}
#' }
MIXC <- function(data = NULL,
                 p,
                 y,
                 cluster,
                 grid,
                 method = c("slope", "intercept"),
                 plot = TRUE,
                 cluster_curves = FALSE,
                 nsims_pi = 10000,
                 CI = TRUE,
                 CI_method = c("naive", "delta"),
                 cl.level = 0.95) {
  # --- Extract from data if provided ---
  callFn <- match.call()
  method <- match.arg(method)
  CI_method <- match.arg(CI_method)

  if (!is.null(data)) {
    if (!all(sapply(c("preds", "y", "cluster"), function(a) as.character(callFn[a])) %in% colnames(data))) {
      stop(paste("Variables", paste0(
        callFn[c("preds", "y", "cluster")],
        collapse = ", "
      ), "were not found in the data.frame."))
    }
    preds <- eval(callFn$preds, data)
    logit <- Logit(preds)
    y <- eval(callFn$y, data)
    cluster <- eval(callFn$cluster, data)
  }
  alpha <- 1 - cl.level

  # --- Base dataframe ---
  df <- data.frame(
    logit_preds = Logit(as.numeric(preds)),
    y           = as.factor(y),
    cluster     = as.factor(cluster)
  )

  # --- Model fitting ---
  args <- list(
    data = df,
    family = "binomial",
    verbose = 0
  )
  args$formula <-
    if (method == "intercept") {
      y ~ rcs(logit_preds, 3) + (1 | cluster)
    } else {
      y ~ rcs(logit_preds, 3) + (rcs(logit_preds, 3) | cluster)
    }
  # Suppress "boundary (singular)" warning: arises from restricted cubic spline basis multicollinearity in GLMM,
  # not from model misspecification, expected in this context.
  fit_model <- suppressMessages(do.call("glmer", args))


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

  avg_cal_data$p_obs_logit <- eta <- predict(fit_model, newdata = avg_cal_data, re.form = NA)
  avg_cal_data$pred_prob <- Ilogit(avg_cal_data$logit_preds)
  avg_cal_data$obs_prob <- Ilogit(avg_cal_data$p_obs_logit)
  avg_cal_data$x <- grid

  # --- Confidence intervals ---
  if (CI) {
    deltaCI <- function(eta, se, alpha) {
      muHat <- Ilogit(eta)
      seMu <- binomial()$mu.eta(eta) * se
      data.frame(
        eta     = eta,
        etaLCL  = eta + qnorm(alpha / 2) * se,
        etaUCL  = eta + qnorm(1 - alpha / 2) * se,
        yHat    = muHat,
        yHatLCL = pmax(0, pmin(muHat + qnorm(alpha / 2) * seMu, 1)),
        yHatUCL = pmax(0, pmin(muHat + qnorm(1 - alpha / 2) * seMu, 1))
      )
    }

    Z <- mm <- model.matrix(~ rcs(logit_preds, 3), avg_cal_data)
    # pvar <- diag(mm %*% Matrix::tcrossprod(vcov(fit_model), mm))
    V <- vcov(fit_model) # dense
    pvar <- rowSums(as.array((mm %*% V) * mm))
    D <- VarCorr(fit_model)$cluster
    rvar <- rowSums((Z %*% D) * Z)
    tvar <- pvar + rvar
    if (CI_method == "delta") {
      deltaResults <- deltaCI(eta, sqrt(tvar), alpha)
      avg_cal_data <- within(
        avg_cal_data,
        {
          p_lower_ci <- deltaResults$yHatLCL
          p_upper_ci <- deltaResults$yHatUCL
          logit_lower_ci <- deltaResults$etaLCL
          logit_upper_ci <- deltaResults$etaUCL
        }
      )
    } else {
      avg_cal_data <- within(
        avg_cal_data,
        {
          logit_lower_ci <- p_obs_logit + qnorm(alpha / 2) * sqrt(tvar)
          logit_upper_ci <- p_obs_logit + qnorm(1 - alpha / 2) * sqrt(tvar)
          p_lower_ci <- Ilogit(logit_lower_ci)
          p_upper_ci <- Ilogit(logit_upper_ci)
        }
      )
    }

    avg_cal_data$cluster <- "average"

    # --- Prediction intervals ---
    # Here, we suppress the following warnings:
    # Warning messages:
    #   1:      The following levels of cluster from newdata
    # -- average -- are not in the model data.
    # Currently, predictions for these values are based only on the
    # fixed coefficients and the observation-level error.
    # 2: executing %dopar% sequentially: no parallel backend registered
    PI <- suppressWarnings(predictInterval(
      merMod = fit_model,
      newdata = avg_cal_data,
      level = cl.level,
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
          aes(x = pred_prob, y = obs_prob, group = cluster),
          linewidth = 0.2, linetype = "solid", show.legend = FALSE
        )
    }
  }

  # --- Patient-level calibration data ---
  patient_cal_data <- data.frame()

  for (current_cluster in unique(df$cluster)) {
    cluster_subset <- df %>% filter(cluster == current_cluster)

    temp_patient_data <- data.frame(
      cluster     = cluster_subset$cluster,
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
