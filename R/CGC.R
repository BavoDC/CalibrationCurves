#' Internal function for the Clustered Grouped Calibration Curve (CGC)
#'
#' Generates a flexible calibration curve using predicted values, observed outcomes,
#' and clustering/grouping information. The function supports two grouping methods:
#' equal-sized groups (`"grouped"`) or interval-based groups (`"interval"`).
#' Optionally, a calibration plot can be produced with cluster-specific curves.
#'
#' @param data Optional data frame containing the columns for `preds`, `y`, and `cluster`.
#'   If provided, `preds`, `y`, and `cluster` should be column names (unquoted).
#'   Default is `NULL`.
#' @param preds A numeric vector of predicted probabilities/values, or a column in `data`.
#' @param y A numeric vector of observed outcomes, or a column in `data`.
#' @param cluster A factor/character vector of cluster identifiers, or a column in `data`.
#' @param ntiles Integer; number of groups (tiles) for calibration. Default is `10`.
#' @param cluster_curves Logical; whether to include cluster-specific calibration
#'   curves in the plot. Default is `FALSE`.
#' @param plot Logical; whether to return a calibration plot. Default is `TRUE`.
#' @param size Numeric; point size for plotted curves. Default is `1`.
#' @param linewidth Numeric; line width for plotted curves. Default is `0.4`.
#' @param univariate Logical; whether to use univariate meta-analysis. Default is `FALSE`.
#' @param method Character; grouping method: `"grouped"` (equal-sized groups) or
#'   `"interval"` (interval-based). Default is `"grouped"`.
#' @param cl.level the confidence level for the calculation of the confidence interval. Default is \code{0.95}.
#'
#' @details
#' - `"grouped"`: predictions are divided into equal-sized bins using quantiles.
#' - `"interval"`: predictions are divided into fixed-width bins across [0, 1].
#'
#' The function performs meta-analysis within each group, with optional univariate
#' analysis. Results are aggregated and plotted as calibration curves.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot_data}{Data frame of meta-analysis calibration estimates.}
#'   \item{trad_grouped}{Data frame with traditional grouped calibration results.}
#'   \item{observed_data}{Data frame with per-observation calibration data.}
#'   \item{cluster_data}{Data frame with cluster-specific calibration summaries.}
#'   \item{plot}{A `ggplot2` object if `plot = TRUE`, otherwise `NULL`.}
#' }
#'
#'
#' @importFrom dplyr group_by summarise ungroup mutate ntile
#' @importFrom meta metaprop
#' @importFrom metafor escalc rma.mv
#' @importFrom ggplot2 ggplot geom_ribbon geom_point geom_line xlab ylab theme_classic theme coord_cartesian scale_x_continuous scale_y_continuous scale_fill_manual scale_color_manual geom_abline geom_errorbar
#' @importFrom ggnewscale new_scale_color
#'
#' @export
CGC <- function(data = NULL,
                preds,
                y,
                cluster,
                cl.level = 0.95,
                ntiles = 10,
                cluster_curves = FALSE,
                plot = TRUE,
                size = 1,
                linewidth = 0.4,
                univariate = FALSE,
                method = "grouped") {
  # --- Extract from data if provided ---
  callFn   = match.call()
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
    preds   = as.numeric(preds),
    y       = as.numeric(y),
    cluster = as.factor(cluster)
  )

  # --- Assign decile/interval groups ---
  if (method == "grouped") {
    df <- df %>%
      group_by(cluster) %>%
      mutate(decile_group = ntile(preds, ntiles))
  } else if (method == "interval") {
    df <- df %>%
      group_by(cluster) %>%
      mutate(decile_group = cut(preds, breaks = seq(0, 1, length.out = ntiles + 1)))
  } else {
    stop("Invalid 'method'. Must be 'grouped' or 'interval'.")
  }

  # --- Cluster-level summaries ---
  deciles_cluster <- df %>%
    group_by(cluster, decile_group) %>%
    mutate(n_tile = n()) %>%
    ungroup() %>%
    mutate(N_tile = n())

  decilesbycluster <- deciles_cluster %>%
    group_by(cluster, decile_group) %>%
    summarise(
      n_tile = n(),
      y_mean = mean(y),
      x_mean = mean(preds),
      .groups = "drop"
    )

  # --- Traditional grouped calibration ---
  deciles_all <- df %>%
    ungroup() %>%
    mutate(decile_group = if (method == "grouped") ntile(preds, ntiles) else cut(preds, breaks = seq(0, 1, length.out = ntiles + 1))) %>%
    group_by(decile_group) %>%
    summarise(
      std_error_y = sd(y) / sqrt(length(y)),
      std_error_x = sd(preds) / sqrt(length(preds)),
      var_x = var(preds),
      var_y = var(y),
      preds = mean(preds),
      y = mean(y),
      .groups = "drop"
    )

  # --- Within-cluster variances ---
  var_150_cluster_decile <- deciles_cluster %>%
    group_by(cluster, decile_group) %>%
    summarise(
      var_x_cluster_tile = var(preds),
      var_y_cluster_tile = var(y),
      preds_150 = mean(preds),
      y_150 = mean(y),
      ntile_150 = n(),
      .groups = "drop"
    )

  # --- Meta-analysis across groups ---
  data_all <- data.frame()
  for (i in unique(var_150_cluster_decile$decile_group)) {
    data_meta <- var_150_cluster_decile %>% filter(decile_group == i)

    if (univariate) {
      x_meta <- metaprop(
        data = data_meta, event = preds_150 * ntile_150,
        n = ntile_150, studlab = cluster,
        method = "Inverse", backtransf = TRUE
      )
      y_meta <- metaprop(
        data = data_meta, event = y_150 * ntile_150,
        n = ntile_150, studlab = cluster,
        method = "Inverse", backtransf = TRUE
      )
      data_meta_curve <- data.frame(
        te_x = x_meta$TE.random, ci_up_x = x_meta$upper.random, ci_low_x = x_meta$lower.random,
        pre_up_x = x_meta$upper.predict, pre_low_x = x_meta$lower.predict,
        te_y = y_meta$TE.random, ci_up_y = y_meta$upper.random, ci_low_y = y_meta$lower.random,
        pre_up_y = y_meta$upper.predict, pre_low_y = y_meta$lower.predict,
        decile_group = i, ntile_plot = sum(data_meta$ntile_150)
      )
    } else {
      preds_escalc <- escalc(
        measure = "PLO", xi = preds_150 * ntile_150,
        ni = ntile_150, data = data_meta
      )
      preds_escalc$group <- "preds"
      y_escalc <- escalc(
        measure = "PLO", xi = y_150 * ntile_150,
        ni = ntile_150, data = data_meta
      )
      y_escalc$group <- "y"
      dat <- rbind(preds_escalc, y_escalc)

      res <- try(rma.mv(yi, vi,
        mods = ~ group - 1,
        random = ~ group | cluster,
        struct = "UN", data = dat,
        control = list(iter.max = 10000, rel.tol = 1e-6)
      ))

      if (inherits(res, "try-error") | nrow(dat) == 2) {
        warning("Meta-analysis failed for decile group ", i, ". Returning NA values.")
        next
      }

      pred_up <- res$b + qt(0.975, df = nrow(y_escalc) - 2) * sqrt(res$tau2 + res$se^2)
      pred_low <- res$b - qt(0.975, df = nrow(y_escalc) - 2) * sqrt(res$tau2 + res$se^2)

      data_meta_curve <- data.frame(
        te_x = res$b[1], ci_up_x = res$ci.ub[1], ci_low_x = res$ci.lb[1],
        pre_up_x = pred_up[1], pre_low_x = pred_low[1],
        te_y = res$b[2], ci_up_y = res$ci.ub[2], ci_low_y = res$ci.lb[2],
        pre_up_y = pred_up[2], pre_low_y = pred_low[2],
        decile_group = i, ntile_plot = sum(data_meta$ntile_150)
      )
    }

    data_all <- rbind(data_all, data_meta_curve)
  }

  data_all <- data_all %>% mutate(across(-c(decile_group, ntile_plot), Ilogit))

  # --- Plotting ---
  curve <- NULL
  if (plot) {
    curve <- ggplot(data_all, aes(x = te_x, y = te_y)) +
      geom_abline(linetype = "dashed", alpha = 0.1) +
      geom_ribbon(aes(ymax = pre_up_y, ymin = pre_low_y, fill = "PI 95%"), alpha = 1) +
      geom_ribbon(aes(ymax = ci_up_y, ymin = ci_low_y, fill = "CI 95%"), alpha = 1) +
      geom_errorbar(aes(x = te_x, y = te_y, xmin = pre_low_x, xmax = pre_up_x),
        alpha = 0.2, width = 0, lwd = linewidth, lty = "dashed"
      ) +
      scale_fill_manual(name = "Heterogeneity", values = c("cornflowerblue", "lightblue")) +
      geom_point(data = deciles_all, aes(x = preds, y = y, color = "Traditional grouped"), size = size) +
      geom_line(data = deciles_all, aes(x = preds, y = y, color = "Traditional grouped"), linewidth = linewidth) +
      geom_point(data = data_all, aes(color = paste0("CGC-C(", method, ")")), size = size * 1.3) +
      geom_line(data = data_all, aes(color = paste0("CGC-C(", method, ")")), linewidth = linewidth) +
      scale_color_manual(name = "Methodology", values = c("black", "gold")) +
      xlab("Estimated probability") +
      ylab("Observed proportion") +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_classic(base_size = 8, base_family = "serif") +
      theme(
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (cluster_curves) {
      curve <- curve +
        geom_line(
          data = decilesbycluster, aes(x_mean, y_mean, group = cluster),
          lwd = linewidth, alpha = 0.5, show.legend = FALSE
        ) +
        geom_point(
          data = decilesbycluster, aes(x_mean, y_mean, group = cluster),
          alpha = 0.8, size = size / 3, show.legend = FALSE
        )
    }
  }

  # --- Return ---
  list(
    plot_data = data_all,
    trad_grouped = deciles_all,
    observed_data = deciles_cluster,
    cluster_data = decilesbycluster,
    plot = curve
  )
}
