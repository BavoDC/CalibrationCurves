<<<<<<< Updated upstream
# New name
CGC_old <- function(data = NULL, preds, y, cluster, ntiles = 10, cluster_curves = F, plot = TRUE, size = 1, linewidth = 0.4, univariate = FALSE, method = "grouped") {
  #' Clustered grouped Calibration Curve (CGC)
  #'
  #' Generates a flexible calibration curve using a vector of predictions, outcomes, and clusters/centers.
  #' The function supports two grouping methods: equal-sized groups ("grouped") and interval-based groups ("interval").
  #' The curve can optionally be plotted with various customizable options.
  #'
  #' @param data Optional data frame containing the columns for `preds`, `y`, and `cluster`. Default is `NULL`.
  #' @param preds A numeric vector of predicted probabilities or values. Ignored if `data` is provided.
  #' @param y A numeric vector of actual outcomes. Ignored if `data` is provided.
  #' @param cluster A factor or character vector indicating the cluster/center for each observation. Ignored if `data` is provided.
  #' @param ntiles An integer specifying the number of groups (tiles) to divide the data into. Default is 10.
  #' @param plot Logical; whether to generate and display a calibration plot. Default is `TRUE`.
  #' @param size Numeric; size of points on the plot. Default is 1.
  #' @param linewidth Numeric; width of the lines in the plot. Default is 0.4.
  #' @param univariate Logical; whether to perform univariate meta-analysis. Default is `FALSE`.
  #' @param method Character; grouping method for calibration. Options are `"grouped"` (equal-sized groups) or `"interval"` (interval-based groups). Default is `"grouped"`.
  #' @param cluster_curves Logical; whether to include cluster-specific calibration curves on the plot. Default is `FALSE`.
  #'
  #' @return A list with the following elements:
  #'   \describe{
  #'     \item{cluster_grouped}{Data frame with cluster-grouped calibration results.}
  #'     \item{trad_grouped}{Data frame with traditional grouped calibration results.}
  #'     \item{bypatient}{Data frame with calibration data for each patient.}
  #'     \item{bycluster}{Data frame with calibration data for each cluster.}
  #'     \item{plot}{The ggplot object for the calibration curve, if `plot = TRUE`. Otherwise, `NULL`.}
  #'   }
  #'
  #' @examples
  #' preds <- runif(1000)
  #' y <- rbinom(1000, 1, preds)
  #' cluster <- rep(1:10, each = 100)
  #' curve_data <- CGC(preds = preds, y = y, cluster = cluster, ntiles = 10, plot = TRUE)
  #'
  #' @importFrom dplyr group_by summarise ungroup mutate ntile
  #' @importFrom meta metaprop
  #' @importFrom metafor escalc rma.mv
  #' @importFrom ggplot2 ggplot geom_ribbon geom_point geom_line xlab ylab theme_classic theme coord_cartesian scale_x_continuous scale_y_continuous scale_fill_manual scale_color_manual
  #' @importFrom plotrix std.error
  #' @importFrom ggnewscale new_scale_color
  #'
  #' @export
  if (!is.null(data)) {
    preds <- data[[deparse(substitute(preds))]]
    y <- data[[deparse(substitute(y))]]
    cluster <- data[[deparse(substitute(cluster))]]
  }
  deciles <- data.frame(preds = as.numeric(preds), y = as.numeric(y), cluster = as.factor(cluster))

  if (method == "grouped") {
    deciles_cluster <- deciles %>%
      group_by(cluster) %>%
      mutate(decile_group = ntile(preds, ntiles))

    decilesbycluster <- deciles_cluster %>%
      group_by(cluster, decile_group) %>%
      summarise(
        n_tile = sum(n()),
        y_mean = mean(y),
        x_mean = mean(preds)
      )


    deciles_cluster <- deciles_cluster %>%
      group_by(cluster, decile_group) %>%
      mutate(
        n_tile = sum(n()),
        y_mean = mean(y)
      ) %>%
      ungroup(cluster) %>%
      mutate(N_tile = sum(n()))

    deciles_all <- deciles %>%
      mutate(decile_group = ntile(preds, ntiles)) %>%
      group_by(decile_group) %>%
      summarise(
        std_error_y = plotrix::std.error(y),
        std_error_x = plotrix::std.error(preds),
        var_x = var(preds),
        var_y = var(y),
        preds = mean(preds),
        y = mean(y),
        .groups = "drop"
      )

    var_15_cluster <- deciles_cluster %>%
      group_by(cluster) %>%
      dplyr::summarise(
        var_x_cluster = var(preds),
        var_y_cluster = var(y),
        .groups = "drop"
      )

    var_150_cluster_decile <- deciles_cluster %>%
      group_by(cluster, decile_group) %>%
      dplyr::summarise(
        var_x_cluster_tile = var(preds),
        var_y_cluster_tile = var(y),
        preds_150 = mean(preds),
        y_150 = mean(y),
        ntile_150 = n(),
        .groups = "drop"
      )
  }

  if (method == "interval") {
    deciles <- deciles %>%
      mutate(decile_group = cut(preds, breaks = seq(0, 1, length.out = ntiles + 1)))

    deciles_cluster <- deciles %>%
      group_by(cluster, decile_group) %>%
      mutate(n_tile = sum(n())) %>%
      ungroup(cluster) %>%
      mutate(N_tile = sum(n()))

    decilesbycluster <- deciles_cluster %>%
      group_by(cluster, decile_group) %>%
      summarise(
        n_tile = sum(n()),
        y_mean = mean(y),
        x_mean = mean(preds)
      )

    deciles_all <- deciles %>%
      group_by(decile_group) %>%
      summarise(
        std_error_y = plotrix::std.error(y),
        std_error_x = plotrix::std.error(preds),
        var_x = var(preds),
        var_y = var(y),
        preds = mean(preds),
        y = mean(y),
        .groups = "drop"
      )

    var_15_cluster <- deciles_cluster %>%
      group_by(cluster) %>%
      dplyr::summarise(
        var_x_cluster = var(preds),
        var_y_cluster = var(y),
        .groups = "drop"
      )

    var_150_cluster_decile <- deciles_cluster %>%
      group_by(cluster, decile_group) %>%
      dplyr::summarise(
        var_x_cluster_tile = var(preds),
        var_y_cluster_tile = var(y),
        preds_150 = mean(preds),
        y_150 = mean(y),
        ntile_150 = n(),
        .groups = "drop"
      )
  }

  data_all <- data.frame()

  for (i in unique(var_150_cluster_decile$decile_group)) {
    data_meta <- var_150_cluster_decile %>%
      filter(decile_group == i)

    if (univariate) {
      x_meta <- meta::metaprop(data = data_meta, event = preds_150 * ntile_150, n = ntile_150, studlab = cluster, method = "Inverse", backtransf = TRUE)
      y_meta <- meta::metaprop(data = data_meta, event = y_150 * ntile_150, n = ntile_150, studlab = cluster, method = "Inverse", backtransf = TRUE)
      data_meta_curve <- data.frame(
        te_x = x_meta$TE.random, ci_up_x = x_meta$upper.random, ci_low_x = x_meta$lower.random, pre_up_x = x_meta$upper.predict,
        pre_low_x = x_meta$lower.predict, te_y = y_meta$TE.random, ci_up_y = y_meta$upper.random, ci_low_y = y_meta$lower.random,
        pre_up_y = y_meta$upper.predict, pre_low_y = y_meta$lower.predict, decile_group = i, ntile_plot = sum(data_meta$ntile_150)
      )
    } else {
      preds_escalc <- metafor::escalc(measure = "PLO", xi = preds_150 * ntile_150, ni = ntile_150, data = data_meta)
      preds_escalc$group <- "preds"
      y_escalc <- metafor::escalc(measure = "PLO", xi = y_150 * ntile_150, ni = ntile_150, data = data_meta)
      y_escalc$group <- "y"
      dat <- rbind(preds_escalc, y_escalc)
      res <- try(metafor::rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | cluster, struct = "UN", data = dat, control = list(iter.max = 10000, rel.tol = 1e-6)))
      if (inherits(res, "try-error") | nrow(dat) == 2) {
        warning("Meta-analysis failed for decile group ", i, ". Returning NA values. Set to 0")
        data_meta_curve <- data.frame(
          te_x = dat$preds_150[1], ci_up_x = 0, ci_low_x = 0, pre_up_x = 0,
          pre_low_x = 0, te_y = dat$y_150[1], ci_up_y = 0, ci_low_y = 0,
          pre_up_y = 0, pre_low_y = 0, decile_group = i, ntile_plot = sum(data_meta$ntile_150)
        )
        next
      }
      pred_up <- res$b + qt((1.95) / 2, df = nrow(y_escalc) - 2) * sqrt(res$tau2 + res$se^2)
      pred_low <- res$b - qt((1.95) / 2, df = nrow(y_escalc) - 2) * sqrt(res$tau2 + res$se^2)
      data_meta_curve <- data.frame(
        te_x = res$b[1], ci_up_x = res$ci.ub[1], ci_low_x = res$ci.lb[1], pre_up_x = pred_up[1],
        pre_low_x = pred_low[1], te_y = res$b[2], ci_up_y = res$ci.ub[2], ci_low_y = res$ci.lb[2],
        pre_up_y = pred_up[2], pre_low_y = pred_low[2], decile_group = i, ntile_plot = sum(data_meta$ntile_150)
      )
    }

    data_all <- rbind(data_all, data_meta_curve)
  }

  data_all <- data_all %>%
    mutate(across(-c(decile_group, ntile_plot), ~ plogis(.)))

  if (plot) {
    curve <- ggplot(data_all, aes(x = te_x, y = te_y)) +
      geom_abline(linetype = "dashed", alpha = 0.1) +
      xlab("Estimated probability") +
      ylab("Observed proportion") +
      theme_classic(base_size = 8, base_family = "serif") +
      theme(legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      geom_ribbon(aes(ymax = pre_up_y, ymin = pre_low_y, fill = "PI 95%"), alpha = 1) +
      # geom_ribbon(aes(xmax = pre_up_x, xmin = pre_low_x, fill = "PI 95%"), alpha = 1) +
      geom_ribbon(aes(ymax = ci_up_y, ymin = ci_low_y, fill = "CI 95%"), alpha = 1) +
      # geom_ribbon(aes(xmax = ci_up_x, xmin = ci_low_x, fill = "CI 95%"), alpha = 1) +
      scale_fill_manual(name = "Heterogeneity", values = c("cornflowerblue", "lightblue")) +
      # geom_rect(data = data_all[ntiles, ], aes(xmin = te_x, ymin = te_y, xmax = pre_up_x, ymax = pre_up_y), fill = "lightblue", alpha = 1) +
      # geom_rect(data = data_all[ntiles, ], aes(xmin = te_x, ymin = te_y, xmax = ci_up_x, ymax = ci_up_y), fill = "cornflowerblue", alpha = 1) +
      # geom_rect(data = data_all[1, ], aes(xmin = te_x, ymin = te_y, xmax = pre_low_x, ymax = pre_low_y), fill = "lightblue", alpha = 1) +
      # geom_rect(data = data_all[1, ], aes(xmin = te_x, ymin = te_y, xmax = ci_low_x, ymax = ci_low_y), fill = "cornflowerblue", alpha = 1)+
      geom_errorbar(data = data_all, aes(
        x = te_x, y = te_y,
        xmin = pre_low_x, xmax = pre_up_x
      ), alpha = 0.2, width = 0, lwd = linewidth, lty = "dashed")
    if (cluster_curves) {
      curve <- curve + geom_line(data = decilesbycluster, aes(x_mean, y_mean, color = cluster), lwd = linewidth, alpha = 0.5, show.legend = F) +
        geom_point(data = decilesbycluster, aes(x_mean, y_mean, color = cluster), alpha = 0.5, size = size / 3, alpha = 0.8, show.legend = F)
    }
    name_leg <- ifelse(method == "interval", "interval", "grouped")

    curve <- curve +
      ggnewscale::new_scale_color() +
      geom_point(data = deciles_all, aes(x = preds, y = y, color = "Traditional grouped"), size = size) +
      geom_line(data = deciles_all, aes(x = preds, y = y, color = "Traditional grouped"), linewidth = linewidth) +
      geom_point(data = data_all, aes(color = paste0("CGC-C(", name_leg, ")")), size = size * 1.3) +
      geom_line(data = data_all, aes(color = paste0("CGC-C(", name_leg, ")")), linewidth = linewidth) +
      scale_color_manual(name = "Methodology", values = c("black", "gold")) +
      xlab("Estimated probability") +
      ylab("Observed proportion")
  } else {
    curve <- NULL
  }

  data <- list(plot_data = data_all, trad_grouped = deciles_all, observed_data = deciles_cluster, cluster_data = decilesbycluster, plot = curve)
  return(data)
}

=======
>>>>>>> Stashed changes
#' Clustered Grouped Calibration Curve (CGC)
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
#' @examples
#' set.seed(123)
#' preds <- runif(1000)
#' y <- rbinom(1000, 1, preds)
#' cluster <- rep(1:10, each = 100)
#' res <- CGC(preds = preds, y = y, cluster = cluster, ntiles = 10, plot = TRUE)
#' res$plot
#'
#' @importFrom dplyr group_by summarise ungroup mutate ntile
#' @importFrom meta metaprop
#' @importFrom metafor escalc rma.mv
#' @importFrom ggplot2 ggplot geom_ribbon geom_point geom_line xlab ylab theme_classic theme coord_cartesian scale_x_continuous scale_y_continuous scale_fill_manual scale_color_manual geom_abline geom_errorbar
#' @importFrom plotrix std.error
#' @importFrom ggnewscale new_scale_color
#'
#' @export
CGC <- function(data = NULL,
                preds,
                y,
                cluster,
                ntiles = 10,
                cluster_curves = FALSE,
                plot = TRUE,
                size = 1,
                linewidth = 0.4,
                univariate = FALSE,
                method = "grouped") {
  # --- Extract from data if provided ---
  if (!is.null(data)) {
    preds <- data[[deparse(substitute(preds))]]
    y <- data[[deparse(substitute(y))]]
    cluster <- data[[deparse(substitute(cluster))]]
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
      std_error_y = plotrix::std.error(y),
      std_error_x = plotrix::std.error(preds),
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
      x_meta <- meta::metaprop(
        data = data_meta, event = preds_150 * ntile_150,
        n = ntile_150, studlab = cluster,
        method = "Inverse", backtransf = TRUE
      )
      y_meta <- meta::metaprop(
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
      preds_escalc <- metafor::escalc(
        measure = "PLO", xi = preds_150 * ntile_150,
        ni = ntile_150, data = data_meta
      )
      preds_escalc$group <- "preds"
      y_escalc <- metafor::escalc(
        measure = "PLO", xi = y_150 * ntile_150,
        ni = ntile_150, data = data_meta
      )
      y_escalc$group <- "y"
      dat <- rbind(preds_escalc, y_escalc)

      res <- try(metafor::rma.mv(yi, vi,
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

  data_all <- data_all %>% mutate(across(-c(decile_group, ntile_plot), plogis))

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
