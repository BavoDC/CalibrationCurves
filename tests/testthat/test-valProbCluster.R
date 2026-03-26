# These tests use the bundled clustered data and a fitted mixed model.
# Most are marked skip_on_cran() due to runtime.

local_cluster_preds <- function(env = parent.frame()) {
  data("clustertraindata", package = "CalibrationCurves", envir = env)
  data("clustertestdata", package = "CalibrationCurves", envir = env)
  mFit <- lme4::glmer(y ~ x1 + x2 + x3 + x5 + (1 | cluster),
                       data = env$clustertraindata, family = "binomial")
  preds   <- predict(mFit, env$clustertestdata, type = "response", re.form = NA)
  y       <- env$clustertestdata$y
  cluster <- env$clustertestdata$cluster
  list(p = unname(preds), y = y, cluster = cluster)
}

# --------------------------------------------------------------------------
# valProbCluster wrapper
# --------------------------------------------------------------------------

test_that("valProbCluster returns correct structure with MIXC", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 50)
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_named(res, c("call", "approach", "cl.level", "grid", "ggPlot", "results"),
               ignore.order = TRUE)
  expect_equal(res$approach, "MIXC")
})

test_that("valProbCluster returns correct structure with CGC", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_equal(res$approach, "CGC")
})

test_that("valProbCluster returns correct structure with MAC2", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 50)
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_equal(res$approach, "MAC2")
})

test_that("valProbCluster errors on non-binary y", {
  d <- local_cluster_preds()
  expect_error(
    valProbCluster(p = d$p, y = d$y + 0.5, cluster = d$cluster, plot = FALSE),
    "binary outcome"
  )
})

test_that("valProbCluster errors on single cluster", {
  d <- local_cluster_preds()
  expect_error(
    valProbCluster(p = d$p, y = d$y, cluster = rep(1, length(d$y)), plot = FALSE),
    "at least two"
  )
})

test_that("valProbCluster produces a ggplot when plot = TRUE", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "MIXC", grid_l = 50)
  )
  expect_s3_class(res$ggPlot, "ggplot")
})

test_that("valProbCluster grid has correct length", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 75)
  )
  expect_length(res$grid, 75)
})

# --------------------------------------------------------------------------
# MIXC-specific arguments
# --------------------------------------------------------------------------

test_that("MIXC with method = 'intercept' errors due to dimension mismatch (known issue)", {
  skip_on_cran()
  d <- local_cluster_preds()
  # Known bug: Z %*% D non-conformable for intercept-only model in MIXC
  expect_error(
    suppressWarnings(
      valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                     plot = FALSE, approach = "MIXC", grid_l = 50,
                     method = "intercept")
    ),
    "non-conformable"
  )
})

test_that("MIXC with method = 'slope'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 50,
                   method = "slope")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true(!is.null(res$results$model))
})

test_that("MIXC with method = 'slope' returns model with random slopes", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 30,
                   method = "slope")
  )
  # Random effects should have > 1 term (intercept + slope)
  vc <- lme4::VarCorr(res$results$model)$cluster
  expect_true(nrow(vc) > 1)
})

test_that("MIXC with CI_method = 'delta'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 50,
                   CI_method = "delta")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true("p_lower_ci" %in% names(res$results$plot_data))
  expect_true("p_upper_ci" %in% names(res$results$plot_data))
})

test_that("MIXC with CI_method = 'naive'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 50,
                   CI_method = "naive")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true("p_lower_ci" %in% names(res$results$plot_data))
})

test_that("MIXC with cluster_curves = TRUE produces plot with cluster lines", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "MIXC", grid_l = 50,
                   cluster_curves = TRUE)
  )
  expect_s3_class(res$ggPlot, "ggplot")
  # Should have more layers than without cluster curves
  res2 <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "MIXC", grid_l = 50,
                   cluster_curves = FALSE)
  )
  expect_true(length(res$ggPlot$layers) > length(res2$ggPlot$layers))
})

test_that("MIXC results contain expected components", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 50)
  )
  expect_true("model" %in% names(res$results))
  expect_true("cluster_data" %in% names(res$results))
  expect_true("plot_data" %in% names(res$results))
  expect_true("observed_data" %in% names(res$results))
  expect_s4_class(res$results$model, "glmerMod")
})

# --------------------------------------------------------------------------
# CGC-specific arguments
# --------------------------------------------------------------------------

test_that("CGC with method = 'grouped'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC",
                   method = "grouped")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true("plot_data" %in% names(res$results))
})

test_that("CGC with method = 'interval'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC",
                   method = "interval")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
})

test_that("CGC with different ntiles", {
  skip_on_cran()
  d <- local_cluster_preds()
  res5 <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC", ntiles = 5)
  )
  res10 <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC", ntiles = 10)
  )
  # Different ntiles → different number of groups in result
  expect_true(nrow(res5$results$plot_data) <= nrow(res10$results$plot_data))
})

test_that("CGC warns when ntiles < 5", {
  skip_on_cran()
  d <- local_cluster_preds()
  expect_warning(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC", ntiles = 3),
    "too low"
  )
})

test_that("CGC with univariate = TRUE", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC",
                   univariate = TRUE, ntiles = 5)
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true(nrow(res$results$plot_data) > 0)
})

test_that("CGC with cluster_curves = TRUE", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "CGC",
                   cluster_curves = TRUE)
  )
  expect_s3_class(res$ggPlot, "ggplot")
})

test_that("CGC results contain expected components", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "CGC")
  )
  expect_true("plot_data" %in% names(res$results))
  expect_true("trad_grouped" %in% names(res$results))
  expect_true("cluster_data" %in% names(res$results))
  expect_true("observed_data" %in% names(res$results))
})

# --------------------------------------------------------------------------
# MAC2-specific arguments
# --------------------------------------------------------------------------

test_that("MAC2 with method_choice = 'log'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 30,
                   method_choice = "log")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_true(nrow(res$results$plot_data) > 0)
})

test_that("MAC2 with method_choice = 'loess'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 30,
                   method_choice = "loess")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
})

test_that("MAC2 with method_choice = 'splines'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 30,
                   method_choice = "splines")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
})

test_that("MAC2 different method_choice values give different results", {
  skip_on_cran()
  d <- local_cluster_preds()
  res_log <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 20,
                   method_choice = "log")
  )
  res_spl <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 20,
                   method_choice = "splines")
  )
  expect_false(identical(res_log$results$plot_data$y,
                         res_spl$results$plot_data$y))
})

test_that("MAC2 with cluster_curves = TRUE", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "MAC2", grid_l = 30,
                   cluster_curves = TRUE)
  )
  expect_s3_class(res$ggPlot, "ggplot")
})

test_that("MAC2 with transf = 'identity'", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 30,
                   transf = "identity")
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
})

test_that("MAC2 results contain expected components", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MAC2", grid_l = 30)
  )
  expect_true("cluster_data" %in% names(res$results))
  expect_true("plot_data" %in% names(res$results))
  # plot_data should have CI/PI columns
  pd <- res$results$plot_data
  expect_true(all(c("y", "upper", "lower", "up_pre", "low_pre") %in% names(pd)))
})

# --------------------------------------------------------------------------
# Shared behavior
# --------------------------------------------------------------------------

test_that("valProbCluster respects cl.level argument", {
  skip_on_cran()
  d <- local_cluster_preds()
  res90 <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 30,
                   cl.level = 0.90)
  )
  res99 <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 30,
                   cl.level = 0.99)
  )
  expect_equal(res90$cl.level, 0.90)
  expect_equal(res99$cl.level, 0.99)
  # Wider CI at 99%
  ci_width_90 <- with(res90$results$plot_data,
                      mean(p_upper_ci - p_lower_ci, na.rm = TRUE))
  ci_width_99 <- with(res99$results$plot_data,
                      mean(p_upper_ci - p_lower_ci, na.rm = TRUE))
  expect_true(ci_width_99 > ci_width_90)
})

test_that("valProbCluster warns when clusters with single outcome are removed", {
  skip_on_cran()
  d <- local_cluster_preds()
  # Create a cluster that only has y = 0
  d$y[d$cluster == d$cluster[1]] <- 0L
  expect_warning(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "MIXC", grid_l = 30),
    "removed"
  )
})

# --------------------------------------------------------------------------
# Default (combined MAC2 + MIXC) approach
# --------------------------------------------------------------------------

test_that("valProbCluster default approach returns correct structure", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = TRUE, approach = "default", grid_l = 50)
  )
  expect_s3_class(res, "ClusteredCalibrationCurve")
  expect_equal(res$approach, "default")

 # results should contain both overall (MAC2) and clusters (MIXC)
  expect_true(all(c("overall", "clusters") %in% names(res$results)))

  # MAC2 overall curve data
  expect_true("plot_data" %in% names(res$results$overall))
  mac2_cols <- c("x", "y", "upper", "lower", "up_pre", "low_pre")
  expect_true(all(mac2_cols %in% names(res$results$overall$plot_data)))

  # MIXC cluster-specific data
  expect_true("cluster_data" %in% names(res$results$clusters))
  mixc_cols <- c("cluster", "pred_prob", "obs_prob")
  expect_true(all(mixc_cols %in% names(res$results$clusters$cluster_data)))

  # ggPlot should be present
  expect_s3_class(res$ggPlot, "ggplot")
})

test_that("valProbCluster default approach works with plot = FALSE", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, approach = "default", grid_l = 50)
  )
  expect_null(res$ggPlot)
  expect_true(all(c("overall", "clusters") %in% names(res$results)))
})

test_that("valProbCluster uses default approach when approach is not specified", {
  skip_on_cran()
  d <- local_cluster_preds()
  res <- suppressWarnings(
    valProbCluster(p = d$p, y = d$y, cluster = d$cluster,
                   plot = FALSE, grid_l = 50)
  )
  expect_equal(res$approach, "default")
})
