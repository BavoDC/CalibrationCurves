local_binary_preds <- function(env = parent.frame()) {
  data("traindata", package = "CalibrationCurves", envir = env)
  data("testdata", package = "CalibrationCurves", envir = env)
  fit <- glm(y ~ ., data = env$traindata, family = binomial)
  p <- predict(fit, newdata = env$testdata, type = "response")
  y <- env$testdata$y
  list(p = unname(p), y = y)
}

test_that("val.prob.ci.2 returns correct structure", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  expect_s3_class(res, "CalibrationCurve")
  expect_named(res, c("call", "stats", "cl.level", "Calibration",
                       "Cindex", "warningMessages", "CalibrationCurves"),
               ignore.order = TRUE)
})

test_that("val.prob.ci.2 stats vector has correct names and types", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  expected_names <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq", "D:p",
                      "U", "U:Chi-sq", "U:p", "Q", "Brier",
                      "Intercept", "Slope", "Emax", "Brier scaled",
                      "Log-loss", "Eavg", "ECI")
  expect_named(res$stats, expected_names)
  expect_type(res$stats, "double")
})

test_that("val.prob.ci.2 calibration and Cindex CIs have 3 elements", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  expect_length(res$Calibration$Intercept, 3)
  expect_length(res$Calibration$Slope, 3)
  expect_length(res$Cindex, 3)
})

test_that("val.prob.ci.2 C-statistic and Brier are in valid ranges", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  expect_true(res$stats["C (ROC)"] >= 0 && res$stats["C (ROC)"] <= 1)
  expect_true(res$stats["Brier"] >= 0 && res$stats["Brier"] <= 1)
})

test_that("val.prob.ci.2 works with smooth = 'rcs'", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE, smooth = "rcs")
  expect_s3_class(res, "CalibrationCurve")
})

test_that("val.prob.ci.2 works with smooth = 'none'", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE, smooth = "none")
  expect_s3_class(res, "CalibrationCurve")
})

test_that("val.prob.ci.2 errors on non-binary y", {
  d <- local_binary_preds()
  expect_error(val.prob.ci.2(d$p, d$y + 0.5, pl = FALSE),
               "binary outcome")
})

test_that("val.prob.ci.2 errors on mismatched lengths", {
  d <- local_binary_preds()
  expect_error(val.prob.ci.2(d$p[1:10], d$y, pl = FALSE),
               "lengths")
})

test_that("val.prob.ci.2 errors on probabilities outside [0, 1]", {
  d <- local_binary_preds()
  bad_p <- d$p
  bad_p[1] <- 1.5
  expect_error(val.prob.ci.2(bad_p, d$y, pl = FALSE))
})

test_that("val.prob.ci.2 warns when allowPerfectPredictions = TRUE and p contains 0/1", {
  d <- local_binary_preds()
  p_perf <- d$p
  p_perf[1] <- 0
  p_perf[2] <- 1
  expect_warning(
    val.prob.ci.2(p_perf, d$y, pl = FALSE, allowPerfectPredictions = TRUE),
    "replaced"
  )
})

test_that("val.prob.ci.2 errors on p = 0 with allowPerfectPredictions = FALSE", {
  d <- local_binary_preds()
  p_perf <- d$p
  p_perf[1] <- 0
  expect_error(
    val.prob.ci.2(p_perf, d$y, pl = FALSE, allowPerfectPredictions = FALSE),
    "Probabilities can not be >= 1 or <= 0"
  )
})

test_that("changing cl.level changes CI widths", {
  d <- local_binary_preds()
  res95 <- val.prob.ci.2(d$p, d$y, pl = FALSE, cl.level = 0.95)
  res80 <- val.prob.ci.2(d$p, d$y, pl = FALSE, cl.level = 0.80)

  width95 <- res95$Cindex[3] - res95$Cindex[1]
  width80 <- res80$Cindex[3] - res80$Cindex[1]
  expect_true(width95 > width80)
})

test_that("different method.ci values give same point estimate", {
  d <- local_binary_preds()
  res_pepe   <- val.prob.ci.2(d$p, d$y, pl = FALSE, method.ci = "pepe")
  res_delong <- val.prob.ci.2(d$p, d$y, pl = FALSE, method.ci = "delong")

  expect_equal(unname(res_pepe$stats["C (ROC)"]),
               unname(res_delong$stats["C (ROC)"]))
  expect_false(identical(res_pepe$Cindex, res_delong$Cindex))
})

test_that("val.prob.ci.2 with logistic calibration curve", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = TRUE, logistic.cal = TRUE)
  expect_s3_class(res, "CalibrationCurve")
})

test_that("val.prob.ci.2 stats match between pl = TRUE and pl = FALSE", {
  d <- local_binary_preds()
  res_pl  <- val.prob.ci.2(d$p, d$y, pl = TRUE)
  res_npl <- val.prob.ci.2(d$p, d$y, pl = FALSE)
  expect_equal(res_pl$stats, res_npl$stats)
})
