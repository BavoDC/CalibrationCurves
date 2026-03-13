local_binary_preds <- function(env = parent.frame()) {
  data("traindata", package = "CalibrationCurves", envir = env)
  data("testdata", package = "CalibrationCurves", envir = env)
  fit <- glm(y ~ ., data = env$traindata, family = binomial)
  p <- predict(fit, newdata = env$testdata, type = "response")
  y <- env$testdata$y
  list(p = unname(p), y = y)
}

test_that("valProbggplot returns correct structure", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = FALSE)

  expect_s3_class(res, "ggplotCalibrationCurve")
  expect_named(res, c("call", "ggPlot", "stats", "cl.level", "Calibration",
                       "Cindex", "warningMessages", "CalibrationCurves"),
               ignore.order = TRUE)
})

test_that("valProbggplot ggPlot slot is NULL when pl = FALSE", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = FALSE)
  expect_null(res$ggPlot)
})

test_that("valProbggplot ggPlot slot is a ggplot object when pl = TRUE", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = TRUE)
  expect_s3_class(res$ggPlot, "ggplot")
})

test_that("valProbggplot stats match val.prob.ci.2 stats", {
  d <- local_binary_preds()
  res_gg   <- valProbggplot(d$p, d$y, pl = FALSE)
  res_base <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  expect_equal(unname(res_gg$stats), unname(res_base$stats))
})

test_that("valProbggplot works with smooth = 'rcs'", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = FALSE, smooth = "rcs")
  expect_s3_class(res, "ggplotCalibrationCurve")
})

test_that("valProbggplot works with smooth = 'none'", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = FALSE, smooth = "none")
  expect_s3_class(res, "ggplotCalibrationCurve")
})

test_that("valProbggplot errors on non-binary y", {
  d <- local_binary_preds()
  expect_error(valProbggplot(d$p, d$y + 0.5, pl = FALSE),
               "binary outcome")
})

test_that("valProbggplot errors on mismatched lengths", {
  d <- local_binary_preds()
  expect_error(valProbggplot(d$p[1:10], d$y, pl = FALSE),
               "lengths")
})

test_that("valProbggplot C-statistic and Brier are in valid ranges", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = FALSE)

  expect_true(res$stats["C (ROC)"] >= 0 && res$stats["C (ROC)"] <= 1)
  expect_true(res$stats["Brier"] >= 0 && res$stats["Brier"] <= 1)
})

test_that("valProbggplot stats match between pl = TRUE and pl = FALSE", {
  d <- local_binary_preds()
  res_pl  <- valProbggplot(d$p, d$y, pl = TRUE)
  res_npl <- valProbggplot(d$p, d$y, pl = FALSE)
  expect_equal(res_pl$stats, res_npl$stats)
})
