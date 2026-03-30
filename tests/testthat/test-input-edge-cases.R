# Edge case and regression tests

test_that("constant predictions (all identical p) returns valid CalibrationCurve", {
  set.seed(6174)
  y <- rbinom(100, 1, 0.3)
  p <- rep(0.3, 100)

  # Constant predictions trigger a special code path (returns early)
  res <- val.prob.ci.2(p, y, pl = FALSE)
  expect_s3_class(res, "CalibrationCurve")
  # C-statistic should be 0.5 for uninformative model
  expect_equal(unname(res$stats["C (ROC)"]), 0.5)
})

test_that("small sample (n = 50) works", {
  set.seed(2847)
  y <- rbinom(50, 1, 0.4)
  p <- runif(50, 0.1, 0.9)

  res <- val.prob.ci.2(p, y, pl = FALSE)
  expect_s3_class(res, "CalibrationCurve")
  expect_length(res$stats, 18)
})

test_that("NAs in weights produce a warning and are handled", {
  set.seed(5012)
  y <- rbinom(100, 1, 0.4)
  p <- runif(100, 0.1, 0.9)
  w <- rep(1, 100)
  w[c(1, 5, 10)] <- NA

  expect_warning(
    res <- val.prob.ci.2(p, y, pl = FALSE, weights = w),
    "missing"
  )
  expect_s3_class(res, "CalibrationCurve")
})

test_that("val.prob.ci.2 with logistic calibration curve", {
  set.seed(3309)
  y <- rbinom(200, 1, 0.5)
  p <- runif(200, 0.05, 0.95)

  res <- val.prob.ci.2(p, y, pl = FALSE, logistic.cal = TRUE)
  expect_s3_class(res, "CalibrationCurve")
})

test_that("valProbggplot with logistic calibration curve", {
  set.seed(3309)
  y <- rbinom(200, 1, 0.5)
  p <- runif(200, 0.05, 0.95)

  res <- valProbggplot(p, y, pl = FALSE, logistic.cal = TRUE)
  expect_s3_class(res, "ggplotCalibrationCurve")
})

test_that("val.prob.ci.2 errors on all-zero predictions with allowPerfectPredictions", {
  y <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
  p <- rep(0, 10)

  expect_error(
    val.prob.ci.2(p, y, pl = FALSE, allowPerfectPredictions = TRUE),
    "deterministic"
  )
})

test_that("val.prob.ci.2 nr.knots validation", {
  set.seed(8816)
  y <- rbinom(200, 1, 0.5)
  p <- runif(200, 0.05, 0.95)

  expect_error(val.prob.ci.2(p, y, pl = FALSE, smooth = "rcs", nr.knots = 6))
  expect_error(val.prob.ci.2(p, y, pl = FALSE, smooth = "rcs", nr.knots = 2))
  # nr.knots = 3 should work
  res <- val.prob.ci.2(p, y, pl = FALSE, smooth = "rcs", nr.knots = 3)
  expect_s3_class(res, "CalibrationCurve")
})

test_that("val.prob.ci.2 dostats = FALSE still returns stats", {
  set.seed(1543)
  y <- rbinom(200, 1, 0.5)
  p <- runif(200, 0.05, 0.95)

  res <- val.prob.ci.2(p, y, pl = FALSE, dostats = FALSE)
  expect_s3_class(res, "CalibrationCurve")
  expect_true(length(res$stats) >= 15)
})

test_that("val.prob.ci.2 dostats with specific stats", {
  set.seed(1543)
  y <- rbinom(200, 1, 0.5)
  p <- runif(200, 0.05, 0.95)

  res <- val.prob.ci.2(p, y, pl = FALSE,
                       dostats = c("C (ROC)", "Intercept", "Slope"))
  expect_s3_class(res, "CalibrationCurve")
})
