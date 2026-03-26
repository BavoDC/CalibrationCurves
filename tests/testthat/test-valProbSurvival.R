# Note: coxph must be called with Surv() (not survival::Surv()) after
# library(survival) is loaded, because riskRegression::Score() does not
# handle namespace-prefixed formulas correctly.

test_that("valProbSurvival returns correct structure", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)

  res <- valProbSurvival(sFit, testDataSurvival, plotCal = "none")

  expect_s3_class(res, "SurvivalCalibrationCurve")
  expect_named(res, c("call", "stats", "alpha", "Calibration",
                       "CalibrationCurves"),
               ignore.order = TRUE)
})

test_that("valProbSurvival stats contain expected components", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)

  res <- valProbSurvival(sFit, testDataSurvival, plotCal = "none")

  expect_true("Calibration" %in% names(res$stats))
  expect_true("Concordance" %in% names(res$stats))
  expect_true("TimeDependentAUC" %in% names(res$stats))
})

test_that("valProbSurvival errors on non-coxph fit", {
  data("testDataSurvival", package = "CalibrationCurves")

  bad_fit <- lm(ryear ~ csize + cnode, data = testDataSurvival)

  expect_error(
    valProbSurvival(bad_fit, testDataSurvival, plotCal = "none"),
    "coxph"
  )
})

test_that("valProbSurvival works with plotCal = 'base'", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)

  expect_no_error(
    valProbSurvival(sFit, testDataSurvival, plotCal = "base")
  )
})

test_that("valProbSurvival works with plotCal = 'ggplot'", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)

  expect_no_error(
    valProbSurvival(sFit, testDataSurvival, plotCal = "ggplot")
  )
})

test_that("valProbSurvival respects timeHorizon argument", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)

  res3 <- valProbSurvival(sFit, testDataSurvival, plotCal = "none", timeHorizon = 3)
  res5 <- valProbSurvival(sFit, testDataSurvival, plotCal = "none", timeHorizon = 5)

  # Different time horizons should produce different calibration stats
  expect_false(identical(res3$stats, res5$stats))
})
