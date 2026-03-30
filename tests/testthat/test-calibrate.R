# calibrationCurve() dispatches to the correct underlying function

test_that("calibrationCurve dispatches to valProbggplot for binary glm", {
  data("traindata", package = "CalibrationCurves")
  data("testdata",  package = "CalibrationCurves")
  glmFit <- glm(y ~ ., data = traindata, family = binomial)
  res <- calibrationCurve(glmFit, newdata = testdata, y = "y")
  expect_s3_class(res, "ggplotCalibrationCurve")
})

test_that("calibrationCurve dispatches to genCalCurve for Poisson glm", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata",  package = "CalibrationCurves")
  glmFit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  res <- calibrationCurve(glmFit, newdata = poissontestdata, y = "Y", plot = FALSE)
  expect_s3_class(res, "GeneralizedCalibrationCurve")
})

test_that("calibrationCurve dispatches to valProbSurvival for coxph", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival",  package = "CalibrationCurves")
  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)
  res <- calibrationCurve(sFit, newdata = testDataSurvival,
                          timeHorizon = 5, plotCal = "none")
  expect_s3_class(res, "SurvivalCalibrationCurve")
})

test_that("calibrationCurve errors on unsupported model", {
  m <- list(class = "unknown_model_xyz")
  class(m) <- "unknown_model_xyz"
  expect_error(calibrationCurve(m, newdata = data.frame(x = 1)))
})
