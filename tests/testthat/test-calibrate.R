# calibrate() dispatches to the correct underlying function

test_that("calibrate dispatches to valProbggplot for binary glm", {
  data("traindata", package = "CalibrationCurves")
  data("testdata",  package = "CalibrationCurves")
  glmFit <- glm(y ~ ., data = traindata, family = binomial)
  res <- calibrate(glmFit, newdata = testdata, y = "y")
  expect_s3_class(res, "ggplotCalibrationCurve")
})

test_that("calibrate dispatches to genCalCurve for Poisson glm", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata",  package = "CalibrationCurves")
  glmFit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  res <- calibrate(glmFit, newdata = poissontestdata, y = "Y", plot = FALSE)
  expect_s3_class(res, "GeneralizedCalibrationCurve")
})

test_that("calibrate dispatches to valProbSurvival for coxph", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival",  package = "CalibrationCurves")
  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)
  res <- calibrate(sFit, newdata = testDataSurvival,
                   timeHorizon = 5, plotCal = "none")
  expect_s3_class(res, "SurvivalCalibrationCurve")
})

test_that("calibrate errors on unsupported model", {
  m <- list(class = "unknown_model_xyz")
  class(m) <- "unknown_model_xyz"
  expect_error(calibrate(m, newdata = data.frame(x = 1)))
})
