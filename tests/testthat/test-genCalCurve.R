test_that("genCalCurve returns correct structure for Poisson family", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  res <- genCalCurve(yOOS, yHat, family = "poisson", plot = FALSE)

  expect_s3_class(res, "GeneralizedCalibrationCurve")
  expect_named(res, c("call", "stats", "cl.level", "Calibration",
                       "warningMessages", "CalibrationCurves"),
               ignore.order = TRUE)
})

test_that("genCalCurve stats contain calibration intercept and slope", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  res <- genCalCurve(yOOS, yHat, family = "poisson", plot = FALSE)

  expect_true("Calibration intercept" %in% names(res$stats))
  expect_true("Calibration slope" %in% names(res$stats))
})

test_that("genCalCurve works with binomial family", {
  data("traindata", package = "CalibrationCurves")
  data("testdata", package = "CalibrationCurves")

  fit <- glm(y ~ ., data = traindata, family = binomial)
  yHat <- predict(fit, newdata = testdata, type = "response")
  yOOS <- testdata$y

  res <- genCalCurve(yOOS, yHat, family = "binomial", plot = FALSE)
  expect_s3_class(res, "GeneralizedCalibrationCurve")
})

test_that("genCalCurve errors on invalid posStats", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  expect_error(
    genCalCurve(yOOS, yHat, family = "poisson", plot = FALSE, posStats = c(1, 2, 3)),
    "Length"
  )
})

test_that("genCalCurve works with Smooth = TRUE and GLMCal = FALSE", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  res <- genCalCurve(yOOS, yHat, family = "poisson", plot = FALSE,
                     Smooth = TRUE, GLMCal = FALSE)
  expect_s3_class(res, "GeneralizedCalibrationCurve")
})

test_that("genCalCurve produces a plot without error", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  expect_no_error(genCalCurve(yOOS, yHat, family = "poisson", plot = TRUE))
})

test_that("genCalCurve with Smooth and pointwise confidence limits", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")

  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")
  yOOS <- poissontestdata$Y

  res <- genCalCurve(yOOS, yHat, family = "poisson", plot = FALSE,
                     Smooth = TRUE, confLimitsSmooth = "pointwise")
  expect_s3_class(res, "GeneralizedCalibrationCurve")
})
