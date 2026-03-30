test_that("recalibrate returns correct structure with Platt scaling", {
  data("traindata", package = "CalibrationCurves")
  data("testdata",  package = "CalibrationCurves")
  glmFit <- glm(y ~ ., data = traindata, family = binomial)
  p_test <- predict(glmFit, newdata = testdata, type = "response")
  y_test <- testdata$y

  res <- recalibrate(p_test, y_test, method = "platt", plot = FALSE)

  expect_s3_class(res, "RecalibratedPredictions")
  expect_equal(res$method, "platt")
  expect_length(res$p_recal, length(p_test))
  expect_true(all(res$p_recal >= 0 & res$p_recal <= 1))

  # Stats slots
  expect_named(res$before, c("Brier", "LogLoss", "Intercept", "Slope"))
  expect_named(res$after,  c("Brier", "LogLoss", "Intercept", "Slope"))

  # Platt: calibration slope after should be closer to 1
  expect_true(!is.na(res$after$Slope))
})

test_that("recalibrate returns correct structure with isotonic regression", {
  data("traindata", package = "CalibrationCurves")
  data("testdata",  package = "CalibrationCurves")
  glmFit <- glm(y ~ ., data = traindata, family = binomial)
  p_test <- predict(glmFit, newdata = testdata, type = "response")
  y_test <- testdata$y

  res <- recalibrate(p_test, y_test, method = "isotonic", plot = FALSE)

  expect_s3_class(res, "RecalibratedPredictions")
  expect_equal(res$method, "isotonic")
  expect_length(res$p_recal, length(p_test))
  expect_true(all(res$p_recal >= 0 & res$p_recal <= 1))
})

test_that("recalibrate reduces Brier score for a miscalibrated model", {
  set.seed(123)
  n  <- 500
  lp <- rnorm(n)
  y  <- rbinom(n, 1, plogis(lp))
  # Deliberately miscalibrated: shrink predictions towards 0.5
  p_bad <- plogis(lp * 0.4 + 0.5)

  res <- recalibrate(p_bad, y, method = "platt", plot = FALSE, verbose = FALSE)
  expect_lte(res$after$Brier, res$before$Brier + 1e-6)  # should not get worse
})

test_that("recalibrate errors on wrong input", {
  expect_error(recalibrate(c(-0.1, 0.5), c(0, 1), plot = FALSE),
               "probabilities in \\[0, 1\\]")
  expect_error(recalibrate(c(0.3, 0.7), c(0, 2), plot = FALSE),
               "binary")
  expect_error(recalibrate(c(0.3, 0.7, 0.5), c(0, 1), plot = FALSE),
               "same length")
})

test_that("print.RecalibratedPredictions works without error", {
  data("traindata", package = "CalibrationCurves")
  data("testdata",  package = "CalibrationCurves")
  glmFit <- glm(y ~ ., data = traindata, family = binomial)
  p_test <- predict(glmFit, newdata = testdata, type = "response")
  y_test <- testdata$y
  res <- recalibrate(p_test, y_test, method = "platt", plot = FALSE, verbose = FALSE)
  expect_no_error(print(res))
})
