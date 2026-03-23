local_binary_preds <- function(env = parent.frame()) {
  data("traindata", package = "CalibrationCurves", envir = env)
  data("testdata", package = "CalibrationCurves", envir = env)
  fit <- glm(y ~ ., data = env$traindata, family = binomial)
  p <- predict(fit, newdata = env$testdata, type = "response")
  y <- env$testdata$y
  list(p = unname(p), y = y)
}

test_that("print.CalibrationCurve works and returns invisibly", {
  d <- local_binary_preds()
  res <- val.prob.ci.2(d$p, d$y, pl = FALSE)

  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl("Call:", out)))
  expect_true(any(grepl("confidence interval", out)))
})

test_that("print.ggplotCalibrationCurve works and returns invisibly", {
  d <- local_binary_preds()
  res <- valProbggplot(d$p, d$y, pl = TRUE)

  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl("Call:", out)))
})

test_that("print.GeneralizedCalibrationCurve works and returns invisibly", {
  data("poissontraindata", package = "CalibrationCurves")
  data("poissontestdata", package = "CalibrationCurves")
  fit <- glm(Y ~ ., data = poissontraindata, family = poisson)
  yHat <- predict(fit, newdata = poissontestdata, type = "response")

  res <- genCalCurve(poissontestdata$Y, yHat, family = "poisson", plot = FALSE)

  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl("Call:", out)))
  expect_true(any(grepl("confidence interval", out)))
})

test_that("print.SurvivalCalibrationCurve works and returns invisibly", {
  data("trainDataSurvival", package = "CalibrationCurves")
  data("testDataSurvival", package = "CalibrationCurves")

  sFit <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
                data = trainDataSurvival, x = TRUE, y = TRUE)
  res <- valProbSurvival(sFit, testDataSurvival, plotCal = "none")

  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl("Call:", out)))
  expect_true(any(grepl("Calibration performance", out)))
})

test_that("print.ClusteredCalibrationCurve works without error", {
  skip_on_cran()
  data("clustertraindata", package = "CalibrationCurves")
  data("clustertestdata", package = "CalibrationCurves")

  mFit <- lme4::glmer(y ~ x1 + x2 + x3 + x5 + (1 | cluster),
                       data = clustertraindata, family = "binomial")
  preds <- predict(mFit, clustertestdata, type = "response", re.form = NA)

  res <- suppressWarnings(
    valProbCluster(p = preds, y = clustertestdata$y,
                   cluster = clustertestdata$cluster,
                   plot = TRUE, approach = "MIXC", grid_l = 50)
  )

  out <- capture.output(ret <- print(res))
  expect_true(any(grepl("Call:", out)))
})
