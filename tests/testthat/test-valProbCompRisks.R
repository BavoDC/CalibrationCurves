# Tests for valProbCompRisks
# Note: CSC/FGR require multi-state data; we simulate a minimal example.

test_that("valProbCompRisks errors on non-CSC/FGR fit", {
  data("testDataSurvival", package = "CalibrationCurves")
  bad_fit <- lm(ryear ~ csize + cnode, data = testDataSurvival)
  expect_error(
    valProbCompRisks(bad_fit, testDataSurvival, plotCal = "none"),
    "CSC or FGR"
  )
})

test_that("valProbCompRisks works with CSC model", {
  skip_if_not_installed("riskRegression")
  library(riskRegression)
  library(survival)

  set.seed(123)
  n_obs <- 400L
  all_df <- data.frame(
    x1     = rnorm(n_obs, 55, 10),
    time   = pmin(rexp(n_obs, 0.1), 10),
    status = sample(0:2, n_obs, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  )
  trainD <- all_df[seq_len(250L), ]
  testD  <- all_df[-seq_len(250L), ]

  # Must assign x1 from trainD to the calling frame: riskRegression uses
  # parent.frame() to resolve formula variables when calling coxph internally.
  # Having x1 here with length == nrow(trainD) avoids 'variable lengths differ'.
  x1 <- trainD$x1
  cscFit <- CSC(Hist(time, status) ~ x1, data = trainD)

  res <- valProbCompRisks(cscFit, testD, cause = 1,
                          timeHorizon = 5, plotCal = "none")

  expect_s3_class(res, "CompRisksCalibrationCurve")
  expect_equal(res$cause, 1)
  expect_equal(res$timeHorizon, 5)
  expect_true(!is.null(res$stats$Calibration$InTheLarge))
})

test_that("print.CompRisksCalibrationCurve works", {
  skip_if_not_installed("riskRegression")
  library(riskRegression)
  library(survival)

  set.seed(99)
  n_obs <- 300L
  all_df <- data.frame(
    x1     = rnorm(n_obs, 60, 10),
    time   = pmin(rexp(n_obs, 0.15), 8),
    status = sample(0:2, n_obs, replace = TRUE, prob = c(0.4, 0.4, 0.2))
  )

  x1     <- all_df[seq_len(200L), "x1"]
  cscFit <- CSC(Hist(time, status) ~ x1, data = all_df[seq_len(200L), ])
  res    <- valProbCompRisks(cscFit, all_df[201:300, ], cause = 1,
                             timeHorizon = 4, plotCal = "none")
  expect_no_error(print(res))
})
