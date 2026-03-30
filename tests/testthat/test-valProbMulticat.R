test_that("valProbMulticat works with multinom fit", {
  skip_if_not_installed("nnet")
  library(nnet)
  set.seed(42)
  n    <- 300
  x    <- rnorm(n)
  prbs <- cbind(exp(0.3 * x), exp(-0.2 * x), rep(1, n))
  prbs <- prbs / rowSums(prbs)
  y    <- apply(prbs, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  df   <- data.frame(y = factor(y), x = x)
  idx  <- seq_len(200)
  trainD <- df[idx, ]
  testD  <- df[-idx, ]

  mFit <- multinom(y ~ x, data = trainD, trace = FALSE)
  res  <- valProbMulticat(mFit, valdata = testD, plot = FALSE, plotCal = "none")

  expect_s3_class(res, "MulticlassCalibrationCurve")
  expect_equal(res$type, "multinomial")
  expect_true("BrierScore" %in% names(res$stats))
  expect_true("LogLoss"    %in% names(res$stats))
  expect_true(!is.null(res$Calibration$Slopes))
})

test_that("valProbMulticat works with polr fit", {
  skip_if_not_installed("MASS")
  library(MASS)

  # Generate synthetic 3-level ordinal data (polr requires >= 3 levels)
  set.seed(7)
  n  <- 400
  df <- local({
    x  <- rnorm(n)
    lp <- 0.8 * x
    y  <- cut(lp + rnorm(n, 0, 1.2),
              breaks = c(-Inf, -1, 0.5, Inf),
              labels = c("low", "med", "high"),
              ordered_result = TRUE)
    data.frame(y = y, x = x)
  })

  oFit <- polr(y ~ x, data = df[1:300, ], Hess = TRUE)
  res  <- valProbMulticat(oFit, valdata = df[301:400, ], plot = FALSE, plotCal = "none")

  expect_s3_class(res, "MulticlassCalibrationCurve")
  expect_equal(res$type, "ordinal")
  expect_true(res$stats$BrierScore["Overall"] >= 0)
  expect_true(res$stats$LogLoss >= 0)
})

test_that("valProbMulticat works with probability matrix input", {
  skip_if_not_installed("nnet")
  library(nnet)
  set.seed(7)
  n    <- 200
  x    <- rnorm(n)
  prbs <- cbind(exp(0.5 * x), exp(-0.5 * x), rep(1, n))
  prbs <- prbs / rowSums(prbs)
  y    <- apply(prbs, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  colnames(prbs) <- c("A", "B", "C")

  res <- valProbMulticat(prbs, y = factor(y), plot = FALSE, plotCal = "none")
  expect_s3_class(res, "MulticlassCalibrationCurve")
})

test_that("valProbMulticat errors on wrong input", {
  expect_error(valProbMulticat(matrix(1:4, 2, 2), y = NULL),
               "'y' must be provided")
})

test_that("print.MulticlassCalibrationCurve works", {
  skip_if_not_installed("nnet")
  library(nnet)
  set.seed(1)
  n    <- 200
  x    <- rnorm(n)
  # Use 3 classes so predict() returns a matrix (2-class returns a vector)
  prbs <- cbind(exp(0.3 * x), exp(-0.3 * x), rep(1, n))
  prbs <- prbs / rowSums(prbs)
  y    <- apply(prbs, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  df   <- data.frame(y = factor(y), x = x)
  mFit <- multinom(y ~ x, data = df[1:150, ], trace = FALSE)
  res  <- valProbMulticat(mFit, valdata = df[151:200, ], plot = FALSE, plotCal = "none")
  expect_no_error(print(res))
})
