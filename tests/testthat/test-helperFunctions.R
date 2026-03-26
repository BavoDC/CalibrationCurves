# Access internal functions
Ilogit_fn <- CalibrationCurves:::Ilogit
Logit_fn  <- CalibrationCurves:::Logit
loess_as  <- CalibrationCurves:::loess.as

test_that("Logit and Ilogit are inverse of each other", {
  p_vals <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  expect_equal(Ilogit_fn(Logit_fn(p_vals)), p_vals, tolerance = 1e-12)
})

test_that("Logit(0.5) = 0", {
  expect_equal(Logit_fn(0.5), 0)
})

test_that("Ilogit(0) = 0.5", {
  expect_equal(Ilogit_fn(0), 0.5)
})

test_that("Logit handles boundaries", {
  expect_equal(Logit_fn(0), -Inf)
  expect_equal(Logit_fn(1), Inf)
})

test_that("Ilogit handles extreme values", {
  expect_equal(Ilogit_fn(-Inf), 0)
  expect_equal(Ilogit_fn(Inf), 1)
})

test_that("LibraryM loads packages without error", {
  expect_no_error(LibraryM(stats))
  expect_invisible(LibraryM(stats))
})

test_that("loess.as returns a loess object", {
  set.seed(7391)
  x <- seq(0, 2 * pi, length.out = 100)
  y <- sin(x) + rnorm(100, sd = 0.2)

  fit <- loess_as(x, y)
  expect_s3_class(fit, "loess")
})

test_that("loess.as errors on missing values in x", {
  x <- c(1, 2, NA, 4)
  y <- c(1, 2, 3, 4)
  expect_error(loess_as(x, y), "missing values")
})

test_that("loess.as errors on missing values in y", {
  x <- c(1, 2, 3, 4)
  y <- c(1, NA, 3, 4)
  expect_error(loess_as(x, y), "missing values")
})

test_that("loess.as respects user.span", {
  set.seed(7391)
  x <- seq(0, 2 * pi, length.out = 100)
  y <- sin(x) + rnorm(100, sd = 0.2)

  fit <- loess_as(x, y, user.span = 0.5)
  expect_s3_class(fit, "loess")
  expect_equal(fit$pars$span, 0.5)
})

test_that("%<=% and %{}% are exported functions", {
  expect_true(is.function(`%<=%`))
  expect_true(is.function(`%{}%`))
})
