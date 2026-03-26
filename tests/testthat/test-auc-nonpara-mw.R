# auc.nonpara.mw is not exported, access via :::
auc_fn <- CalibrationCurves:::auc.nonpara.mw
getAUCmw_fn <- CalibrationCurves:::getAUCmw

test_that("auc.nonpara.mw returns a length-3 numeric vector", {
  set.seed(4821)
  x <- rnorm(50, mean = 1)
  y <- rnorm(50, mean = 0)

  for (m in c("newcombe", "pepe", "delong", "jackknife")) {
    res <- auc_fn(x, y, method = m)
    expect_length(res, 3)
    expect_type(res, "double")
  }
})

test_that("all methods give the same AUC point estimate", {
  set.seed(4821)
  x <- rnorm(50, mean = 1)
  y <- rnorm(50, mean = 0)

  estimates <- vapply(
    c("newcombe", "pepe", "delong", "jackknife"),
    function(m) auc_fn(x, y, method = m)[1],
    numeric(1)
  )
  # All point estimates should be equal (values, ignore names)
  expect_equal(unname(estimates), rep(unname(estimates[1]), length(estimates)),
               tolerance = 1e-10)
})

test_that("CI ordering: lower <= point <= upper", {
  set.seed(4821)
  x <- rnorm(50, mean = 1)
  y <- rnorm(50, mean = 0)

  for (m in c("newcombe", "pepe", "delong", "jackknife")) {
    res <- auc_fn(x, y, method = m)
    expect_true(res[2] <= res[1], label = paste(m, "lower <= point"))
    expect_true(res[1] <= res[3], label = paste(m, "point <= upper"))
  }
})

test_that("perfect discrimination gives AUC = 1", {
  x <- c(5, 6, 7, 8)
  y <- c(1, 2, 3, 4)

  auc <- getAUCmw_fn(x, y)
  expect_equal(auc, 1)
})

test_that("no discrimination gives AUC = 0.5", {
  x <- c(1, 2, 3, 4)
  y <- c(1, 2, 3, 4)

  auc <- getAUCmw_fn(x, y)
  expect_equal(auc, 0.5)
})

test_that("getAUCmw with known values", {
  # x always > y → AUC = 1
  expect_equal(getAUCmw_fn(c(3, 4, 5), c(1, 2)), 1)

  # x always < y → AUC = 0
  expect_equal(getAUCmw_fn(c(1, 2), c(3, 4, 5)), 0)
})
