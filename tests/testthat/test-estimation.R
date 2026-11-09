test_that("skellam.mle() does not change options(warn)  [E1]", {
  old <- options(warn = 0); on.exit(options(old))
  set.seed(42); x <- rpois(200, 6) - rpois(200, 4)
  skellam.mle(x)
  expect_equal(getOption("warn"), 0)
})

test_that("skellam.reg() does not change options(warn), even if the fit fails  [E1]", {
  old <- options(warn = 0); on.exit(options(old))
  set.seed(42); y <- rpois(200, 6) - rpois(200, 4)
  try(suppressWarnings(skellam.reg(y, rnorm(200))), silent = TRUE)
  expect_equal(getOption("warn"), 0)
})

test_that("skellam.mle recovers the rates and keeps its return shape", {
  set.seed(11)
  x <- rpois(5000, 10) - rpois(5000, 6)
  f <- skellam.mle(x)
  expect_named(f, c("iters", "loglik", "param"))
  expect_named(f$param, c("mu 1", "mu 2"))
  expect_equal(unname(f$param[1] - f$param[2]), mean(x), tolerance = 1e-8)   # MLE matches the mean
  expect_lt(abs(f$param[[1]] - 10), 3 * sqrt(2 * 16 / 5000) * 4)
  expect_equal(f$loglik, sum(dskellam(x, f$param[[1]], f$param[[2]], log = TRUE)), tolerance = 1e-12)
  # the maximum: nudging lambda2 either way (with lambda1 tracking) lowers the log-likelihood
  ll <- function(d) sum(dskellam(x, f$param[[1]] + d, f$param[[2]] + d, log = TRUE))
  expect_lt(ll(0.01), f$loglik); expect_lt(ll(-0.01), f$loglik)
})

test_that("skellam.mle returns non-negative rates for under-dispersed samples  [E2]", {
  f <- skellam.mle(c(rep(2, 45), rep(3, 5)))     # var 0.09 < mean 2.1
  expect_true(all(f$param >= 0))
  expect_equal(unname(f$param), c(2.1, 0))
})

test_that("skellam.mle handles degenerate samples with the boundary estimate  [E3]", {
  expect_equal(unname(skellam.mle(rep(0, 50))$param), c(0, 0))
  expect_equal(skellam.mle(rep(0, 50))$loglik, 0)
  expect_equal(unname(skellam.mle(rep(3, 50))$param), c(3, 0))
  expect_equal(unname(skellam.mle(rep(-3, 50))$param), c(0, 3))
})

test_that("skellam.mle rejects invalid input with clear messages  [E3]", {
  expect_error(skellam.mle(c(1, NA, 2)), "missing values")
  expect_error(skellam.mle(c(1.5, 2)), "integers")
  expect_error(skellam.mle(3), "two observations")
  expect_error(skellam.mle("a"), "numeric")
})

test_that("skellam.mle works at very large and very small rates", {
  set.seed(9)
  x <- rpois(500, 1e5) - rpois(500, 1e5 - 50)
  f <- skellam.mle(x)
  expect_true(all(is.finite(f$param)) && all(f$param > 0))
  expect_equal(unname(f$param[1] - f$param[2]), mean(x), tolerance = 1e-8)
  x <- rpois(2000, 0.05) - rpois(2000, 0.03)
  expect_true(all(is.finite(skellam.mle(x)$param)))
})

# ---- skellam.reg -------------------------------------------------------------

sim_reg <- function(n, b1, b2, seed) {
  set.seed(seed)
  x <- rnorm(n)
  y <- rpois(n, exp(b1[1] + b1[2] * x)) - rpois(n, exp(b2[1] + b2[2] * x))
  list(x = x, y = y)
}

test_that("skellam.reg: param1 is the rate of the added counts, param2 of the subtracted  [E7]", {
  d <- sim_reg(3000, c(1.5, 0.5), c(0.2, -0.3), seed = 1)
  r <- skellam.reg(d$y, d$x)
  expect_lt(max(abs(r$param1[, 1] - c(1.5, 0.5)) / r$param1[, 2]), 4)
  expect_lt(max(abs(r$param2[, 1] - c(0.2, -0.3)) / r$param2[, 2]), 4)
})

test_that("skellam.reg keeps its return shape", {
  d <- sim_reg(300, c(1, 0.3), c(0.5, 0), seed = 2)
  r <- skellam.reg(d$y, d$x)
  expect_named(r, c("loglik", "param1", "param2"))
  expect_identical(colnames(r$param1), c("Estimate", "Std. Error", "Wald value", "p-value"))
  expect_identical(rownames(r$param1), c("(Intercept)", "x"))
  b1 <- r$param1[, 1]; b2 <- r$param2[, 1]
  expect_equal(r$loglik, sum(dskellam(d$y, exp(b1[1] + b1[2] * d$x), exp(b2[1] + b2[2] * d$x), log = TRUE)),
               tolerance = 1e-10)
})

test_that("skellam.reg works when the rates are in the hundreds  [E5]", {
  d <- sim_reg(300, c(6, 0.2), c(5.9, 0.2), seed = 7)
  expect_silent(r <- skellam.reg(d$y, d$x))
  expect_true(is.finite(r$loglik))
  expect_lt(abs(r$param1[1, 1] - 6), 0.3)
  expect_lt(abs(r$param2[1, 1] - 5.9), 0.3)
})

test_that("skellam.reg is reproducible: identical results on repeated calls  [E6]", {
  d <- sim_reg(500, c(1.8, 0), c(1.4, 0), seed = 42)
  set.seed(1); a <- skellam.reg(d$y, d$x)
  set.seed(2); b <- skellam.reg(d$y, d$x)
  expect_identical(a, b)
})

test_that("skellam.reg removes incomplete observations with a warning and handles factors", {
  set.seed(3)
  dd <- data.frame(a = rnorm(400), g = factor(sample(c("u", "v", "w"), 400, TRUE)))
  y <- rpois(400, exp(1 + 0.3 * dd$a + (dd$g == "v"))) - rpois(400, 2)
  r <- skellam.reg(y, dd)
  expect_identical(rownames(r$param1), c("(Intercept)", "a", "gv", "gw"))
  dd$a[3] <- NA; y[5] <- NA
  expect_warning(skellam.reg(y, dd), "2 observation")
})

test_that("skellam.reg rejects invalid input", {
  expect_error(skellam.reg(c(1.5, 2, 3, 4, 5, 6), 1:6), "integers")
  expect_error(skellam.reg(1:5, 1:6), "same number")
  expect_error(skellam.reg(1:4, 1:4), "too few")
})
