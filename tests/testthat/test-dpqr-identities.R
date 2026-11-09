# Properties that already hold in 0.2.4; these tests keep them from regressing.

test_that("pmf sums to 1 for moderate rates", {
  for (L in list(c(1, 1), c(3, 2), c(0.1, 5), c(50, 50), c(1e3, 1e3), c(1e4, 1e4))) {
    x <- support_grid(L[1], L[2])
    expect_equal(sum(dskellam(x, L[1], L[2])), 1, tolerance = 1e-13)
  }
})

test_that("pmf agrees with the exact Poisson-product sum", {
  for (L in list(c(0.01, 0.01), c(3, 2), c(100, 80), c(1e3, 1e3), c(50, 0.001))) {
    s <- sqrt(sum(L)); m <- L[1] - L[2]
    x <- unique(round(m + s * c(-8, -3, 0, 3, 8)))
    expect_equal(dskellam(x, L[1], L[2], log = TRUE),
                 vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]), tolerance = 1e-12)
  }
})

test_that("cdf equals the cumulative sum of the pmf", {
  for (L in list(c(3, 2), c(0.1, 5), c(1e3, 1e3))) {
    x <- support_grid(L[1], L[2])
    expect_equal(pskellam(x, L[1], L[2]), cumsum(dskellam(x, L[1], L[2])), tolerance = 1e-12)
  }
})

test_that("lower and upper tails add up to 1", {
  x <- -10:10
  expect_equal(pskellam(x, 3, 2) + pskellam(x, 3, 2, lower.tail = FALSE), rep(1, 21), tolerance = 1e-15)
})

test_that("log scales agree with the natural scale", {
  x <- -10:10
  expect_equal(dskellam(x, 3, 2, log = TRUE), log(dskellam(x, 3, 2)), tolerance = 1e-14)
  expect_equal(pskellam(x, 3, 2, log.p = TRUE), log(pskellam(x, 3, 2)), tolerance = 1e-14)
})

test_that("lower-tail quantile inverts the cdf", {
  for (L in list(c(3, 2), c(1e3, 1e3), c(1e4, 1), c(1, 0))) {
    x <- support_grid(L[1], L[2], w = 6)
    if (L[2] == 0) x <- x[x >= 0]
    expect_identical(qskellam(pskellam(x, L[1], L[2]), L[1], L[2]), as.numeric(x))
  }
})

test_that("a zero rate reduces to the Poisson", {
  expect_equal(dskellam(0:10, 5, 0), dpois(0:10, 5))
  expect_equal(dskellam(-(0:10), 0, 5), dpois(0:10, 5))
})

test_that("rskellam has the right mean and variance", {
  set.seed(20261108)
  z <- rskellam(1e5, 7, 3)
  expect_lt(abs(mean(z) - 4) / sqrt(10 / 1e5), 5)
  expect_lt(abs(var(z) - 10) / (10 * sqrt(2 / 1e5)), 5)
})
