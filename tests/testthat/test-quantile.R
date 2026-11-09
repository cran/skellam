# For a discrete distribution the quantile is defined by an inequality:
#   lower tail: smallest integer x with P(X <= x) >= p
#   upper tail: smallest integer x with P(X >  x) <= p
# so feeding a cdf value back in must return the same x ("round trip").

roundtrip <- function(x, l1, l2, lower.tail = TRUE, log.p = FALSE) {
  qskellam(pskellam(x, l1, l2, lower.tail = lower.tail, log.p = log.p),
           l1, l2, lower.tail = lower.tail, log.p = log.p)
}

test_that("round trip holds on both tails and both scales  [N3]", {
  for (L in list(c(3, 2), c(0.5, 4), c(1e3, 1e3), c(1e4, 1))) {
    x <- support_grid(L[1], L[2], w = 5)
    for (lt in c(TRUE, FALSE)) for (lg in c(FALSE, TRUE)) {
      expect_identical(roundtrip(x, L[1], L[2], lt, lg), as.numeric(x),
                       label = sprintf("lambda=(%g,%g) lower.tail=%s log.p=%s", L[1], L[2], lt, lg))
    }
  }
})

test_that("upper and lower tails agree: q(p, upper) == q(1 - p, lower)", {
  p <- c(0.01, 0.1, 0.3, 0.77)
  expect_identical(qskellam(p, 3, 2, lower.tail = FALSE), qskellam(1 - p, 3, 2))
})

test_that("p = 0 and p = 1 give the ends of the support  [N2]", {
  expect_identical(qskellam(c(0, 0.5, 1), 3, 2), c(-Inf, 1, Inf))
  expect_identical(qskellam(c(0, 0.5, 1), 3, 2, lower.tail = FALSE), c(Inf, 1, -Inf))
  expect_identical(qskellam(c(-Inf, log(0.5), 0), 3, 2, log.p = TRUE), c(-Inf, 1, Inf))
  expect_identical(qskellam(c(0, 1), 3, 0), qpois(c(0, 1), 3))       # X >= 0
  expect_identical(qskellam(c(0, 1), 0, 3), c(-Inf, 0))               # X <= 0
  expect_identical(qskellam(c(0, 0.5, 1), 0, 0), c(0, 0, 0))          # X = 0
})

test_that("a zero rate gives the (reflected) Poisson quantile  [N4]", {
  x <- -12:0
  expect_identical(qskellam(pskellam(x, 0, 1), 0, 1), as.numeric(x))
  expect_identical(qskellam(pskellam(x, 0, 1, lower.tail = FALSE), 0, 1, lower.tail = FALSE), as.numeric(x))
  x <- 0:12
  expect_identical(qskellam(pskellam(x, 2, 0), 2, 0), qpois(ppois(x, 2), 2))
})

test_that("invalid input gives NA/NaN, not errors or Inf  [N12]", {
  expect_warning(q <- qskellam(1.1, 3, 2), "NaNs produced");  expect_true(is.nan(q))
  expect_warning(q <- qskellam(-0.1, 3, 2), "NaNs produced"); expect_true(is.nan(q))
  expect_warning(q <- qskellam(0.5, 3, 2, log.p = TRUE), "NaNs produced"); expect_true(is.nan(q))
  expect_warning(q <- qskellam(0.5, -1, 2), "NaNs produced"); expect_true(is.nan(q))
  expect_identical(qskellam(NA, 3, 2), NA_real_)
  expect_identical(qskellam(c(0.5, NA), 3, 2), c(1, NA))
  expect_identical(qskellam(numeric(0), 3, 2), numeric(0))
})

test_that("recycling warning and names  [N13, A10]", {
  expect_warning(qskellam(c(.1, .5, .9), 1:2, 1), "multiple")
  expect_named(qskellam(c(a = 0.5), 3, 2), "a")
})

test_that("options(verbose) does not change the result  [N14]", {
  old <- options(verbose = TRUE); on.exit(options(old))
  v <- qskellam(c(0, 0.5, 1), 3, 0)
  options(verbose = FALSE)
  expect_identical(v, qskellam(c(0, 0.5, 1), 3, 0))
  expect_identical(v, c(0, 3, Inf))
})
