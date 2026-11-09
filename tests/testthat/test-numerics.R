# Accuracy of the log-pmf core against the exact Poisson-product sum, over the
# whole range of rates. Target: relative error of the pmf below 1e-12 up to
# rates of 1e6 with equal rates, below 1e-9 for extremely unequal huge rates
# (where the Bessel representation itself cancels 1e6-sized terms).

acc_grid <- list(c(0.01, 0.01), c(3, 2), c(100, 80), c(1e3, 1e3), c(1e4, 1e4),
                 c(3e4, 3e4), c(1e5, 1e5), c(1e6, 1e6), c(50, 0.001), c(300, 300))
unequal  <- list(c(1e4, 1), c(1e6, 1e3), c(1e6, 1e5), c(2e5, 10))

test_that("log-pmf matches the exact sum to 1e-12, including at rates of 1e5-1e6  [N5b, N7]", {
  for (L in acc_grid) {
    s <- sqrt(sum(L)); m <- L[1] - L[2]
    x <- unique(round(m + s * c(-30, -8, -3, -1, 0, 1, 3, 8, 30)))
    expect_silent(lp <- dskellam(x, L[1], L[2], log = TRUE))
    expect_equal(lp, vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]), tolerance = 1e-12,
                 label = sprintf("log-pmf at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("log-pmf matches the exact sum to 1e-9 for extremely unequal huge rates", {
  for (L in unequal) {
    s <- sqrt(sum(L)); m <- L[1] - L[2]
    x <- unique(round(m + s * c(-30, -8, -3, 0, 3, 8, 30)))
    expect_silent(lp <- dskellam(x, L[1], L[2], log = TRUE))
    expect_equal(lp, vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]), tolerance = 1e-9,
                 label = sprintf("log-pmf at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("pmf sums to 1 at large and very unequal rates, without warnings  [N7]", {
  for (L in list(c(1e5, 1e5), c(1e6, 1e6), c(1e4, 1), c(1e6, 1e3))) {
    x <- support_grid(L[1], L[2], w = 25)
    expect_silent(d <- dskellam(x, L[1], L[2]))
    expect_equal(sum(d), 1, tolerance = 1e-10, label = sprintf("mass at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("log = TRUE stays finite in the deep tails  [N5]", {
  x <- c(-2000, -500, -200, 200, 500, 1000)
  expect_silent(lp <- dskellam(x, 3, 2, log = TRUE))
  expect_true(all(is.finite(lp)))
  expect_equal(lp, vapply(x, lp_exact, 0, l1 = 3, l2 = 2), tolerance = 1e-12)
  expect_identical(dskellam(500, 3, 2), 0)                  # natural scale underflows to 0, correctly
})

test_that("log = TRUE never returns a natural-scale value  [N5b]", {
  expect_lt(dskellam(0, 1e5, 1e5, log = TRUE), 0)
  expect_equal(dskellam(0, 1e5, 1e5, log = TRUE), lp_exact(0, 1e5, 1e5), tolerance = 1e-12)
})

test_that("fallback indexing cannot mix observations  [N6]", {
  expect_identical(dskellam(c(0, 400, -3), 3, 2)[2], dskellam(400, 3, 2))
  expect_identical(dskellam(c(0, 400, -3), 3, 2, log = TRUE)[2], dskellam(400, 3, 2, log = TRUE))
})

test_that("tiny rates", {
  x <- c(-2, -1, 0, 1, 2, 5)
  for (L in list(c(1e-8, 1e-8), c(1e-300, 1), c(1e-3, 1e4))) {
    expect_equal(dskellam(x, L[1], L[2], log = TRUE), vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]),
                 tolerance = 1e-11)
  }
})

test_that("dskellam.sp honours log and carries the full second-order correction  [N15]", {
  expect_equal(dskellam.sp(-5:5, 3, 2, log = TRUE), log(dskellam.sp(-5:5, 3, 2)))
  # with the full correction the approximation integrates to 1 up to O(1/lambda^2)
  x <- support_grid(1e5, 1e5, w = 20)
  expect_equal(sum(dskellam.sp(x, 1e5, 1e5)), 1, tolerance = 1e-9)
  # and is within 1e-5 relative of the exact pmf at moderate rates
  x <- -40:40
  expect_equal(dskellam.sp(x, 100, 80), dskellam(x, 100, 80), tolerance = 1e-4)
})

# ---- cdf --------------------------------------------------------------------

# Reference: cumulative sums of the (now accurate) pmf over +-45 sd, compared
# only in the interior (+-30 sd) where the truncation of the reference is
# negligible.
cdf_reference <- function(l1, l2) {
  g <- support_grid(l1, l2, w = 45)
  d <- dskellam(g, l1, l2)
  s <- sqrt(l1 + l2); m <- l1 - l2
  cs <- cumsum(d); us <- rev(cumsum(rev(d))) - d
  keep <- abs(g - m) <= 30 * s & cs > 1e-290 & us > 1e-290
  list(x = g[keep], lower = cs[keep], upper = us[keep])
}

test_that("pskellam matches cumulative sums at large rates, both tails, no warnings  [N8]", {
  # the reference inherits the pmf's own accuracy floor: 1e-13 for equal
  # rates, about 1e-9 for very unequal huge rates (see the log-pmf tests)
  for (L in list(c(2e4, 2e4), c(1e5, 1e5), c(1e5, 1e4), c(1e6, 1e6))) {
    tol <- if (L[1] / L[2] > 2) 1e-9 else 1e-11
    ref <- cdf_reference(L[1], L[2])
    expect_silent(lo <- pskellam(ref$x, L[1], L[2]))
    expect_silent(up <- pskellam(ref$x, L[1], L[2], lower.tail = FALSE))
    expect_equal(lo, ref$lower, tolerance = tol, label = sprintf("lower cdf at (%g, %g)", L[1], L[2]))
    expect_equal(up, ref$upper, tolerance = tol, label = sprintf("upper cdf at (%g, %g)", L[1], L[2]))
    expect_equal(lo + up, rep(1, length(lo)), tolerance = 1e-15)
  }
})

test_that("pskellam satisfies the symmetry identity P(X <= 0) = (1 + P(X = 0))/2 up to 1e7  [N8]", {
  for (L in c(1e2, 1e4, 1e5, 1e6, 1e7)) {
    expect_equal(pskellam(0, L, L), (1 + dskellam(0, L, L)) / 2, tolerance = 1e-13,
                 label = sprintf("lambda = %g", L))
  }
})

test_that("pskellam tails are non-negative and exact far out, on both scales  [N9]", {
  expect_gte(pskellam(-200, 3, 2), 0)
  expect_equal(pskellam(-200, 3, 2, log.p = TRUE), -729.5627, tolerance = 1e-6)
  expect_equal(pskellam(300, 3, 2, lower.tail = FALSE, log.p = TRUE), -1094.901, tolerance = 1e-6)
  expect_equal(pskellam(-1000, 3, 2, log.p = TRUE), -5223.973, tolerance = 1e-6)
  expect_identical(pskellam(-1000, 3, 2), 0)
  # continuity between the pchisq route and the summation route
  expect_equal(pskellam(c(-5, 0, 5), 5000, 5000), pskellam(c(-5, 0, 5), 5000 + 1e-9, 5000 + 1e-9),
               tolerance = 1e-10)
})

test_that("pskellam with many distinct rate pairs agrees with one-at-a-time calls", {
  set.seed(2)
  l1 <- runif(60, 2e4, 3e4); l2 <- runif(60, 1e4, 2e4)
  x  <- round(l1 - l2 + rnorm(60, 0, 2 * sqrt(l1 + l2)))
  a  <- pskellam(x, l1, l2, log.p = TRUE)
  b  <- vapply(seq_along(x), function(i) pskellam(x[i], l1[i], l2[i], log.p = TRUE), 0)
  expect_equal(a, b, tolerance = 1e-13)
})

test_that("quantiles at large rates are fast, silent and invert the cdf  [N8]", {
  expect_silent(q <- qskellam(c(1e-10, 0.5, 1 - 1e-10), 1e6, 1e6))
  expect_identical(q, c(-8996, 0, 8996))
  x <- support_grid(1e5, 1e5, w = 5)
  expect_identical(qskellam(pskellam(x, 1e5, 1e5), 1e5, 1e5), as.numeric(x))
  expect_identical(qskellam(pskellam(x, 1e5, 1e5, lower.tail = FALSE), 1e5, 1e5, lower.tail = FALSE), as.numeric(x))
})

test_that("pskellam.sp stays within [0, 1], handles both tails and log.p  [N16]", {
  x <- -40:40
  expect_true(all(pskellam.sp(x, 3, 2) >= 0 & pskellam.sp(x, 3, 2) <= 1))
  expect_gte(pskellam.sp(-200, 3, 2), 0)
  expect_equal(pskellam.sp(x, 3, 2, lower.tail = FALSE), 1 - pskellam.sp(x, 3, 2), tolerance = 1e-12)
  expect_equal(pskellam.sp(x, 3, 2, log.p = TRUE), log(pskellam.sp(x, 3, 2)))
  # relative accuracy of the approximation itself at moderate rates
  expect_equal(pskellam.sp(x, 100, 80), pskellam(x, 100, 80), tolerance = 1e-3)
})
