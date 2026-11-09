## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4,
                      out.width = "100%")
library(skellam)
set.seed(20261108)

## ----pmf-exact----------------------------------------------------------------
pmf_by_sum <- function(x, l1, l2) {
  k <- max(0, -x):(abs(x) + 200 + 20 * (l1 + l2))
  sum(dpois(k + x, l1) * dpois(k, l2))
}
x <- -10:15
for (L in list(c(1.5, 0.5), c(12, 8))) {
  ref <- vapply(x, pmf_by_sum, 0, l1 = L[1], l2 = L[2])
  cat(sprintf("lambda = (%g, %g): largest relative difference %.1e\n",
              L[1], L[2], max(abs(dskellam(x, L[1], L[2]) / ref - 1))))
}

## ----pmf-mass-----------------------------------------------------------------
for (L in list(c(1.5, 0.5), c(12, 8), c(1e4, 1e4), c(1e6, 1e6))) {
  s <- sqrt(sum(L)); m <- L[1] - L[2]
  x <- seq(floor(m - 30 * s), ceiling(m + 30 * s))
  cat(sprintf("lambda = (%g, %g): 1 - sum = %.1e\n", L[1], L[2],
              1 - sum(dskellam(x, L[1], L[2]))))
}

## ----pmf-poisson--------------------------------------------------------------
all.equal(dskellam(0:10, 5, 0), dpois(0:10, 5))
all.equal(dskellam(-(0:10), 0, 5), dpois(0:10, 5))

## ----cdf----------------------------------------------------------------------
l1 <- 12; l2 <- 8
x <- -40:50
max(abs(pskellam(x, l1, l2) - cumsum(dskellam(x, l1, l2))))
max(abs(pskellam(x, l1, l2) + pskellam(x, l1, l2, lower.tail = FALSE) - 1))

## ----quantile-----------------------------------------------------------------
x <- -15:25
all(qskellam(pskellam(x, l1, l2), l1, l2) == x)
all(qskellam(pskellam(x, l1, l2, lower.tail = FALSE), l1, l2, lower.tail = FALSE) == x)
all(qskellam(pskellam(x, l1, l2, log.p = TRUE), l1, l2, log.p = TRUE) == x)
qskellam(c(0, 0.025, 0.5, 0.975, 1), l1, l2)

## ----random-------------------------------------------------------------------
N <- 20000
lambda1 <- 1.5; lambda2 <- 0.5
z <- rskellam(N, lambda1, lambda2)
c(mean = mean(z), expected = lambda1 - lambda2)
c(variance = var(z), expected = lambda1 + lambda2)

lo <- qskellam(0.001, lambda1, lambda2); hi <- qskellam(0.999, lambda1, lambda2)
cells <- lo:hi
prob <- dskellam(cells, lambda1, lambda2)
prob[1] <- pskellam(lo, lambda1, lambda2)                                 # all x <= lo
prob[length(prob)] <- pskellam(hi - 1, lambda1, lambda2, lower.tail = FALSE)  # all x >= hi
obs <- table(factor(pmin(pmax(z, lo), hi), levels = cells))
chisq.test(obs, p = prob)

## ----random-plot, fig.cap = "Observed frequencies of 20,000 draws (bars) and the expected frequencies from dskellam (points)."----
plot(cells, as.vector(obs), type = "h", lwd = 6, col = "grey70",
     xlab = "x", ylab = "Frequency", main = "")
points(cells, N * prob, pch = 19, col = "#0072B2")
legend("topright", c("rskellam", "N * dskellam"), lwd = c(6, NA), pch = c(NA, 19),
       col = c("grey70", "#0072B2"), bty = "n")

## ----mle----------------------------------------------------------------------
x <- rpois(2000, 12) - rpois(2000, 8)
fit <- skellam.mle(x)
fit$param

## ----mle-boundary-------------------------------------------------------------
skellam.mle(rep(3, 50))$param
skellam.mle(c(rep(2, 45), rep(3, 5)))$param

## ----reg----------------------------------------------------------------------
n <- 1000
covariate <- rnorm(n)
y <- rpois(n, exp(1 + 0.5 * covariate)) - rpois(n, exp(-0.5 + 0.3 * covariate))
reg <- skellam.reg(y, covariate)
round(reg$param1, 3)   # true values 1 and 0.5
round(reg$param2, 3)   # true values -0.5 and 0.3

