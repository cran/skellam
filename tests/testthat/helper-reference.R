# Exact reference for the Skellam log-pmf:
#   P(X = x) = sum_k P(Y1 = k + x) P(Y2 = k),  Y1 ~ Pois(l1), Y2 ~ Pois(l2),
# summed on the log scale (log-sum-exp) so nothing underflows. Slow but
# accurate to about 1e-13; the yardstick for the accuracy tests.
lp_exact <- function(x, l1, l2) {
  if (l1 == 0) return(dpois(-x, l2, log = TRUE))
  if (l2 == 0) return(dpois(x, l1, log = TRUE))
  k <- max(0, -x):ceiling(max(l1, l2) + 60 * sqrt(max(l1, l2)) + abs(x) + 100)
  m <- dpois(k + x, l1, log = TRUE) + dpois(k, l2, log = TRUE)
  mm <- max(m)
  if (!is.finite(mm)) return(-Inf)
  mm + log(sum(exp(m - mm)))
}

# Integer grid covering +-w standard deviations around the mean.
support_grid <- function(l1, l2, w = 20) {
  s <- sqrt(l1 + l2); m <- l1 - l2
  seq(floor(m - w * s), ceiling(m + w * s))
}
