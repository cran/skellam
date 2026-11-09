## --------------------------------------------------------------------------
## pskellam.r -- Skellam cumulative distribution function
##
## Public functions:
##   pskellam()      P(X <= q) or P(X > q): non-central chi-square identity
##                   for lambda1 + lambda2 <= 1e4, log-space summation of the
##                   pmf above that and wherever pchisq() underflows
##
## Internal functions:
##   .skellam_logcdf_sum()  log P(X <= x) and log P(X > x) by summing the
##                          log-pmf from x toward the nearer tail
##   .logaddexp(), .log1mexp()  stable log-space arithmetic
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
pskellam <- function(q, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE){
 # CDF of Skellam distribution (difference of Poissons)
    lower.tail <- isTRUE(as.logical(lower.tail[1]))
    log.p <- isTRUE(as.logical(log.p[1]))
    nms <- names(q)
    nq <- length(q)
    a <- .skellam_recycle(q, lambda1, lambda2)
    x <- as.numeric(a[[1]]); lambda1 <- as.numeric(a[[2]]); lambda2 <- as.numeric(a[[3]])
    len <- length(x)
    ret <- rep(NaN, length.out = len)
    bad <- .skellam_bad_rate(lambda1) | .skellam_bad_rate(lambda2)
    if (any(bad)) warning("NaNs produced")
    miss <- !bad & (is.na(x) | is.na(lambda1) | is.na(lambda2))
    ret[miss] <- (x + lambda1 + lambda2)[miss]
    ok <- !bad & !miss
   # CDF is a step function, so convert to integer values without warning
   # (the small tolerance keeps 2.9999999999996 at 3, as ppois does)
    fin <- ok & is.finite(x)
    x[fin] <- floor(x[fin] + 1e-7 * pmax(1, abs(x[fin])))
   # q = +-Inf: the cdf is 0 or 1
    p0 <- if (log.p) -Inf else 0
    p1 <- if (log.p) 0 else 1
    ninf <- ok & x == -Inf
    pinf <- ok & x == Inf
    ret[ninf] <- if (lower.tail) p0 else p1
    ret[pinf] <- if (lower.tail) p1 else p0
    ok <- ok & is.finite(x)
   # a zero rate reduces to a (reflected) Poisson; handled exactly (N10)
   #   lambda2 = 0: X = Y1,  P(X <= x) = P(Y1 <= x)
   #   lambda1 = 0: X = -Y2, P(X <= x) = P(Y2 >= -x) = P(Y2 > -x - 1)
    z2 <- ok & lambda2 == 0
    z1 <- ok & lambda1 == 0 & !z2
    ret[z2] <- stats::ppois(x[z2], lambda1[z2], lower.tail = lower.tail, log.p = log.p)
    ret[z1] <- stats::ppois(-x[z1] - 1, lambda2[z1], lower.tail = !lower.tail, log.p = log.p)
   # both rates positive. Johnson's identity through pchisq() is accurate to
   # about 1e-12 while the non-centrality 2*lambda is moderate; above
   # lambda1 + lambda2 = 1e4 it slows down (cost grows with the non-centrality),
   # loses accuracy, and beyond ~5e5 stops converging and returns a wrong value
   # (N8). It also returns 0 in tails where the answer is ~1e-16 and, on the
   # natural scale, cannot represent anything below 1e-308. Those cases go to a
   # direct log-space summation of the pmf, which is exact to ~1e-13.
    nz  <- ok & !z1 & !z2
    chi <- nz & (lambda1 + lambda2 <= 1e4)
    neg <- chi & (x <  0)
    pos <- chi & (x >= 0)
    if (lower.tail) {
        ret[neg] <- stats::pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],log.p=log.p)
        ret[pos] <- stats::pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],lower.tail=FALSE,log.p=log.p)
    } else {
        ret[neg] <- stats::pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],lower.tail=FALSE,log.p=log.p)
        ret[pos] <- stats::pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],log.p=log.p)
    }
    tiny <- if (log.p) -700 else 1e-300
    sum_ <- nz & (!chi | !is.finite(ret) | ret < tiny)
    if (any(sum_)) {
        lc <- .skellam_logcdf_sum(x[sum_], lambda1[sum_], lambda2[sum_])
        v  <- if (lower.tail) lc$lower else lc$upper
        ret[sum_] <- if (log.p) v else exp(v)
    }
    if (length(nms) && nq == len) names(ret) <- nms
    ret
}


## log(exp(a) + exp(b)), elementwise, safe for -Inf.
.logaddexp <- function(a, b) {
  m <- pmax(a, b)
  d <- -abs(a - b)
  out <- m + log1p(exp(d))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

## log(1 - exp(a)) for a <= 0, accurate for a near 0 and near -Inf.
.log1mexp <- function(a) {
  ifelse(a > -log(2), log(-expm1(a)), log1p(-exp(a)))
}

## Cumulative log-sum-exp along a vector (forward), in blocks of 512: each
## block is exponentiated relative to its own maximum, so nothing overflows
## and deep-tail blocks keep their relative accuracy.
.cumlogsumexp <- function(lp) {
  n <- length(lp)
  out <- numeric(n)
  total <- -Inf
  B <- 512L
  for (start in seq.int(1L, n, by = B)) {
    i  <- start:min(start + B - 1L, n)
    mx <- max(lp[i])
    cs <- cumsum(exp(lp[i] - mx))
    out[i] <- mx + log(exp(total - mx) + cs)
    total <- out[i[length(i)]]
  }
  out
}

## log P(X <= x) and log P(X > x), both rates positive, x integer.
## The tail that holds the smaller mass is summed directly, term by term from
## x outward, so that deep tails keep full relative accuracy; the other tail is
## its complement. Terms more than W = 9 sd + 30 steps from x are below 1e-16
## of the total and are dropped.
## With a single (lambda1, lambda2) pair (the usual call) the pmf is evaluated
## once on a grid covering all queries and accumulated; with many distinct
## pairs the terms are evaluated in blocks of 128 per element until each
## element has converged.
.skellam_logcdf_sum <- function(x, l1, l2) {
  n  <- length(x)
  m  <- l1 - l2
  W  <- ceiling(9 * sqrt(l1 + l2) + 30)
  lowdir <- x < m                    # sum the lower tail directly, else the upper
  direct <- rep(NA_real_, n)
  one_pair <- all(l1 == l1[1]) && all(l2 == l2[1])
  if (one_pair) {
    if (any(lowdir)) {
      xs <- x[lowdir]
      g  <- (min(xs) - W[1]):max(xs)
      cl <- .cumlogsumexp(.skellam_logpmf(g, rep(l1[1], length(g)), rep(l2[1], length(g))))
      direct[lowdir] <- cl[xs - g[1] + 1]
    }
    if (any(!lowdir)) {
      xs <- x[!lowdir]
      g  <- (max(xs) + 1 + W[1]):(min(xs) + 1)          # descending, from the far tail inward
      cl <- .cumlogsumexp(.skellam_logpmf(g, rep(l1[1], length(g)), rep(l2[1], length(g))))
      direct[!lowdir] <- cl[g[1] - xs]                  # position of x + 1 in g
    }
  } else {
    B <- 128L
    S <- rep(-Inf, n)
    act <- seq_len(n)
    k0 <- 0L
    repeat {
      ks  <- k0 + 0:(B - 1L)
      d   <- ifelse(lowdir[act], -1, 1)
      pos <- outer(ks, d) + rep(x[act], each = B) + rep(d > 0, each = B)
      lp  <- matrix(.skellam_logpmf(as.vector(pos), rep(l1[act], each = B), rep(l2[act], each = B)), B)
      mx  <- apply(lp, 2, max)
      chunk <- mx + log(colSums(exp(lp - rep(mx, each = B))))
      S[act] <- .logaddexp(S[act], chunk)
      done <- (lp[B, ] < S[act] - 40) | (k0 + B > W[act])
      act <- act[!done]
      k0  <- k0 + B
      if (!length(act)) break
    }
    direct <- S
  }
  other <- .log1mexp(direct)
  list(lower = ifelse(lowdir, direct, other),
       upper = ifelse(lowdir, other, direct))
}
