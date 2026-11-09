## --------------------------------------------------------------------------
## bessel.R -- log of the scaled modified Bessel function, for all orders
##
## Internal functions (not exported):
##   .log_besselI_scaled()  log(I_nu(z)) - z for integer nu >= 0 and z > 0,
##                          finite for every input, accurate to ~1e-13
##   .skellam_logpmf()      log P(X = x) for valid, equal-length inputs
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

## Why this file exists
##
## The Skellam pmf is  P(X = x) = e^{-l1-l2} (l1/l2)^{x/2} I_|x|(2 sqrt(l1 l2)).
## R's besselI(z, nu, expon.scaled = TRUE) is accurate in the middle of the
## distribution but fails in two places that matter for likelihoods:
##   * z > 1e5 (rates above about 5e4): it returns 0 without any message;
##   * large order relative to z (deep tails): the scaled value underflows to
##     0, and R warns "precision lost in result".
## log(0) = -Inf then poisons any log-likelihood. The function below computes
## log(I_nu(z)) - z directly in the three regimes where besselI cannot:
##
##   nu >= 250              : the uniform asymptotic (Debye) expansion,
##                            DLMF 10.41.3 with coefficients u_1..u_5;
##   nu < 250, z > 1000     : the Debye values at nu = 250 and 251, then the
##                            three-term recurrence run DOWNWARD in nu, which is
##                            the numerically stable direction for I_nu;
##   nu < 250, z <= 1000    : besselI, and the power series of I_nu (summed on
##                            the log scale) where besselI underflows.
##
## The Debye expansion with six terms has relative error about
## |u_6(p)| / nu^6, below 1e-15 for nu >= 250 at every t.

## Debye coefficients u_k(p), DLMF 10.41.10 (p = 1 / sqrt(1 + t^2)).
.debye_u <- function(p) {
  p2 <- p * p
  list(
    u1 = p * (3 - 5 * p2) / 24,
    u2 = p2 * (81 - 462 * p2 + 385 * p2^2) / 1152,
    u3 = p * p2 * (30375 - 369603 * p2 + 765765 * p2^2 - 425425 * p2^3) / 414720,
    u4 = p2^2 * (4465125 - 94121676 * p2 + 349922430 * p2^2 - 446185740 * p2^3 +
                   185910725 * p2^4) / 39813120,
    u5 = p * p2^2 * (1519035525 - 49286948607 * p2 + 284499769554 * p2^2 -
                       614135872350 * p2^3 + 566098157625 * p2^4 -
                       188699385875 * p2^5) / 6688604160
  )
}

## log(I_nu(nu t)) - nu t  by the Debye expansion, for nu > 0.
## eta(t) - t = 1/(t + sqrt(1+t^2)) - asinh(1/t)   (both terms evaluated stably)
.log_besselI_debye <- function(z, nu) {
  t  <- z / nu
  s  <- sqrt(1 + t * t)
  g  <- 1 / (t + s) - asinh(1 / t)
  p  <- 1 / s
  u  <- .debye_u(p)
  nu * g - 0.5 * log(2 * pi * nu) - 0.5 * log(s) +
    log1p((((u$u5 / nu + u$u4) / nu + u$u3) / nu + u$u2) / nu^2 + u$u1 / nu)
}

## log(I_nu(z)) - z  by the power series, for z^2/4 not large compared to nu:
## I_nu(z) = (z/2)^nu / nu! * sum_k (z^2/4)^k / (k! (nu+1)_k)
.log_besselI_series <- function(z, nu) {
  q    <- z * z / 4
  term <- rep(1, length(z))
  sum  <- term
  for (k in 1:200) {
    term <- term * q / (k * (nu + k))
    sum  <- sum + term
    if (all(term < 1e-17 * sum)) break
  }
  nu * log(z / 2) - lgamma(nu + 1) - z + log(sum)
}

## Large z, nu < 250: start the recurrence at nu = 251, 250 (Debye) and step
## down, I_{k-1} = I_{k+1} + (2k/z) I_k. All terms are positive, so there is
## no cancellation, and the scaled values stay well inside double range.
.log_besselI_recur <- function(z, nu) {
  top  <- 250L
  ikp1 <- exp(.log_besselI_debye(z, top + 1))
  ik   <- exp(.log_besselI_debye(z, top))
  out  <- rep(NA_real_, length(z))
  for (k in top:1) {
    ikm1 <- ikp1 + (2 * k / z) * ik
    hit  <- nu == k - 1
    if (any(hit)) out[hit] <- ikm1[hit]
    ikp1 <- ik
    ik   <- ikm1
  }
  log(out)
}

## log(I_nu(z)) - z  for z > 0 and integer nu >= 0, vectorised, always finite.
##   nu >= 250            : Debye, whatever z (relative error below 1e-15)
##   nu <  250, z <= 1000 : besselI; the power series if besselI underflowed
##   nu <  250, z >  1000 : downward recurrence from the Debye values at 250
## besselI costs O(max(nu, z)) per element, so it is only used where that is
## cheap; everything else is O(1) per element (plus 250 vectorised steps).
.log_besselI_scaled <- function(z, nu) {
  out <- rep(NA_real_, length(z))
  deb <- nu >= 250
  rec <- !deb & z > 1000
  bes <- !deb & !rec
  if (any(deb)) out[deb] <- .log_besselI_debye(z[deb], nu[deb])
  if (any(rec)) out[rec] <- .log_besselI_recur(z[rec], nu[rec])
  if (any(bes)) {
    v  <- suppressWarnings(besselI(z[bes], nu[bes], expon.scaled = TRUE))
    ok <- is.finite(v) & v > 1e-290
    idx <- which(bes)
    out[idx[ok]] <- log(v[ok])
    if (any(!ok)) out[idx[!ok]] <- .log_besselI_series(z[idx[!ok]], nu[idx[!ok]])
  }
  out
}

## log P(X = x) for equal-length vectors with l1 > 0, l2 > 0, integer x.
## -l1 - l2 + z = -(sqrt(l1) - sqrt(l2))^2 avoids cancellation when l1 ~ l2.
.skellam_logpmf <- function(x, l1, l2) {
  z <- 2 * sqrt(l1) * sqrt(l2)
  .log_besselI_scaled(z, abs(x)) - (sqrt(l1) - sqrt(l2))^2 + 0.5 * x * (log(l1) - log(l2))
}
