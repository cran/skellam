## --------------------------------------------------------------------------
## skellam.reg.R -- Skellam regression
##
## Public functions:
##   skellam.reg()   log-linear regression on both rates
##
## Original author: Michail Tsagris
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' Skellam Regression
#'
#' Fits a regression model assuming a Skellam distribution for the response variable.
#'
#' @param y A vector of integers (positive or negative)
#' @param x A matrix, vector or data.frame of covariates
#'
#' @details
#' The response is modelled as \eqn{y_i = Y_{1i} - Y_{2i}}{y = Y1 - Y2} with
#' \eqn{Y_{1i} \sim \mathrm{Poisson}(\lambda_{1i})}{Y1 ~ Poisson(lambda1)} and
#' \eqn{Y_{2i} \sim \mathrm{Poisson}(\lambda_{2i})}{Y2 ~ Poisson(lambda2)}, so that
#' \code{param1} describes the rate of the counts that are added and
#' \code{param2} the rate of the counts that are subtracted, as in
#' \code{dskellam(y, lambda1, lambda2)}. Both rates use a log link,
#' \eqn{\log\lambda_{ki} = x_i^\top\beta_k}{log(lambda_k) = x'beta_k}, with the same
#' covariates (an intercept is added).
#'
#' The log-likelihood is computed with \code{dskellam(log = TRUE)} and
#' maximized by \code{\link[stats]{optim}} (BFGS) with the analytic gradient,
#' starting deterministically from the moment estimates of an intercept-only
#' model, so repeated calls give identical results. Standard errors come from
#' the inverse of the observed information (\code{\link[stats]{optimHess}}).
#' Observations with a missing response or covariate are removed with a
#' warning. A warning is also given if the optimizer does not converge or the
#' information matrix is not invertible (standard errors are then \code{NA}).
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{loglik}{Maximized log-likelihood value}
#'   \item{param1}{Matrix for \eqn{\lambda_1} parameters:
#'     \itemize{
#'       \item Column 1: Estimated coefficients
#'       \item Column 2: Standard errors
#'       \item Column 3: Wald z statistics (coef/se)
#'       \item Column 4: p-values (Wald test)
#'     }
#'   }
#'   \item{param2}{Matrix for \eqn{\lambda_2} parameters (same structure as param1)}
#' }
#'
#' @references
#' \itemize{
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between
#'   two Poisson variates belonging to different populations. \emph{Journal of the
#'   Royal Statistical Society, Series A} \bold{109}(3), 296.
#'
#'   \item Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of
#'   the difference between two Poisson variates. \emph{Statistica Neerlandica}
#'   \bold{16}(1), 17-23.
#'
#'   \item Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football outcomes:
#'   using the Skellam's distribution for the goal difference. \emph{IMA Journal of
#'   Management Mathematics} \bold{20}(2), 133-145.
#' }
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(100)
#' y1 <- rpois(100, exp(1 + 1 * x))
#' y2 <- rpois(100, exp(-1 + 1 * x))
#' y <- y1 - y2
#' skellam.reg(y, x)   # param1 near (1, 1), param2 near (-1, 1)
#'
#' @keywords models regression
#' @author Michail Tsagris
#' @export
skellam.reg <- function(y, x) {

  X <- stats::model.matrix(~., stats::model.frame(~., data.frame(x), na.action = stats::na.pass))
  if (length(y) != nrow(X)) stop("'y' and 'x' must have the same number of observations")
  if (!is.numeric(y)) stop("'y' must be a numeric vector of integers")
  # remove incomplete observations, as lm() and glm() do, but say so
  ok <- !is.na(y) & stats::complete.cases(X)
  if (!all(ok)) {
    warning(sprintf("%d observation(s) with missing values removed", sum(!ok)), call. = FALSE)
    y <- y[ok]; X <- X[ok, , drop = FALSE]
  }
  if (any(!is.finite(y)) || any(.skellam_nonint(y)))
    stop("'y' must contain finite integers (differences of counts)")
  y <- round(y)
  n <- length(y)
  p <- ncol(X)
  if (n <= 2 * p) stop("too few observations for the number of parameters")

  # minus log-likelihood and its gradient in (beta1, beta2), log links.
  # d log f / d log lambda1 = (z/2) I'_nu(z)/I_nu(z) - lambda1 + y/2,
  # d log f / d log lambda2 = (z/2) I'_nu(z)/I_nu(z) - lambda2 - y/2,
  # with nu = |y|, z = 2 sqrt(lambda1 lambda2), I'_nu/I_nu = I_{nu+1}/I_nu + nu/z.
  rates <- function(pa) {
    list(l1 = exp(drop(X %*% pa[1:p])), l2 = exp(drop(X %*% pa[-(1:p)])))
  }
  # the line search can try steps where a rate overflows to Inf or
  # underflows to 0; report those points as "infinitely bad" instead of
  # evaluating the density there
  usable <- function(r) all(is.finite(r$l1) & is.finite(r$l2) & r$l1 > 0 & r$l2 > 0)
  nll <- function(pa) {
    r <- rates(pa)
    if (!usable(r)) return(1e300)
    v <- -sum(dskellam(y, r$l1, r$l2, log = TRUE))
    if (is.finite(v)) v else 1e300
  }
  grad <- function(pa) {
    r  <- rates(pa)
    if (!usable(r)) return(rep(0, 2 * p))
    nu <- abs(y)
    z  <- 2 * sqrt(r$l1) * sqrt(r$l2)
    ratio <- exp(.log_besselI_scaled(z, nu + 1) - .log_besselI_scaled(z, nu))
    common <- 0.5 * z * ratio + 0.5 * nu
    g1 <- common - r$l1 + 0.5 * y
    g2 <- common - r$l2 - 0.5 * y
    -c(crossprod(X, g1), crossprod(X, g2))
  }

  # deterministic start: intercept-only moment estimates, slopes 0
  m  <- mean(y)
  l2 <- max((stats::var(y) - m) / 2, 0.5)
  l1 <- max(l2 + m, 0.5)
  start <- c(log(l1), rep(0, p - 1), log(l2), rep(0, p - 1))

  mod <- stats::optim(start, nll, grad, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-12))
  if (mod$convergence != 0)
    warning("optimizer did not converge (code ", mod$convergence, "); estimates may be unreliable",
            call. = FALSE)

  b1 <- mod$par[1:p]    ;    b2 <- mod$par[ -c(1:p) ]
  H <- stats::optimHess(mod$par, nll, grad)
  V <- tryCatch(solve(H), error = function(e) NULL)
  s <- if (is.null(V)) rep(NA_real_, 2 * p) else diag(V)
  if (is.null(V) || any(!is.finite(s)) || any(s <= 0)) {
    warning("the information matrix is not positive definite; standard errors are not available",
            call. = FALSE)
    s[!is.finite(s) | s <= 0] <- NA_real_
  }
  s1 <- sqrt(s[1:p])    ;    s2 <- sqrt(s[ -c(1:p) ])
  param1 <- cbind(b1, s1, b1 / s1, stats::pchisq( (b1 / s1)^2, 1, lower.tail = FALSE) )
  param2 <- cbind(b2, s2, b2 / s2, stats::pchisq( (b2 / s2)^2, 1, lower.tail = FALSE) )
  rownames(param1) <- rownames(param2) <- colnames(X)
  colnames(param1) <- colnames(param2) <- c("Estimate", "Std. Error", "Wald value", "p-value")

  list(loglik = -mod$value, param1 = param1, param2 = param2)
}
