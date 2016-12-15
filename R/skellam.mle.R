#' @export
skellam.mle <- function(x) {
  
  n <- length(x)
  sx <- sum(x)
  sx2 <- sx / 2
  mx <- sx / n
  theta <- stats::var(x)/2  - mx/2
  
  skel <- function(theta) {
    a <- 2 * sqrt(theta^2 + theta * mx)
    n * 2 * theta + sx - sx2 * log(1 + mx / theta) -
    sum( log( besselI(a, x, expon.scaled = TRUE) ) ) - n * a
  }
  
  options(warn = -1)
  mod <- stats::nlm(skel, theta) 
  param <- c( mod$estimate + mx, mod$estimate )
  names(param) <- c("mu 1", "mu 2")
  list(iters = mod$iterations, loglik = -mod$minimum, param = param)

}  






