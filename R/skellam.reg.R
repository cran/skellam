#' @export
skellam.reg <- function(y, x) {
  
  n <- length(y)
  x <- stats::model.matrix( ~., data.frame(x) )
  p <- dim(x)[2]
 
  skelreg <- function(pa) {
    b1 <- pa[1:p]   ;   b2 <- pa[ -c(1:p) ]
    a1 <- x %*% b1      ;     a2 <- x %*% b2
    lam1 <- exp(a1)     ;     lam2 <- exp(a2)
    a <- 2 * sqrt(lam1 * lam2)
    sum(lam1 + lam2)  + 0.5 * sum(y * (a1 - a2) ) - sum( log( besselI(a, y) ) )
  }
  
  options(warn = -1)
  mod <- stats::nlm(skelreg, stats::rnorm(2 * p), iterlim = 5000 ) 
  mod <- stats::nlm(skelreg, mod$estimate, iterlim = 5000 ) 
  mod <- stats::optim(mod$estimate, skelreg, hessian = TRUE, control = list(maxit = 5000) )
  b1 <- mod$par[1:p]    ;    b2 <- mod$par[ -c(1:p) ]
  s <- diag( solve(mod$hessian) ) 
  s1 <- sqrt(s[1:p])    ;    s2 <- sqrt(s[ -c(1:p) ])
  param1 <- cbind(b1, s1, b1 / s1, stats::pchisq( (b1 / s1)^2, 1, lower.tail = FALSE) )
  param2 <- cbind(b2, s2, b2 / s2, stats::pchisq( (b2 / s2)^2, 1, lower.tail = FALSE) )
  rownames(param1) <- rownames(param2) <- colnames(x)
  colnames(param1) <- colnames(param2) <- c("Estimate", "Std. Error", "Wald value", "p-value")

  list(loglik = -mod$value, param1 = param1, param2 = param2)
}  

