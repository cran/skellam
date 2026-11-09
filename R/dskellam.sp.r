## --------------------------------------------------------------------------
## dskellam.sp.r -- Saddlepoint approximation to the Skellam pmf
##
## Public functions:
##   dskellam.sp()   second-order saddlepoint density
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
dskellam.sp <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # second-order saddlepoint approximation to the Skellam PMF (Butler 2007, eq. 1.8)
    log <- isTRUE(as.logical(log[1]))
    s  <- log(0.5*(x+sqrt(x^2+4*lambda1*lambda2))/lambda1)  # the saddlepoint, K'(s) = x
    K  <- lambda1*(exp(s)-1)+lambda2*(exp(-s)-1)             # CGF(s)
    K2 <- lambda1*exp(s)+lambda2*exp(-s)                     # K''(s)
    K3 <- lambda1*exp(s)-lambda2*exp(-s)                     # K'''(s);  K''''(s) equals K''(s)
   # correction  kappa4/8 - 5 kappa3^2/24,  kappa3 = K3/K2^1.5,  kappa4 = 1/K2
    corr <- (1 - 5/3*(K3/K2)^2)/(8*K2)
    ret <- K - x*s - 0.5*log(2*pi*K2) + log1p(corr)
    if (log) ret else exp(ret)
}
