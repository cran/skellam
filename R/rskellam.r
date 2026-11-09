## --------------------------------------------------------------------------
## rskellam.r -- Skellam random variates
##
## Public functions:
##   rskellam()      draws as differences of two Poisson variates
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
rskellam <- function(n, lambda1, lambda2=lambda1){
 # Skellam random variables
    if (length(n)>1) n <- length(n)
    lambda1 <- rep(lambda1,length.out=n)
    lambda2 <- rep(lambda2,length.out=n)
    oops <- !(is.finite(lambda1)&(lambda1>=0)&is.finite(lambda2)&(lambda2>=0))
    if(any(oops)) warning("NAs produced")
   # integer result, like rpois (N17); NA for invalid rates, as rpois does.
   # The draws (and the random-number stream) are unchanged.
    ret <- rep(NA_integer_,length.out=n)
    n <- n-sum(oops)
    ret[!oops] <- stats::rpois(n,lambda1[!oops])-stats::rpois(n,lambda2[!oops])
    ret
}
