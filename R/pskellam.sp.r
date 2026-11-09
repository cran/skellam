## --------------------------------------------------------------------------
## pskellam.sp.r -- Saddlepoint approximation to the Skellam cdf
##
## Public functions:
##   pskellam.sp()   Lugannani-Rice cdf with Butler's second continuity correction
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
pskellam.sp <- function(q, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE) {
 # Luganni-Rice saddlepoint CDF with Butler's 2nd continuity correction
    lower.tail <- isTRUE(as.logical(lower.tail[1]))
    log.p <- isTRUE(as.logical(log.p[1]))
    if (lower.tail) {
        xm <- -floor(q)-0.5                                     # continuity corrected x
       # distribution specific variables
        s <- log(0.5*(xm+sqrt(xm^2+4*lambda2*lambda1))/lambda2) # the saddlepoint
        K <- lambda2*(exp(s)-1)+lambda1*(exp(-s)-1)             # CGF(s)
        K2 <- lambda2*exp(s)+lambda1*exp(-s)                    # CGF''(s)
       # code depending on distribution only through previous variables
        u2 <- 2*sinh(0.5*s)*sqrt(K2)
        w2 <- sign(s)*sqrt(2*(s*xm-K))
        ret <- stats::pnorm(-w2)-stats::dnorm(w2)*(1/w2-1/u2)
       # avoid numeric problems near the removable discontinuity
        xe <- (xm+(lambda1-lambda2))/sqrt(lambda1+lambda2)
        g1 <- (lambda1-lambda2)/(lambda1+lambda2)^1.5
        ew <- abs(xe) < 1e-4
        ret[ew] <- (stats::pnorm(-xe)+stats::dnorm(xe)*g1/6*(1-xe^2))[ew]
       # the approximation is not guaranteed to stay inside [0, 1] far in the tails
        ret <- pmin(pmax(ret, 0), 1)
    } else ret <- pskellam.sp(-q-1,lambda2,lambda1)
    if (log.p) log(ret) else ret
}
