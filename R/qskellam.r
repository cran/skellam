## --------------------------------------------------------------------------
## qskellam.r -- Skellam quantile function
##
## Public functions:
##   qskellam()      smallest integer x with P(X <= x) >= p (Cornish-Fisher start + search)
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
qskellam <- function(p, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE){
 # inverse CDF of Skellam distriubition (difference of Poissons)
    lower.tail <- isTRUE(as.logical(lower.tail[1]))
    log.p <- isTRUE(as.logical(log.p[1]))
    nms <- names(p)
    np <- length(p)
   # make all args the same length; warn when lengths are not multiples (N13)
    a <- .skellam_recycle(p, lambda1, lambda2)
    p <- as.numeric(a[[1]]); lambda1 <- as.numeric(a[[2]]); lambda2 <- as.numeric(a[[3]])
    len <- length(p)
    ret <- rep(NaN, length.out = len)
   # validate (N12): missing -> NA; invalid rate or probability -> NaN + warning
    bad <- .skellam_bad_rate(lambda1) | .skellam_bad_rate(lambda2)
    miss <- !bad & (is.na(p) | is.na(lambda1) | is.na(lambda2))
    ret[miss] <- (p + lambda1 + lambda2)[miss]
    pp <- if (log.p) exp(p) else p                  # probability on the natural scale
    badp <- !bad & !miss & (pp < 0 | pp > 1)
    if (any(bad | badp)) warning("NaNs produced")
    ok <- !bad & !miss & !badp
   # ends of the support (N2): X >= 0 when lambda2 = 0, X <= 0 when lambda1 = 0.
   # "towards the left end": p = 0 on the lower tail, p = 1 on the upper tail.
   # (match qpois in assuming that p within 2.25e-16 of 1 is actually 1)
    is0 <- pp == 0
    is1 <- pp + 1.01*.Machine$double.eps >= 1
    lend <- ok & (if (lower.tail) is0 else is1)
    rend <- ok & (if (lower.tail) is1 else is0)
    ret[lend] <- ifelse(lambda2[lend] == 0, 0, -Inf)
    ret[rend] <- ifelse(lambda1[rend] == 0, 0,  Inf)
    ok <- ok & !lend & !rend
   # lambda2 = 0: X = Y1 is Poisson, and qpois is exact (and faster).
   # lambda1 = 0 goes through the general search below: reflecting qpois
   # turned P(X <= x) into P(Y >= y) and was off by one at the atoms (N4).
    z2 <- ok & lambda2 == 0
    ret[z2] <- stats::qpois(p[z2], lambda1[z2], lower.tail = lower.tail, log.p = log.p)
    nz <- ok & !z2
    p <- p[nz]
    lambda1 <- lambda1[nz]
    lambda2 <- lambda2[nz]
   # Cornish-Fisher approximations
    z <- stats::qnorm(p,lower.tail=lower.tail,log.p=log.p)
    mu <- lambda1-lambda2
    vr <- lambda1+lambda2
    sg <- sqrt(vr)
    c0 <- mu+z*sg
    c1 <- (z^2-1)*mu/vr/6
    c2 <- -( c1*mu - 2*lambda1*lambda2*(z^2-3)/vr )*z/12/vr/sg
   # test and linear search (slow if p extreme or lambda1+lambda2 small)
    q0 <- round(c0+c1+c2)
    p0 <- pskellam(q0,lambda1,lambda2,lower.tail=lower.tail,log.p=log.p)
   # match qpois in assuming that a value within 1.4e-14 of p actually equals p.
   # The nudge must make the stopping condition easier to meet: shrink p for
   # the lower tail (F(x) >= p), enlarge it for the upper tail (P(X > x) <= p),
   # and on the log scale add log(1 -/+ eps) rather than multiply (N3).
    fz <- 64*.Machine$double.eps
    if (log.p) {
        p <- p + if (lower.tail) log1p(-fz) else log1p(fz)
    } else {
        p <- p * if (lower.tail) (1 - fz) else (1 + fz)
    }
    if (lower.tail){    # smallest x such that F(x) >= p (consider floor for greater efficiency?)
        up <- up1 <- p0 < p
        while(any(up1)){
            q0[up1] <- q0[up1]+1
            up1[up1] <- pskellam(q0[up1],lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) < p[up1]
        }
        down1 <- (!up)&(p0>p)   # oops: \qskellam(p,lambda,0) == qpois(p,lambda)+1 usually
        down1[down1] <- (pskellam(q0[down1]-1,lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1])
        while(any(down1)){
            q0[down1] <- q0[down1]-1
            down1[down1] <- (pskellam(q0[down1]-1,lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1])
        }
    } else {        # largest x such that F(x,lower.tail=FALSE) <= p (consider ceiling for greater efficiency?)
        down <- down1 <- p0 > p
        while(any(down1)){
            q0[down1] <- q0[down1]+1
            down1[down1] <- pskellam(q0[down1],lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1]
        }
        up1 <- (!down)&(p0<p)
        up1[up1] <- !(pskellam(q0[up1]-1,lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) > p[up1])
        while(any(up1)){
            q0[up1] <- q0[up1]-1
            up1[up1] <- !(pskellam(q0[up1]-1,lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) > p[up1])
        }
    }
    ret[nz] <- q0
    if (length(nms) && np == len) names(ret) <- nms
    ret
}
