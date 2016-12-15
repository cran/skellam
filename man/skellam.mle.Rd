\name{skellam.mle}
\alias{skellam.mle}
\title{ MLE of the Skellam distribution}
\description{
MLE of the Skellam distribution.
}
\usage{
skellam.mle(x)
}
\arguments{
\item{x}{
A vector of integers, positive or negative.
}
}
\details{
Instead of having to maximise the log-likelihood with respect 
to the two parameters, \eqn{\lambda_1} and \eqn{\lambda_2}, we maximise with respect to 
\eqn{\lambda_2} and then \eqn{\lambda_1 = \lambda_2 + \bar{x}}. This makes it faster. 
The command "nlm" is used to optimise the log-likelihood as it proved to be faster than the "optimise". 
}

\value{
A list including:
\item{iters}{
The number of iterations required by "nlm".
}  
\item{loglik}{
The maximised log-likelihood value.
}
\item{param}{
The estimated parameters, \eqn{\hat{\lambda}_1} and \eqn{\hat{\lambda}_2}. 
}  
}

\references{
  Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
   Cambridge University Press, Cambridge & New York, p.17.

  Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{chi^2} distributions.
   \emph{Biometrika} \bold{46}, 352--362.

  Johnson, N. L.; Kotz, S.; Kemp, A. W. (1993) \emph{Univariate Discrete Distributions}, 2nd ed.,
   John Wiley and Sons, New York, pp.190-192.

  Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates
   belonging to different populations. \emph{Journal of the Royal Statistical Society, series A} \bold{109}/3, 26.

  Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of the difference between two Poisson variates.
   \emph{Statistica Neerlandica} \bold{16}/1, 17-23.

   Abdulhamid, A. A.; Maha, A. O. (2010) On The Poisson Difference Distribution Inference and Applications.
   \emph{BULLETIN of the Malaysian Mathematical Sciences Society}, \bold{33}/1, 17--45.

  Wikipedia. \emph{Skellam distribution} \url{http://en.wikipedia.org/wiki/Skellam_distribution}
}
\examples{
require('skellam')

x1 <- rpois(1000, 10)
x2 <- rpois(1000, 6)
x <- x1 - x2
skellam.mle(x)

x1 <- rpois(10000, 10)
x2 <- rpois(10000, 6)
x <- x1 - x2
skellam.mle(x)
}
\keyword{ Maximum likelihood estimation }
\author{ Michail Tsagris }
