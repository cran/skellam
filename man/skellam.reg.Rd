\name{skellam.reg}
\alias{skellam.reg}
\title{ Regression assuming a Skellam distribution }
\description{
Regression assuming a Skellam distribution.
}
\usage{
skellam.reg(y, x)
}
\arguments{
\item{y}{
A vector of integers, positive or negative.
}
\item{x}{
A matrix, a vector or a data.frame with the covariates.
}
}
\details{
We use the exponential link function to ensure that the both \eqn{\lambda_s} are positive. 
The command \code{\link[stats]{nlm}} does the main job.
}

\value{
A list including:
\item{loglik}{
The maximised log-likelihood value.
}  
\item{param1}{
The estimated regression coefficients of \eqn{\lambda_1}. This is matrix, with the first 
column being the estimated regression coefficients. The second column is their relevant standard error.
The third column is the t value (coef/se(coef)) and the final column is the p-value of the Wald test.
}
\item{param2}{
The estimated regression coefficients of \eqn{\lambda_2}. This is matrix, with the first 
column being the estimated regression coefficients. The second column is their relevant standard error.
The third column is the t value (coef/se(coef)) and the final column is the p-value of the Wald test.
}  
}

\references{
  Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates
   belonging to different populations. \emph{Journal of the Royal Statistical Society, series A} \bold{109}/3, 26.

  Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of the difference between two Poisson variates.
   \emph{Statistica Neerlandica} \bold{16}/1, 17-23.

   \emph{Karlis and Ntzoufras IMA 2009 presentation} \url{http://stat-athens.aueb.gr/~jbn/papers/files/20_Karlis_Ntzoufras_2009_IMA_presentation_handouts_v01.pdf}
}
\examples{
require('skellam')

set.seed(0)

x <- rnorm(1000)
y1 <- rpois(1000, exp(1 + 1 * x) )
y2 <- rpois(1000 , exp(-1 + 1 * x) )
y <- y2 - y1
skellam.reg(y, x)
}
\keyword{ Skellam regression }
\author{ Michail Tsagris }
