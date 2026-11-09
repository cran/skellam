# skellam

[![CRAN status](https://www.r-pkg.org/badges/version/skellam)](https://cran.r-project.org/package=skellam)
[![R-CMD-check](https://github.com/monty-se/skellam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/monty-se/skellam/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>= 2)](https://img.shields.io/badge/License-GPL%20(%3E%3D%202)-blue.svg)](https://www.r-project.org/Licenses/)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/skellam)](https://cran.r-project.org/package=skellam)

The Skellam distribution is the distribution of the difference of two independent
Poisson counts: home goals minus away goals, up-ticks minus down-ticks, spikes in a
stimulus window minus spikes in a baseline window. **skellam** provides its density,
distribution function, quantile function and random generation, saddlepoint
approximations, maximum likelihood estimation, and a regression model in which both
rates depend on covariates.

The distribution functions are accurate to a relative error below 1e-12 for rates up
to at least 1e7, and remain finite on the log scale far into the tails, so they can be
used safely inside likelihoods. They follow the conventions of R's own `dpois`,
`ppois`, `qpois` and `rpois`.

Documentation, the changelog and the validation article are on the package website:
<https://monty-se.github.io/skellam/>.

## Installation

```r
install.packages("skellam")              # CRAN
# install.packages("remotes")
remotes::install_github("monty-se/skellam")   # development version
```

## Usage

```r
library(skellam)

# X = Y1 - Y2 with Y1 ~ Poisson(12), Y2 ~ Poisson(8)
dskellam(-2:2, 12, 8)                    # P(X = x)
pskellam(0, 12, 8)                       # P(X <= 0)
qskellam(c(0.025, 0.975), 12, 8)         # central 95% interval
set.seed(1)
rskellam(5, 12, 8)                       # random draws (integers)

# Log scale stays finite where the probability itself underflows
dskellam(500, 3, 2, log = TRUE)
pskellam(-1000, 3, 2, log.p = TRUE)

# Large rates
pskellam(0, 1e6, 1e6)
```

### Estimation

```r
set.seed(1)
x <- rpois(2000, 12) - rpois(2000, 8)
skellam.mle(x)$param                     # close to (12, 8)
```

### Regression

For `y = Y1 - Y2`, both rates follow a log-linear model in the covariates. `param1`
describes the rate of the counts that are added and `param2` the rate of the counts
that are subtracted.

```r
set.seed(123)
n <- 1000
covariate <- rnorm(n)
y <- rpois(n, exp(1 + 0.5 * covariate)) - rpois(n, exp(-0.5 + 0.3 * covariate))
fit <- skellam.reg(y, covariate)
fit$param1                               # intercept ~ 1,    slope ~ 0.5
fit$param2                               # intercept ~ -0.5, slope ~ 0.3
```

> **Note for users of version 0.2.4 or earlier:** in those versions `skellam.reg()`
> returned the two parameter blocks the wrong way round (`param1` estimated the rate
> of the subtracted counts). Results obtained with them should be re-run. See
> [NEWS](NEWS.md) for all changes.

## Background

If $Y_1 \sim \mathrm{Poisson}(\lambda_1)$ and $Y_2 \sim \mathrm{Poisson}(\lambda_2)$
are independent, $X = Y_1 - Y_2$ has probability mass function

$$
P(X = x) = e^{-(\lambda_1+\lambda_2)} \left(\frac{\lambda_1}{\lambda_2}\right)^{x/2}
I_{|x|}\left(2\sqrt{\lambda_1\lambda_2}\right), \qquad x \in \mathbb{Z},
$$

where $I_\nu$ is the modified Bessel function of the first kind, with mean
$\lambda_1 - \lambda_2$ and variance $\lambda_1 + \lambda_2$. The help page
`?skellam` describes the numerical methods; the vignette
`vignette("validate_skellam", package = "skellam")` checks every function against
the definition.

## References

- Skellam, J. G. (1946). The frequency distribution of the difference between two
  Poisson variates belonging to different populations. *Journal of the Royal
  Statistical Society, Series A*, 109(3), 296.
- Johnson, N. L. (1959). On an extension of the connection between Poisson and
  χ² distributions. *Biometrika*, 46, 352–362.
- Karlis, D., & Ntzoufras, I. (2009). Bayesian modelling of football outcomes: using
  the Skellam's distribution for the goal difference. *IMA Journal of Management
  Mathematics*, 20(2), 133–145.

## Authors

Jerry W. Lewis, Patrick E. Brown and Michail Tsagris wrote the original package.
Montasser Ghachem (maintainer) and Oguz Ersan made the corrections and extensions
from version 0.2.5 on. Please report bugs at
<https://github.com/monty-se/skellam/issues>, and cite the package with
`citation("skellam")`.

## License

GPL (>= 2)
