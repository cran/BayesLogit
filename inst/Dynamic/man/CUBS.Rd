% package: BayesLogit
\name{CUBS}

\alias{CUBS.R}
\alias{CUBS.C}
%\alias{FFBS.R}
%\alias{FFBS.C}

\title{Conjugate Update Backward Sample}

\description{

  Draw an approximation of the states in a dynamic generalized linear model
  using conjugate updating and a backwards sampling.

  THIS IS A BETA VERSION.  AN ERROR CAN CRASH R.  BE CAREFUL.

}

\details{

  In a generalized linear model, the observation \eqn{y_i} is related,
  non-linearly, to some parameter \eqn{\psi_i}, which is modeled linearly:
  \eqn{x_i \beta}.  In a dynamic model, \eqn{beta} is allowed to change in time
  so that \eqn{\psi_i = x_i \beta_i} where \eqn{beta_i} is governed by some
  pre-specified dynamics.  For the purposes of this routine we assume that
  \eqn{beta_i \sim AR(1)}.  One may approximately draw \eqn{\beta | y} by using
  linear Bayes with a conjugate updating procedure and the using backward
  sampling formulas.

  The parameter \code{n} has different interpretations depending on the type of
  observation.  For "binom", \code{n} is the number of trials at each
  observation.  For "nbinom", \code{n} is the over-dispersion parameter at each
  response.  For "norm", \code{n} is the observation variance.

}

\usage{

CUBS.C(z, X, n, mu, phi, W, m0, C0, obs=c("binom", "nbinom", "norm"), eps.rel=1e-8, max.iter=100)
CUBS.R(z, X, n, mu, phi, W, m0, C0, obs=c("binom", "nbinom", "norm"), max.iter=100)

}

\arguments{
  \item{z}{An N dimensional vector; \eqn{z_i} is the response at \eqn{x_i}.}
  \item{X}{An N x P dimensional design matrix.  X must be a _matrix_.}
  \item{n}{An N dimensional vector; the interpretation of \eqn{n_i} changes with obs.}
  \item{mu}{A P.b dim. vector; the mean of beta.}
  \item{phi}{A P.b dim. vector; the autocorrelation of beta.}
  \item{W}{A P.b x P.b dim matrix; the innovation variance of beta.}
  \item{m0}{A P.b dimensional prior mean.}
  \item{C0}{A P.b x P.b dimensional prior precision.  C0 must be a _matrix_.}
  \item{obs}{The type of response.}
  \item{eps.rel}{The relative error condition used when calculating the conjugate update parameters.}
  \item{max.iter}{The maximum number of iterations used when calculating the conjugate update parameters.}
}	

\value{

  \item{alpha}{A P.a array; the posterior sample of the static regression
  coefficient.  When there is no static component this returns 0.}

  \item{beta}{A (N+1) x P.b array; the posterior sample of the dynamic regression
  coefficients.}

  \item{log.dens}{The log density (up to a shift) of the sample.}

}

\examples{

  data(rain)

  y = rain$y
  X = rain$X
  n = rain$n
  N = length(y+1)

  mu  = 0.0
  phi = 1.0
  W   = 0.1
  m0  = -1.0
  C0  = matrix(1.0)

  out = CUBS.C(y, X, n, mu, phi, W, m0, C0, "binom")		

  samp = 1
  beta = matrix(0, nrow=samp, ncol=N+1)

  for (i in 1:samp) {
    ## if (i%%100==0) {cat("Iteration", i, "\n")}
    out = CUBS.C(y, X, n, mu, phi, W, m0, C0, "binom")
    beta[i,] = out$beta
  }	

  qt = apply(beta, 2, function(x){quantile(x, c(0.25, 0.5, 0.75))})

  ## plot (qt[2,], col=1, type="l")
  ## lines(qt[1,], col=3)
  ## lines(qt[3,], col=3)
  
}

\references{

  Mike West, Jeff Harrison, and Helio S. Migon.  Dynamics Generalized Linear
  Models and Bayesian Forecasting.  Journal of the American Statistical
  Association, 1985.

}

\keyword{logit}
\keyword{regression}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.
