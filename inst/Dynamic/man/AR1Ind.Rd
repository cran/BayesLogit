% BayesLogit
\name{draw.ar1.ind}

\alias{draw.ar1.ind}
\alias{draw.mu.ar1.ind}
\alias{draw.phi.ar1.ind}
\alias{draw.W.ar1.ind}

\title{AR1 Tools}

\description{

  Generate posterior draw of parameters from a collection of
  independent AR(1) processes.

}

\details{

  These routines work for a collection of independent AR(1)'s, i.e. a
  VAR(1), 
  \deqn{x_t = \mu + \phi \odot (x_{t-1} - \mu) + \omega_t, \, \omega_t \sim N(0, W),}
  where \eqn{\mu} and \eqn{\phi} are arrays
  and \eqn{W} is diagonal.

  We assume the posterior draw for \eqn{\phi} must lie in \eqn{(0,1)}
  and thus return the previous value of \eqn{\phi} when that is not
  the case.

}

\usage{

  draw.W.ar1.ind(x, mu, phi, a.0, b.0)
  draw.phi.ar1.ind(x, mu, W, phi.m0, phi.P0, phi.prev)
  draw.mu.ar1.ind(x, phi, W, mu.m0, mu.P0)

}

\arguments{

  \item{x}{P x N+1 matrix.}
  \item{mu}{P dim vector}
  \item{phi}{P dim vector}
  \item{W}{A P dim vector.}
  \item{a.0}{A P dim vector corresponding to IGa(a.0/2, b.0/2) prior for W.}
  \item{b.0}{A P dim vector corresponding to IGa(a.0/2, b.0/2) prior for W.}
  \item{mu.m0}{A P dim vector corresponding to prior mean.}
  \item{mu.P0}{A P dim vector corresponding to component precisions.}
  \item{phi.m0}{A P dim vector corresponding to prior mean.}
  \item{phi.P0}{A P dim vector corresponding to component precisions.}
  \item{phi.prev}{Previous phi.}

}

\value{

  Return posterior draw, in each case a one-dimensional array.

}

\examples{

  N   = 100
  P.b = 2
  mu  = rep(0.3, P.b)
  phi = rep(0.95, P.b)
  W   = rep(0.1, P.b)
  m0  = rep(0, P.b)
  C0  = rep(1, P.b)

  beta0  = m0 + sqrt(C0) * rnorm(P.b, 0, 1)

  beta = matrix(0, nrow=P.b, ncol=N+1)
  beta[,1] = beta0
  
  for (i in 2:(N+1)) {
    beta[,i] = mu + phi * (beta[,i-1] - mu) + sqrt(W) * rnorm(P.b, 0, 1);
  }   

  W.vec = diag(W)
  mu  = draw.mu.ar1.ind (beta, phi, W.vec, rep(0, P.b), rep(0, P.b))
  phi = draw.phi.ar1.ind(beta, mu , W.vec, rep(0.9, P.b), rep(0, P.b), 0.5)
  W   = draw.W.ar1.ind  (beta, mu, phi, rep(2, P.b), rep(2, P.b))

}

% \seealso{
%   \code{\link{compute.mixture}}
% }

% \references{

% }

\keyword{draw.indicators}
\keyword{normal mixture}
