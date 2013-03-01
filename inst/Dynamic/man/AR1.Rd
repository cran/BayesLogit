% BayesLogit
\name{ar1.llh}

\alias{ar1.llh.C}

\title{AR1 Tools}

\description{

  Calculate log-likelihood of a VAR(1).

}

\details{

  Calculate log-density for a VAR(1) of the form 

  \deqn{\beta_t = \mu + \phi \odot (\beta_{t-1} - \mu) + \omega_t, \, \omega_t \sim N(0, W),} 

  where \eqn{\mu} and \eqn{\phi} are arrays and \eqn{W} is a
  matrix and \deqn{(\alpha, \beta_0) \sim N(m0, C0).}

}

\usage{

  ar1.llh.C(beta, mu, phi, W, m0, C0, alpha=NULL)

}

\arguments{

  \item{beta}{A P.b x (N+1) matrix.}
  \item{alpha}{A P.a dim vector.}
  \item{mu}{A P.b dim vector.}
  \item{phi}{A P.b dim vector.}
  \item{W}{A P.b x P.b matrix.}
  \item{m0}{A P dim vector}
  \item{C0}{A P x P dim matrix.}	

}

\value{

  Return the log-density, up to a shift by a function of C0.

}

\examples{

  N   = 100
  P.a = 1
  P.b = 2
  P   = P.a + P.b
  mu  = rep(0.3, P.b)
  phi = rep(0.95, P.b)
  W   = diag(0.1, P.b)
  m0  = rep(0, P)
  C0  = diag(1, P)

  ab  = m0 + t(chol(C0)) \%*\% rnorm(P, 0, 1)
  alpha = ab[1:P.a]
  beta0 = ab[(P.a+1):P.b]

  beta = matrix(0, nrow=P.b, ncol=N+1)
  beta[,1] = beta0
  
  U = t(chol(W))
  for (i in 2:(N+1)) {
    beta[,i] = mu + phi * (beta[,i-1] - mu) + (U \%*\% rnorm(P.b, 0, 1));
  }   

  llh = ar1.llh.C(beta, mu, phi, W, m0, C0, alpha=alpha)
  llh

}

% \seealso{
%   \code{\link{compute.mixture}}
% }

% \references{

% }

\keyword{draw.indicators}
\keyword{normal mixture}
