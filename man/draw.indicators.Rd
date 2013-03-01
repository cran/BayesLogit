% package: BayesLogit
\name{draw.indicators}

\alias{draw.indicators}
\alias{draw.indicators.C}
\alias{draw.indicators.R}

\title{Draw Indicators}

\description{

  Draw indicator variables from a normal mixture.  

}

\details{

  Suppose \eqn{e_i \sim N(m_{\gamma_i}, v_{\gamma_i})}.  Then, given
  \eqn{e_i}, draw \eqn{\gamma_i} from \eqn{P(\gamma_i = j | e_i)}
  for each \eqn{i}.

}

\usage{

draw.indicators(res, nmix)
draw.indicators.C(res, nmix)
draw.indicators.R(res, nmix)

}

\arguments{

  \item{res}{A one-dimensional array of residuals.}

  \item{nmix}{A list representing the mixture of normals with
  components \code{m}, \code{v}, \code{p}, representing the arrays of
  means, variances, weights of the normal mixture.}
	
}

\value{

  Returns a one-dimensional array of mixture component identifiers.

}

\examples{

## Approximate - log[Ga] using mixture of normals.
sp = 9
nm = compute.mixture(sp, "log.gamma")
nc = length(nm$m)

r  = sample.int(nc, 10, replace=TRUE, prob=nm$p)
e  = rnorm(10, nm$m[r], sqrt(nm$v[r]))

r.post = draw.indicators(e, nm)

## Approximate type III logistic using mixture of normals.
sp = 9
nm = compute.mixture(sp, "logistic.iii")
nc = length(nm$m)

r  = sample.int(nc, 10, replace=TRUE, prob=nm$p)
e  = rnorm(10, nm$m[r], sqrt(nm$v[r]))

r.post = draw.indicators(e, nm)

}

\seealso{
  \code{\link{compute.mixture}}
}

% \references{

% }

\keyword{draw.indicators}
\keyword{normal mixture}
