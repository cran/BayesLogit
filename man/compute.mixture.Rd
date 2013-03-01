% package: BayesLogit
\name{compute.mixture}

\alias{compute.mixture}
\alias{compute.mixture.lg}
\alias{compute.mixture.l3}

\title{Compute Mixture}

\description{

  Compute the means, variances, and probabilities to approximate a
  -log(Ga) distribution or a type III logistic distribution.

}

\details{

  Fruhwirth-Schnatter et al. use tables of discrete mixtures of
  normals to approximate entire families of distributions.  These
  approximations are useful for data augmentation techniques that lead
  to convenient posterior simulation for certain classes of
  generalized linear models.  In particular, approximations to the
  -log(Ga) distribution and the type III logistic distribution are
  useful for negative binomial regression and for binomial logistic
  regression respectively.

  This function generates the means, variances, and probabilities that
  approximate either a -log(Ga) or type III logistic distribution.

  The code that produces these values can be found in the functions
  \code{compute.mixture.lg} and \code{compute.mixture.l3}, which
  cannot be called directly.  \code{compute.mixture.lg} is an R
  translation of MATLAB code found in the \code{bayesf} package.
  \code{compute.mixture.l3} is taken from the R package
  \pkg{binomlogit}.  The function \code{compute.mixture.l3} returns a
  list without an array of means, since they are identically zero.

}

\usage{

compute.mixture(shape, type=c("log.gamma", "logistic.iii"))

}

\arguments{

  \item{shape}{The shape parameter.}

  \item{type}{The type of distribution to approximate, either a
  -log(Ga(shape,1)) distribution or a type III logistic(shape)
  distribution.}
	
}

\value{

  Returns a list with components \code{m}, \code{v}, \code{p}
  representing the mean, variance, and probability of the specific
  discrete mixture of normals.

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
  \code{\link{draw.indicators}}
}

\references{

  Agnes Fussl.  binomlogit: Efficient MCMC for Binomial Logit Models (2012).

  Sylvia Fruhwirth-Schnatter.  bayesf.  URL
  <http://statmath.wu.ac.at/~fruehwirth/monographie/> (2007).

  Sylvia Fruhwirth-Schnatter, Rudolf Fruhwirth, Leonhard Held, Havard
  Rue.  Improved auxiliary mixture sampling for hierarchical models of
  non-Gaussian data.  Statistics and Computing (2009).

  Sylvia Fruhwirth-Schnatter and Rudolf Fruhwirth. Data Augmentation
  and MCMC for Binary and Multinomial Logit Models.  Statistical
  Modelling and Regression Structures (2010), Springer-Verlag.

  Agnes Fussl, Sylvia Fruhwirth-Schnatter and Rudolf Fruhwirth.
  Efficient {MCMC} for Binomial Logit Models (2013).

}

\keyword{normal mixture}

%  Sylvia Fruhwirth-Schnatter. Finite Mixture and Markov Switch Models
%  (2007), Springer.