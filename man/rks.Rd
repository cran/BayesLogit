% package: BayesLogit
\name{rks}

\alias{rks}

\title{The Kolmogorov-Smirnov distribution}

\description{
  Generate a random variate from the Kolmogorov-Smirnov distribution.

  This is not directly related to the Polya-Gamma technique, but it is
  a nice example of using an alternating sum to generate a random
  variate.  

}

\details{

  The density function of the KS distribution is 

  \deqn{ f(x) = 8 \sum_{i=1}^\infty (-1)^{n+1} n^2 x e^{-2 n^2 x^2}. }

  We follow Devroye (1986) p. 161 to generate random draws from KS.

}

\usage{

rks(N=1)

}

\arguments{
  \item{N}{The number of random variates to generate.}
}

\examples{

X = rks(1000)

}

\references{
  L. Devroye.  Non-Uniform Random Variate Generation, 1986.
}

\keyword{rks}
