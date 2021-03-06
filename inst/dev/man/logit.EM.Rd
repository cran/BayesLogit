% package: BayesLogit
\name{logit.EM}

\alias{logit.EM}
\alias{logit.EM.R}

\title{Logistic Regression Expectation Maximization}

\description{
  Expectation maximization for logistic regression.
}

\details{

  Logistic regression is a classification mechanism.  Given the binary
  data \eqn{\{y_i\}} and the p-dimensional predictor variables
  \eqn{\{x_i\}}, one wants to forecast whether a future data point y*
  observed at the predictor x* will be zero or one.  Logistic
  regression stipulates that the statistical model for observing a
  success=1 or failure=0 is governed by

  \deqn{ P(y^* = 1 | x^*, \beta) = (1 + \exp(-x^* \beta))^{-1}. }

  Instead of representing data as a collection of binary outcomes, one
  may record the average response \eqn{y_i} at each unique \eqn{x_i}
  given a total number of \eqn{n_i} observations at \eqn{x_i}.  We
  follow this method of encoding data.

  A non-informative prior is used.

}

  % Polson and Scott suggest placing a Jeffrey's Beta prior
  % Be(1/2,1/2) on

  % \deqn{ m(\beta) := P(y_0 = 1 | x_0, \beta) = (1 + \exp(-x_0 \beta))^{-1}, }

  % which generates a Z-distribution prior for \eqn{\beta},

  % \deqn{ p(\beta) = \exp(0.5 x_0 \beta) / (1 + \exp(0.5 x_0 \beta)). }

  % One may interpret this as "prior" data where the average response at
  % \eqn{x_0} is \eqn{1/2} based upon a "single" observation.  The
  % default value of \eqn{x_0=mean(x), x=\{x_i\}}.

\usage{
logit.EM(y, X, n=rep(1,length(y)), tol=1e-9, max.iter=100)
}

\arguments{
  \item{y}{An N dimensional vector; \eqn{y_i} is the average response at \eqn{x_i}.}
  \item{X}{An N x P dimensional design matrix; \eqn{x_i} is the ith row.}
  \item{n}{An N dimensional vector; n_i is the number of observations at each \eqn{x_i}.}
  \item{tol}{Threshold at which algorithm stops.}
  \item{max.iter}{Maximum number of iterations.}
}

  % \item{y.prior}{Average response at \code{x.prior}.}
  % \item{x.prior}{Prior predictor variable.}
  % \item{n.prior}{Number of observations at \code{x.prior}.}

\value{

  \item{beta}{The posterior mode.}

  \item{iter}{The number of iterations.}

}

\examples{

## From UCI Machine Learning Repository.
data(spambase);

## A subset of the data.
sbase = spambase[seq(1,nrow(spambase),10),];

X = model.matrix(is.spam ~ word.freq.free + word.freq.1999, data=sbase);
y = sbase$is.spam;

## Run logistic regression.
output = logit.EM(y, X);

}

\seealso{
  \code{\link{rpg}}, \code{\link{logit}}, \code{\link{mlogit}}
}

\references{
  Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian inference for
  logistic models using Polya-Gamma latent variables.
  \url{http://arxiv.org/abs/1205.0310}

  Nicholas G. Poslon and James G. Scott.  Default Bayesian analysis for
  multi-way tables: a data-augmentation approach.
  \url{http://arxiv.org/pdf/1109.4180}
}

\keyword{logit}
\keyword{regression}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.
