% package: BayesLogit
\name{mlogit}

\alias{mlogit}
\alias{mlogit.R}

\title{Bayesian Multinomial Logistic Regression}

\description{
  Run a Bayesian multinomial logistic regression.
}

\details{

  Multinomial logistic regression is a classifiction mechanism.  Given
  the multinomial data \eqn{\{y_i\}} with J categories and the
  p-dimensional predictor variables \eqn{\{x_i\}}, one wants to
  forecast whether a future data point y* at the predictor x*.
  Multinomial Logistic regression stiuplates that the statistical
  model for observing a draw category j after rolling the multinomial
  die \eqn{n^*=1} time is governed by

  \deqn{ P(y^* = j | x^*, \beta, n^*=1) = e^{x^* \beta_j} / \sum_{k=1}^J e^{x^* \beta_k}. }

  Instead of representing data as the total number of responses in
  each category, one may record the average number of responses in
  each category and the total number of responses \eqn{n_i} at
  \eqn{x_i}.  We follow this method of encoding data.

  We assume that \eqn{\beta_J = 0} for purposes of identification!

  You may use mlogit for binary logistic regression with a normal prior.
}

\usage{

mlogit(y, X, n=rep(1,nrow(as.matrix(y))),
       m.0=array(0, dim=c(ncol(X), ncol(y))),
       P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
       samp=1000, burn=500)

}

% mlogit.R(y, X, n=rep(1,nrow(as.matrix(y))),
%          m.0=array(0, dim=c(ncol(X), ncol(y))),
%          P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
%          samp=1000, burn=500, verbose=1000)

\arguments{
  \item{y}{An N x J-1 dimensional matrix; 
  	   \eqn{y_{ij}} is the average response for category j at \eqn{x_i}.}
  \item{X}{An N x P dimensional design matrix; \eqn{x_i} is the ith row.}
  \item{n}{An N dimensional vector; \eqn{n_i} is the total number of observations at each \eqn{x_i}.}
  \item{m.0}{A P x J-1 matrix with the \eqn{beta_j}'s prior means.}
  \item{P.0}{A P x P x J-1 array of matrices with the \eqn{beta_j}'s prior precisions.}
  \item{samp}{The number of MCMC iterations saved.}
  \item{burn}{The number of MCMC iterations discarded.}
  % \item{verbose}{Print message every \code{verbose} number of iterations.}
}

\value{

  \code{mlogit} returns a list.

  \item{beta}{A samp x P x J-1 array; the posterior sample of the regression
  coefficients.}

  \item{w}{A samp x N' x J-1 array; the posterior sample of the latent variable.
  WARNING: N' may be less than N if data is combined.}

  \item{y}{The response matrix--different than input if data is combined.}

  \item{X}{The design matrix--different than input if data is combined.}

  \item{n}{The number of samples at each observation--different than input if
  data is combined.}

}

\examples{

## Use the iris dataset.
data(iris)
N = nrow(iris)
P = ncol(iris)
J = nlevels(iris$Species)

X     = model.matrix(Species ~ ., data=iris);
y.all = model.matrix(~ Species - 1, data=iris);
y     = y.all[,-J];

out = mlogit(y, X, samp=1000, burn=100);

}

\seealso{
\code{\link{rpg}}, \code{\link{logit.EM}}, \code{\link{logit}}
}

\references{
  Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian inference for
  logistic models using Polya-Gamma latent variables.
  \url{http://arxiv.org/abs/1205.0310}
}

\keyword{logit}
\keyword{regression}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.
