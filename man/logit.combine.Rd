% package: BayesLogit
\name{logit.combine}

\alias{logit.combine}

\title{Collapse Data for Binomial Logistic Regression}

\description{
  Collapse data for binomial logistic regression.
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
  given a total number of \eqn{n_i} observations at \eqn{x_i}.

  Thus, when a predictor is repeated the two reponses may be collapsed
  into a single observation representing multiple trials.  This
  function collapses data in this way.

}

\usage{

logit.combine(y, X, n=rep(1,length(y)))

}

\arguments{
  \item{y}{An N dimensional vector; \eqn{y_i} is the average response at \eqn{x_i}.}
  \item{X}{An N x P dimensional design matrix; \eqn{x_i} is the ith row.}
  \item{n}{An N dimensional vector; n_i is the number of observations at each \eqn{x_i}.}
}

\value{

  \code{logit.combine} returns a list.

  \item{y}{The new response.}

  \item{X}{The new design matrix.}

  \item{n}{The number of samples at each revised observation.}

}

\examples{

## From UCI Machine Learning Repository.
data(spambase);

## A subset of the data.
sbase = spambase[seq(1,nrow(spambase),10),];

X = model.matrix(is.spam ~ word.freq.free + word.freq.1999, data=sbase);
y = sbase$is.spam;

## Actually unnecessary as logit.EM automatically tries to compress.
new.data = logit.combine(y, X)
mode.spam = logit.EM(new.data$y, new.data$X, new.data$n)
mode.spam

}

\seealso{
\code{\link{logit}}, \code{\link{logit.EM}}, \code{\link{mlogit}}
}

% \references{

% }

\keyword{logit}
\keyword{regression}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.