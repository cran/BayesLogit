% TAKEN FROM UCI MACHINE LEARNING REPOSITORY

\name{spambase}
\alias{spambase}
\title{Spambase Data}

\description{

The \code{spambase} data has 57 real valued explanatory variables which
characterize the contents of an email and and one binary response
variable indicating if the email is spam.  There are 4601 observations.

}

\format{
  A data frame: the first column is a binary response variable
  indicating if the email is spam.  The remaining 57 columns are real
  valued explanatory variables.
}

\details{

Of the 57 explanatory variables, 48 describe word frequency, 6 describe
character frequency, and 3 describe sequences of capital letters.

  \describe{
    \item{word.freq.<word>}{A continuous explanatory variable describing
      the frequency with which the word <word> appears; measured in percent.}
    \item{char.freq.<char>}{A continuous explanatory variable describing
      the frequency with which the character <char> appears; measured in
      percent.}
    \item{capital.run.length.<oper>}{A statistic involving the
      length of consecutive capital letters.}
  }

Use \code{names} to see the specific words, characters, or statistics for
each respective class of variable.
}
% \source{}

\references{
Mark Hopkins, Erik Reeber, George Forman, and Jaap Suermondt of
Hewlett-Packard Labs (1999).  Spambase Data Set.
\url{http://archive.ics.uci.edu/ml/datasets/Spambase}

Frank, A. & Asuncion, A. (2010). UCI Machine Learning Repository
[http://archive.ics.uci.edu/ml]. Irvine, CA: University of California,
School of Information and Computer Science.
}
\keyword{datasets}
