################################################################################

## Copyright 2012 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit.

## BayesLogit is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or any later version.

## BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## BayesLogit.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
                            ## POSTERIOR BY GIBBS ##
################################################################################

## Bayesian logistic regression
##------------------------------------------------------------------------------
logit.R <- function(y, X, n=rep(1, length(y)),
                          y.prior=0.5, x.prior=colMeans(as.matrix(X)), n.prior=1.0,
                          samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, avg response
  ## n: n by 1 vector, # of obs at distinct x

  ## Combine data.
  new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
  y = new.data$y;
  X = new.data$X;
  n = new.data$n;
  n.prior = 0.0;

  ## X = as.matrix(X);

  p = ncol(X)
  N = nrow(X)

  alpha = (y-1/2)*n

  Z = colSums(X*alpha)
  ## PsiToBeta = solve(t(X) %*% X) %*% t(X);

  w = rep(0,N)
  ## w = w.known;
  beta = rep(0.0, p)

  output <- list(w = matrix(nrow=samp, ncol=N),
                 beta = matrix(nrow=samp, ncol=p)
                 )

  c_k = (1:200-1/2)^2 * pi^2 * 4;

  ## Sample
  for ( j in 1:(samp+burn) )
  {

    ## draw w
    psi = drop(X%*%beta)
    ## Sum of gamma: poor approximation when psi is large!  Causes crash.
    ## Devroye is faster anyway.
    w = rpg.devroye(N, n, psi);

    ## # draw beta - Gibbs sampling
    ## ups <- t(X) %*% (X*w)
    ## for ( i in 1:p )
    ## {
    ##   mi = Z[i]/ups[i,i] - crossprod(ups[i,-i],beta[-i])/ups[i,i]
    ##   vi = 1/ups[i,i]
    ##   beta[i] = rnorm(1,mean=mi,sd=sqrt(vi))
    ## }

    # draw beta - Joint Sample.
    PC = t(X) %*% (X * w);
    ## S = solve(PC); ## chol2inv works better for larger P?
    S = chol2inv(chol(PC));
    m = S %*% as.vector(Z);
    beta = m + t(chol(S)) %*% rnorm(p);

    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,] <- w
        output$beta[j-burn,] <- beta
    }

    if (j %% verbose == 0) { print(paste("Iteration", j)); }
  }

  ## Add new data to output.
  output$"y" = y;
  output$"X" = X;
  output$"n" = n;

  output
} ## logit.gibbs.R

## Bayesian logistic regression - Normal Prior
##------------------------------------------------------------------------------
logit.gibbs.np.R <- function(y, X, n=rep(1, length(y)),
                          P.0 = diag(0.0, ncol(X)),
                          samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, avg response
  ## n: n by 1 vector, # of obs at distinct x

  ## DO NOT USE DEFAULT PRIOR
  y.prior=0.5;
  x.prior=colMeans(as.matrix(X));
  n.prior=0.0;
  
  ## Combine data.
  new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
  y = new.data$y;
  X = new.data$X;
  n = new.data$n;
  n.prior = 0.0;

  ## X = as.matrix(X);

  p = ncol(X)
  N = nrow(X)

  alpha = (y-1/2)*n

  Z = colSums(X*alpha)

  w = rep(0,N)
  ## w = w.known;
  beta = rep(0.0, p)

  output <- list(w = matrix(nrow=samp, ncol=N),
                 beta = matrix(nrow=samp, ncol=p)
                 )

  ## Sample
  for ( j in 1:(samp+burn) )
  {

    ## draw w
    psi = drop(X%*%beta)
    for ( i in 1:N )
    {
        w[i] = rpg.devroye(1, n[i], psi[i]);
    }

    # draw beta - Joint Sample.
    PC = t(X) %*% (X * w) + P.0;
    ## S = solve(PC); ## chol2inv works better for larger P?
    S = chol2inv(chol(PC));
    m = S %*% as.vector(Z);
    beta = m + t(chol(S)) %*% rnorm(p);

    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,] <- w
        output$beta[j-burn,] <- beta
    }

    if (j %% verbose == 0) { print(paste("Iteration", j)); }
  }

  ## Add new data to output.
  output$"y" = y;
  output$"X" = X;
  output$"n" = n;

  output
} ## logit.gibbs.np.R

################################################################################
                                 ## TESTING ##
################################################################################

#data = read.table("orings.dat",header=TRUE)
#attach(data)
#failure = 2*failure-1
## x = c(53,56,57,63,66,67,68,69, 70,72,73, 75,76,78,79,80,81)
## y = c( 1, 1, 1, 0, 0, 0, 0, 0,3/4, 0, 0,1/2, 0, 0, 0, 0, 0)
## n = c( 1, 1, 1, 1, 1, 3, 1, 1,  4, 1, 1,  2, 2, 1, 1, 1, 1)
## ans = logit.MCMC(100000,cbind(1,x),y,n)

## hist(ans$beta[,1])
## hist(ans$beta[,2])

## mean(ans$beta[,1])
## mean(ans$beta[,2])
