## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


pigauss <- function(x, mu, lambda=1.0)
{
    Z = 1.0 / mu;
    b = sqrt(lambda / x) * (x * Z - 1);
    a = -1.0 * sqrt(lambda / x) * (x * Z + 1);
    y = exp(pnorm(b, log.p=TRUE)) + exp(2 * lambda * Z + pnorm(a, log.p=TRUE));
    # y2 = 2 * pnorm(-1.0 / sqrt(x));
    y
}

rtigauss <- function(Z, R=Inf)
{
    Z = abs(Z);
    mu = 1/Z;
    X = R + 1;
    if (mu > R) {
        alpha = 0.0;
        while (runif(1) > alpha) {
            ## X = R + 1
            ## while (X > R) {
            ##     X = 1.0 / rgamma(1, 0.5, rate=0.5);
            ## }
            E = rexp(2)
            while ( E[1]^2 > 2 * E[2] / R) {
                E = rexp(2)
            }
            X = R / (1 + R*E[1])^2
            alpha = exp(-0.5 * Z^2 * X);
        }
    }
    else {
        while (X > R) {
            lambda = 1.0;
            Y = rnorm(1)^2;
            X = mu + 0.5 * mu^2 / lambda * Y -
                0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
            if ( runif(1) > mu / (mu + X) ) {
                X = mu^2 / X;
            }
        }
    }
    X;
}

rtigauss.2 <- function(Z, R=Inf)
{
  mu = 1/Z;
  t  = R;

  tx = (t-mu)^2 / (mu^2 * t);
  tn = sqrt(tx);
  X  = 0;
  
  if (mu > R) {
    Z = rtnorm(1, mean=0, sd=1, lower=tn);
    Y = Z^2;
    X = mu + 0.5 * mu^2 * Y +
        0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
    X = mu^2 / X;
  }
  else {
    Z = rtnorm(1, mean=0, sd=1, lower=-1*tn);
    Y = Z^2;
    X2 = mu + 0.5 * mu^2 * Y +
         0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
    X1 = mu + 0.5 * mu^2 * Y -
         0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
    X = X1
    if ( Z < tn ) {
      if ( runif(1) > mu / (mu + X1) ) {
        X = X2;
      }
    }
    else {
      ## if ( runif(1) > mu / (mu + X1) ) {
      ##   X = -0.01;
      ## }
    }
  }
  
  X
}

rtigauss.3 <- function(Z, R=Inf)
{
  mu = 1/Z;
  t  = R;

  tx = (t-mu)^2 / (mu^2 * t);
  tl = min(tx, mu^2 / tx);
  tn = sqrt(tx);
  X  = 0;
  
  if (mu > R) {
    Z = rtnorm(1, mean=0, sd=1, lower=tn);
    Y = Z^2;
    X = mu + 0.5 * mu^2 * Y +
        0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
    X = mu^2 / X;
  }
  else {
    
    U = runif(1);

    ## Take the left side.
    if (U < pigauss(tl, mu) / pigauss(t,mu)) {

      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 * Y -
          0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
      
    }
    ## Take the right side.
    else {
      Z = rtnorm(1, mean=0, sd=1, lower=-1*tn, upper=tn);
      Y = Z^2;
      X = mu + 0.5 * mu^2 * Y +
          0.5 * mu * sqrt(4 * mu * Y + (mu * Y)^2);
      
    }
    
  }
  
  X
}


N = 10000;
Z = 1.0;
R = 1.0;

samp.1 = rep(0, N);
samp.2 = rep(0, N);
samp.3 = rep(0, N);

for (i in 1:N) {
  samp.1[i] = rtigauss(Z, R);
  samp.2[i] = rtigauss.2(Z, R);
  samp.3[i] = rtigauss.3(Z, R);
}

par(mfrow=c(1,3));
hist(samp.1, breaks=20)
hist(samp.2, breaks=20)
hist(samp.3, breaks=20)
