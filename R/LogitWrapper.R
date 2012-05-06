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

## dyn.load("Logit.so");

################################################################################
                               ## POLYAGAMMA ##
################################################################################

## Draw PG(a, b)
##------------------------------------------------------------------------------
rpg.gamma <- function(num=1, n=1, z=0.0, trunc=200)
{
    ## Check Parameters.
    if (sum(n<=0)!=0) {
        print("a must be greater than zero.");
        return(NA);
    }
    if (trunc < 1) {
        print("trunc must be > 0.");
        return(NA);
    }

    x = rep(0, num);

    if (length(n) != num) { n = array(n, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_gamma", x, n, z, as.integer(num), as.integer(trunc), PACKAGE="BayesLogit");

    OUT[[1]]
}

## Draw PG(n, z) where n is a natural number.
##------------------------------------------------------------------------------
rpg.devroye <- function(num=1, n=1, z=0.0)
{
    ## Check Parameters.
    if (sum(n<=0)!=0) {
        print("n must be greater than zero.");
        return(NA);
    }

    x = rep(0, num);

    if (length(n) != num) { n = array(n, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_devroye", x, as.integer(n), z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}

################################################################################
                                 ## Utility ##
################################################################################

## Check parameters to prevent an obvious error.
##------------------------------------------------------------------------------
check.parameters <- function(y, n, y.prior, x.prior, n.prior, R.X, C.X, samp, burn)
{
    ok = rep(TRUE, 8);
    ok[1] = all(y >= 0);
    ok[2] = all(n > 0);
    ok[3] = (y.prior >= 0);
    ok[4] = (n.prior >= 0);
    ok[5] = (length(y) == length(n) && length(y) == R.X);
    ok[6] = (length(x.prior) == C.X);
    ok[7] = (samp > 0);
    ok[8] = (burn >=0);
    ok[9] = all(y <= 1);

    if (!ok[1]) print("y must be >= 0.");
    if (!ok[9]) print("y is a proportion; it must be <= 1.");
    if (!ok[2]) print("n must be > 0.");
    if (!ok[3]) print(paste("y.prior must be >= 0: y.prior =", y.prior));
    if (!ok[4]) print(paste("n.prior must be >= 0. n.prior =", n.prior));
    if (!ok[5]) print(paste("Dimensions do not conform for y, X, and n.",
                            "len(y) =", length(y),
                            "dim(x) =", R.X, C.X,
                            "len(n) =", length(n)));
    if (!ok[6]) print(paste("x prior does not conform to data: x.prior =", x.prior));
    if (!ok[7]) print("samp must be > 0.");
    if (!ok[8]) print("burn must be >=0.");

    ok = all(ok)
}

## Combine
##------------------------------------------------------------------------------
logit.combine <- function(y, X, n=rep(1,length(y)),
                          y.prior=0.5, x.prior=colMeans(as.matrix(X)), n.prior=1.0)
{
    X = as.matrix(X);

    N = dim(X)[1];
    P = dim(X)[2];

    ok = check.parameters(y, n, y.prior, x.prior, n.prior, N, P, 1, 0);
    if (!ok) return(-1);

    ## Our combine_data function, written in C, uses t(X).
    tX = t(X);

    ## Add prior as extra data point.
    if (n.prior > 0) {
        y  = c(y, y.prior);
        tX = cbind(tX, x.prior);
        n  = c(n, n.prior);
        N  = N +1;
    }
    ## Don't need prior now.
    n.prior = 0.0

    OUT = .C("combine",
             as.double(y), as.double(tX), as.double(n),
             as.double(y.prior), as.double(x.prior), as.double(n.prior),
             as.integer(N), as.integer(P),
             PACKAGE="BayesLogit");

    N = OUT[[7]];

    y  = array(as.numeric(OUT[[1]]), dim=c(N));
    tX = array(as.numeric(OUT[[2]]), dim=c(P, N));
    n  = array(as.numeric(OUT[[3]]), dim=c(N));

    list("y"=as.numeric(y), "X"=t(tX), "n"=as.numeric(n));
}

################################################################################
                           ## POSTERIOR INFERENCE ##
################################################################################

## Posterior by Gibbs
##------------------------------------------------------------------------------
logit <- function(y, X, n=rep(1,length(y)),
                        y.prior=0.5, x.prior=colMeans(as.matrix(X)), n.prior=1.0,
                        samp=1000, burn=500)
{
    ## In the event X is one dimensional.
    X = as.matrix(X);

    ## Combine data.  We do this so that the auxiliary variable matches the
    ## data.
    new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;
    n.prior = 0.0;

    ## Check that the data and priors are okay.
    N = dim(X)[1];
    P = dim(X)[2];

    ok = check.parameters(y, n, y.prior, x.prior, n.prior, N, P, samp, burn);
    if (!ok) return(-1)

    ## Initialize output.
    output = list();

    ## w    = array(known.w, dim=c(N, samp));
    ## beta = array(known.beta, dim=c(P  , samp));
    w    = array(0.0, dim=c(N, samp));
    beta = array(0.0, dim=c(P  , samp));

    ## Our Logit function, written in C, uses t(X).
    tX = t(X);

    OUT = .C("gibbs",
             w, beta,
             as.double(y), as.double(tX), as.double(n),
             as.double(y.prior), as.double(x.prior), as.double(n.prior),
             as.integer(N), as.integer(P),
             as.integer(samp), as.integer(burn),
             PACKAGE="BayesLogit");

    N = OUT[[9]];

    tempw = array( as.numeric(OUT[[1]]), dim=c(N, samp) );

    output = list("w"=t(tempw), "beta"=t(OUT[[2]]), "y"=y, "X"=X, "n"=n);

    output
}

## Posterior mode by EM
##------------------------------------------------------------------------------
logit.EM <- function(y, X, n=rep(1, length(y)),
                     y.prior=0.5, x.prior=colMeans(as.matrix(X)), n.prior=1.0,
                     tol=1e-9, max.iter=100)
{

    ## In the event X is one dimensional.
    X = as.matrix(X);

    ## Check that the data and priors are okay.
    N = dim(X)[1];
    P = dim(X)[2];

    ok = check.parameters(y, n, y.prior, x.prior, n.prior, N, P, 1, 0);
    if (!ok) return(-1);

    ## Initialize output.
    beta = array(0, P);

    ## Our Logit function, written in C, uses t(X).
    tX = t(X);

    OUT = .C("EM",
             beta,
             as.double(y), as.double(tX), as.double(n),
             as.double(y.prior), as.double(x.prior), as.double(n.prior),
             as.integer(N), as.integer(P),
             as.double(tol), as.integer(max.iter),
             PACKAGE="BayesLogit");

    list("beta"=OUT[[1]], "iter"=OUT[[11]]);
}

################################################################################
                             ## Multinomial Case ##
################################################################################

## Check parameters to prevent an obvious error.
##------------------------------------------------------------------------------
mult.check.parameters <- function(y, X, n, m.0, P.0, samp, burn)
{
    ok = rep(TRUE, 6);
    ok[1] = all(y >= 0);
    ok[2] = all(n > 0);
    ok[3] = (nrow(y) == length(n) && nrow(y) == nrow(X));
    ok[4] = (samp > 0);
    ok[5] = (burn >=0);
    ok[6] = all(rowSums(y) <= 1);
    ok[7] = (ncol(y)==ncol(m.0) && ncol(X)==nrow(m.0));
    ok[8] = (ncol(X)==dim(P.0)[1] && ncol(X)==dim(P.0)[2] && ncol(y)==dim(P.0)[3]);

    if (!ok[1]) print("y must be >= 0.");
    if (!ok[6]) print("y[i,] are proportions and must sum <= 1.");
    if (!ok[2]) print("n must be > 0.");
    if (!ok[3]) print(paste("Dimensions do not conform for y, X, and n.",
                            "dim(y) =", nrow(y), ncol(y),
                            "dim(x) =", nrow(X), ncol(X),
                            "len(n) =", length(n)));
    if (!ok[4]) print("samp must be > 0.");
    if (!ok[5]) print("burn must be >=0.");
    if (!ok[7]) print("m.0 does not conform.");
    if (!ok[8]) print("P.0 does not conform.");

    ok = all(ok)
}

## Combine for multinomial logit.
##------------------------------------------------------------------------------

mlogit.combine <- function(y, X, n=rep(1,nrow(as.matrix(y))))
{
    X = as.matrix(X);
    y = as.matrix(y);
    
    N = dim(X)[1];
    P = dim(X)[2];
    J = dim(y)[2]+1;

    m.0=array(0, dim=c(ncol(X), ncol(y)));
    P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y)));
    
    ok = mult.check.parameters(y, X, n, m.0, P.0, 1, 0);
    if (!ok) return(NA);

    ## Our combine_data function, written in C, uses t(X), t(y).
    ty = t(y);
    tX = t(X);

    OUT = .C("mult_combine",
             as.double(ty), as.double(tX), as.double(n),
             as.integer(N), as.integer(P), as.integer(J),
             PACKAGE="BayesLogit");

    N = OUT[[4]];

    ty = array(as.numeric(OUT[[1]]), dim=c(J-1, N));
    tX = array(as.numeric(OUT[[2]]), dim=c(P, N));
    n  = array(as.numeric(OUT[[3]]), dim=c(N));

    list("y"=t(ty), "X"=t(tX), "n"=as.numeric(n));
}

## Posterior for multinomial logistic regression
##------------------------------------------------------------------------------
mlogit <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                   m.0=array(0, dim=c(ncol(X), ncol(y))),
                   P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
                   samp=1000, burn=500)
{
    ## In the event y or X is one dimensional.
    X = as.matrix(X);
    y = as.matrix(y);

    ## Combine data.  We do this so that the auxiliary variable matches the
    ## data.
    new.data = mlogit.combine(y, X, n);
    if (!is.list(new.data)) return(NA);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;

    N = dim(X)[1];
    P = dim(X)[2];
    J = dim(y)[2]+1;

    ## Check that the data and priors are okay.
    ok = mult.check.parameters(y, X, n, m.0, P.0, samp, burn);
    if (!ok) return(NA)

    ## Initialize output.
    output = list();

    ## w    = array(known.w, dim=c(N, samp));
    ## beta = array(known.beta, dim=c(P  , samp));
    w    = array(0.0, dim=c(N, J-1, samp));
    beta = array(0.0, dim=c(P, J-1, samp));

    ## Our Logit function, written in C, uses t(X), t(y).
    tX = t(X);
    ty = t(y);

    OUT = .C("mult_gibbs",
             w, beta,
             as.double(ty), as.double(tX), as.double(n),
             as.double(m.0), as.double(P.0),
             as.integer(N), as.integer(P), as.integer(J),
             as.integer(samp), as.integer(burn),
             PACKAGE="BayesLogit");

    N = OUT[[8]];

    ## Transpose for standard output.
    w    = array(0, dim=c(samp, N, J-1));
    beta = array(0, dim=c(samp, P, J-1));
    for (i in 1:samp) {
      w[i,,]    = OUT[[1]][,,i]
      beta[i,,] = OUT[[2]][,,i]
    }
    
    output = list("w"=w, "beta"=beta, "y"=y, "X"=X, "n"=n);

    output
}



