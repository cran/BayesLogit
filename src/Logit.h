////////////////////////////////////////////////////////////////////////////////

// Copyright 2012 Nick Polson, James Scott, and Jesse Windle.

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or any later version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <http://www.gnu.org/licenses/>.

////////////////////////////////////////////////////////////////////////////////

// This class possess the functions needed to implement Polson and Scott's
// "Default Bayesian Logistic Regression."  See <http://arxiv.org/pdf/1109.4180>
// for details for the algorithm.

#ifndef __LOGIT__
#define __LOGIT__

#include "Matrix.h"
#include "RNG.h"
#include "PolyaGamma.h"
#include "Normal.h"
#include <stdexcept>
#include <list>
#include <time.h>
#include <algorithm>
#include <stdio.h>

// I need this to interrupt the Gibbs sampler from R.
#ifdef USE_R
#include <R_ext/Utils.h>
#endif

using std::list;

////////////////////////////////////////////////////////////////////////////////
				  // LOGIT //
////////////////////////////////////////////////////////////////////////////////

class Logit{

  // Variables.
  uint P;
  uint N;

  // Sufficient Statistics
  Matrix tX;
  Matrix n;
  Matrix Z;
  Matrix y; // Not a sufficient statistic, but need when combining data.

  // Random variates.
  Normal mvnorm;
  PolyaGamma pg;

public:

  // Constructors.
  Logit();
  Logit(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data,
	double y_prior, const Matrix& x_prior, double n_prior);

  // Utilities.
  void set_data(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data,
		double y_prior, const Matrix& x_prior, double  n_prior);

  bool data_conforms(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data,
		     const Matrix& x_prior);

  void get_data(Matrix& y_data, Matrix& tX_data, Matrix& n_data);

  uint get_N() { return N; }
  uint get_P() { return P; }

  // Gibbs sampling -- default
  inline void draw_w   (MF w, MF psi, RNG& r);
  inline void draw_beta(MF beta, MF w, MF beta_prev, RNG& r);
  void gibbs(Matrix& w, Matrix& beta, int samp, int burn, RNG& r);

  // Exepectation Maximization.
  int EM(Matrix& beta, double tol, int max_iter);

}; // Logit

////////////////////////////////////////////////////////////////////////////////
			       // CONSTRUCTORS //
////////////////////////////////////////////////////////////////////////////////

Logit::Logit()
{
  // We do not want to call this constructor.
  throw std::invalid_argument("Logit: default constructor called.");
} // Logit

Logit::Logit(const Matrix& y_data , const Matrix& tX_data, const Matrix& n_data ,
	           double  y_prior, const Matrix& x_prior,       double  n_prior)
  : mvnorm(tX_data.rows())
{
  set_data(y_data, tX_data, n_data, y_prior, x_prior, n_prior);
} // Logit

////////////////////////////////////////////////////////////////////////////////
				// UTILITIES //
////////////////////////////////////////////////////////////////////////////////

bool Logit::data_conforms(const Matrix& y_data, const Matrix& tX_data,
			  const Matrix& n_data, const Matrix& x_prior)
{
  bool ok = true;
  bool check[3];

  ok *= check[0] = y_data.area()  == tX_data.cols();
  ok *= check[1] = y_data.area()  == n_data.area();
  ok *= check[2] = x_prior.area() == tX_data.rows();

  for(int i = 0; i < 3; i++)
    if (!check[i]) Rprintf("Problem with check %i .\n", i);

  return ok;
}

void Logit::set_data(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data,
		     double y_prior, const Matrix& x_prior, double n_prior)
{

  // Check that the data is valid.
  if (!data_conforms(y_data, tX_data, n_data, x_prior)) {
    throw std::invalid_argument("set_data: data does not conform.");
  }

  Matrix tX_temp(tX_data);

  P = tX_data.rows();
  N = tX_data.cols();

  // Push everything into a list.
  list<double> ylist;
  list<Matrix> xlist;
  list<double> nlist;

  // Our data should not have n_data(i)=0.
  for(uint i=0; i < N; ++i){
      ylist.push_back(y_data(i));
      xlist.push_back(tX_temp.col(i));
      nlist.push_back(n_data(i));
  }

  // Presumably we could have a flat prior, i.e. n_prior=0.
  if (n_prior > 0) {
    N = N + 1; // Add one more "data" point.
    ylist.push_back(y_prior);
    xlist.push_back(x_prior);
    nlist.push_back(n_prior);
  }

  // Merge data.
  list<double>::iterator y_i;
  list<Matrix>::iterator x_i;
  list<double>::iterator n_i;

  list<double>::iterator y_j;
  list<Matrix>::iterator x_j;
  list<double>::iterator n_j;

  x_i = xlist.begin();
  y_i = ylist.begin();
  n_i = nlist.begin();

  while(x_i != xlist.end()){
    x_j = x_i; y_j = y_i; n_j = n_i;
    ++x_j; ++y_j; ++n_j;
    while(x_j != xlist.end()){
      if (*x_i == *x_j) {
  	double sum = *n_i + *n_j;
   	*y_i = (*n_i/sum) * *y_i + (*n_j/sum) * *y_j;
  	*n_i = sum;
  	// Increment THEN erase!
	// Actually, send pointer to erase, then increment, then erase.
  	xlist.erase(x_j++);
  	ylist.erase(y_j++);
  	nlist.erase(n_j++);
      }
      else {
  	++x_j; ++y_j; ++n_j;
      }

    }
    ++x_i; ++y_i; ++n_i;
  }

  uint old_N = N; // Record to see if we have changed data.

  // Set up y, tX, n.
  N = xlist.size();

  // cout << "Old N: " << old_N << " N: " << N << "\n";
  // Warning...
  if (old_N != N) {
    Rprintf("Warning: data was combined!\n");
    Rprintf("N: %i, P: %i \n", N, P);
  }

  // Matrix y(N);
  y.resize(N);
  tX.resize(P, N);
  n.resize(N);

  x_i = xlist.begin();
  y_i = ylist.begin();
  n_i = nlist.begin();
  for(uint i = 0; i < N; ++i){
    y(i)      = *y_i++;
    tX.col(i) = *x_i++;
    n(i)      = *n_i++;
  }

  // cout << "y:\n" << y;
  // cout << "tX:\n" << tX;
  // cout << "n:\n" << n;

  // Set up Z and alpha.
  Matrix alpha(N);
  Z.resize(P);
  for(uint i = 0; i < N; ++i)
    alpha(i) = n(i) * (y(i) - 0.5);
  gemm(Z, tX, alpha);

}

void Logit::get_data(Matrix& y_data, Matrix& tX_data, Matrix& n_data)
{
  y_data  = y;
  tX_data = tX;
  n_data  = n;
}

////////////////////////////////////////////////////////////////////////////////
			    // POSTERIOR BY GIBBS //
////////////////////////////////////////////////////////////////////////////////

inline void Logit::draw_w(MF w, MF psi, RNG& r)
{
  for(int i = 0; i < (int)N; ++i) {
    // w(i) = pg.draw(n(i), psi(i), r);
    w(i) = pg.draw((int)n(i), psi(i), r);
  }
}

inline void Logit::draw_beta(MF beta, MF w, MF beta_prev, RNG& r)
{
    Matrix tXOmega(tX);
    prodonrow(tXOmega, w);
    Matrix tXOmX(tXOmega, tX, 'N', 'T');

    // Joint draw.
    mvnorm.set_from_likelihood(Z, tXOmX);
    mvnorm.draw(beta[0], r);

    // // Gibbs sampling.
    // beta.copy(beta_prev);
    // for(uint i = 0; i < P; ++i){
    //   double m_i = Z(i);
    //   for(uint j = 0; j < P; ++j){
    // 	m_i -= j != i ? tXOmX(i,j) * beta(j) : 0.0;
    //   }
    //   m_i /= tXOmX(i,i);
    //   beta(i) = r.norm(m_i, sqrt(1/tXOmX(i,i)));
    // }
}

// #define GIBBS_CORE(CIDX, PIDX)
//   draw_w   (w[CIDX]   ,  psi, r);
//   draw_beta(beta[CIDX], w[CIDX], beta[PIDX], r);
//   gemm(psi, X, beta[CIDX], 'T');
//   R_CheckUserInterrupt(void);

// Gibbs sampling -- Default Logit.
void Logit::gibbs(Matrix& w, Matrix& beta, int samp, int burn, RNG& r)
{
  uint M = (uint)samp;

  w.resize(N, 1, M);
  beta.resize(P, 1, M);
  // beta.fill(-1.0);

  Matrix psi(N);
  gemm(psi, tX, beta[0], 'T');
  // psi.fill(0.0);

  // Keep track of time.
  clock_t start, end;

  start = clock();
  // Burn-in - Take an extra for first sample of MCMC.
  for(int m = 0; m < burn+1; ++m){
    draw_w   (w[0]   ,  psi, r);
    draw_beta(beta[0], w[0], beta[0], r);
    gemm(psi, tX, beta[0], 'T');
    // In case we are using R.
    #ifdef USE_R
    if (m%1==0) R_CheckUserInterrupt();
    #endif
  }
  end = clock();

  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  Rprintf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  Rprintf("Expect approx. %g sec. for %i samples.\n", total_time * samp / burn, samp);

  start = clock();
  // Sample - Already took one sample from burn-in.
  for(int m = 1; m < samp; ++m){
    draw_w   (w[m]   ,  psi, r);
    draw_beta(beta[m], w[m], beta[m-1], r);
    gemm(psi, tX, beta[m], 'T');
    // In case we are using R.
    #ifdef USE_R
    if (m%1==0) R_CheckUserInterrupt();
    #endif
  }
  end = clock();

  total_time = (double)(end - start) / CLOCKS_PER_SEC;
  Rprintf("Sampling complete: %g sec. for %i iterations.\n", total_time, samp);

}

////////////////////////////////////////////////////////////////////////////////
			   // POSTERIOR MODE BY EM //
////////////////////////////////////////////////////////////////////////////////

// Exepectation Maximization.
int Logit::EM(Matrix& beta, double tol, int max_iter)
{
  Matrix psi(N);
  Matrix w(N);
  double dist = tol + 1.0;

  // Proper size.
  beta.resize(P);

  int  iter = 0;
  while (dist > tol && iter < max_iter) {

    // w: posterior mean
    gemm(psi, tX, beta, 'T', 'N');
    for (int i = 0; i < (int)N; ++i) {
      double hpsi = psi(i) * 0.5;
      // n * tanh(psi/2) / (psi/2) * 0.5
      if ( hpsi < 0.01 ) {
	w(i) = n(i) / cosh(hpsi)
	  * (1 + hpsi*hpsi / 6.0 + pow(hpsi, 4.0) / 120.0 + pow(hpsi, 6) / 5040.0)
	  * 0.25 ;
      }
      else
	w(i) = n(i) * tanh(hpsi) / hpsi * 0.25;
    }

    // beta: posterior mode
    Matrix old_beta(beta);

    Matrix tXOmega(tX);
    prodonrow(tXOmega, w);
    Matrix tXOmX(tXOmega, tX, 'N', 'T');
    beta.clone(Z); symsolve(tXOmX, beta);

    // Check how much we improved.
    // Matrix diff = beta - old_beta;
    // dist = sqrt( dot(diff, diff) );
    Matrix diff = fabs(beta - old_beta);
    dist = *std::max_element(diff.begin(), diff.end());

    ++iter;
  }

  return iter;
}

////////////////////////////////////////////////////////////////////////////////
			       // END OF CODE //
////////////////////////////////////////////////////////////////////////////////

#endif

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

/*

  The biggest bottleneck in our algorithm appears to be drawing from the Poyla
  Gamma distribution.  There may not be anyway around this.  I tried sampling
  all the gamma variables at once, but this didn't seem to improve things.  In
  fact, it seemed to make things take longer.  Basically, it seems that you the
  overhead to create an array to to sample a bunch of gammas at "one time" does
  not overcome any optential speedup from preenting repeatedly calling the
  function.  Of course it may be that the compiler is optimizing things and
  inlining things.

  Basically, it appears to just be a problem to sample a Poyla Gamma.

  UPDATE: You can sample a PG more quickly following the method of Devroye.

  ---

  The algorithm specifies drawing a block of omegas and then a block of psi's.
  But due to the independence of each we can sample them in an alternatating
  fashion.  But this doesn't work.  Why?  Too many degrees of freedom when using
  psi?  We could parallelize this... if this worked.

 */
