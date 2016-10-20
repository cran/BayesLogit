#include "DynExpFamMH.h"
#include <iostream>

using std::pow;
using std::fabs;
using std::sqrt;
using std::log;
using std::exp;

void bidc_to_midc(int* mstart, int* mnum, const int bstart, const int bnum, const int sob)
{
  // Assume blocks are sob x sob.
  // bstart: block index start.
  // sob: size of block.
  // num: number of blocks.
  *mstart = sob * bstart;
  *mnum   = sob * bnum;
}

void make_A(MatrixXd& A, const VectorXd& phi, const unsigned int num)
{
  // Assume A is 0's.

  unsigned int bsize = phi.size();
  unsigned int N = num * bsize;

  A.resize(N, N);
  for (unsigned int i=0; i<N; i++) {
    A(i,i) = 1.0;
    if (i < N - bsize) 
      A(i+bsize,i) = -1.0 * phi(i % bsize);
  }
}

// void log_logit_likelihood(const double* y, const double* ntrials, llh_view& llh, int block_start, int num_blocks)
// {
//   int block_end = block_start + num_blocks;

//   if (block_end > (int)llh.N) Rprintf( "Warning: log_logit_likelihood block_end out of bounds.\n");

//   for (int i = block_start; i < block_end; i++) {
//     double ntrialsd = (double) ntrials[i];
//     llh.psi[i]  = llh.psi_stc[i] + llh.psi_dyn[i];
//     double psi  = llh.psi[i];
//     double epsi = exp(psi);
//     double p    = epsi / (1 + epsi);
//     llh.l0[i]   = psi * y[i] - ntrialsd * log(1.0 + epsi);
//     llh.l1[i]   = y[i] - ntrialsd * p;
//     llh.l2[i]   = -1.0 * ntrialsd * (p - p*p);
//     llh.pl2[i]  = psi * llh.l2[i];
//     // Rprintf("i: %2i, psi: %3.2g, y: %3.2g, l0: %3.2g, l1: %3.2g, l2: %3.2g\n", i, psi, y[i], llh.l0[i], llh.l1[i], llh.l2[i]);
//   }
// }

// void log_logit_likelihood(const double* y, const double* ntrials, llh_view& llh, int block_start)
// {
//   log_logit_likelihood(y, ntrials, llh, block_start, (int)llh.N - block_start);
// }

void log_logit_likelihood(const double* y, const double* ntrials, llh_struct& llh, int block_start, int num_blocks)
{
  int block_end = block_start + num_blocks;

  if (block_end > llh.psi.size()) Rprintf( "Warning: log_logit_likelihood block_end out of bounds.\n");

  for (int i = block_start; i < block_end; i++) {
    double ntrialsd = (double) ntrials[i];
    llh.psi[i]  = llh.psi_stc[i] + llh.psi_dyn[i];
    double psi  = llh.psi[i];
    double epsi = exp(psi);
    double p    = epsi / (1 + epsi);
    llh.l0[i]   = psi * y[i] - ntrialsd * log(1.0 + epsi);
    llh.l1[i]   = y[i] - ntrialsd * p;
    llh.l2[i]   = -1.0 * ntrialsd * (p - p*p);
    llh.pl2[i]  = psi * llh.l2[i];
    // Rprintf("i: %2i, psi: %3.2g, y: %3.2g, l0: %3.2g, l1: %3.2g, l2: %3.2g\n", i, psi, y[i], llh.l0[i], llh.l1[i], llh.l2[i]);
  }
}

void log_logit_likelihood(const double* y, const double* ntrials, llh_struct& llh, int block_start)
{
  log_logit_likelihood(y, ntrials, llh, block_start, llh.psi.size() - block_start);
}

void log_nbinom_p_likelihood(const double* y, const double* d, llh_struct& llh, int block_start, int num_blocks)
{
  int block_end = block_start + num_blocks;

  if (block_end > llh.psi.size()) Rprintf( "Warning: log_logit_likelihood block_end out of bounds.\n");

  for (int i = block_start; i < block_end; i++) {
    double bexpon = y[i] + d[i];
    llh.psi[i]  = llh.psi_stc[i] + llh.psi_dyn[i];
    double psi  = llh.psi[i];
    double epsi = exp(psi);
    double p    = epsi / (1 + epsi);
    llh.l0[i]   = psi * y[i] - bexpon * log(1.0 + epsi);
    llh.l1[i]   = y[i] - bexpon * p;
    llh.l2[i]   = -1.0 * bexpon * (p - p*p);
    llh.pl2[i]  = psi * llh.l2[i];
    // Rprintf("i: %2i, psi: %3.2g, y: %3.2g, l0: %3.2g, l1: %3.2g, l2: %3.2g\n", i, psi, y[i], llh.l0[i], llh.l1[i], llh.l2[i]);
  }
}

void log_nbinom_p_likelihood(const double* y, const double* d, llh_struct& llh, int block_start)
{
  log_nbinom_p_likelihood(y, d, llh, block_start, llh.psi.size() - block_start);
}

void log_nbinom_mu_likelihood(const double* y, const double* d, llh_struct& llh, int block_start, int num_blocks)
{
  int block_end = block_start + num_blocks;

  if (block_end > llh.psi.size()) Rprintf( "Warning: log_logit_likelihood block_end out of bounds.\n");

  for (int i = block_start; i < block_end; i++) {
    double bexpon = y[i] + d[i];
    llh.psi[i]  = llh.psi_stc[i] + llh.psi_dyn[i];
    double psi  = llh.psi[i];
    double epsi = exp(psi);
    double p    = epsi / (d[i] + epsi);
    llh.l0[i]   = psi * y[i] - bexpon * log(d[i] + epsi);
    llh.l1[i]   = y[i] - bexpon * p;
    llh.l2[i]   = -1.0 * bexpon * (p - p*p);
    llh.pl2[i]  = psi * llh.l2[i];
    // Rprintf("i: %2i, psi: %3.2g, y: %3.2g, l0: %3.2g, l1: %3.2g, l2: %3.2g\n", i, psi, y[i], llh.l0[i], llh.l1[i], llh.l2[i]);
  }
}

void log_nbinom_mu_likelihood(const double* y, const double* d, llh_struct& llh, int block_start)
{
  log_nbinom_mu_likelihood(y, d, llh, block_start, llh.psi.size() - block_start);
}

void draw_norm_upper_prec(double* draw_, double* m_, double* U_, const int* n_) 
{
  int n = *n_;

  Map<VectorXd> m_map(m_, n);
  Map<MatrixXd> U_map(U_, n, n);

  Gaussian nbp(n);
  nbp.m = m_map;
  nbp.U = U_map;

  Map<VectorXd> d_map(draw_, n);
  VectorXd d(n);

  RNG r;
  rNorm(d, nbp, r);
  
  d_map = d;
}

void draw_omega(double* omega, double* beta,
		double* psi_dyn, double* psi_stc,
		double* psi, double* l0, double* l1, double* l2, double* pl2,
		double* y, double* tX, double* ntrials, double* offset,
		double* prior_prec, double* phi,
		int* starts,  const int* N_,  const int* B_,  const int* num_starts, 
		int* naccept, bool* just_maximize, int *type)
{
  const int N = *N_;
  const int B = *B_;
  int NB = N * B;

  Map<MatrixXd> omega_m(omega, B, N);
  Map<MatrixXd> beta_m (beta , B, N);

  Map<VectorXd> y_m(y, N);
  Map<MatrixXd> tX_m(tX, B, N);
  Map<VectorXd> ntrials_m(ntrials, N);
  Map<VectorXd> offset_m(offset, N);
   
  MatrixXd Xi(N, NB);
  MatrixXd L (N, NB);

  Map<MatrixXd> prior_prec_m(prior_prec, NB, NB);
  Map<VectorXd> phi_m(phi, B);
  
  MatrixXd Phi(B,B);
  Phi = phi_m.asDiagonal();

  Map<MatrixXd> Phi_m(&Phi(0), B, B);

  make_L (L, phi_m, N);
  make_Xi(Xi, tX_m, L);

  // cout << "Xi:\n" << Xi << "\n";

  Map<MatrixXd> L_m(&L(0), NB, NB);
  Map<MatrixXd> Xi_m(&Xi(0), N, NB);

  Map<Matrix<int, Dynamic, 1> > starts_map(starts, *num_starts, 1);
  MatrixXi starts_mat = starts_map;

  Map<VectorXd> psi_m(psi, N);
  Map<VectorXd> psi_dyn_m(psi_dyn, N);
  Map<VectorXd> psi_stc_m(psi_stc, N);
  Map<VectorXd> l0_m(l0, N);
  Map<VectorXd> l1_m(l1, N);
  Map<VectorXd> l2_m(l2, N);
  Map<VectorXd> pl2_m(pl2, N);

  llh_struct llh(N);
  llh.psi = psi_m;
  llh.psi_dyn = psi_dyn_m;
  llh.psi_stc = psi_stc_m;
  // Or calculate again using psi.
  llh.l0  = l0_m;
  llh.l1  = l1_m;
  llh.l2  = l2_m;
  llh.pl2 = pl2_m;

  RNG r;

  // Overloading confusing ? : expression.
  // log_likelihood function_call = *type == 'L' ? &log_logit_likelihood : &log_nbinom_mu_likelihood;
  log_likelihood function_call = &log_logit_likelihood;
  if (*type == 1) function_call = &log_nbinom_mu_likelihood;

  #ifdef USE_R
  GetRNGstate();
  #endif

  *naccept = 0;
  *naccept = draw_omega(omega_m, beta_m, llh,
			y_m, tX_m, ntrials_m, offset_m,
			Xi_m, L_m, 
			prior_prec_m, Phi_m,
			starts_mat, r, function_call, *just_maximize);

  #ifdef USE_R
  PutRNGstate();
  #endif

  psi_m = llh.psi;
  l0_m  = llh.l0;
  l1_m  = llh.l1;
  l2_m  = llh.l2;
  pl2_m = llh.pl2;
  psi_dyn_m = llh.psi_dyn;
  psi_stc_m = llh.psi_stc;
  
}

void draw_stc_beta(double* beta,
		   double* psi_dyn, double* psi_stc,
		   double* psi, double* l0, double* l1, double* l2, double* pl2,
		   double*y, double* X, double* ntrials, double* offset,
		   double* b0, double* P0,
		   const int* N_, const int* B_,
		   int* naccept, bool* just_maximize, int* type)
{
  int N = *N_;
  int B = *B_;

  Map<VectorXd> beta_m(beta, B);
  Map<VectorXd> psi_m(psi, N);
  Map<VectorXd> psi_dyn_m(psi_dyn, N);
  Map<VectorXd> psi_stc_m(psi_stc, N);
  Map<VectorXd> l0_m(l0, N);
  Map<VectorXd> l1_m(l1, N);
  Map<VectorXd> l2_m(l2, N);
  Map<VectorXd> pl2_m(pl2, N);

  Map<VectorXd> y_m(y, N);
  Map<MatrixXd> X_m(X, N, B);
  Map<VectorXd> ntrials_m(ntrials, N);
  Map<VectorXd> offset_m(offset, N);

  Map<VectorXd> b0_m(b0, B);
  Map<MatrixXd> P0_m(P0, B, B);

  llh_struct llh(N);
  llh.psi = psi_m;
  llh.psi_dyn = psi_dyn_m;
  llh.psi_stc = psi_stc_m;
  llh.l0  = l0_m;
  llh.l1  = l1_m;
  llh.l2  = l2_m;
  llh.pl2 = pl2_m;

  RNG r;

  log_likelihood function_call = &log_logit_likelihood;
  if (*type == 1) function_call = &log_nbinom_mu_likelihood;

  #ifdef USE_R
  GetRNGstate();
  #endif

  *naccept = 0;
  *naccept = draw_stc_beta(beta_m, llh, y_m, X_m, ntrials_m, offset_m, b0_m, P0_m, r, function_call, *just_maximize);

  #ifdef USE_R
  PutRNGstate();
  #endif

  psi_m = llh.psi;
  l0_m  = llh.l0;
  l1_m  = llh.l1;
  l2_m  = llh.l2;
  pl2_m = llh.pl2;
  psi_dyn_m = llh.psi_dyn;
  psi_stc_m = llh.psi_stc;
}

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

  // MatrixXd omega_new(B,N);
  // int right_start = idx_start + num_idxs;
  // if (idx_start > 0) 
  //   omega_new.array().segment(0, idx_start) = omega.array().segment(0, idx_start);
  // if (right_start < NB)
  //   omega_new.array().segment(right_start, NB-right_start) = 
  //     omega.array().segment(right_start, NB-right_start);
