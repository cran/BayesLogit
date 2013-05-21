
#ifndef __DYNEXPFAMMH__
#define __DYNEXPFAMMH__

#include <Eigen/Dense>
#include <stdio.h>
#include "RNG.h"
#include <vector>

#ifdef USE_R
#include <R_ext/Utils.h>
#endif

using std::vector;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::MatrixBase;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Map;
using Eigen::MatrixBase;

// using std::cout;
// using std::cerr;

class llh_view {
public:
  double* psi_dyn;
  double* psi_stc;
  double* psi;
  double* l0;
  double* l1;
  double* l2;
  double* pl2;
  unsigned int N;

  llh_view(const llh_view& l)
  : psi_dyn(l.psi_dyn)
  , psi_stc(l.psi_stc)
  , psi(l.psi)
  , l0(l.l0)
  , l1(l.l1)
  , l2(l.l2)
  , pl2(l.pl2)
  , N(l.N) {};

  llh_view(double* psi_dyn_,
		   double* psi_stc_,
		   double* psi_,
		   double* l0_,
		   double* l1_,
		   double* l2_,
		   double* pl2_,
		   unsigned int N_)
  : psi_dyn(psi_dyn_)
  , psi_stc(psi_stc_)
  , psi(psi_)
  , l0(l0_)
  , l1(l1_)
  , l2(l2_)
  , pl2(pl2_)
  , N(N_) {};

  void copy(llh_view& l) {
	N = l.N;
	size_t amount = N * sizeof(double);
	memcpy(psi_dyn, l.psi_dyn, amount);
	memcpy(psi_stc, l.psi_stc, amount);
	memcpy(psi    , l.psi    , amount);
	memcpy(l0     , l.l0     , amount);
	memcpy(l1     , l.l1     , amount);
	memcpy(l2     , l.l2     , amount);
	memcpy(pl2    , l.pl2    , amount);
  }

  void copy(llh_view& l, unsigned int start) {
	N = l.N;
	unsigned int num = N - start;
	size_t amount = num * sizeof(double);
	memcpy(psi_dyn + start, l.psi_dyn + start, amount);
	memcpy(psi_stc + start, l.psi_stc + start, amount);
	memcpy(psi + start    , l.psi + start   , amount);
	memcpy(l0 + start     , l.l0 + start     , amount);
	memcpy(l1 + start     , l.l1 + start     , amount);
	memcpy(l2 + start     , l.l2 + start     , amount);
	memcpy(pl2 + start    , l.pl2 + start    , amount);
  }
  
};

class llh_struct {
public:
  VectorXd psi;
  VectorXd psi_stc;
  VectorXd psi_dyn;
  VectorXd l0;
  VectorXd l1;
  VectorXd l2;
  VectorXd pl2;

  llh_struct();

  llh_struct(int size)
  : psi(size)
  , psi_stc(size)
  , psi_dyn(size)
  , l0(size)
  , l1(size)
  , l2(size)
  , pl2(size) {};

  llh_struct(const llh_struct& l)
  : psi(l.psi)
  , psi_stc(l.psi_stc)
  , psi_dyn(l.psi_dyn)
  , l0(l.l0)
  , l1(l.l1)
  , l2(l.l2)
  , pl2(l.pl2) {};

  // void resize(int size) {
  //   psi.resize(size);
  //   l0.resize(size);
  //   l1.resize(size);
  //   l2.resize(size);
  //   pl2.resize(size);
  // }

};

class Gaussian {
public:
  VectorXd m;
  MatrixXd P;
  MatrixXd U;

  Gaussian();

  Gaussian(int dim) 
  : m(dim)
  , P(dim, dim)
  , U(dim, dim) {};

  // void resize(int dim) {
  //   m.resize(dim);
  //   P.resize(dim, dim);
  //   L.resize(dim, dim);
  // }

};

void bidc_to_midc(int* mstart, int* mnum, const int bstart, const int bnum, const int sob);

void make_A(MatrixXd& A, const VectorXd& phi, const unsigned int num);
void make_A(MatrixXd& A, const MatrixXd& G, const unsigned int num);

// void make_L(MatrixXd& L, const VectorXd& phi, const unsigned int num);

// void make_Xi(MatrixXd& XiL, const MatrixXd& tX, const MatrixXd& L);

// double rNorm(VectorXd& d, const Gaussian& nbp, RNG& r);

// double dNorm(const VectorXd& x, const Gaussian& nbp, bool use_log=false);

// void log_logit_likelihood(const VectorXd& y, const int n, llh_struct& llh, const int block_start=0);

// void log_logit_likelihood(const VectorXd& y, const int n, llh_struct& llh, 
//			  const int block_start, const int block_end);

//void omega_to_beta(MatrixXd& beta, const MatrixXd& Phi, const MatrixXd& omega, const int block_start);

// void laplace_omega(Gaussian& nbp, const VectorXd& omega, const llh_struct& llh,
// 		   const VectorXd& y, const MatrixXd& Xi, const MatrixXd& prior_prec,
// 		   const int block_start, const int num_blocks);

// bool draw_omega_block(MatrixXd& omega, MatrixXd& beta, llh_struct& llh,
// 		      VectorXd& y, MatrixXd& tX, int ntrials, 
// 		      MatrixXd& Xi, MatrixXd& L,
// 		      MatrixXd& prior_prec, MatrixXd& Phi,
// 		      int block_start, int num_blocks, 
// 		      RNG& r);

void log_logit_likelihood(const double* y, const double* ntrials, llh_view& llh, int total_blocks, int block_start, int num_blocks);
void log_logit_likelihood(const double* y, const double* ntrials, llh_view& llh, int block_start);

void log_logit_likelihood(const double* y, const double* ntrials, llh_struct& llh, int total_blocks, int block_start, int num_blocks);
void log_logit_likelihood(const double* y, const double* ntrials, llh_struct& llh, int block_start);

void log_nbinom_likelihood(const double* y, const double* ntrials, llh_struct& llh, int total_blocks, int block_start, int num_blocks);
void log_nbinom_likelihood(const double* y, const double* ntrials, llh_struct& llh, int block_start);

typedef void (*log_likelihood)(const double* y, const double* p1, llh_struct& llh, int block_start);

extern "C" {

  void draw_omega(double* omega, double* beta,
		  double* psi_dyn, double* psi_stc,
		  double* psi, double* l0, double* l1, double* l2, double* pl2,
		  double* y, double* tX, double* ntrials, double* offset,
		  double* prior_prec, double* Phi,
		  int* starts,  const int* N_,  const int* B_,  const int* num_starts, 
		  int* naccept, bool* just_maximize, int *type);

  void draw_stc_beta(double* beta,
		     double* psi_dyn, double* psi_stc,
		     double* psi, double* l0, double* l1, double* l2, double* pl2,
		     double*y, double* X, double* ntrials, double* offset, 
		     double* b0, double* P0,
		     const int* N_, const int* B_,
		     int* naccept, bool* just_maximize, int *type);

}



//------------------------------------------------------------------------------

template<typename dV>
double dNorm(const MatrixBase<dV>& x, const Gaussian& nbp, bool use_log) 
{
  int size = nbp.m.size();
  VectorXd d(size);

  d = nbp.U * (x - nbp.m);

  double ldens = size * log(SQRT2PI);
  ldens += nbp.U.diagonal().array().log().sum();
  ldens += -0.5 * d.dot(d);

  ldens = use_log ? ldens : exp(ldens);
  return ldens;
}

template <typename dV>
double rNorm(MatrixBase<dV>& d, const Gaussian& nbp, RNG& r)
{
  int size = nbp.m.size();
  d.resize(size);
  r.norm(d, 0, 1);

  // cout << "ep: " << d.transpose() << "\n";

  double ldens = size * log(SQRT2PI);
  ldens += nbp.U.diagonal().array().log().sum();
  ldens += -0.5 * d.dot(d);

  // nbp.P.inverse().llt().matrixL().solveInPlace(d);

  nbp.U.triangularView<Eigen::Upper>().solveInPlace(d);
  d += nbp.m;

  // cout << "m: " << nbp.m.transpose() << "\n";
  // cout << "d: " << d.transpose() << "\n";
  
  return ldens;
}

// template<typename dV1>
// void log_logit_likelihood(const MatrixBase<dV1>& y,
// 			  const MatrixBase<dV1>& ntrials,
// 			  llh_struct& llh,
// 			  const int block_start,
// 			  const int num_blocks)
// {
//   // Assume psi_stc and psi_dyn are set correct.
//   // int nsize = y.size();
  
//   int block_end = block_start + num_blocks;

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

// template<typename dV1>
// void log_logit_likelihood(const MatrixBase<dV1>& y,
// 			  const MatrixBase<dV1>& ntrials,
// 			  llh_struct& llh,
// 			  const int block_start)
// {
//   int nsize = y.size();
//   log_logit_likelihood(y, ntrials, llh, block_start, nsize - block_start);
// }

//------------------------------------------------------------------------------

template<typename dV2, typename dM1, typename dM2>
void dyn_beta_to_psi(double* psi_, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, const MatrixBase<dV2>& offset,
		     const int block_start, const int num_blocks)
{
  int B = beta.rows();
  int N = beta.cols();
  // std::cout << "offset:\n" << offset.transpose() << "\n";

  Map<VectorXd> psi(psi_, N);

  MatrixXd tXbeta_vec;
  tXbeta_vec = tX.block(0, block_start, B, num_blocks).array() *
    beta.block(0, block_start, B, num_blocks).array();

  Map<MatrixXd> tXbeta(&tXbeta_vec(0), B, num_blocks);

  // cout << "tXbeta:\n" << tXbeta << "\n";
  
  psi.segment(block_start, num_blocks) = tXbeta.colwise().sum();
  psi.segment(block_start, num_blocks) += offset.segment(block_start, num_blocks);

}

template<typename dV2, typename dM1, typename dM2>
void dyn_beta_to_psi(double* psi, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, const MatrixBase<dV2>& offset,
		     const int block_start)
{
  int N = tX.cols();

  int num_blocks = N - block_start;

  dyn_beta_to_psi(psi, tX, beta, offset, block_start, num_blocks);
}

template<typename dV1, typename dV2, typename dM1, typename dM2>
void dyn_beta_to_psi(MatrixBase<dV1>& psi, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, const MatrixBase<dV2>& offset,
		     const int block_start, const int num_blocks)
{
  int B = beta.rows();
  // int N = beta.cols();
  // std::cout << "offset:\n" << offset.transpose() << "\n";

  MatrixXd tXbeta_vec;
  tXbeta_vec = tX.block(0, block_start, B, num_blocks).array() *
    beta.block(0, block_start, B, num_blocks).array();

  Map<MatrixXd> tXbeta(&tXbeta_vec(0), B, num_blocks);

  // cout << "tXbeta:\n" << tXbeta << "\n";
  
  psi.segment(block_start, num_blocks) = tXbeta.colwise().sum();
  psi.segment(block_start, num_blocks) += offset.segment(block_start, num_blocks);

}

template<typename dV1, typename dV2, typename dM1, typename dM2>
void dyn_beta_to_psi(MatrixBase<dV1>& psi, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, const MatrixBase<dV2>& offset,
		     const int block_start)
{
  int N = tX.cols();

  int num_blocks = N - block_start;

  dyn_beta_to_psi(psi, tX, beta, offset, block_start, num_blocks);
}

template<typename dV, typename dM1, typename dM2>
void dyn_beta_to_psi(MatrixBase<dV>& psi, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, 
		     const int block_start, const int num_blocks)
{
  int B = tX.rows();
  // int N = tX.cols();

  MatrixXd tXbeta_vec;
  tXbeta_vec = tX.block(0, block_start, B, num_blocks).array() *
    beta.block(0, block_start, B, num_blocks).array();

  Map<MatrixXd> tXbeta(&tXbeta_vec(0), B, num_blocks);

  // cout << "tXbeta:\n" << tXbeta << "\n";

  psi.segment(block_start, num_blocks) = tXbeta.colwise().sum();
}

template<typename dV, typename dM1, typename dM2>
void dyn_beta_to_psi(MatrixBase<dV>& psi, 
		     const MatrixBase<dM1>& tX, const MatrixBase<dM2>& beta, const int block_start)
{
  int N = tX.cols();
  
  int num_blocks = N - block_start;

  dyn_beta_to_psi(psi, tX, beta, block_start, num_blocks);
}

////////////////////////////////////////////////////////////////////////////////

template <typename dV1, typename dV2, typename dM1>
void laplace_stc_beta(Gaussian& nbp, const MatrixBase<dV1>& beta, const llh_struct& llh,
		      const MatrixBase<dV2>& y, const MatrixBase<dM1>& X,
		      const MatrixBase<dV2>& b0, const MatrixBase<dM1>& P0)
{
  // int N = X.rows();
  // int B = X.cols();

  MatrixXd Hl = X.transpose() * (llh.l2.asDiagonal() * X);
  VectorXd bl = X.transpose() * llh.l1 - Hl * beta;
  // VectorXd bl = X.transpose() * (llh.l1 - llh.psi_stc * llh.l2);

  // std::cout << "beta:\n" << beta << "\n";
  // std::cout << "bl:\n" << bl << "\n";
  // std::cout << "Pl:\n" << Pl << "\n";

  nbp.P = P0 - Hl;
  nbp.U = nbp.P.llt().matrixU();

  nbp.m = b0 + bl;
  nbp.P.llt().solveInPlace(nbp.m);

}

template <typename dM, typename dV1, typename dV2>
bool draw_stc_beta(MatrixBase<dV2>& beta, llh_struct& llh,
		   MatrixBase<dV1>& y, MatrixBase<dM>& X, MatrixBase<dV1>& ntrials, 
		   MatrixBase<dV1>& offset, 
		   MatrixBase<dV1>& b0, MatrixBase<dM>& P0,
		   RNG& r, log_likelihood log_like, bool just_maximize=false)
{

  // std::cout << "y:\n" << y.transpose() << "\n";
  // std::cout << "X:\n" << X << "\n";
  // std::cout << "b0:\n" << b0.transpose() << "\n";
  // std::cout << "P0:\n" << P0 << "\n";
  // std::cout << "offset:\n" << offset.transpose() << "\n";
  // std::cout << "ntrials:\n" << ntrials.transpose() << "\n";

  // int N = X.rows();
  int B = X.cols();

  int num_idxs    = B;
  int block_start = 0;

  // -- New --

  Gaussian nbp_new(num_idxs);
  VectorXd beta_new(num_idxs);
  double   lppsl_new = 0.0;

  laplace_stc_beta(nbp_new, beta, llh, y, X, b0, P0);
  
  if (just_maximize) {
    beta_new  = nbp_new.m;
    lppsl_new = dNorm(beta_new, nbp_new, true);
  } else {
    lppsl_new = rNorm(beta_new, nbp_new, r);
  }

  llh_struct llh_new(llh);
  llh_new.psi_stc = X * beta_new;
  llh_new.psi_dyn = offset;
  // log_logit_likelihood(y, ntrials, llh_new, block_start);
  // log_logit_likelihood(&y(0), &ntrials(0), llh_new, block_start);
  log_like(&y(0), &ntrials(0), llh_new, block_start);

  double llike_new  = llh_new.l0.sum();

  double lprior_new = (-0.5) * (beta_new.transpose() * (P0 * beta_new))(0) + beta_new.transpose() * b0;

  // -- Old --

  Gaussian nbp_old(num_idxs);
  Map<VectorXd> beta_old(&beta(0), num_idxs);

  laplace_stc_beta(nbp_old, beta_new, llh_new, y, X, b0, P0);

  double lppsl_old = dNorm(beta_old, nbp_old, true);
  double llike_old = llh.l0.sum();
  double lprior_old = -0.5 * (beta_old.transpose() * (P0 * beta_old))(0) + beta_old.transpose() * b0;

  double lratio = (llike_new + lprior_new) - (llike_old + lprior_old) - (lppsl_new - lppsl_old);
  
  // Rprintf("llike_new: %g, lprior_new: %g, lppsl_new %g\n", llike_new, lprior_new, lppsl_new);
  // Rprintf("llike_old: %g, lprior_old: %g, lppsl_old %g\n", llike_old, lprior_old, lppsl_old);
  // Rprintf("lratio: %g\n", lratio);

  bool accept = log(r.unif()) < lratio;
  accept += just_maximize;

  if (accept) {
    beta   = beta_new;
    // llh    = llh_new;
    llh.psi = llh_new.psi;
    llh.l0  = llh_new.l0;
    llh.l1  = llh_new.l1;
    llh.l2  = llh_new.l2;
    llh.psi_dyn = llh_new.psi_dyn;
    llh.psi_stc = llh_new.psi_stc;
    llh.pl2 = llh_new.pl2;
  }

  return accept;
}

////////////////////////////////////////////////////////////////////////////////

template<typename dM1, typename dM2, typename dM3>
void omega_to_dyn_beta(MatrixBase<dM1>& beta, 
		       MatrixBase<dM2>& Phi, 
		       MatrixBase<dM3>& omega, 
		       const int block_start)
{
  int c = beta.cols();

  int start = block_start;

  if (block_start == 0) {
    beta.col(0) = omega.col(0);
    start++;
  }

  for (int i=start; i < c; i++) {
    beta.col(i) = Phi * beta.col(i-1) + omega.col(i);
  }

}

template <typename dV>
void make_L(MatrixXd& L, const MatrixBase<dV>& phi, const unsigned int num)
{
  unsigned int bsize = phi.size();
  unsigned int N = num * bsize;

  L.resize(N, N); L.fill(0.0);

  //unsigned int initial_offset = 0;
  //unsigned int nrep           = num;
  
  VectorXd phi_n(bsize);

  for (unsigned int i = 0; i < num; i++) {

    unsigned int offset = i*bsize;
    unsigned int J      = N-offset;

    phi_n = phi.array().pow((double) i);

    for (unsigned int j = 0; j < J; j++) {
      L(j+offset,j) = phi_n(j % bsize);
    }

  }
  
}

template <typename dM>
void make_Xi(MatrixXd& Xi, const MatrixBase<dM>& tX, const MatrixXd& L)
{
  unsigned int B  = tX.rows();
  unsigned int N  = tX.cols();
  unsigned int BN = L.rows();

  // Xi.resize(BN, N);
  Xi.resize(N, BN);
  Xi.fill(0.0);

  for (unsigned int j=0; j < N; j++) {
    for (unsigned int i=j; i < N; i++) {
      // Xi.block(j*B, i, B, 1) = L.block(i*B, j*B, B, B).transpose() * tX.col(j);
      // cout << tX.col(j).transpose() << "\n";
      Xi.block(i, j*B, 1, B) = tX.col(j).transpose() * L.block(i*B, j*B, B, B);
    }
  }
  
}

template <typename dV1, typename dV2, typename dM1>
void laplace_omega(Gaussian& nbp, const MatrixBase<dV1>& omega, const llh_struct& llh,
		   const MatrixBase<dV2>& y, const MatrixBase<dM1>& Xi, const MatrixBase<dM1>& prior_prec,
		   const int block_start, const int num_blocks)
{
  int N         = Xi.rows();
  int NB        = Xi.cols();
  int bsize     = NB / N;
  int idx_start = block_start * bsize;
  int num_idxs  = num_blocks * bsize;

  // std::cout << "N: " << N << ", B: " << bsize << ", NB:" << NB 
  //	    << ", idx_start: " << idx_start << ", num_idxs: " << num_idxs << "\n";

  // This doesn't work.
  // MatrixXd XiSub(Xi.block(block_start, idx_start, num_blocks, num_idxs));
  // VectorXd Gtil = XiSub.transpose() * llh.l1.segment(block_start, num_blocks);
  // VectorXd Ftil = XiSub.transpose() * llh.pl2.segment(block_start, num_blocks);

  MatrixXd XiSub(Xi.block(0, idx_start, N, num_idxs));
  VectorXd Gtil = XiSub.transpose() * llh.l1;
  VectorXd pdl2 = llh.psi_dyn.array() * llh.l2.array();
  VectorXd Ftil = XiSub.transpose() * (pdl2);

  // std::cout << "pl2:\n" << llh.pl2.segment(block_start, num_blocks).transpose() << "\n";
  // std::cout << "pdl2:\n" << pdl2.transpose() << "\n";

  // std::cout << "Gtil:\n" << Gtil.transpose() << "\n";
  // std::cout << "Gtil2:\n" << Gtil2.transpose() << "\n";

  // std::cout << "Ftil:\n" << Ftil.transpose() << "\n";
  // std::cout << "Ftil2:\n" << Ftil2.transpose() << "\n";

  // There may be a more clever way to do this.
  // MatrixXd blah(Xi); blah.array() *= llh.l2.array();
  // MatrixXd Htil1_ = Xi.block(0, idx_start, N, num_idxs).transpose() * (blah);
  MatrixXd Htil1_ = Xi.block(0, idx_start, N, num_idxs).transpose() * (llh.l2.asDiagonal() * Xi);

  // O11
  nbp.P  = prior_prec.block(idx_start, idx_start, num_idxs, num_idxs);
  nbp.P += -1.0 * Htil1_.block(0, idx_start, num_idxs, num_idxs);

  MatrixXd O12;
  VectorXd O12w(num_idxs); O12w.fill(0.0);

  // If we are not using the whole thing.
  int O12_cols = NB - num_idxs;

  if (O12_cols > 0) {

    O12.resize(num_idxs, O12_cols);
    VectorXd om2(O12_cols);
    
    int num_left  = idx_start;
    int num_right = NB - (idx_start + num_idxs);
    
    // cout << "num_left: " << num_left << ", num_right: " << num_right << "\n";

    if (block_start > 0) {
      O12.block(0, 0, num_idxs, num_left) = -1.0 * Htil1_.block(0, 0, num_idxs, num_left);
      om2.segment(0, num_left) = omega.segment(0, num_left);
    }
    
    if ( num_right > 0 ) {
      int right_start = NB - num_right;
      O12.block(0, num_left, num_idxs, num_right) = -1.0 * Htil1_.block(0, right_start, num_idxs, num_right); 
      om2.segment(num_left, num_right) = omega.segment(right_start, num_right);
    }

    O12w = O12 * om2;  

    // std::cout << "O12:\n" << O12 << "\n";
    // std::cout << "om2: " << om2.transpose() << "\n";
    // std::cout << "O12w: " << O12w.transpose() << "\n";
  }

  // std::cout << "Htil1_:\n" << Htil1_ << "\n";
  // std::cout << "O12:\n" << O12 << "\n";
  // std::cout << "O11:\n" << O11 << "\n";

  VectorXd btil  = Gtil - Ftil;
  // cout << "btil: " << btil.transpose() << "\n";
  
  // Solve for m.
  
  // nbp.P = O11;
  nbp.U = nbp.P.llt().matrixU();
  nbp.m = btil - O12w;
  
  //nbp.U.triangularView<Eigen::Upper>().solveInPlace(nbp.m);
  //nbp.U.transpose().triangularView<Eigen::Lower>().solveInPlace(nbp.m);
  nbp.P.llt().solveInPlace(nbp.m);

  // std::cout << "btil:\n" << btil << "\n";
  // std::cout << "O12w:\n" << O12w << "\n";

  // cout << "P:\n" << nbp.P << "\n";
  // cout << "U:\n" << nbp.U << "\n";
  // std::cout << "m:\n" << nbp.m.transpose() << "\n";

}

template <typename dM, typename dV>
bool draw_omega_block(MatrixBase<dM>& omega, MatrixBase<dM>& beta, llh_struct& llh,
					  MatrixBase<dV>& y, MatrixBase<dM>& tX, MatrixBase<dV>& ntrials, MatrixBase<dV>& offset,
					  MatrixBase<dM>& Xi, MatrixBase<dM>& L,
					  MatrixBase<dM>& prior_prec, MatrixBase<dM>& Phi,
					  int block_start, int num_blocks, 
					  RNG& r, log_likelihood log_like, bool just_maximize=false)
{
  int N  = Xi.rows();
  int NB = Xi.cols();
  int B  = NB / N;

  int idx_start = block_start * B;
  int num_idxs  = num_blocks * B;

  // -- New --

  Gaussian   nbp_new(num_idxs);
  VectorXd draw_new(num_idxs);
  double   lppsl_new = 0.0;

  Map<VectorXd> omega_vec(&omega(0), NB);

  laplace_omega(nbp_new, omega_vec, llh, y, Xi, prior_prec, block_start, num_blocks);
  
  if (just_maximize) {
    draw_new  = nbp_new.m;
    lppsl_new = dNorm(draw_new, nbp_new, true);
  } else {
    lppsl_new = rNorm(draw_new, nbp_new, r);
  }

  MatrixXd omega_new(omega);
  Map<VectorXd> omega_new_vec(&omega_new(0), NB);

  omega_new_vec.segment(idx_start, num_idxs) = draw_new;
  // omega_new.block(0, block_start, B, num_blocks).array() = draw_new;

  MatrixXd beta_new(beta);
  omega_to_dyn_beta(beta_new, Phi, omega_new, block_start);

  llh_struct llh_new(llh);

  // dyn_beta_to_psi(llh_new.psi_dyn, tX, beta_new, block_start);
  // llh_new.psi_stc = offset;
  // log_logit_likelihood(y, ntrials, llh_new, block_start);

  dyn_beta_to_psi(llh_new.psi_dyn, tX, beta_new, block_start);
  llh_new.psi_stc = offset;
  // log_logit_likelihood(&y(0), &ntrials(0), llh_new, block_start);
  log_like(&y(0), &ntrials(0), llh_new, block_start);

  double llike_new  = llh_new.l0.sum();

  double lprior_new = draw_new.transpose() * 
    prior_prec.block(idx_start, idx_start, num_idxs, num_idxs) * draw_new;
  lprior_new *= -0.5;

  // -- Old --

  Gaussian nbp_old(num_idxs);
  Map<VectorXd> draw_old(&omega(0, block_start), num_idxs);

  laplace_omega(nbp_old, omega_new_vec, llh_new, y, Xi, prior_prec, block_start, num_blocks);

  double lppsl_old = dNorm(draw_old, nbp_old, true);
  double llike_old = llh.l0.sum();
  double lprior_old = draw_old.transpose() * 
    prior_prec.block(idx_start, idx_start, num_idxs, num_idxs) * draw_old;
  lprior_old *= -0.5;

  double lratio = (llike_new + lprior_new) - (llike_old + lprior_old) - (lppsl_new - lppsl_old);
  
  bool accept = log(r.unif()) < lratio;
  accept += just_maximize;

  // cout << "\nblock_start: " << block_start << "\n";
  // cout << "\naccept: " << accept << "\n";
  // cout << "omega:\n" << omega << "\n";
  // cout << "omega_new:\n" << omega_new << "\n";
  
  // cout << "beta:\n" << beta << "\n";
  // cout << "beta new:\n" << beta_new << "\n";

  // cout << "llh.psi:\n" << llh.psi.transpose() << "\n";
  // cout << "llh_new.psi:\n" << llh_new.psi.transpose() << "\n"; 

  // Rprintf("llike_new: %g, llike_old: %g, diff: %g\n", llike_new, llike_old, llike_new - llike_old);
  // Rprintf("lpp_new: %g, lpp_old: %g, diff: %g\n", lppsl_new, lppsl_old, lppsl_new - lppsl_old);
  // Rprintf("lprior_new: %g, lprior_old: %g, diff: %g\n", lprior_new, lprior_old, lprior_new - lprior_old);
  // Rprintf("block_start: %i, lratio: %g\n", block_start, lratio);

  if (accept) {
    omega  = omega_new;
    beta   = beta_new;
    // llh    = llh_new;
    llh.psi = llh_new.psi;
    llh.l0  = llh_new.l0;
    llh.l1  = llh_new.l1;
    llh.l2  = llh_new.l2;
    llh.pl2 = llh_new.pl2;
    llh.psi_dyn = llh_new.psi_dyn;
    llh.psi_stc = llh_new.psi_stc;
  }

  return accept;
}

template<typename dM, typename dV>
int draw_omega(MatrixBase<dM>& omega, MatrixBase<dM>& beta, llh_struct& llh,
			   MatrixBase<dV>& y, MatrixBase<dM>& tX, MatrixBase<dV>& ntrials, MatrixBase<dV>& offset,
			   MatrixBase<dM>& Xi, MatrixBase<dM>& L,
			   MatrixBase<dM>& prior_prec, MatrixBase<dM>& Phi,
			   MatrixXi& starts,
			   RNG& r, log_likelihood log_like, bool just_maximize=false)
{
  int nsize = y.size();
  int ssize = starts.size();
  
  // This doesn't work.
  // Vector<int> the_breaks(ssize+1);
  // the_breaks.segment(0, ssize);

  MatrixXi the_breaks(ssize+1, 1);  

  the_breaks.block(0, 0, ssize, 1) = starts;
  the_breaks(ssize) = nsize;

  int naccept = 0;

  // std::cout << "omega:\n" << omega << "\n";
  // std::cout << "y:\n" << y.transpose() << "\n";
  // std::cout << "tX:\n" << tX << "\n";
  // std::cout << "Xi:\n" << Xi << "\n";
  // std::cout << "L:\n" << L << "\n";
  // std::cout << "prior.prec:\n" << prior_prec << "\n";
  // std::cout << "Phi:\n" << Phi << "\n";

  // std::cout << "psi:\n" << llh.psi.transpose() << "\n";
  // std::cout << "l0:\n"  << llh.l0.transpose()  << "\n";
  // std::cout << "l1:\n"  << llh.l1.transpose()  << "\n";
  // std::cout << "l2:\n"  << llh.l2.transpose()  << "\n";
  // std::cout << "pl2:\n" << llh.pl2.transpose() << "\n";

  for (int i=0; i<ssize; i++) {
    int num_blocks = the_breaks(i+1) - the_breaks(i);
    naccept +=
      draw_omega_block(omega, beta, llh, y, tX, ntrials, 
		       offset, Xi, L, prior_prec, Phi, starts(i), num_blocks, 
		       r, log_like, just_maximize);
  }

  return naccept;
}

#endif

// TODO

// Add offset to draw omega.
// Create wrapper for static beta draw.
// Create block dynamic beta draw? 
//   - would need to get a new laplace approximation.

// APPENDIX

// Having too many or too few template parameters can create problems.

// I had problems when trying to pass a character to an R wrapping function!
