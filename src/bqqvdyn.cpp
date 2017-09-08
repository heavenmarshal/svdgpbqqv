#include<algorithm>
#include<cmath>
#include<R.h>
#include"bqqv.h"
#define LOG2PI 1.83787706641
//#define PI 3.14159265359
using Eigen::LLT;
typedef LLT<MatrixXd> LLTMd;
void bohmanCorr(const VectorXd& phi, const MatrixXd& X, MatrixXd& corrmat)
{
  int n = X.rows();
  int d = X.cols();
  corrmat = MatrixXd::Identity(n,n);
  double tmp, absval;
  for(int i = 0; i != n; ++i)
  {
    for(int j = 0; j!=i; ++j)
    {
      tmp = 1.0;
      for(int k = 0; k!=d; ++k)
      {
	absval = fabs(X(i,k)-X(j,k));
	if(absval > phi(k)) {
	  tmp = 0.0;
	  break;
	}
	tmp *= (1.0-absval/phi(k))*cos(PI*absval/phi(k))+sin(PI*absval/phi(k))/PI;
      }
      corrmat(i,j) = corrmat(j,i) = tmp;
    }
  }
}
double loglikAddDyn(const MatrixXd& phis, const VectorXd& phit,const VectorXd& vtheta,
		    const VectorXd& vsigma2, const VectorXi& levels, const MatrixXd& X,
		    const MatrixXi& Z, const VectorXd& tv, const MatrixXd& resp)
{
  int n = X.rows();
  int L = tv.size();
  int q = levels.size();
  double dn = (double)n, dL = (double)L;
  MatrixXd corrs, corrt, lmats, lmatt;
  VectorXd oneq = VectorXd::Constant(q,1.0);
  
  addCorr(phis, vtheta, vsigma2, levels, X, Z, corrs);
  bohmanCorr(phit, tv, corrt);

  LLTMd lltcorrs(corrs), lltcorrt(corrt);
  if(lltcorrs.info() != Eigen::Success)
    error("Cholesky failure for spatial correlation matrix!");
  if(lltcorrt.info() != Eigen::Success)
    error("Cholesky failure for temporal correlation matrix!");

  lmats = lltcorrs.matrixL();
  lmatt = lltcorrt.matrixL();
  
  VectorXd onen = VectorXd::Constant(n, 1.0);
  VectorXd oneL = VectorXd::Constant(L, 1.0);
  MatrixXd soly = lltcorrs.solve(lltcorrt.solve(resp).transpose()).transpose();
  VectorXd sol1n = lltcorrs.solve(onen);
  VectorXd sol1L = lltcorrt.solve(oneL);

  double quady = (soly.array()*resp.array()).sum();
  double quad1 = sol1n.sum()*sol1L.sum();
  double bily1 = soly.sum();
  
  double rss = quady - bily1 * bily1 / quad1;
  double logdet = lmats.diagonal().array().log().sum() * dn;
  logdet += lmatt.diagonal().array().log().sum() * dL;

  double loglik = -0.5*(dn*dL*LOG2PI+rss)-logdet;
  return loglik;
}

void addDynPred(const MatrixXd& X0, const MatrixXi& Z0, const MatrixXd& X,
		const MatrixXi& Z,  const VectorXd& tv, const MatrixXd& resp,
		const MatrixXd& phis, const VectorXd& phit, const VectorXd& vtheta,
		const VectorXd& vsigma2, const VectorXi& levels,
		MatrixXd& pmean)
{
  int n = X.rows();
  int q = levels.size();
  int L = resp.rows();
  MatrixXd corrs, corrt, crosscorrs;
  VectorXd oneq = VectorXd::Constant(q,1.0);
  addCorr(phis, vtheta, vsigma2, levels, X, Z, corrs);
  bohmanCorr(phit, tv, corrt);
  
  LLTMd lltcorrs(corrs), lltcorrt(corrt);
  if(lltcorrs.info() != Eigen::Success)
    error("Cholesky failure for spatial correlation matrix!");
  if(lltcorrt.info() != Eigen::Success)
    error("Cholesky failure for temporal correlation matrix!");

  VectorXd onen = VectorXd::Constant(n, 1.0);
  VectorXd oneL = VectorXd::Constant(L, 1.0);
  MatrixXd soly = lltcorrs.solve(resp.transpose()).transpose();
  VectorXd sol1n = lltcorrs.solve(onen);
  VectorXd sol1L = lltcorrt.solve(oneL);
  
  double bily1 = sol1L.dot(resp*sol1n);
  double quad1 = sol1n.sum()*sol1L.sum();
  double muhat = bily1/quad1;

  addCrossCorr(phis, vtheta, vsigma2, levels, X, X0, Z, Z0, crosscorrs);
  pmean = soly - sol1n.transpose().replicate(L,1)*muhat;
  pmean = pmean * crosscorrs;
  pmean = muhat + pmean.array();
}

extern "C"
{
  void loglikAddDyn_R(const double* phis_, const double* phit_, const double* vtheta_,
		      const double* vsigma2_, const int* levels_, const double* X_,
		      const int* Z_, const double* tv_, const double* resp_,
		      const int* n_, const int* L_, const int* d_, const int* q_,
		      double* loglik_)
  {
    int n = *n_, L = *L_, d = *d_, q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;
    mapMat phis(phis_,d,q);
    mapVec phit(phit_,1);
    mapVec vtheta(vtheta_,lentht);
    mapVec vsigma2(vsigma2_, q);
    mapVeci levels(levels_, q);
    mapMat X(X_,n,d);
    mapMati Z(Z_,n,q);
    mapVec tv(tv_,L);
    mapMat resp(resp_,L,n);
    double loglik = loglikAddDyn(phis, phit, vtheta, vsigma2, levels,
				 X, Z, tv, resp);
    *loglik_ = loglik;
  }
  void addDynPred_R(const double* X0_, const int* Z0_, const double* X_,
		    const int* Z_, const double* tv_, const double* resp_,
		    const double* phis_, const double* phit_, const double* vtheta_,
		    const double* vsigma2_, const int* levels_, const int* n0_,
		    const int* n_, const int* L_, const int* d_, const int* q_,
		    double* pmean_)
  {
    int n0 = *n0_, n = *n_, L = *L_, d = *d_, q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;
    mapMat X0(X0_, n0, d);
    mapMati Z0(Z0_, n0, q);
    mapMat X(X_, n, d);
    mapMati Z(Z_, n, q);
    mapVec tv(tv_, L);
    mapMat resp(resp_, L, n);
    mapMat phis(phis_, d, q);
    mapVec phit(phit_, 1);
    mapVec vtheta(vtheta_, lentht);
    mapVec vsigma2(vsigma2_, q);
    mapVeci levels(levels_, q);

    MatrixXd pmean;
    addDynPred(X0, Z0, X, Z, tv, resp, phis, phit,
	       vtheta, vsigma2, levels, pmean);
    std::copy(pmean.data(), pmean.data()+n0*L, pmean_);
  }
}
  
