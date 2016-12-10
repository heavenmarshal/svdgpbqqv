#include "bqqv.h"
#include<algorithm>
extern "C"
{
  void loglikProd_R(const double* phi_, const double* vtheta_,
		    const int* levels_, const double* X_,
		    const int* Z_, const double* resp_,
		    const int* n_, const int* d_, const int* q_,
		    double *loglik_)
  {
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;
    mapVec phi(phi_,d);
    mapVec vtheta(vtheta_,lentht);
    mapVeci levels(levels_, q);
    mapMat X(X_,n,d);
    mapMati Z(Z_,n,q);
    mapVec resp(resp_,n);
    double loglik = loglikProd(phi,vtheta,levels,X,Z,resp);
    *loglik_ = loglik;
  }
  void loglikAdd_R(const double* phi_, const double* vtheta_, const double* vsigma2_,
		   const int* levels_, const double* X_, const int* Z_,
		   const double* resp_, const int* n_, const int* d_, const int* q_,
		   double *loglik_)
  {
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;
    mapMat phi(phi_,d,q);
    mapVec vtheta(vtheta_,lentht);
    mapVec vsigma2(vsigma2_, q);
    mapVeci levels(levels_, q);
    mapMat X(X_,n,d);
    mapMati Z(Z_,n,q);
    mapVec resp(resp_,n);
    double loglik = loglikAdd(phi,vtheta,vsigma2,levels,X,Z,resp);
    *loglik_ = loglik;
  }
  void loglikAddHom_R(const double* phi_, const double* vtheta_,
		      const int* levels_, const double* X_,
		      const int* Z_, const double* resp_,
		      const int* n_, const int* d_, const int* q_,
		      double *loglik_)
  {
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;
    mapVec phi(phi_,d);
    mapVec vtheta(vtheta_,lentht);
    mapVeci levels(levels_, q);
    mapMat X(X_,n,d);
    mapMati Z(Z_,n,q);
    mapVec resp(resp_,n);
    double loglik = loglikAddHom(phi,vtheta,levels,X,Z,resp);
    *loglik_ = loglik;
  }

  void prodPredict_R(const double* X0_, const int* Z0_, const double* X_,
		     const int* Z_, const double* resp_, const double* phi_,
		     const double *vtheta_, const int *levels_,
		     const int *n0_, const int* n_, const int *d_, const int* q_,
		     double* pmean_, double* psig2_)
  {
    int n0 = *n0_;
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;

    mapVec phi(phi_,d);
    mapVec vtheta(vtheta_, lentht);
    mapVeci levels(levels_, q);
    mapMat X0(X0_,n0,d);
    mapMat X(X_, n, d);
    mapMati Z0(Z0_, n0, q);
    mapMati Z(Z_, n, q);
    mapVec resp(resp_, n);
    VectorXd pmean, psig2;
    prodPredict(X0, Z0, X, Z, resp, phi, vtheta, levels, pmean, psig2);
    std::copy(pmean.data(), pmean.data()+n0, pmean_);
    std::copy(psig2.data(), psig2.data()+n0, psig2_);
  }
  void addPredict_R(const double* X0_, const int* Z0_, const double* X_,
		    const int* Z_, const double* resp_, const double* phi_,
		    const double* vtheta_, const double* vsigma2_, const int* levels_,
		    const int* n0_, const int* n_, const int* d_, const int* q_,
		    double* pmean_, double* psig2_)
  {
    int n0 = *n0_;
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;

    mapMat phi(phi_,d,q);
    mapVec vtheta(vtheta_, lentht);
    mapVec vsigma2(vsigma2_, q);
    mapVeci levels(levels_, q);
    mapMat X0(X0_,n0,d);
    mapMat X(X_, n, d);
    mapMati Z0(Z0_, n0, q);
    mapMati Z(Z_, n, q);
    mapVec resp(resp_, n);
    VectorXd pmean, psig2;
    addPredict(X0, Z0, X, Z, resp, phi, vtheta, vsigma2, levels, pmean, psig2);
    std::copy(pmean.data(), pmean.data()+n0, pmean_);
    std::copy(psig2.data(), psig2.data()+n0, psig2_);
  }
  void addHomPredict_R(const double* X0_, const int* Z0_, const double* X_,
		       const int* Z_, const double* resp_, const double* phi_,
		       const double *vtheta_, const int *levels_,
		       const int *n0_, const int* n_, const int *d_, const int* q_,
		       double* pmean_, double* psig2_)
  {
    int n0 = *n0_;
    int n = *n_;
    int d = *d_;
    int q = *q_;
    int i, lentht = 0;
    for(i = 0; i != q; ++i)
      lentht += levels_[i]*(levels_[i]-1)/2;

    mapVec phi(phi_,d);
    mapVec vtheta(vtheta_, lentht);
    mapVeci levels(levels_, q);
    mapMat X0(X0_,n0,d);
    mapMat X(X_, n, d);
    mapMati Z0(Z0_, n0, q);
    mapMati Z(Z_, n, q);
    mapVec resp(resp_, n);
    VectorXd pmean, psig2;
    addHomPredict(X0, Z0, X, Z, resp, phi, vtheta, levels, pmean, psig2);
    std::copy(pmean.data(), pmean.data()+n0, pmean_);
    std::copy(psig2.data(), psig2.data()+n0, psig2_);
  }

}
