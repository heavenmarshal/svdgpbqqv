#include"bqqv.h"
#include<cmath>
#include<R.h>
#include<R_ext/Lapack.h>
#define SQ(x) ((x)*(x))
#define LOG2PI 1.83787706641

using Eigen::LLT;
typedef LLT<MatrixXd> LLTMd;
typedef vector<MatrixXd>::const_iterator vmciter;
// note that the levels of the qualitative factor start from 0 to nlevel - 1;

// the input mat must be a symetric matrix
void nugDetermine(const MatrixXd& mat, double& cond, double& anorm)
{
  int n = mat.cols();
  int info;
  char typeNorm[] = {'1','\0'};
  double val;
  double *matcp = new double[n*n];
  double *work = new double[4*n];
  int *iwork = new int[n];
  
  std::copy(mat.data(), mat.data()+n*n,matcp);
  anorm = F77_CALL(dlange)(typeNorm,&n,&n,matcp,&n,work);
  F77_CALL(dgetrf)(&n,&n,matcp,&n, iwork, &info);
  F77_CALL(dgecon)(typeNorm, &n, matcp, &n, &anorm, &val, work, iwork, &info);
  cond = 1.0/val;
  delete matcp;
  delete work;
  delete iwork;
}
void qualCorr(const VectorXd& theta, int nlevel, MatrixXd& corrmat)
{
  MatrixXd umat = MatrixXd::Zero(nlevel,nlevel);
  umat(0,0) = 1.0;
  double tmp, ctheta;
  int idx = 0;
  for(int i = 1; i != nlevel; ++i)
  {
    ctheta = theta(idx++);
    tmp = cos(ctheta);
    umat(0,i) = tmp;
    for(int j=1; j != i; ++j)
    {
      tmp *= tan(ctheta);
      ctheta = theta(idx++);
      tmp *= cos(ctheta);
      umat(j,i) = tmp;
    }
    umat(i,i) = tmp*tan(ctheta);
  }
  corrmat = umat.transpose() * umat;
}
void vecQualCorr(const VectorXd& vtheta, const VectorXi& levels,
		 vector<MatrixXd>& vcorr)
{
  int nfrac = levels.size();
  mapVec mtheta(NULL,0);
  MatrixXd qualcorr;
  int beg = 0, step;
  int nlevel;
  for(int i=0; i!=nfrac; ++i)
  {
    nlevel = levels(i);
    step = nlevel*(nlevel-1)/2;
    new (&mtheta) mapVec(vtheta.data()+beg, step);
    qualCorr(mtheta, nlevel, qualcorr);
    vcorr.push_back(qualcorr);
    beg += step;
  }
}

void quanCorr(const VectorXd& phi, const MatrixXd& X, MatrixXd& corrmat)
{
  int n = X.rows();
  int d = X.cols();
  corrmat = MatrixXd::Identity(n,n);
  double tmp;
  for(int i = 0; i != n; ++i)
  {
    for(int j=0; j!=i; ++j)
    {
      tmp = 0.0;
      for(int k=0; k!=d; ++k)
	tmp -= SQ(X(i,k)-X(j,k))/phi(k);
      corrmat(i,j) = corrmat(j,i) = exp(tmp);
    }
  }
}
void quanCrossCorr(const VectorXd& phi, const MatrixXd& X1, const MatrixXd& X2,
		   MatrixXd& corrmat)
{
  int n1 = X1.rows();
  int n2 = X2.rows();
  int d = X1.cols();
  corrmat = MatrixXd::Zero(n1,n2);
  double tmp;
  for(int i = 0; i != n1; ++i)
  {
    for(int j = 0; j != n2; ++j)
    {
      tmp = 0.0;
      for(int k = 0; k != d; ++k)
	tmp -= SQ(X1(i,k)-X2(j,k))/phi(k);
      corrmat(i,j) = exp(tmp);
    }
  }
}

void prodCorr(const VectorXd& phi, const VectorXd& vtheta,
	      const VectorXi& levels, const MatrixXd& X,
	      const MatrixXi& Z, MatrixXd& corr)
{
  int n = X.rows();
  double tmp;
  vector<MatrixXd> vQualCorr;
  MatrixXd quancorr;

  vecQualCorr(vtheta, levels, vQualCorr);
  quanCorr(phi, X, quancorr);
  corr = quancorr;
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k = 0; it != vQualCorr.end(); ++it, ++k)
  {
    for(int i=0; i!= n; ++i)
      for(int j=0; j!=i; ++j)
      {
	tmp = (*it)(Z(i,k),Z(j,k));
 	corr(i,j) *= tmp;
	corr(j,i) *= tmp;
      }
  }
}

void addCorr(const MatrixXd& phi, const VectorXd& vtheta, const VectorXd& vsigma2,
	     const VectorXi& levels, const MatrixXd& X, const MatrixXi& Z,
	     MatrixXd& corr)
{
  int n = X.rows();
  vector<MatrixXd> vQualCorr;
  MatrixXd quancorr;
  double tmp, sigma2;
  vecQualCorr(vtheta, levels,vQualCorr);
  corr = MatrixXd::Zero(n,n);
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k=0; it != vQualCorr.end(); ++it, ++k)
  {
    quanCorr(phi.col(k), X, quancorr);
    sigma2 = vsigma2(k);
    for(int i=0; i != n; ++i)
    {
      corr(i,i) += sigma2*quancorr(i,i);
      for(int j=0; j!=i; ++j)
      {
	tmp = (*it)(Z(i,k),Z(j,k))*sigma2*quancorr(i,j);
	corr(i,j) += tmp;
	corr(j,i) += tmp;
      }
    }
  }

}
// additive correlation for qualitative factors with homogoneous correlation
// for quantitative variables
void addHomCorr(const VectorXd& phi, const VectorXd& vtheta,
		const VectorXi& levels, const MatrixXd& X,
		const MatrixXi& Z, MatrixXd& corr)
{
  int n = X.rows();
  int q = Z.cols();
  double tmp;
  vector<MatrixXd> vQualCorr;
  MatrixXd quancorr;
  MatrixXd qualcorr = MatrixXd::Zero(n,n);
  
  vecQualCorr(vtheta, levels, vQualCorr);
  quanCorr(phi, X, quancorr);
  corr = quancorr;
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k = 0; it != vQualCorr.end(); ++it, ++k)
    for(int i=0; i!= n; ++i)
      for(int j=0; j!=i; ++j)
      {
	tmp = (*it)(Z(i,k),Z(j,k));
	qualcorr(i,j) += tmp;
      }
  qualcorr /= (double) q;

  for(int i = 0; i != n; ++i)
    for(int j = 0; j != i; ++j)
    {
      tmp = qualcorr(i,j);
      corr(i,j) *= tmp;
      corr(j,i) *= tmp;
    }
}

double loglikProd(const VectorXd& phi, const VectorXd& vtheta,
		  const VectorXi& levels, const MatrixXd& X,
		  const MatrixXi& Z, const VectorXd& resp,
		  double condthres, double& nug)
{
  MatrixXd corr;
  int n = X.rows();
  double dn = (double)n;
  double cond, anorm, ethres;
  MatrixXd lmat;
  prodCorr(phi, vtheta, levels, X, Z, corr);
  nugDetermine(corr, cond, anorm);

  if(log(cond) > condthres)
  {
    ethres = exp(condthres);
    nug = anorm*(cond-ethres)/cond/(ethres-1.0);
    corr += nug * MatrixXd::Identity(n,n);
  }
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");
  lmat = lltcorr.matrixL();
  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quady = resp.dot(soly);
  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);
  double rss = quady - bily1*bily1/quad1;

  double sigma2 = rss/dn;
  double logdet = lmat.diagonal().array().log().sum();

  double loglik = -0.5*dn*(LOG2PI+log(sigma2)+1.0)-logdet;
  return loglik;
}

double loglikAdd(const MatrixXd& phi, const VectorXd& vtheta, const VectorXd& vsigma2,
		 const VectorXi& levels, const MatrixXd& X, const MatrixXi& Z,
		 const VectorXd& resp, double condthres, double& nug)
{
  int n = X.rows();
  double dn = (double)n;
  double cond, anorm, ethres;
  MatrixXd corr, lmat;
  addCorr(phi, vtheta, vsigma2, levels, X, Z, corr);
  nugDetermine(corr, cond, anorm);
  if(log(cond) > condthres)
  {
    ethres = exp(condthres);
    nug = anorm*(cond-ethres)/cond/(ethres-1.0);
    corr += nug * MatrixXd::Identity(n,n);
  }
  
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");
  lmat = lltcorr.matrixL();

  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quady = resp.dot(soly);
  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);

  double rss = quady - bily1*bily1/quad1;
  double logdet = lmat.diagonal().array().log().sum();
  double loglik = -0.5*(dn*LOG2PI+rss)-logdet;
  return loglik;
}

double loglikAddHom(const VectorXd& phi, const VectorXd& vtheta,
		    const VectorXi& levels, const MatrixXd& X,
		    const MatrixXi& Z, const VectorXd& resp,
		    double condthres, double& nug)
{
  MatrixXd corr;
  int n = X.rows();
  double cond, anorm, ethres;
  double dn = (double)n;
  MatrixXd lmat;
  addHomCorr(phi, vtheta, levels, X, Z, corr);
  nugDetermine(corr, cond, anorm);
    
  if(log(cond) > condthres)
  {
    ethres = exp(condthres);
    nug = anorm*(cond-ethres)/cond/(ethres-1.0);
    corr += nug * MatrixXd::Identity(n,n);
  }
  
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");
  lmat = lltcorr.matrixL();
  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quady = resp.dot(soly);
  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);
  double rss = quady - bily1*bily1/quad1;

  double sigma2 = rss/dn;
  double logdet = lmat.diagonal().array().log().sum();

  double loglik = -0.5*dn*(LOG2PI+log(sigma2)+1.0)-logdet;
  return loglik;
}

void prodCrossCorr(const VectorXd& phi, const VectorXd& vtheta,
		   const VectorXi& levels, const MatrixXd& X1,
		   const MatrixXd& X2, const MatrixXi& Z1,
		   const MatrixXi& Z2, MatrixXd& corr)
{
  int n1 = X1.rows();
  int n2 = X2.rows();
  vector<MatrixXd> vQualCorr;
  vecQualCorr(vtheta, levels, vQualCorr);
  quanCrossCorr(phi, X1, X2, corr);
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k = 0; it != vQualCorr.end(); ++it, ++k)
    for(int i = 0; i != n1; ++i)
      for(int j = 0; j != n2; ++j)
	corr(i,j) *= (*it)(Z1(i,k),Z2(j,k));
}

void addCrossCorr(const MatrixXd& phi, const VectorXd& vtheta, const VectorXd& vsigma2,
		  const VectorXi& levels, const MatrixXd& X1, const MatrixXd& X2,
		  const MatrixXi& Z1, const MatrixXi& Z2, MatrixXd& corr)
{
  int n1 = X1.rows();
  int n2 = X2.rows();
  vector<MatrixXd> vQualCorr;
  MatrixXd quancorr;
  double sigma2;
  vecQualCorr(vtheta, levels, vQualCorr);
  corr = MatrixXd::Zero(n1,n2);
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k = 0; it != vQualCorr.end(); ++it, ++k)
  {
    quanCrossCorr(phi.col(k), X1, X2, quancorr);
    sigma2 = vsigma2(k);
    for(int i = 0; i != n1; ++i)
      for(int j = 0; j != n2; ++j)
	corr(i,j) += (*it)(Z1(i,k), Z2(j,k)) * sigma2 * quancorr(i,j);
  }
}
void addHomCrossCorr(const VectorXd& phi, const VectorXd& vtheta,
		     const VectorXi& levels, const MatrixXd& X1,
		     const MatrixXd& X2, const MatrixXi& Z1,
		     const MatrixXi& Z2, MatrixXd& corr)
{
  int n1 = X1.rows();
  int n2 = X2.rows();
  int q = Z1.cols();
  vector<MatrixXd> vQualCorr;
  MatrixXd qualcrosscorr = MatrixXd::Zero(n1,n2);
  vecQualCorr(vtheta, levels, vQualCorr);
  quanCrossCorr(phi, X1, X2, corr);
  vmciter it;
  int k;
  for(it = vQualCorr.begin(), k = 0; it != vQualCorr.end(); ++it, ++k)
    for(int i = 0; i != n1; ++i)
      for(int j = 0; j != n2; ++j)
	qualcrosscorr(i,j) += (*it)(Z1(i,k),Z2(j,k));
  qualcrosscorr /= (double) q;
  corr.noalias() = corr.cwiseProduct(qualcrosscorr);
}

void prodPredict(const MatrixXd& X0, const MatrixXi& Z0, const MatrixXd& X,
		 const MatrixXi& Z,  const VectorXd& resp,
		 const VectorXd& phi, const VectorXd& vtheta,
		 const VectorXi& levels, double nug,
		 VectorXd& pmean, VectorXd& psig2)
{
  int n = X.rows();
  double dn = (double)n;
  MatrixXd corr, crosscorr;
  prodCorr(phi, vtheta, levels, X, Z, corr);
  if(nug>0.0)
    corr += nug*MatrixXd::Identity(n,n);
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");

  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quady = resp.dot(soly);
  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);
  double muhat = bily1/quad1;
  double sigma2 = quady - bily1*bily1/quad1;
  sigma2 /= dn;
  prodCrossCorr(phi, vtheta, levels, X0, X, Z0, Z, crosscorr);
  pmean = soly - muhat * sol1;
  pmean = crosscorr * pmean;
  pmean = pmean.array() + muhat;
  MatrixXd solcross = lltcorr.solve(crosscorr.transpose());
  crosscorr = crosscorr.array() * solcross.transpose().array();
  psig2 = crosscorr.colwise().sum();
  psig2 = 1.0 - psig2.array();
  psig2 = psig2 * sigma2;
}

void addPredict(const MatrixXd& X0, const MatrixXi& Z0, const MatrixXd& X,
		const MatrixXi& Z,  const VectorXd& resp, const MatrixXd& phi,
		const VectorXd& vtheta, const VectorXd& vsigma2,
		const VectorXi& levels, double nug,
		VectorXd& pmean, VectorXd& psig2)
{
  int n = X.rows();
  MatrixXd corr, crosscorr;
  addCorr(phi, vtheta, vsigma2, levels, X, Z, corr);
  if(nug > 0.0)
    corr += nug*MatrixXd::Identity(n,n);
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");

  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);
  double muhat = bily1/quad1;
  double sigma2 = vsigma2.sum();
  addCrossCorr(phi, vtheta, vsigma2, levels, X0, X, Z0, Z, crosscorr);
  pmean = soly - muhat * sol1;
  pmean = crosscorr * pmean;
  pmean = pmean.array() + muhat;
  MatrixXd solcross = lltcorr.solve(crosscorr.transpose());
  crosscorr = crosscorr.array() * solcross.transpose().array();
  psig2 = crosscorr.colwise().sum();
  psig2 = sigma2 - psig2.array();
}

void addHomPredict(const MatrixXd& X0, const MatrixXi& Z0, const MatrixXd& X,
		   const MatrixXi& Z,  const VectorXd& resp,
		   const VectorXd& phi, const VectorXd& vtheta,
		   const VectorXi& levels, double nug,
		   VectorXd& pmean, VectorXd& psig2)
{
  int n = X.rows();
  double dn = (double)n;
  MatrixXd corr, crosscorr;
  addHomCorr(phi, vtheta, levels, X, Z, corr);
  if(nug > 0.0)
    corr += nug*MatrixXd::Identity(n,n);
  LLTMd lltcorr(corr);
  if(lltcorr.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix!");

  VectorXd one = VectorXd::Constant(n,1.0);
  VectorXd soly = lltcorr.solve(resp);
  VectorXd sol1 = lltcorr.solve(one);

  double quady = resp.dot(soly);
  double quad1 = one.dot(sol1);
  double bily1 = resp.dot(sol1);
  double muhat = bily1/quad1;
  double sigma2 = quady - bily1*bily1/quad1;
  sigma2 /= dn;
  addHomCrossCorr(phi, vtheta, levels, X0, X, Z0, Z, crosscorr);
  pmean = soly - muhat * sol1;
  pmean = crosscorr * pmean;
  pmean = pmean.array() + muhat;
  MatrixXd solcross = lltcorr.solve(crosscorr.transpose());
  crosscorr = crosscorr.array() * solcross.transpose().array();
  psig2 = crosscorr.colwise().sum();
  psig2 = 1.0 - psig2.array();
  psig2 = psig2 * sigma2;
}
