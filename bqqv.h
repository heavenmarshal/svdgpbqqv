#ifndef BQQV_H
#define BQQV_H
#include<Eigen/Core>
#include<Eigen/Dense>
#include<vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::Map;
typedef Map<const MatrixXd> mapMat;
typedef Map<const VectorXd> mapVec;
typedef Map<const MatrixXi> mapMati;
typedef Map<const VectorXi> mapVeci;
using std::vector;

void qualCorr(const VectorXd&, int, MatrixXd&);
void vecQualCorr(const VectorXd&, const VectorXi&, vector<MatrixXd>&);
void quanCorr(const VectorXd&, const MatrixXd&, MatrixXd&);
void prodCorr(const VectorXd&, const VectorXd&, const VectorXi&,
	      const MatrixXd&, const MatrixXi&, MatrixXd&);
void addCorr(const MatrixXd&, const VectorXd&, const VectorXd&,
	     const VectorXi&, const MatrixXd&, const MatrixXi&,
	     MatrixXd&);
void addHomCorr(const VectorXd&,const VectorXd&,
		const VectorXi&, const MatrixXd&,
		const MatrixXi&, MatrixXd&);
double loglikProd(const VectorXd&, const VectorXd&,
		  const VectorXi&, const MatrixXd&,
		  const MatrixXi&, const VectorXd&,
		  double, double&);
double loglikAdd(const MatrixXd&, const VectorXd&, const VectorXd&,
		 const VectorXi&, const MatrixXd&, const MatrixXi&,
		 const VectorXd&, double, double&);
double loglikAddHom(const VectorXd&, const VectorXd&,
		    const VectorXi&, const MatrixXd&,
		    const MatrixXi&, const VectorXd&,
		    double, double&);
void prodCrossCorr(const VectorXd&, const VectorXd&,
		   const VectorXi&, const MatrixXd&,
		   const MatrixXd&, const MatrixXi&,
		   const MatrixXi&, MatrixXd&);
void addCrossCorr(const MatrixXd&, const VectorXd&, const VectorXd&,
		  const VectorXi&, const MatrixXd&, const MatrixXd&,
		  const MatrixXi&, const MatrixXi&, MatrixXd&);
void addHomCrossCorr(const VectorXd&, const VectorXd&,
		     const VectorXi&, const MatrixXd&,
		     const MatrixXd&, const MatrixXi&,
		     const MatrixXi&, MatrixXd&);
void prodPredict(const MatrixXd&, const MatrixXi&, const MatrixXd&,
		 const MatrixXi&, const VectorXd&,
		 const VectorXd&, const VectorXd&,
		 const VectorXi&, double, VectorXd&, VectorXd&);

void addPredict(const MatrixXd&, const MatrixXi&, const MatrixXd&,
		const MatrixXi&, const VectorXd&, const MatrixXd&,
		const VectorXd&, const VectorXd&,
		const VectorXi&, double, VectorXd&, VectorXd&);
void addHomPredict(const MatrixXd&, const MatrixXi&, const MatrixXd&,
		   const MatrixXi&,  const VectorXd&,
		   const VectorXd&, const VectorXd&,
		   const VectorXi&, double, VectorXd&, VectorXd&);
#endif
