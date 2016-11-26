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
double loglikProd(const VectorXd&, const VectorXd&,
		  const VectorXi&, const MatrixXd&,
		  const MatrixXi&, const VectorXd&);
double loglikAdd(const MatrixXd&, const VectorXd&, const VectorXd&,
		 const VectorXi&, const MatrixXd&, const MatrixXi&,
		 const VectorXd&);
void prodCrossCorr(const VectorXd&, const VectorXd&,
		   const VectorXi&, const MatrixXd&,
		   const MatrixXd&, const MatrixXi&,
		   const MatrixXi&, MatrixXd&);
void addCrossCorr(const MatrixXd&, const VectorXd&, const VectorXd&,
		  const VectorXi&, const MatrixXd&, const MatrixXd&,
		  const MatrixXi&, const MatrixXi&, MatrixXd&);

void prodPredict(const MatrixXd&, const MatrixXi&, const MatrixXd&,
		 const MatrixXi&, const VectorXd&,
		 const VectorXd&, const VectorXd&,
		 const VectorXi&, VectorXd&, VectorXd&);

void addPredict(const MatrixXd&, const MatrixXi&, const MatrixXd&,
		const MatrixXi&, const VectorXd&, const MatrixXd&,
		const VectorXd&, const VectorXd&,
		const VectorXi&, VectorXd&, VectorXd&);

#endif
