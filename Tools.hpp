
#ifndef TOOLS_HPP
#define TOOLS_HPP

	#include <Eigen>
	
	using namespace Eigen;

	int bijection(int i, int j, int nx);
	void charge(int n, int np, int me, int recouvr, int& i1, int& im);
	void relativeNormL2(const VectorXd& v, const VectorXd& vref, double& norm, double& normRef);


#endif
