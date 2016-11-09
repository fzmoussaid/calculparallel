
#ifndef SOLVERCG_HPP
#define SOLVERCG_HPP

	#include <Eigen>
	#include "Tools.hpp"
	
	using namespace Eigen;

	class SolverCG
	{
		protected:
			double _alpha, _beta, _gamma, _eps;
			int _nx, _ny;
		public:
			SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny);
			void matmulA(const VectorXd& X, VectorXd& Y) const;
			int gradConj(VectorXd& X, const VectorXd& B, int Niter);	
	};


#endif

