
#ifndef SOLVERCG_HPP
#define SOLVERCG_HPP

	#include <Eigen>
	#include "Tools.hpp"
	#include "Parameters.hpp"
	#include <iostream>

	using namespace Eigen;
	using namespace std;

	class SolverCG
	{
		protected:
			double _alpha, _beta, _gamma, _eps;
			int _nx, _ny;
		public:
			SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny);
			void matmulA(const VectorXd& X, VectorXd& Y) const;
			void matmulArobin_decentre(const VectorXd& X, VectorXd& Y, int proc, int np) const;
			int gradConj(VectorXd& X, const VectorXd& B, int Niter, int proc, int np) const;
			int bicgstab(VectorXd& X, const VectorXd& B, int Niter, int proc, int np) const;
	};


#endif
