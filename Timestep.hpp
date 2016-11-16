
#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

	#include <Eigen>
	#include "Parameters.hpp"
	#include "Tools.hpp"
	#include "Functions.hpp"
	#include "SolverCG.hpp"
	#include <mpi.h>
	
	void timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int nx, int ny, int recouvr, int me, int np, int i1, int im);

#endif
