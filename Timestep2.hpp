
#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

	#include <Eigen>
	#include "Parameters.hpp"
	#include "Tools.hpp"
	#include "Functions.hpp"
	#include "SolverCG2.hpp"
	#include <mpi.h>
	#include <iostream>
	#include <fstream>
	
	using namespace std;
	
	int timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int recouvr, int me, int np, int i1, int im);

#endif
