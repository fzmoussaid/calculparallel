#include <Eigen>
#include "Tools.hpp"
#include "SolverCG.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
//~ #include <mpi.h>

#define N 3

using namespace std;
using namespace Eigen;


int main(int argc, char** argv)
{
	int me, np, i1, im;
	double t, tp;
	
	VectorXd X(N*N),Y;
	double alpha(1.0), beta(1.0), gamma(1.0), eps(1e-8);
	X = 4.2*VectorXd::Ones(N*N);
	
	SolverCG var(alpha, beta, gamma, eps, N, N);
	
	var.matmulA(X,Y);
	
	X.setZero();
	cout << "nb iter : " << var.gradConj(X, Y, N*N) << endl;
	
	for(int i(0); i < N*N; ++i) {
		cout << X(i) << endl;
	}
	
	//~ MPI_Init (&argc, &argv);
	//~ tp = MPI_Wtime();
	//~ 
	//~ MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//~ MPI_Comm_size(MPI_COMM_WORLD, &np);
	//~ 
	//~ int nx(10), ny(10);
	//~ charge(ny, np, me, i1, im);
	//~ 
	//~ 
	//~ 
	//~ 
	//~ tp = MPI_Wtime() - tp;
	//~ MPI_Finalize();
	//~ 
	//~ cout << "Elapsed time : " << tp << " s" << endl;
	
	return 0;
}
