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
	double alpha(1.0), beta(1.0), gamma(1.0);
	X.setOnes();
	
	SolverCG var(alpha, beta, gamma, N, N);
	
	var.matmulA(X,Y);
	
	for(int i(0); i < N*N; ++i) {
		cout << Y(i) << endl;
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
