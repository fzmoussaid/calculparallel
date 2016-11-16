#include <mpi.h>
#include <iostream>
#include <fstream>
#include <Eigen>
#include "Parameters.hpp"
#include "Tools.hpp"
#include "Functions.hpp"
#include "SolverCG.hpp"
#include "Timestep.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	int me, np, i1, im, recouvr(3);
	double alpha = 1./dt + 2.*D*(1./(dx*dx) + 1./(dy*dy)), beta = -D/(dx*dx), gamma = -D/(dy*dy);

	MPI_Init (&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	charge(ny, np, me, recouvr, i1, im);

	VectorXd U(nl);

	SolverCG solver(alpha,beta,gamma,eps,nx,ny);
	cout << "eps : " << eps << endl;
	for (int time=0; time<10; time++)
	{
		//! construction de RHS
		if(me == 0)
		{
			cout << "iter " << time << endl;
		}
		timeStep(solver, U, eps, beta, gamma, time*dt, nIterMax, nx, im-i1+1, recouvr, me, np, i1, im);
	}


	ofstream file_sol("sol/Sol" + to_string(me) + ".dat"), file_sol_exact("sol/Sol_exacte" + to_string(me) + ".dat");
	double norm_diff=0., norm_exact=0.;
	for (int j=0; j<(im-i1+1); j++)
	{
		for (int i=0; i<nx; i++)
		{
			double xi = (i+1)*dx;
			double yi = (i1+j+1)*dy;
			file_sol << xi << " " << yi << " " << U(bijection(i,j,nx)) << endl;

			if (p==2)
			{
				file_sol_exact << xi << " " << yi << " " << sin(xi)+cos(yi) << endl;
				norm_diff = max(abs(U(bijection(i,j,nx)) - (sin(xi)+cos(yi))),norm_diff);
				norm_exact = max(abs(sin(xi)+cos(yi)),norm_exact);
			}
			else if (p==1)
			{
				file_sol_exact << xi << " " << yi << " " << xi*(1-xi)*yi*(1-yi) << endl;
				norm_diff = max(abs(U(bijection(i,j,nx)) - (xi*(1-xi)*yi*(1-yi))),norm_diff);
				norm_exact = max(abs(xi*(1-xi)*yi*(1-yi)),norm_exact);
			}
		}
		file_sol << endl;
		file_sol_exact << endl;
	}
	if ((p==1)||(p==2))
	{
		cout << "erreur relative en norme infinie : " << norm_diff/norm_exact << endl;
	}

	file_sol.close();
	file_sol_exact.close();

	MPI_Finalize();



	return 0;
}
