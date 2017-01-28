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
	int me, np, i1, im, recouvr(4);
	double alpha = 1./dt + 2.*D*(1./(dx*dx) + 1./(dy*dy)), beta = -D/(dx*dx), gamma = -D/(dy*dy);

	MPI_Init (&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	charge(ny, np, me, recouvr, i1, im);

	VectorXd U(nx*(im-i1+1));
	//~ U.setZero();
	U.setOnes();

	SolverCG solver(alpha,beta,gamma,eps,nx,im-i1+1);
	for (int time=0; time<1; time++)
	{
		//! construction de RHS
		if(me == 0)
		{
			cout << "iter " << time << endl;
			cout << "Niter Schwartz : " << timeStep(solver, U, 1e-12, beta, gamma, time*dt, nIterMax, recouvr, me, np, i1, im) << endl;
		} else {
			timeStep(solver, U, 1e-12, beta, gamma, time*dt, nIterMax, recouvr, me, np, i1, im);
		}
	}

	ofstream file_sol("sol/Sol" + to_string(me) + ".dat"), file_sol_exact("sol/Sol_exacte" + to_string(me) + ".dat");
	double norm_diff=0., norm_exact=0., x;
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
				x = sin(xi)+cos(yi);
				norm_exact += x*x;
				x = U(bijection(i,j,nx)) - x;
				norm_diff += x*x;
			}
			else if (p==1)
			{
				file_sol_exact << xi << " " << yi << " " << xi*(1-xi)*yi*(1-yi) << endl;
				x = xi*(1-xi)*yi*(1-yi);
				norm_exact += x*x;
				x = U(bijection(i,j,nx)) - x;
				norm_diff += x*x;
			}
		}
		file_sol << endl;
		file_sol_exact << endl;
	}
	if ((p==1)||(p==2))
	{
		MPI_Allreduce(&norm_diff, &norm_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&norm_exact, &norm_exact, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if(me == 0)
		{
			cout << "Erreur relative en norme L2 : " << norm_diff/norm_exact << endl;
		}
	}

	file_sol.close();
	file_sol_exact.close();

	MPI_Finalize();



	return 0;
}
