#include <mpi.h>
#include <iostream>
#include <fstream>
#include <Eigen>
#include "Parameters.hpp"
#include "Tools.hpp"
#include "Functions.hpp"
#include "SolverCG2.hpp"
#include "Timestep2.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	int me, np, i1, im, recouvr(2);
	//~ double alpha = 1./dt + 2.*D*(1./(dx*dx) + 1./(dy*dy)), beta = -D/(dx*dx), gamma = -D/(dy*dy);
	double alpha = 2., beta = 2., gamma = 2.;

	MPI_Init (&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	charge(ny, np, me, recouvr, i1, im);
	
	double nyLocal(im-i1+1);
	VectorXd U(nx*(nyLocal+(me > 0)+(me < np-1)));
	
	//~ U.setZero();
	U.setOnes();
	
	SolverCG solver(alpha,beta,gamma,eps,nx,nyLocal,me,np);
	for (int time=0; time<1; time++)
	{
		//! construction de RHS
		if(me == 0)
		{
			cout << "iter " << time << endl;
		}
		cout << "Niter Schwartz : " << timeStep(solver, U, 1e-12, beta, gamma, time*dt, nIterMax, recouvr, me, np, i1, im) << endl;
	}

	ofstream file_sol("sol/Sol" + to_string(me) + ".dat"), file_sol_exact("sol/Sol_exacte" + to_string(me) + ".dat");
	double norm_diff=0., norm_exact=0.;
	
	int decal(me > 0); // On décale d'une ligne en présence d'une ligne fantome
	double yi, xi;
	for (int j = decal; j < nyLocal+decal; j++)
	{
		yi = (i1+j+1-decal)*dy;
		for (int i = 0; i < nx; i++)
		{
			xi = (i+1)*dx;
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
