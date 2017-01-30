#include <mpi.h>
#include <iostream>
#include <fstream>
#include <Eigen>			// On utilise la bibliothèque Eigen pour la gestion des tableaux
#include "Parameters.hpp"	// Variables globales : (tailles du domaine, choix du solveur, choix des conditions de bords (Dirichlet ou Robin), ...)
#include "Tools.hpp"
#include "Functions.hpp"	// Choix des conditions de bords (fonctions f,g,h)
#include "SolverCG.hpp"		// Solveurs
#include "Timestep.hpp"		// Iteration de l'algorithme de Schwarz

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	int me; // numéro du proc courant
	int np; // nombre de procs
	int i1, im; // intervalle de la charge
	int recouvr = 6; // taille du recouvrement
	
	// Constantes dans la matrice A :
	double alpha = 1./dt + 2.*D*(1./(dx*dx) + 1./(dy*dy));
	double beta = -D/(dx*dx);
	double gamma = -D/(dy*dy);


	// Initialisation de MPI :
	MPI_Init (&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	charge(ny, np, me, recouvr, i1, im); // On répartit la charge suivant l'axe y
	int nyLocal(im-i1+1); // Charge du proc courant

	// Déclaration et initialisation du vecteur solution;
	VectorXd U(nx*nyLocal);
	U.setZero();
	
	// Déclaration du solveur, on lui passe les constantes de la matrice, le critère d'arret eps des solveurs, les tailles locales, me et np
	SolverCG solver(alpha, beta, gamma, eps, nx, nyLocal, me, np);
	
	// Boucle en temps :
	for (int time = 0; time < 1; time++)
	{
		if(me == 0)
		{
			cout << "Iteration en temps numero " << time << endl;
		}
		timeStep(solver, U, 1e-12, beta, gamma, time*dt, nIterMax, recouvr, me, np, i1, im); // Iterations de Schartz
	}

	// Ecriture des fichiers de sortie et calcul de l'erreur L2 :
	
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
