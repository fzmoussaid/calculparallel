
#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

	#include <Eigen>
	#include "Parameters.hpp"
	#include "Tools.hpp"
	#include "Functions.hpp"
	#include "SolverCG.hpp"
	#include <mpi.h>
	#include <iostream>
	#include <fstream>
	
	using namespace std;
	
	/**
	 *  Communication avec <proc> :
	 * 		sendbuf : Vecteur à envoyer
	 * 		recvbuf : Vecteur reçu
	 * */
	void commMPI(void *sendbuf, void *recvbuf, int proc, int tag1, int tag2);
	
	/**
	 * Gère les communication d'une itération de Swcharz additif
	 * 		U : solution calculée localement
	 * 		Uhaut : ligne du haut reçue
	 * 		Ubas : ligne du bas reçue
	 *		nyLocal : charge
	 * 		recouvr : taille du recouvrement
	 * */
	void commSchwarz(VectorXd& U, VectorXd& Uhaut, VectorXd& Ubas, int nyLocal, int recouvr, int me, int np);
	
	/**
	 * Construit le second membre :
	 * 		U : solution calculée localement
	 * 		Uhaut : ligne du haut reçue
	 * 		Ubas : ligne du bas reçue
	 * 		nyLocal : charge
	 * 		recouvr : taille du recouvrement
	 * */
	void computeRHS(VectorXd& Rhs, const VectorXd& Un, const VectorXd& Ubas, const VectorXd& Uhaut, double gamma, double beta, double t, int nyLocal, int me, int np, int i1);
	
	/**
	 * Iterations de Schwarz additif :
	 * 		var : solveur
	 * 		Un : solution calculée au temps précédent
	 * 		eps : critère d'arret des itérations de Schwarz
	 * 		beta, gamma : constantes dans la matrice A
	 * 		t : temps
	 * 		Niter : nombre maximal d'itérations de Schwarz
	 * 		recouvre : taille du recouvrement
	 * */
	void timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int recouvr, int me, int np, int i1, int im);

#endif
