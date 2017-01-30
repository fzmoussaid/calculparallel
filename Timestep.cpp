
#include "Timestep.hpp"

/**
 *  Communication avec <proc> :
 * 		sendbuf : Vecteur à envoyer
 * 		recvbuf : Vecteur reçu
 * */
void commMPI(void *sendbuf, void *recvbuf, int proc, int tag1, int tag2)
{
	MPI_Status status;
	MPI_Sendrecv(sendbuf, nx, MPI_DOUBLE, proc, tag1,
			recvbuf, nx, MPI_DOUBLE, proc, tag2,
			MPI_COMM_WORLD, &status);
}

/**
 * Gère les communication d'une itération de Swcharz additif
 * 		U : solution calculée localement
 * 		Uhaut : ligne du haut reçue
 * 		Ubas : ligne du bas reçue
 *		nyLocal : charge
 * 		recouvr : taille du recouvrement
 * */
void commSchwarz(VectorXd& U, VectorXd& Uhaut, VectorXd& Ubas, int nyLocal, int recouvr, int me, int np)
{
	void *sendbuf = 0;
	VectorXd vecSend;
	//---------------- Communication : reçoit le bord bas de me+1 et envoie son bord haut à me+1  -----------------------
	if(me < np-1) {
		if(BC == 0) // Condition de Dirichlet
		{
			sendbuf = &U.data()[(nyLocal-(recouvr+1))*nx];
		}
		else if(BC == 1) // Condition de Robin
		{
			vecSend.resize(nx);
			for (int i=0; i<nx; i++)
			{
				vecSend(i) = U((nyLocal-(recouvr+1))*nx+i) - (b/(a*dy+b))*U((nyLocal-(recouvr+1))*nx+i+nx); // calcul de a*u + b*du/dn
			}
			sendbuf = vecSend.data();
		}
		commMPI(sendbuf, Uhaut.data(), me+1, 200+me, 300+me+1);
	}
	//---------------- Communication : reçoit le bord haut de me-1 et envoie son bord bas à me-1  -----------------------
	if(me > 0) {
		if(BC == 0) // Condition de Dirichlet
		{
			sendbuf = &U.data()[(recouvr)*nx];
		}
		else if(BC == 1) // Condition de Robin
		{
			vecSend.resize(nx);
			for (int i=0; i<nx; i++)
			{
				vecSend(i) = U(recouvr*nx+i) - (b/(a*dy+b))*U(recouvr*nx+i-nx); // calcul de a*u + b*du/dn
			}
			sendbuf = vecSend.data();
		}
		commMPI(sendbuf, Ubas.data(), me-1, 300+me, 200+me-1);
	}
}

/**
 * Construit le second membre :
 * 		U : solution calculée localement
 * 		Uhaut : ligne du haut reçue
 * 		Ubas : ligne du bas reçue
 * 		nyLocal : charge
 * 		recouvr : taille du recouvrement
 * */
void computeRHS(VectorXd& Rhs, const VectorXd& Un, const VectorXd& Ubas, const VectorXd& Uhaut, double gamma, double beta, double t, int nyLocal, int me, int np, int i1)
{
	Rhs.setZero();
	//---------------- Contribution du bord haut -----------------------------------------------------------------------
	if(me < np-1) {
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,nyLocal-1,nx)) = -gamma*Uhaut(i);
		}
	} else {
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
		}
	}
	
	//---------------- Contribution du bord bas -------------------------------------------------------------------------
	if(me > 0) {
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
		}
	} else {
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
		}
	}
	

	//---------------- Contribution du terme source ---------------------------------------------------------------------
	for (int j=0; j<nyLocal; j++)
	{
		for (int i=1; i<nx-1; i++)
		{
			Rhs(bijection(i,j,nx)) += functionF((i+1)*dx,(i1+j+1)*dy,t,p);
		}
		Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(i1+j+1)*dy,t,p);
		Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(i1+j+1)*dy,t,p);
	}

	//---------------- Schéma d'Euler explicite -------------------------------------------------------------------------
	Rhs += Un/dt;
}

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
void timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int recouvr, int me, int np, int i1, int im)
{
	double norm, normRef, normHaut, normHautRef, normBas, normBasRef;
	int nyLocal(im-i1+1);
	VectorXd Ubas(nx), Uhaut(nx), UbasNext(nx), UhautNext(nx), Unext(nx*nyLocal), Rhs(nx*nyLocal);
	
	int iter(0); // Compteur d'itérations
	
	if(np == 1) // Cas à un seul proc, on effectue une résolution standard
	{
		Rhs.setZero();
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
		}

		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
		}

		for (int j=0; j<ny; j++)
		{
			for (int i=1; i<nx-1; i++)
			{
				Rhs(bijection(i,j,nx)) += functionF((i+1)*dx,(i1+j+1)*dy,t,p);
			}
			Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(i1+j+1)*dy,t,p);
			Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(i1+j+1)*dy,t,p);
		}

		Rhs += Un/dt;

		if(SOLVER == 0)
			var.gradConj(Unext, Rhs, nIterMax);
		else if(SOLVER == 1)
			var.bicgstab(Unext, Rhs, nIterMax);

		Un = Unext;
		iter = 1;
	}
	else // Plus d'un proc, on effectue l'algorithme de Schwarz additif
	{
		//---------------- Communications des conditions de bords initiales -------------------------------------------------------------

		commSchwarz(Un, Uhaut, Ubas, nyLocal, recouvr, me, np);
		
		do
		{
			//---------------- Construction du second membre ----------------------------------------------------------------------------
			
			computeRHS(Rhs, Un, Ubas, Uhaut, gamma, beta, t, nyLocal, me, np, i1);
			
			//---------------- Résolution du système linéaire (0 = Gradient conjugué; 1 = BICGSTAB) -------------------------------------

			if(SOLVER == 0)
				var.gradConj(Unext, Rhs, nIterMax);
			else if(SOLVER == 1)
				var.bicgstab(Unext, Rhs, nIterMax);

			//---------------- Communications des conditions de bords ---------------------------------------------------------------------

			commSchwarz(Unext, UhautNext, UbasNext, nyLocal, recouvr, me, np);
			
			//---------------- Calcul de la norme des conditions de bords -----------------------------------------------------------------

			norm = 0.;
			normRef = 0.;
			if(me < np-1) {
				relativeNormL2(UhautNext, Uhaut, normHaut, normHautRef);
				norm += normHaut;
				normRef += normHautRef;
			}
			if(me > 0) {
				relativeNormL2(UbasNext, Ubas, normBas, normBasRef);
				norm += normBas;
				normRef += normBasRef;
			}

			MPI_Allreduce(&norm, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&normRef, &normRef, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			Un = Unext;
			Uhaut = UhautNext;
			Ubas = UbasNext;
			iter++;
		} while(norm > (eps*normRef) && iter < Niter);
	}
}
