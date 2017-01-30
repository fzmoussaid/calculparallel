#include "SolverCG.hpp"

SolverCG::SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny, int me, int np)
	:	_alpha(alpha), _beta(beta), _gamma(gamma), _eps(eps), _nx(nx), _ny(ny), _me(me), _np(np)
{ }


/**
 * Algorithme du gradient conjugué
 * 
 * X : Vecteur en sortie solution de AX=B
 * B : Vecteur du second membre
 * Niter : Nombre maximum d'itérations
 * 
 * */
void SolverCG::gradConj(VectorXd& X, const VectorXd& B, int Niter) const
{
	int n(X.size());
	int iter(1); // Compteur d'itérations
	
	VectorXd R(n), W(n), D(n);
	double nr, nrb, a, b;

	X.setZero();
	
	// R = AX, en fonction du choix de condition aux bords :
	if (BC==0)
		matmulA(X, R);
	else if (BC==1)
		matmulArobin_decentre(X,R);
		
	R -= B; // Au final : R = AX - B
	D = R;
	nr = R.norm();  // ||r||
	nrb = B.norm(); // ||b||

	while(nr > _eps*nrb && iter < Niter) // Critère d'arret : ||r|| / ||b|| <= epsilon
	{
		// W = AD, en fonction du choix de condition aux bords :
		if (BC==0)
			matmulA(D, W);
		else if (BC==1)
			matmulArobin_decentre(D,W);

		a = D.dot(R)/D.dot(W);

		X -= a*D;

		b = 1.0/(nr*nr);
		R -= a*W;

		nr = R.norm();
		b *= nr*nr;

		D = R + b*D;
		iter++;
	}
}


/**
 * Algorithme bicgstab sans préconditionnement
 * 
 * X : Vecteur en sortie solution de AX=B
 * B : Vecteur du second membre
 * Niter : Nombre maximum d'itérations
 * 
 * */
void SolverCG::bicgstab(VectorXd& X, const VectorXd& B, int Niter) const
{
	int n(X.size());
	int iter(1); // Compteur d'itérations
	
	VectorXd R(n), R0(n), V(n), P(n), T(n);
	double nr, nrb, rho, alpha, omega, beta;

	X.setZero();
	
	// R = AX, en fonction du choix de condition aux bords :
	if (BC==0)
		matmulA(X, R);
	else if (BC==1)
		matmulArobin_decentre(X,R);
		
	
	R = B - R; // Au final : R = B - AX
	R0 = R;
	rho = alpha = omega = 1.0;
	P.setZero();
	V.setZero();
	
	nr = R.norm();  // ||r||
	nrb = B.norm(); // ||b||

	while(nr > _eps*nrb && iter < Niter) // Critère d'arret : ||r|| / ||b|| <= epsilon
	{
		beta = (rho*omega);
		rho = R0.dot(R);
		beta = (alpha*rho)/beta;
		
		P = R + beta*(P - omega*V);
		
		// V = AP, en fonction du choix de condition aux bords :
		if (BC==0)
			matmulA(P, V);
		else if (BC==1)
			matmulArobin_decentre(P,V);
		
		alpha = rho/(R0.dot(V));
		X += alpha*P;
		
		R -= alpha*V;
		nr = R.norm();
		
		// Si le critère d'arret est vérifié on s'arrête ici, sinon :
		if(nr > _eps*nrb) {
			
			// T = AR, en fonction du choix de condition aux bords :
			if (BC==0)
				matmulA(R, T);
			else if (BC==1)
				matmulArobin_decentre(R,T);
			
			omega = T.dot(R)/T.dot(T);
			X += omega*R;
			
			R -= omega*T;
			nr = R.norm();
		}
		iter++;
	}
}

/**
 * Produit matriciel pour les conditions de Dirichlet aux bords
 * 
 * X : Vecteur d'entrée
 * Y : Vecteur résultat du produit A*X
 * 
 * */
void SolverCG::matmulA(const VectorXd& X, VectorXd& Y) const
{
	Y.resize(X.size());

	// 1er bloc
	Y(bijection(0,0,_nx)) = _alpha*X(bijection(0,0,_nx)) + _beta*X(bijection(1,0,_nx)) + _gamma*X(bijection(0,1,_nx));
	for(int i = 1; i < _nx-1; i++)
	{
		Y(bijection(i,0,_nx)) = _beta*X(bijection(i-1,0,_nx)) + _alpha*X(bijection(i,0,_nx)) + _beta*X(bijection(i+1,0,_nx)) + _gamma*X(bijection(i,1,_nx));
	}
	Y(bijection(_nx-1,0,_nx)) = _beta*X(bijection(_nx-2,0,_nx)) + _alpha*X(bijection(_nx-1,0,_nx)) + _gamma*X(bijection(_nx-1,1,_nx));


	// Milieu Matrice
	for(int j = 1; j < _ny-1; j++)
	{
		Y(bijection(0,j,_nx)) = _gamma*X(bijection(0,j-1,_nx)) + _alpha*X(bijection(0,j,_nx)) + _beta*X(bijection(1,j,_nx)) + _gamma*X(bijection(0,j+1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,j,_nx)) = _gamma*X(bijection(i,j-1,_nx)) + _beta*X(bijection(i-1,j,_nx)) + _alpha*X(bijection(i,j,_nx)) + _beta*X(bijection(i+1,j,_nx)) + _gamma*X(bijection(i,j+1,_nx));
        }
		Y(bijection(_nx-1,j,_nx)) = _gamma*X(bijection(_nx-1,j-1,_nx)) + _beta*X(bijection(_nx-2,j,_nx)) + _alpha*X(bijection(_nx-1,j,_nx)) + _gamma*X(bijection(_nx-1,j+1,_nx));
	}

	// Dernier bloc
	Y(bijection(0,_ny-1,_nx)) = _gamma*X(bijection(0,_ny-2,_nx)) + _alpha*X(bijection(0,_ny-1,_nx)) + _beta*X(bijection(1,_ny-1,_nx));
	for(int i = 1; i < _nx-1; i++)
	{
		Y(bijection(i,_ny-1,_nx)) = _gamma*X(bijection(i,_ny-2,_nx)) + _beta*X(bijection(i-1,_ny-1,_nx)) + _alpha*X(bijection(i,_ny-1,_nx)) + _beta*X(bijection(i+1,_ny-1,_nx));
	}
	Y(bijection(_nx-1,_ny-1,_nx)) = _gamma*X(bijection(_nx-1,_ny-2,_nx)) + _beta*X(bijection(_nx-2,_ny-1,_nx)) + _alpha*X(bijection(_nx-1,_ny-1,_nx));

}

/**
 * Produit matriciel pour les conditions de Robin aux bords avec un schéma de difference finies décentré amont pour le calcul de la dérivée
 * 
 * X : Vecteur d'entrée
 * Y : Vecteur résultat du produit A*X
 * 
 * */
void SolverCG::matmulArobin_decentre(const VectorXd& X, VectorXd& Y) const
{
	Y.resize(X.size());

	if(_me == 0) // Les bords haut et bas du domaine entier présentent toujours des conditions de Dirichlet :
	{
		// 1er bloc
		Y(bijection(0,0,_nx)) = _alpha*X(bijection(0,0,_nx)) + _beta*X(bijection(1,0,_nx)) + _gamma*X(bijection(0,1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,0,_nx)) = _beta*X(bijection(i-1,0,_nx)) + _alpha*X(bijection(i,0,_nx)) + _beta*X(bijection(i+1,0,_nx)) + _gamma*X(bijection(i,1,_nx));
		}
		Y(bijection(_nx-1,0,_nx)) = _beta*X(bijection(_nx-2,0,_nx)) + _alpha*X(bijection(_nx-1,0,_nx)) + _gamma*X(bijection(_nx-1,1,_nx));
	}
	else
	{
		// Pour les conditions de Robin, seul le facteur devant la diagonale de la matrice, et le second membre changent :
		double alpha_robin = _alpha - (b*D)/((b+a*dy)*(dy*dy));

		// 1er bloc
		Y(bijection(0,0,_nx)) = alpha_robin*X(bijection(0,0,_nx)) + _beta*X(bijection(1,0,_nx)) + _gamma*X(bijection(0,1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,0,_nx)) = _beta*X(bijection(i-1,0,_nx)) + alpha_robin*X(bijection(i,0,_nx)) + _beta*X(bijection(i+1,0,_nx)) + _gamma*X(bijection(i,1,_nx));
		}
		Y(bijection(_nx-1,0,_nx)) = _beta*X(bijection(_nx-2,0,_nx)) + alpha_robin*X(bijection(_nx-1,0,_nx)) + _gamma*X(bijection(_nx-1,1,_nx));
	}



	// Milieu Matrice
	for(int j = 1; j < _ny-1; j++)
	{
		Y(bijection(0,j,_nx)) = _gamma*X(bijection(0,j-1,_nx)) + _alpha*X(bijection(0,j,_nx)) + _beta*X(bijection(1,j,_nx)) + _gamma*X(bijection(0,j+1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,j,_nx)) = _gamma*X(bijection(i,j-1,_nx)) + _beta*X(bijection(i-1,j,_nx)) + _alpha*X(bijection(i,j,_nx)) + _beta*X(bijection(i+1,j,_nx)) + _gamma*X(bijection(i,j+1,_nx));
        }
		Y(bijection(_nx-1,j,_nx)) = _gamma*X(bijection(_nx-1,j-1,_nx)) + _beta*X(bijection(_nx-2,j,_nx)) + _alpha*X(bijection(_nx-1,j,_nx)) + _gamma*X(bijection(_nx-1,j+1,_nx));
	}

	if(_me == _np-1) // Les bords haut et bas du domaine entier présentent toujours des conditions de Dirichlet :
	{
		// Dernier bloc
		Y(bijection(0,_ny-1,_nx)) = _gamma*X(bijection(0,_ny-2,_nx)) + _alpha*X(bijection(0,_ny-1,_nx)) + _beta*X(bijection(1,_ny-1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,_ny-1,_nx)) = _gamma*X(bijection(i,_ny-2,_nx)) + _beta*X(bijection(i-1,_ny-1,_nx)) + _alpha*X(bijection(i,_ny-1,_nx)) + _beta*X(bijection(i+1,_ny-1,_nx));
		}
		Y(bijection(_nx-1,_ny-1,_nx)) = _gamma*X(bijection(_nx-1,_ny-2,_nx)) + _beta*X(bijection(_nx-2,_ny-1,_nx)) + _alpha*X(bijection(_nx-1,_ny-1,_nx));
	}
	else
	{
		// Pour les conditions de Robin, seul le facteur devant la diagonale de la matrice, et le second membre changent :
		double alpha_robin = _alpha - (b*D)/((b+a*dy)*(dy*dy));
		
		// Dernier bloc
		Y(bijection(0,_ny-1,_nx)) = _gamma*X(bijection(0,_ny-2,_nx)) + alpha_robin*X(bijection(0,_ny-1,_nx)) + _beta*X(bijection(1,_ny-1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,_ny-1,_nx)) = _gamma*X(bijection(i,_ny-2,_nx)) + _beta*X(bijection(i-1,_ny-1,_nx)) + alpha_robin*X(bijection(i,_ny-1,_nx)) + _beta*X(bijection(i+1,_ny-1,_nx));
		}
		Y(bijection(_nx-1,_ny-1,_nx)) = _gamma*X(bijection(_nx-1,_ny-2,_nx)) + _beta*X(bijection(_nx-2,_ny-1,_nx)) + alpha_robin*X(bijection(_nx-1,_ny-1,_nx));
	}
}
