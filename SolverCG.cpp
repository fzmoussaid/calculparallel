#include "SolverCG.hpp"

SolverCG::SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny)
	:	_alpha(alpha), _beta(beta), _gamma(gamma), _eps(eps), _nx(nx), _ny(ny)
{ }

int SolverCG::gradConj(VectorXd& X, const VectorXd& B, int Niter, int proc, int np) const
{
	int n(X.size()), iter(1);
	VectorXd R(n), W(n), D(n);
	double nr, a, b;

	X.setZero();
	if (BC==0)
		matmulA(X, R);
	else if (BC==1)
		matmulArobin_decentre(X,R,proc,np);
	else if (BC==2)
		matmulArobin_centre(X,R,proc,np);
	R -= B;
	D = R;
	nr = R.norm();

	while((nr>_eps) && (iter<Niter))
	{
		if (BC==0)
			matmulA(D, W);
		else if (BC==1)
			matmulArobin_decentre(D,W,proc,np);
		else if (BC==2)
			matmulArobin_centre(D,W,proc,np);

		a = D.dot(R)/D.dot(W);

		X -= a*D;

		b = 1.0/(nr*nr);
		R -= a*W;

		nr = R.norm();
		b *= nr*nr;

		D = R + b*D;
		iter++;
	}
	//~ cout << "iter " << iter << " residu " << nr << endl;
	return iter;
}


int SolverCG::bicgstab(VectorXd& X, const VectorXd& B, int Niter, int proc, int np) const
{
	int n(X.size()), iter(1);
	VectorXd R(n), R0(n), V(n), P(n), T(n);
	double nr, rho, alpha, omega, beta;

	X.setZero();
	if (BC==0)
		matmulA(X, R);
	else if (BC==1)
		matmulArobin_decentre(X,R,proc,np);
	else if (BC==2)
		matmulArobin_centre(X,R,proc,np);
	R = B - R;
	R0 = R;
	rho = alpha = omega = 1.0;
	P.setZero();
	V.setZero();
	nr = R.norm();

	while(nr > _eps && iter < Niter)
	{
		beta = (rho*omega);
		rho = R0.dot(R);
		beta = (alpha*rho)/beta;
		
		P = R + beta*(P - omega*V);
		if (BC==0)
			matmulA(P, V);
		else if (BC==1)
			matmulArobin_decentre(P,V,proc,np);
		else if (BC==2)
			matmulArobin_centre(P,V,proc,np);
		
		alpha = rho/(R0.dot(V));
		X += alpha*P;
		
		R -= alpha*V;
		nr = R.norm();
		
		if(nr > _eps) {
			if (BC==0)
				matmulA(R, T);
			else if (BC==1)
				matmulArobin_decentre(R,T,proc,np);
			else if (BC==2)
				matmulArobin_centre(R,T,proc,np);
			
			omega = T.dot(R)/T.dot(T);
			X += omega*R;
			
			R -= omega*T;
			nr = R.norm();
		}
		iter++;
	}

	return iter;
}

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

void SolverCG::matmulArobin_decentre(const VectorXd& X, VectorXd& Y, int proc, int np) const
{
	Y.resize(X.size());

	if (proc==0)
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

	if (proc == np-1)
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

void SolverCG::matmulArobin_centre(const VectorXd& X, VectorXd& Y, int proc, int np) const
{
	Y.resize(X.size());

	if (proc==0)
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
		double alpha_robin = _alpha + (b*D)/(2.*a*dy*dy*dy);

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

	if (proc == np-1)
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
		double alpha_robin = _alpha + (b*D)/(2.*a*dy*dy*dy);
		// Dernier bloc
		Y(bijection(0,_ny-1,_nx)) = _gamma*X(bijection(0,_ny-2,_nx)) + alpha_robin*X(bijection(0,_ny-1,_nx)) + _beta*X(bijection(1,_ny-1,_nx));
		for(int i = 1; i < _nx-1; i++)
		{
			Y(bijection(i,_ny-1,_nx)) = _gamma*X(bijection(i,_ny-2,_nx)) + _beta*X(bijection(i-1,_ny-1,_nx)) + alpha_robin*X(bijection(i,_ny-1,_nx)) + _beta*X(bijection(i+1,_ny-1,_nx));
		}
		Y(bijection(_nx-1,_ny-1,_nx)) = _gamma*X(bijection(_nx-1,_ny-2,_nx)) + _beta*X(bijection(_nx-2,_ny-1,_nx)) + alpha_robin*X(bijection(_nx-1,_ny-1,_nx));
	}


}
