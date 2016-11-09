
#include "SolverCG.hpp"

SolverCG::SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny)
	:	_alpha(alpha), _beta(beta), _gamma(gamma), _eps(eps), _nx(nx), _ny(ny)
{ }

int SolverCG::gradConj(VectorXd& X, const VectorXd& B, int Niter)
{
	int n(X.size()), iter(1);
	VectorXd R(n), W(n), D(n);
	double nr, a, b;
	
	X.setZero();
	matmulA(X, R);
	R -= B;
	D = R;
	nr = R.norm();
	
	while(nr > _eps && iter < Niter)
	{
		matmulA(D, W);
		a = D.dot(R)/D.dot(W);
		
		X -= a*D;
		
		b = 1.0/(nr*nr);
		R -= a*W;
		
		nr = R.norm();
		b *= nr*nr;
		
		D = R + b*D;
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
