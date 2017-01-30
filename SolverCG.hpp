
#ifndef SOLVERCG_HPP
#define SOLVERCG_HPP

	#include <Eigen>
	#include "Tools.hpp"
	#include "Parameters.hpp"
	#include <iostream>

	using namespace Eigen;
	using namespace std;

	class SolverCG
	{
		protected:
			double _alpha, _beta, _gamma, _eps;
			int _nx, _ny, _me, _np;
		public:
			SolverCG(double alpha, double beta, double gamma, double eps, int nx, int ny, int me, int np);
			
			/**
			 * Produit matriciel pour les conditions de Dirichlet aux bords
			 * 
			 * X : Vecteur d'entrée
			 * Y : Vecteur résultat du produit A*X
			 * 
			 * */
			void matmulA(const VectorXd& X, VectorXd& Y) const;
			
			
			/**
			 * Produit matriciel pour les conditions de Robin aux bords avec un schéma de difference finies décentré amont pour le calcul de la dérivée
			 * 
			 * X : Vecteur d'entrée
			 * Y : Vecteur résultat du produit A*X
			 * 
			 * */
			void matmulArobin_decentre(const VectorXd& X, VectorXd& Y) const;
			
			/**
			 * Algorithme du gradient conjugué
			 * 
			 * X : Vecteur en sortie solution de AX=B
			 * B : Vecteur du second membre
			 * Niter : Nombre maximum d'itérations
			 * 
			 * */
			void gradConj(VectorXd& X, const VectorXd& B, int Niter) const;
			
			
			/**
			 * Algorithme bicgstab sans préconditionnement
			 * 
			 * X : Vecteur en sortie solution de AX=B
			 * B : Vecteur du second membre
			 * Niter : Nombre maximum d'itérations
			 * 
			 * */
			void bicgstab(VectorXd& X, const VectorXd& B, int Niter) const;
	};


#endif
