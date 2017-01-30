
#ifndef TOOLS_HPP
#define TOOLS_HPP

	#include <Eigen>
	
	using namespace Eigen;
	
	/**
	 * Bijection entre coordonnées 2D et coordonnées 1D
	 * pour un stockage ligne par ligne
	 * */
	int bijection(int i, int j, int nx);
	
	/**
	 * Calcul de la charge en tenant compte du recouvrement
	 * */
	void charge(int n, int np, int me, int recouvr, int& i1, int& im);
	
	/**
	 * Calcul de l'erreur relative en norme L2
	 * 		v : Vecteur calculé
	 * 		vref : vecteur de référence
	 * En sortie : 
	 * 		norm : erreur relative
	 * 		normRef : norme L2 de référence
	 * */
	void relativeNormL2(const VectorXd& v, const VectorXd& vref, double& norm, double& normRef);


#endif
