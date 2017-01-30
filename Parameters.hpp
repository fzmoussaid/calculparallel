
#ifndef PARAMETER_HPP
#define PARAMETER_HPP

const double eps = 1.e-12;					// Critère d'arret du solveur
const int nx = 40, ny = 40;						// Nombre de points suivant chaque direction
const double lx = 1., ly = 1.;					// Tailles du domaine
const double dt = 10.;						// Pas de temps
const double dx = lx/(nx+1), dy = ly/(ny+1);	// Tailles des mailles 
const int nIterMax = 10000;					// Nombre maximum d'itérations du solveur
const double pi = 3.1415926535897932384;	// Constante pi
const int D = 1;							// Coefficient devant le Laplacien dans l'équation de la chaleur (diffusivité thermique)

// Choix des fonctions conditions aux bords :
//	p = 1 : f = 2(y − y^2 + x − x^2), g = 0, h = 0
//	p = 2 : f = sin(x) + cos(y), g = sin(x) + cos(y), h = sin(x) + cos(y)
//	p = 3 : f = exp(-(x-lx/2)^2)*exp(-(y-ly/2)^2)*cos(pi*t/2), g = 0, h = 1
const int p = 2;


// Choix du type de condition de bord dans l'algorithme de Schwarz additif :
//	BC = 0 : Dirichlet
//	BC = 1 : Robin de coefficients a et b : a*u + b*du/dn
const int BC = 1;

// Coefficients de Robin :
const double a = 1.0;
const double b = 1.0;


// Choix du solveur :
//	SOLVER = 1 : Gradient conjugué
//	SOLVER = 2 : Bicgstab
const int SOLVER = 1;

#endif
