
#ifndef PARAMETER_HPP
#define PARAMETER_HPP

const double eps = 1.e-12;
const int nx=40, ny=40;
const double lx=1., ly=1., dt=10.;
const double dx=lx/(nx+1), dy=ly/(ny+1);
const int nl = nx*ny;
const int nIterMax = 10000;
const double pi = 3.1415926535897932384;
const int p=1;
const int D = 1;
const int BC = 2;
const int SOLVER = 0;
const double a = 1.0;
const double b = 1.0;

#endif
