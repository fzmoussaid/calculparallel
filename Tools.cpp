
#include "Tools.hpp"


int bijection(int i, int j, int nx)
{
	return i + j*nx;
}

void charge(int n, int np, int me, int recouvr, int& i1, int& im)
{
	int r = n%np;
	if (me <= r-1)
	{
		  i1 = me*(n/np) + me;
		  im = i1 + n/np;
	}
	else
	{
		  i1 = me*(n/np) + r;
		  im = i1 + n/np  -1;
	}
	
	
	int rc = recouvr>>1;
	if(me > 0)
	{
		i1 -= rc;
	}
	if(me < np-1)
	{
		im += rc+recouvr%2;
	}
}

void relativeNormL2(const VectorXd& v, const VectorXd& vref, double& norm, double& normRef)
{
	double x;
	norm = 0.0;
	normRef = 0.0;
	for(int i(0); i < v.size(); ++i) {
		x = v(i) - vref(i);
		norm += x*x;
		
		x = vref(i);
		normRef += x*x;
	}
}
