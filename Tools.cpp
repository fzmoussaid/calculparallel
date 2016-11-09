
#include "Tools.hpp"


int bijection(int i, int j, int nx)
{
	return i + j*nx;
}

void charge(int n, int np, int me, int& i1, int& im)
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
}
