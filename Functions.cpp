double functionF(double x, double y, double t, int n)
{
  double f = 0.;
  if (n==1)
  {
    f = 2*(y-y*y+x-x*x);
  }
  else if (n==2)
  {
    f = sin(x) + cos(y);
  }
  else if (n==3)
  {
    f = exp(-(x-lx/2)*(x-lx/2))*exp(-(y-ly/2)*(y-ly/2))*cos(pi/2*t);
  }
}

double functionG(double x, double y, double t, int n)
{
  double g = 0.;
  if (n==1)
  {
       g = 0.;
  }
  else if (n==2)
  {
       g = sin(x) + cos(y);
  }
}

double functionH(double x, double y, double t, int n)
{
  double h = 0.;
  if (n==1)
  {
       h = 0;
  }
  else if (n==2)
  {
       h = sin(x) + cos(y);
  }
  else if (n==3)
  {
       h = 1.;
  }
}
