
#include "Timestep.hpp"

using namespace std;

// ny = im -i1 + 1
void timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int nx, int ny, int recouvr, int me, int np, int i1, int im)
{
	double nr, maxnr, x;
	VectorXd Ubas(nx), Uhaut(nx), Unext(nx*ny), Rhs(nx*ny);
	MPI_Status status;
<<<<<<< HEAD
	
	Rhs.setZero();	
	
	int rc(recouvr>>1);
=======

	Rhs.setZero();

>>>>>>> 072a6480bfb4d31caf01de9228b6769bf8c9d9e5
	if(me == 0)
	{
		if(np > 1)
		{
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,ny-1,nx)) = -gamma*Un(bijection(i,ny-1,nx));
				Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
			}
		}
		else
		{
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
			}
		}
	}
	if(me == np-1)
	{
		if(np > 1)
		{
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,ny-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
				Rhs(bijection(i,0,nx)) = -gamma*Un(bijection(i,0,nx));
			}
		}
		else
		{
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,ny-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
			}
		}
		
	}
	else if(me > 0)
	{
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,ny-1,nx)) = -gamma*Un(bijection(i,ny-1,nx));
			Rhs(bijection(i,0,nx)) = -gamma*Un(bijection(i,0,nx));
		}
	}

	for (int j=0; j<ny; j++)
	{
		for (int i=1; i<nx-1; i++)
		{
			Rhs(bijection(i,j,nx)) += dt*functionF((i+1)*dx,(j+1)*dy,t,p);
		}
		Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(j+1)*dy,t,p);
		Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(j+1)*dy,t,p);
	}

	Rhs += Un/dt;
<<<<<<< HEAD
	
	cout << "Gradconj " << var.gradConj(Unext, Rhs, Niter) << endl;
	
=======

>>>>>>> 072a6480bfb4d31caf01de9228b6769bf8c9d9e5
	nr = 0.0;

	if(me == 0)
	{
		for(int i(0); i < nx; i++)
		{
			x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
			nr += x*x;
		}
	}
	if(me == np-1)
	{
		for(int i(0); i < nx; i++)
		{
			x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
			nr += x*x;
		}
	}
	else if(me > 0)
	{
		for(int i(0); i < nx; i++)
		{
			x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
			nr += x*x;

			x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
			nr += x*x;
		}
	}

	MPI_Allreduce(&nr, &maxnr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
<<<<<<< HEAD
	
	if(me == 0)
		cout << "nr = " << maxnr << endl;
	
	Un = Unext;
	
	
	while(maxnr > eps*0.01)
	//~ for(int k=0; k < 5; ++k)
	{
		Rhs.setZero();
		if(me == 0)
		{
			MPI_Sendrecv(&Un.data()[(ny-rc+recouvr%2)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					MPI_COMM_WORLD, &status);
		
		
			for(int i(0); i < nx; i++)
			{
				//~ double xi = (i+1)*dx;
				//~ double yi = (i1+ny)*dy;
				//~ Uhaut(i) = (xi*(1-xi)*yi*(1-yi))
				Rhs(bijection(i,ny-1,nx)) = -gamma*Uhaut(i);
				Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
			}
		}
		else if(me == np-1)
		{
			MPI_Sendrecv(&Un.data()[rc*nx], nx, MPI_DOUBLE, me-1, 300+me,
					Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					MPI_COMM_WORLD, &status);
		
		
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,ny-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
				Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
			}
		}
		else
		{
			MPI_Sendrecv(&Un.data()[rc*nx], nx, MPI_DOUBLE, me-1, 300+me,
					Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					MPI_COMM_WORLD, &status);
					
			MPI_Sendrecv(&Un.data()[(ny-rc+recouvr%2)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					MPI_COMM_WORLD, &status);
		
		
			for(int i(0); i < nx; i++)
			{
				Rhs(bijection(i,ny-1,nx)) = -gamma*Uhaut(i);
				Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
			}
		}
		
		for (int j=0; j<ny; j++)
		{
			for (int i=1; i<nx-1; i++)
			{
				Rhs(bijection(i,j,nx)) += dt*functionF((i+1)*dx,(j+1)*dy,t,p);
			}
			Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(j+1)*dy,t,p);
			Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(j+1)*dy,t,p);
		}

		Rhs += Un/dt;
		
		cout << "Gradconj " << var.gradConj(Unext, Rhs, Niter) << endl;
		
		nr = 0.0;
		
		if(me == 0)
		{
			for(int i(0); i < nx; i++)
			{				
				x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
				nr += x*x;
			}
		}
		else if(me == np-1)
		{
			for(int i(0); i < nx; i++)
			{
				x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
				nr += x*x;
			}
		}
		else
		{
			for(int i(0); i < nx; i++)
			{
				x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
				nr += x*x;
				
				x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
				nr += x*x;
			}
		}
		
        MPI_Allreduce(&nr, &maxnr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
		if(me == 0)
			cout << "nr = " << maxnr << endl;
        
        Un = Unext;
	}
=======
	Un = Unext;


	 while(maxnr > eps)
	 {
		 Rhs.setZero();
		 if(me == 0)
		 {
			 MPI_Sendrecv(&Un.data()[(ny-recouvr)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					 Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					 MPI_COMM_WORLD, &status);


			 for(int i(0); i < nx; i++)
			 {
				 Rhs(bijection(i,ny-1,nx)) = -gamma*Uhaut(i);
				 Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
			 }
		 }
		 if(me == np-1)
		 {
			 MPI_Sendrecv(&Un.data()[recouvr*nx], nx, MPI_DOUBLE, me-1, 300+me,
					 Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					 MPI_COMM_WORLD, &status);


			 for(int i(0); i < nx; i++)
			 {
				 Rhs(bijection(i,ny-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
				 Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
			 }
		 }
		 else
		 {
			 MPI_Sendrecv(&Un.data()[recouvr*nx], nx, MPI_DOUBLE, me-1, 300+me,
					 Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					 MPI_COMM_WORLD, &status);

			 MPI_Sendrecv(&Un.data()[(ny-recouvr)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					 Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					 MPI_COMM_WORLD, &status);


			 for(int i(0); i < nx; i++)
			 {
				 Rhs(bijection(i,ny-1,nx)) = -gamma*Uhaut(i);
				 Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
			 }
		 }

		 for (int j=0; j<ny; j++)
		 {
			 for (int i=1; i<nx-1; i++)
			 {
				 Rhs(bijection(i,j,nx)) += dt*functionF((i+1)*dx,(j+1)*dy,t,p);
			 }
			 Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(j+1)*dy,t,p);
			 Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(j+1)*dy,t,p);
		 }

		 Rhs += Un/dt;

		 var.gradConj(Unext, Rhs, Niter);

		 nr = 0.0;

		 if(me == 0)
		 {
			 for(int i(0); i < nx; i++)
			 {
				 x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
				 nr += x*x;
			 }
		 }
		 else if(me == np-1)
		 {
			 for(int i(0); i < nx; i++)
			 {
				 x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
				 nr += x*x;
			 }
		 }
		 else
		 {
			 for(int i(0); i < nx; i++)
			 {
				 x = Unext(bijection(i,0,nx)) - Un(bijection(i,0,nx));
				 nr += x*x;

				 x = Unext(bijection(i,ny-1,nx)) - Un(bijection(i,ny-1,nx));
				 nr += x*x;
			 }
		 }

         MPI_Allreduce(&nr, &maxnr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	 }
>>>>>>> 072a6480bfb4d31caf01de9228b6769bf8c9d9e5
}
