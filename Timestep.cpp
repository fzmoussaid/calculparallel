
#include "Timestep.hpp"

// nyLocal = im -i1 + 1
int timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int nx, int ny, int recouvr, int me, int np, int i1, int im)
{
	double nr, nr2, maxnr, maxnr2, x;
	int nyLocal(im-i1+1);
	VectorXd Ubas(nx), Uhaut(nx), UbasNext(nx), UhautNext(nx), Unext(nx*nyLocal), Rhs(nx*nyLocal);
	MPI_Status status;

	Rhs.setZero();
	Ubas.setZero();
	Uhaut.setZero();
	UbasNext.setZero();
	UhautNext.setZero();

	int rc(recouvr>>1);

	int i(0);

	if(np == 1)
	{
		Rhs.setZero();
		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
		}

		for(int i(0); i < nx; i++)
		{
			Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
		}

		for (int j=0; j<ny; j++)
		{
			for (int i=1; i<nx-1; i++)
			{
				Rhs(bijection(i,j,nx)) += functionF((i+1)*dx,(i1+j+1)*dy,t,p);
			}
			Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(i1+j+1)*dy,t,p);
			Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(i1+j+1)*dy,t,p);
		}

		Rhs += Un/dt;

		var.gradConj(Unext, Rhs, Niter);

		Un = Unext;
	}
	else
	{
		do
		{
			Rhs.setZero();
			if(me == 0)
			{
				for(int i(0); i < nx; i++)
				{
					Rhs(bijection(i,nyLocal-1,nx)) = -gamma*UhautNext(i);
					Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
				}
			}
			else if(me == np-1)
			{
				for(int i(0); i < nx; i++)
				{
					Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
					Rhs(bijection(i,0,nx)) = -gamma*UbasNext(i);
				}
			}
			else
			{
				for(int i(0); i < nx; i++)
				{
					Rhs(bijection(i,nyLocal-1,nx)) = -gamma*UhautNext(i);
					Rhs(bijection(i,0,nx)) = -gamma*UbasNext(i);
				}
			}

			for (int j=0; j<nyLocal; j++)
			{
				for (int i=1; i<nx-1; i++)
				{
					Rhs(bijection(i,j,nx)) += functionF((i+1)*dx,(i1+j+1)*dy,t,p);
				}
				Rhs(bijection(0,j,nx)) -=  beta*functionH(0.,(i1+j+1)*dy,t,p);
				Rhs(bijection(nx-1,j,nx)) -= beta*functionH(lx,(i1+j+1)*dy,t,p);
			}

			Rhs += Un/dt;

			var.gradConj(Unext, Rhs, Niter);

			if(me == 0)
			{
				MPI_Sendrecv(&Unext.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
						UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}
			else if(me == np-1)
			{
				MPI_Sendrecv(&Unext.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
						UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);
			}
			else
			{
				MPI_Sendrecv(&Unext.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
						UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);

				MPI_Sendrecv(&Unext.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
						UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}

			nr = 0.0;
			nr2 = 0.0;

			if(me == 0)
			{
				for(int i(0); i < nx; i++)
				{
					x = UhautNext(i) - Uhaut(i);
					nr += x*x;

					x = Uhaut(i);
					nr2 += x*x;
				}
			}
			else if(me == np-1)
			{
				for(int i(0); i < nx; i++)
				{
					x = UbasNext(i) - Ubas(i);
					nr += x*x;

					x = Ubas(i);
					nr2 += x*x;
				}
			}
			else
			{
				for(int i(0); i < nx; i++)
				{
					x = UhautNext(i) - Uhaut(i);
					nr += x*x;

					x = UbasNext(i) - Ubas(i);
					nr += x*x;

					x = Uhaut(i);
					nr2 += x*x;

					x = Ubas(i);
					nr2 += x*x;
				}
			}

			MPI_Allreduce(&nr, &maxnr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&nr2, &maxnr2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			Un = Unext;
			Uhaut = UhautNext;
			Ubas = UbasNext;
			i++;
		} while(maxnr > (eps*maxnr2));
	}


	//~ ofstream file_sol("sol/Sol" + to_string(me) + ".dat");
	//~ if(me > 0) {
		//~ for (int i=0; i<nx; i++)
		//~ {
			//~ double xi = (i+1)*dx;
			//~ double yi = (i1)*dy;
			//~ file_sol << xi << " " << yi << " " << Ubas(i) << endl;
		//~ }
		//~ file_sol << endl;
	//~ }
	//~ for (int j=0; j<(im-i1+1); j++)
	//~ {
		//~ for (int i=0; i<nx; i++)
		//~ {
			//~ double xi = (i+1)*dx;
			//~ double yi = (i1+j+1)*dy;
			//~ file_sol << xi << " " << yi << " " << Un(bijection(i,j,nx)) << endl;
		//~ }
		//~ file_sol << endl;
	//~ }
	//~ if(me < np-1) {
		//~ for (int i=0; i<nx; i++)
		//~ {
			//~ double xi = (i+1)*dx;
			//~ double yi = (im+2)*dy;
			//~ file_sol << xi << " " << yi << " " << Uhaut(i) << endl;
		//~ }
		//~ file_sol << endl;
	//~ }
//~
	//~ file_sol.close();

	return i;
}
