
#include "Timestep.hpp"

// nyLocal = im -i1 + 1
int timeStep(const SolverCG& var, VectorXd& Un, double eps, double beta, double gamma, double t, int Niter, int nx, int ny, int recouvr, int me, int np, int i1, int im)
{
	double nr, nr2, maxnr, maxnr2, x;
	int nyLocal(im-i1+1);
	VectorXd Ubas(nx), Ubas_robin(nx), Uhaut(nx), Uhaut_robin(nx), UbasNext(nx), UbasNext_robin(nx), UhautNext(nx), UhautNext_robin(nx), Unext(nx*nyLocal), Rhs(nx*nyLocal);
	MPI_Status status;

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

		var.gradConj(Unext, Rhs, nIterMax, me, np);

		Un = Unext;
	}
	else
	{
		if(me == 0)
		{
			if (BC==0)
			{
				MPI_Sendrecv(&Un.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
						Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}
			else if (BC==1)
			{
				for (int i=0; i<nx; i++)
				{
					Ubas_robin(i) = a*Un((nyLocal-(recouvr+1))*nx+i) + b*(Un((nyLocal-(recouvr+1))*nx+i) - Un((nyLocal-(recouvr+1))*nx+i+nx))/dy;
				}

				MPI_Sendrecv(Ubas_robin.data(), nx, MPI_DOUBLE, me+1, 200+me,
						Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}
		}
		else if(me == np-1)
		{
			if (BC==0)
			{
				MPI_Sendrecv(&Un.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
						Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);
			}
			else if (BC==1)
			{
				for (int i=0; i<nx; i++)
				{
					Uhaut_robin(i) = a*Un(recouvr*nx+i) + b*(Un(recouvr*nx+i) - Un(recouvr*nx+i-nx))/dy;
				}

				MPI_Sendrecv(Uhaut_robin.data(), nx, MPI_DOUBLE, me-1, 300+me,
						Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);
			}
		}
		else
		{
			if (BC==0)
			{
				MPI_Sendrecv(&Un.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
						Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);

				MPI_Sendrecv(&Un.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
						Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}
			else if (BC==1)
			{
				for (int i=0; i<nx; i++)
				{
					Ubas_robin(i) = a*Un((nyLocal-(recouvr+1))*nx+i) + b*(Un((nyLocal-(recouvr+1))*nx+i) - Un((nyLocal-(recouvr+1))*nx+i+nx))/dy;
					Uhaut_robin(i) = a*Un(recouvr*nx+i) + b*(Un(recouvr*nx+i) - Un(recouvr*nx+i-nx))/dy;
				}

				MPI_Sendrecv(Uhaut_robin.data(), nx, MPI_DOUBLE, me-1, 300+me,
						Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
						MPI_COMM_WORLD, &status);

				MPI_Sendrecv(Ubas_robin.data(), nx, MPI_DOUBLE, me+1, 200+me,
						Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
						MPI_COMM_WORLD, &status);
			}
		}
		do
		{
			Rhs.setZero();
			if(me == 0)
			{
				for(int i(0); i < nx; i++)
				{
					if (BC==0)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*Uhaut(i);
						Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
					}
					else if (BC==1)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*Uhaut(i)*(dy/(b+a*dy));
						Rhs(bijection(i,0,nx)) = -gamma*functionG((i+1)*dx,0.,t,p);
					}
				}
			}
			else if(me == np-1)
			{
				if (BC==0)
				{
					for(int i(0); i < nx; i++)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
						Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
					}
				}
				else if (BC==1)
				{
					for(int i(0); i < nx; i++)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*functionG((i+1)*dx,ly,t,p);
						Rhs(bijection(i,0,nx)) = -gamma*Ubas(i)*dy/(b+a*dy);
					}
				}
			}
			else
			{
				if (BC==0)
				{
					for(int i(0); i < nx; i++)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*Uhaut(i);
						Rhs(bijection(i,0,nx)) = -gamma*Ubas(i);
					}
				}
				else if (BC==1)
				{
					for(int i(0); i < nx; i++)
					{
						Rhs(bijection(i,nyLocal-1,nx)) = -gamma*Uhaut(i)*dy/(b+a*dy);
						Rhs(bijection(i,0,nx)) = -gamma*Ubas(i)*dy/(b+a*dy);
					}
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

			var.gradConj(Unext, Rhs, nIterMax, me, np);

			if(me == 0)
			{
				if (BC==0)
				{
					MPI_Sendrecv(&Unext.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
							UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
							MPI_COMM_WORLD, &status);
				}
				else if (BC==1)
				{
					for (int i=0; i<nx; i++)
					{
						UbasNext_robin(i) = a*Unext((nyLocal-(recouvr+1))*nx+i) + b*(Unext((nyLocal-(recouvr+1))*nx+i) - Unext((nyLocal-(recouvr+1))*nx+i+nx))/dy;
					}

					MPI_Sendrecv(UbasNext_robin.data(), nx, MPI_DOUBLE, me+1, 200+me,
							UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
							MPI_COMM_WORLD, &status);
				}
			}
			else if(me == np-1)
			{
				if (BC==0)
				{
					MPI_Sendrecv(&Unext.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
							UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
							MPI_COMM_WORLD, &status);
				}
				else if (BC==1)
				{
					for (int i=0; i<nx; i++)
					{
						UhautNext_robin(i) = a*Unext(recouvr*nx+i) + b*(Unext(recouvr*nx+i) - Unext(recouvr*nx+i-nx))/dy;
					}

					MPI_Sendrecv(UhautNext_robin.data(), nx, MPI_DOUBLE, me-1, 300+me,
							UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
							MPI_COMM_WORLD, &status);
				}
			}
			else
			{
				if (BC==0)
				{
					MPI_Sendrecv(&Unext.data()[(recouvr)*nx], nx, MPI_DOUBLE, me-1, 300+me,
							UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
							MPI_COMM_WORLD, &status);

					MPI_Sendrecv(&Unext.data()[(nyLocal-(recouvr+1))*nx], nx, MPI_DOUBLE, me+1, 200+me,
							UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
							MPI_COMM_WORLD, &status);
				}
				else if (BC==1)
				{
					for (int i=0; i<nx; i++)
					{
						UbasNext_robin(i) = a*Unext((nyLocal-(recouvr+1))*nx+i) + b*(Unext((nyLocal-(recouvr+1))*nx+i) - Unext((nyLocal-(recouvr+1))*nx+i+nx))/dy;
						UhautNext_robin(i) = a*Unext(recouvr*nx+i) + b*(Unext(recouvr*nx+i) - Unext(recouvr*nx+i-nx))/dy;
					}

					MPI_Sendrecv(UhautNext_robin.data(), nx, MPI_DOUBLE, me-1, 300+me,
							UbasNext.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
							MPI_COMM_WORLD, &status);

					MPI_Sendrecv(UbasNext_robin.data(), nx, MPI_DOUBLE, me+1, 200+me,
							UhautNext.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
							MPI_COMM_WORLD, &status);
				}
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
		} while(maxnr > (eps*maxnr2) && i < Niter);
	}

	return i;
}
