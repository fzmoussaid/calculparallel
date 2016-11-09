
void timeStep(const SolverCG& var, VectorXd& Un, double eps, int nx, int ny, int recouvr, int me, int np, int i1, int im)
{
	double nr, maxnr;
	VectorXd Ubas(nx), Uhaut(nx), Unext(nx*ny);
	MPI_Status status;
	
	do
	{
		if(me == 0)
		{
			MPI_Sendrecv(&Un.data()[(ny-recouvr)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					MPI_COMM_WORLD, &status);
		}
		else if(me == np-1)
		{
			MPI_Sendrecv(&Un.data()[recouvr*nx], nx, MPI_DOUBLE, me-1, 300+me,
					Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Sendrecv(&Un.data()[recouvr*nx], nx, MPI_DOUBLE, me-1, 300+me,
					Ubas.data(), nx, MPI_DOUBLE, me-1, 200+me-1,
					MPI_COMM_WORLD, &status);
					
			MPI_Sendrecv(&Un.data()[(ny-recouvr)*nx], nx, MPI_DOUBLE, me+1, 200+me,
					Uhaut.data(), nx, MPI_DOUBLE, me+1, 300+me+1,
					MPI_COMM_WORLD, &status);
		}
		
		
                
        MPI_Allreduce(&nr, &maxnr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	} while(maxnr > eps);
}
