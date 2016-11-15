//#include <mpi.h>
#include <iostream>
#include <fstream>
#include <Eigen>
#include "Parameters.cpp"
#include "Tools.cpp"
#include "Functions.cpp"
#include "SolverCG.cpp"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	//int me, np, i1, im, p=1;
	double alpha = 1/dt + 2*D*(1/(dx*dx) + 1/(dy*dy)), beta = -D/(dx*dx), gamma = -D/(dy*dy);

	//MPI_Init (&argc, &argv);

	//MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//MPI_Comm_size(MPI_COMM_WORLD, &np);

	//charge(ny, np, me, i1, im);

  VectorXd U(nl),U0(nl),X,RHS(nl),RHS0(nl);
  U0.setZero();

  SolverCG solver(alpha,beta,gamma,eps,nx,ny);
	cout << "eps : " << eps << endl;
  for (int time=0; time<1; time++)
  {
     //! construction de RHS
     RHS0.setZero();
     for (int j=0; j<ny; j++)
     {
       for (int i=1; i<nx-1; i++)
       {
         RHS0(bijection(i,j,nx)) = functionF(i*dx,j*dy,time*dt,p);
       }
       RHS0(bijection(0,j,nx)) -=  beta*functionH(0.,j*dy,time*dt,p);
       RHS0(bijection(nx-1,j,nx)) -= beta*functionH(lx,j*dy,time*dt,p);
		 }
		 for (int i=0; i<nx; i++)
		 {
			 RHS0(bijection(i,0,nx)) -= gamma*functionG(i*dx,0.,time*dt,p);
			 RHS0(bijection(i,ny-1,nx)) -= gamma*functionG(i*dx,ly,time*dt,p);
		 }

		 RHS = RHS0 + U0/dt;
		 //! gradient conjugue adapte
		 cout << "iter : " << solver.gradConj(U,RHS,nIterMax) << endl;
		 U0 = U;

	 }

  ofstream file_sol("Sol.dat"), file_sol_exact("Sol_exacte.dat");
	double norm_diff=0., norm_exact=0.;
	for (int j=0; j<ny; j++)
	{
		for (int i=0; i<nx; i++)
		{
			double xi = i*dx;
			double yi = j*dy;
			file_sol << xi << " " << yi << " " << U(bijection(i,j,nx)) << endl;

			if (p==2)
			{
				file_sol_exact << xi << " " << yi << " " << sin(xi)+cos(yi) << endl;
				norm_diff = max(abs(U(bijection(i,j,nx)) - (sin(xi)+cos(yi))),norm_diff);
				norm_exact = max(abs(sin(xi)+cos(yi)),norm_exact);
			}
			else if (p==1)
			{
				file_sol_exact << xi << " " << yi << " " << xi*(1-xi)*yi*(1-yi) << endl;
				norm_diff = max(abs(U(bijection(i,j,nx)) - (xi*(1-xi)*yi*(1-yi))),norm_diff);
				norm_exact = max(abs(xi*(1-xi)*yi*(1-yi)),norm_exact);
			}
		}
		file_sol << endl;
		file_sol_exact << endl;
	}
	if ((p==1)||(p==2))
	{
		cout << "erreur relative en norme infinie : " << norm_diff/norm_exact << endl;
	}

	file_sol.close();
	file_sol_exact.close();
	//MPI_Finalize();



	return 0;
}


/*
  !write(*,*) "Me ",Me," connait ",U
  call Rename(Me,name)
  open(unit=11,file='Sol_exacte.dat',form='formatted',status='unknown',action='write')
  open(unit=12+Me,file=name,form='formatted',status='unknown',action='write')
  norme_diff = 0.
  norme_exacte = 0.
  do i=i1,iN
     do j=1,Nx
        xi = dx*j
        yi = dy*i
        write(12+Me,*) xi,yi,U((i-1)*Nx+j)
        write(12+Me,*)
        if (p==2) then
           norme_diff = max(abs(U((i-1)*Nx+j)-(sin(xi)+cos(yi))),norme_diff)
           norme_exacte =  max(abs(sin(xi)+cos(yi)),norme_exacte)
        else if (p==1) then
           norme_diff = max(abs(U((i-1)*Nx+j)-(xi*(1-xi)*yi*(1-yi))),norme_diff)
           norme_exacte = max(xi*(1-xi)*yi*(1-yi),norme_exacte)
        end if
     end do
  end do
  if (Me==0) then
     do i=1,Ny
        do j=1,Nx
           xi = dx*j
           yi = dy*i
           if (p==2) then
              write(11,*) xi,yi,sin(xi)+cos(yi)
           else if (p==1) then
              write(11,*) xi,yi,xi*(1-xi)*yi*(1-yi)
           end if
        end do
        write(11,*)
     end do
  end if

  if ((p==1).or.(p==2)) then
     call MPI_ALLREDUCE(norme_diff,norme_diff_tot,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,statinfo)
     call MPI_ALLREDUCE(norme_exacte,norme_exacte_tot,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,statinfo)
     if (Me==0) then
        write (*,*) "Erreur relative en norme infinie : ",norme_diff_tot/norme_exacte_tot
     end if
  end if
  !stop
  close(11)
  close(12+Me)
  deallocate(U,U0,RHS,RHS0)

  t2 = MPI_WTIME()
  call MPI_ALLREDUCE(t2-t1,t_total,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,statinfo)
  if (Me==0) then
     write(*,*)"Temps total (s) : ",t_total
  end if

  call MPI_FINALIZE(statinfo)

end program
*/
