//#include <mpi.h>
#include <Eigen>
#include "Parameters.cpp"
#include "Tools.hpp"
#include "Functions.cpp"
#include "SolverCG.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	//int me, np, i1, im, p=1;
	double alpha = 1/dt + 2*D*(1/(dx*dx) + 1/(dy*dy)), beta = -D/(dx*dx), gamma = -D/(dy*dy);
  double xi, yi, t_total,norme_diff,norme_exacte,norme_diff_tot,norme_exacte_tot;
  string name;
  int iter=0;

	//MPI_Init (&argc, &argv);

	//MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//MPI_Comm_size(MPI_COMM_WORLD, &np);

	//charge(ny, np, me, i1, im);

  VectorXd U(nl),U0(nl),X,RHS(nl),RHS0(nl);
  U0.setZero();

  SolverCG solver(alpha,beta,gamma,eps,nx,ny);

  for (int time=1; time<10; time++)
  {
     //! construction de RHS
     RHS0.setZero();
     for (int j=0; j<ny; j++)
     {
       for (int i=0; i<nx; i++)
       {
         RHS0(bijectionectionection(i,j,nx)) = f((i+1)*dx,(j+1)*dy,(time-1)*dt,p);
       }
       RHS0(bijectionectionection(i,j,nx)) -=  beta*h(0.,(i+1)*dy,(time-1)*dt,p);
       RHS0(bijectionectionection(0,j,nx)) -= beta*h(lx,(i+1)*dy,(time-1)*dt,p);
       for (int i=0; i<nx; i++)
       {
         if (j==0)
         {
           RHS0(i) -= gamma*g((j+1)*dx,0.,(time-1)*dt,p);
         }
         else if (j==ny-1)
         {
           RHS0(bijectionectionection(i,j,nx)) -= gamma*g((j+1)*dx,ly,(time-1)*dt,p);
         }
       }
       RHS(bijectionectionection(i,j,nx)) = RHS0(bijectionectionection(i,j,nx)) + U0(bijectionectionection(i,j,nx))/dt;
     }

     //! gradient conjugue adapte
     iter = solver.gradConj(U,RHS,nIterMax);
     U0 = U
}

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
