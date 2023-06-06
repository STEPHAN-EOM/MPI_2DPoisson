/**
 * @file: main.c
 * @author: Chanho Eom
 * @date: 8-Apr-2023
 * */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"
#include "function.h"
#include "decomp1d.h"

#define maxit 10000

int main(int argc, char **argv)
{
double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
int nx, ny;
int myid, nprocs;

int nbrleft, nbrright;
int s, e, it;
double glob_diff;
double glob_grid_diff;
double ldiff;
double t1, t2;
double tol=1.0E-11;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &myid);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

MPI_Comm comm;

int dims[1], periods[1];
dims[0] = nprocs;
periods[0] = 0;


if( myid == 0 ){
if(argc > 2){
fprintf(stderr,"---->Usage: mpirun -np <nproc> %s <nx>\n",argv[0]);
fprintf(stderr,"---->(for this code nx=ny)\n");
MPI_Abort(MPI_COMM_WORLD, 1);}
if( argc == 2 ){
nx = atoi(argv[1]);}
if( argc == 1 ){
nx=15;}

if( nx > maxn-2 ){
fprintf(stderr,"grid size too large\n");
exit(1);}

}


MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
printf("(myid: %d) nx = %d\n", myid, nx);

ny = nx;

init_full_grids(a, b, f);

// Question.1 => MPI_Cart_create/shift by C.Eom
MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &comm);
MPI_Cart_shift(comm, 0, 1, &nbrleft, &nbrright);

// Question.2 => 1D Decomposition by C.Eom
decomp_1d(nx, myid, nprocs, &s, &e);


printf("(myid: %d) nx: %d s: %d; e: %d; nbrleft: %d; nbrright: %d\n",myid, nx , s, e, nbrleft, nbrright);

// Original Given Initial Condition
//onedinit_basic(a, b, f, nx, ny, s, e);

// Question.2 => Poisson Initial Condition by C.Eom
init_basic_1d(a, b, f, nx, ny, s, e);

printf("==================== Initial Condition of grid ====================\n");
print_in_order(a, MPI_COMM_WORLD, nx);
if( nprocs == 1  ){
print_grid_to_file("grid", a,  nx, ny);
print_full_grid(a, nx);
}

t1 = MPI_Wtime();

glob_diff = 1000;
for(it=0; it<maxit; it++){

exchang3(a, ny, s, e, MPI_COMM_WORLD, nbrleft, nbrright);
sweep1d(a, f, nx, s, e, b);

exchang3(b, nx, s, e, MPI_COMM_WORLD, nbrleft, nbrright);
sweep1d(b, f, nx, s, e, a);

ldiff = griddiff(a, b, nx, s, e);
MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(myid==0 && it%10==0){
printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
}

if( glob_diff < tol ){      
if(myid==0){  	
printf("iterative solve converged\n");
      }      
break;
    } 
}
    
t2=MPI_Wtime();  
printf("DONE! (it: %d)\n",it);


//MPI_Barrier(MPI_COMM_WORLD);

if( myid == 0 ){    
if( it == maxit ){      
fprintf(stderr,"Failed to converge\n");    
}    
printf("Run took %lf s\n",t2-t1);  
}


// Question.2 => Define Analytic solution by C.Eom
double solution[maxn][maxn];
solution_grid(solution, nx, ny);

if (myid == 0){
printf("\n================== Show me the analytic solution ==================\n");
print_full_grid(solution, nx);}

// Question.2 => Compare the results from finite difference method and analytic solution
ldiff = griddiff(solution, a, nx, s, e);
MPI_Allreduce(&ldiff, &glob_grid_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (myid == 0){
printf("\nThe global difference to the analytic solution on a grid size of %d is %f\n", nx, glob_grid_diff);
}

print_in_order(a, MPI_COMM_WORLD, nx);  
if( nprocs == 1  ){    
print_grid_to_file("grid", a,  nx, ny);    
print_full_grid(a, nx);  
}

// Question.3 => Gather the grid from the processors onto rank 0
gather_grid(a, nx, nprocs, myid, MPI_COMM_WORLD);
if (myid == 0){
printf("\n================== The gathered grid onto rank 0 ==================\n");
print_full_grid(a, nx);
printf("\n================== Analytic Solution ====================\n");
print_full_grid(solution, nx);
}

// Question.3 => Write the grid of each processor onto the file
write_grid(a, nx, myid, s, e);
if (myid == 0){
char filename[50] = "NumericalApproximation_1Dpoission.txt";
print_grid_to_file(filename, a, nx, ny);
}
 

MPI_Finalize();
  
return 0;
}

