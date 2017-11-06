#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

// ***  Solution of Laplace's Equation.
// ***
// ***  Uxx + Uyy = 0
// ***  0 <= x <= pi, 0 <= y <= pi
// ***  U(x,pi) = sin(x), U(x,0) = U(0,y) = U(pi,y) = 0
// ***
// ***  then U(x,y) = (sinh(y)*sin(x)) / sinh(pi)
// ***
// ***  Should converge with
// ***                 tol = 0.001 and M = 22  in  60 iterations.
// ***   and with tol = 0.001 and M = 102 in 200 iterations.
// ***   and with tol = 0.001 and M = 502 in 980 iterations.
// *** 

#define N 502
#define MAX(a,b)  ( ( (a)>(b) ) ? (a) : (b) )
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

double x[N][N], xnew[N][N], solution[N][N];

double calcerror(double g[][N], int iter);


void *main(int argc, char *argv[]){
	double tol=0.001, h, omega, error;
	double pi = (double)4.0*atan((double)1.0);
	int iter=0, i, j;
	int myid, numprocs, rc, ierr;
	double total_start;
	double total_time = 0.0;

	ierr = MPI_Init(NULL, NULL); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Status status;
	printf ("Process %d of %d is alive\n", myid, numprocs);
	
	total_start = omp_get_wtime();
	h = M_PI/(double)(N-1);

	for(i=0; i<N; i++)
		x[i][N-1] = sin((double)i*h);
	

	for(i=0; i<N; i++)
		for(j=0; j<N-1; j++)
			x[i][j] = (double)j*h*x[i][N-1];

	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			solution[i][j] = sinh((double)j*h) * sin((double)i*h)/sinh(M_PI);
	
	omega = 2.0/(1.0+sin(M_PI/(double)(N-1)));
   
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			xnew[i][j] = x[i][j];
		}
	}

	error = calcerror(x, iter);

	while(error >= tol){

		//red for id =0
		if(myid==0){
			for(i=1; i<N-1; i++){
				for(j=1; j<N-1; j++){
					if((i+j)%2==0){
						xnew[i][j] = x[i][j]+0.25*omega*(x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1] - (4*x[i][j]));
					}
				}
			}
		}

		// black for id==1
		if(myid==1){
			for(i=1; i<N-1; i++){
				for(j=1; j<N-1; j++){
					if((i+j)%2==1){
						xnew[i][j] = x[i][j]+0.25*omega*(x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1] - (4*x[i][j]));
					}
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//MPI Bcast red to black
		for(i=0; i<N; i++){
				for(j=0; j<N; j++){
					MPI_Bcast(&i
					MPI_Bcast(&j
					MPI_Bcast(&i

//MPI Bcast black to red	
		

		iter++;

		if (fmod(iter, 20) == 0)
		  error = calcerror(x, iter);
		
	}
	total_time = omp_get_wtime() - total_start;
	printf("Omega = %0.20f\n", omega);
	printf("Convergence in %d iterations for %dx%d grid with tolerance %f.\n", iter, N, N, tol);
	printf("Total time to convergence: %f seconds\n", total_time);

	return 0;
}

double calcerror(double g[][N], int iter){

	int i,j;
	double error = 0.0;

	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			error = MAX(error, fabs(solution[i][j] - x[i][j]));

	printf("On iteration %d error= %f\n",iter, error);
	return error;
}
