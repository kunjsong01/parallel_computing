#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// ***  Solution of Laplace's Equation.
// ***
// ***  Uxx + Uyy = 0
// ***  0 <= x <= pi, 0 <= y <= pi
// ***  U(x,pi) = sin(x), U(x,0) = U(0,y) = U(pi,y) = 0
// ***
// ***  then U(x,y) = (sinh(y)*sin(x)) / sinh(pi)
// ***
// ***  Should converge with: 
// ***   tol = 0.001 and N = 22  in  60 iterations.
// ***   and with tol = 0.001 and N = 102 in 200 iterations.
// ***   and with tol = 0.001 and N = 502 in 980 iterations.
// *** 

// now we try to do dynamic mem allocation to facilitate the test of various N values
//#define N 502
#define MAX(a,b)  ( ( (a)>(b) ) ? (a) : (b) )

// To paralise the code, x[][], xnew[][] and solution[][] should not be global
//double x[N][N], xnew[N][N], solution[N][N];

double calcerror(double **g, int iter, double **s, int N);
void matrix_initialise(double **x_matrix, double **solution_matrix, double h, int N);
void boundary_condition(double **x_matrix, double **xnew_matrix, int N);

int main(int argc, char *argv[]){
	int N;
	double tol=0.001, h, omega, error;
	double pi = (double)4.0*atan((double)1.0);
	int iter=0, i, j;
	double **x, **x_ptr; 
	double **xnew, **xnew_ptr;
	double **solution, **solution_ptr;

	// variables for omp timmer
	double calcerror_start; 
	double calcerror_time = 0.0;

	/* argument parsing, making N configurable to investigate the
	 * correlation between problem size and iterations to converge
	 */
	if ( argc == 2 ) {
		N = atoi(argv[1]);
		printf("Running with N = %d \n", N);
	}
	else{
		printf("Incorrect number of arguments supplied. Please give the problem size(integer) \n");
	}

	// dynamic memory allocation for x, xnew and solution
	x = malloc(N * sizeof(double *));
	xnew = malloc(N * sizeof(double *));
	solution = malloc(N * sizeof(double *));
	for (i=0; i < N; i++) {
		x[i] = (double*)malloc(N * sizeof(double));
		xnew[i] = (double*)malloc(N * sizeof(double));
		solution[i] = (double*)malloc(N * sizeof(double));
	}

	// calculate constant values, h and omega 
	h = M_PI/(double)(N-1);
	omega = 2.0/(1.0+sin(M_PI/(double)(N-1)));

	// set up pointers to allow x[][], xnew[][] and solution[][] access in calcerror
	x_ptr = (double **)malloc((unsigned) N * sizeof(double));
	for (i=0; i < N; i++) {
		x_ptr[i] = x[i];
	}

	xnew_ptr = (double **)malloc((unsigned) N * sizeof(double));
	for (i=0; i < N; i++) {
		xnew_ptr[i] = xnew[i];
	}

	solution_ptr = (double **)malloc((unsigned) N * sizeof(double));
	for (i=0; i < N; i++) {
		solution_ptr[i] = solution[i];
	}
   
	// initialise x, y
	matrix_initialise(x_ptr, solution_ptr, h, N);

	// apply boundary conditions
	boundary_condition(x_ptr, xnew_ptr, N);
	
	// start time
    calcerror_start = omp_get_wtime();

	// calculate the initial error
	error = calcerror(x_ptr, iter, solution_ptr, N);

	while(error >= tol){

		// Increase N will affect the loading of xnew into cache (L1/L2/L3)
		for(i=1; i<N-1; i++)
			for(j=1; j<N-1; j++){
				xnew[i][j] = x[i][j]+0.25*omega*(xnew[i-1][j] + xnew[i][j-1] + x[i+1][j] + x[i][j+1] - (4*x[i][j]));
			}

		for(i=1; i<N-1; i++)
			for(j=1; j<N-1; j++)
				x[i][j] = xnew[i][j];

		iter++;

		if (fmod(iter, 20) == 0)
			error = calcerror(x_ptr, iter, solution_ptr, N);
	}
	printf("Omega = %0.20f\n", omega);
	printf("Convergence in %d iterations for %dx%d grid with tolerance %f.\n", iter, N, N, tol);

	calcerror_time = calcerror_time + (omp_get_wtime() - calcerror_start);

	printf("Total iteration to converge = %d \n", iter);
	printf("Total time elapsed to converge = %f \n", calcerror_time);

	return 0;
}

void matrix_initialise(double **x_matrix, double **solution_matrix, double h, int N){
	int i, j;
	for(i=0; i<N; i++)
		x_matrix[i][N-1] = sin((double)i*h);
	for(i=0; i<N; i++)
		for(j=0; j<N-1; j++)
			x_matrix[i][j] = (double)j*h*x_matrix[i][N-1];
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			solution_matrix[i][j] = sinh((double)j*h) * sin((double)i*h)/sinh(M_PI);
}

void boundary_condition(double **x_matrix, double **xnew_matrix, int N){
	int i,j;
	// apply boundary conditions, copy the box boundaires from x to xnew matrix
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			xnew_matrix[0][j] = x_matrix[0][j];
			xnew_matrix[N-1][j] = x_matrix[N-1][j];
		}
		xnew_matrix[i][0] = x_matrix[i][0];
		xnew_matrix[i][N-1] = x_matrix[i][N-1];
	}
}

double calcerror(double **g, int iter, double **s, int N){
	int i,j;
	double error = 0.0;

	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			error = MAX(error, fabs(s[i][j] - g[i][j]));

	printf("On iteration %d error= %f\n",iter, error);
	return error;
}
