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
// ***  Should converge with
// ***   tol = 0.001 and N = 22  in  60 iterations.
// ***   and with tol = 0.001 and N = 102 in 200 iterations.
// ***   and with tol = 0.001 and N = 502 in 980 iterations.
// *** 

#define N 502
#define MAX(a,b)  ( ( (a)>(b) ) ? (a) : (b) )

// To paralise the code, x[][], xnew[][] and solution[][] should not be global
//double x[N][N], xnew[N][N], solution[N][N];

double calcerror(double **g, int iter, double **s);
void matrix_initialise(double **x_matrix, double **xnew_matrix, double **solution_matrix, double h);
void update_points(double **x_matrix, double **xnew_matrix, double omega); 

int main(int argc, char *argv[]){
	double tol=0.001, h, omega, error;
	double pi = (double)4.0*atan((double)1.0);
	int iter=0, i, j;
	double x[N][N], **x_ptr; 
	double xnew[N][N], **xnew_ptr;
	double solution[N][N], **solution_ptr;

	// variables for omp timmer
	double calcerror_start; 
	double calcerror_time = 0.0;

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
   
	// initialise x, xnew, y
	matrix_initialise(x_ptr, xnew_ptr, solution_ptr, h);

	// boundary conditions, copy the box boundaires from x to xnew array
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			xnew[0][j] = x[0][j];
			xnew[N-1][j] = x[N-1][j];
		}
		xnew[i][0] = x[i][0];
		xnew[i][N-1] = x[i][N-1];
	}

	// start time
    calcerror_start = omp_get_wtime();

	// calculate the initial error
	error = calcerror(x_ptr, iter, solution_ptr);

	while(error >= tol){
		// update all the points in x if error is not acceptable
		// Increase N will affect the loading of xnew into cache (L1/L2/L3)
		update_points(x_ptr, xnew_ptr, omega);

		iter++;

		if (fmod(iter, 20) == 0)
			error = calcerror(x_ptr, iter, solution_ptr);
	}
	printf("Omega = %0.20f\n", omega);
	printf("Convergence in %d iterations for %dx%d grid with tolerance %f.\n", iter, N, N, tol);

	calcerror_time = calcerror_time + (omp_get_wtime() - calcerror_start);

	printf("Total iteration to converge = %d \n", iter);
	printf("Total time elapsed to converge = %f \n", calcerror_time);

	return 0;
}

void matrix_initialise(double **x_matrix, double **xnew_matrix, double **solution_matrix, double h){
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

void update_points(double **x_matrix, double **xnew_matrix, double omega){

	int i, j;

	// TO-DO: wavefront array access
	
	// update points in x, put it in the xnew and copy back to x 
	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++){
			xnew_matrix[i][j] = x_matrix[i][j]+0.25*omega*(xnew_matrix[i-1][j] + xnew_matrix[i][j-1] + x_matrix[i+1][j] + x_matrix[i][j+1] - (4*x_matrix[i][j]));
		}

	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			x_matrix[i][j] = xnew_matrix[i][j];
}

double calcerror(double **g, int iter, double **s){

	int i,j;
	double error = 0.0;

	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			error = MAX(error, fabs(s[i][j] - g[i][j]));

	printf("On iteration %d error= %f\n",iter, error);
	return error;
}