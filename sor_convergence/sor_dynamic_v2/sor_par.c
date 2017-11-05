#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/*
 * This is the parallel version of sor algorithm to solve the Laplace Equation using wavefront 
 * parallelism, i.e. accessing the array in diagonal strips to obviate the dependency when calculating
 * U(k+1).
 * The main strategy is to:
 * 		- Separate the iterative computation blocks from the main routine to
 * 		  make them easier for parallelisation
 * 		- Apply OpenMP parallelisation techniques in the iteration space and data space 
 * The mathematical problem it's trying to solve is explained below. 
 *
 *  Solution of Laplace's Equation.
 *
 *  Uxx + Uyy = 0
 *  0 <= x <= pi, 0 <= y <= pi
 *  U(x,pi) = sin(x), U(x,0) = U(0,y) = U(pi,y) = 0
 *
 *  then U(x,y) = (sinh(y)*sin(x)) / sinh(pi)
 *
 *  Should converge with: 
 *   tol = 0.001 and N = 22  in  60 iterations.
 *   and with tol = 0.001 and N = 102 in 200 iterations.
 *   and with tol = 0.001 and N = 502 in 980 iterations.
 * 
 *  Author: Kunjian Song
 *  Date: 04/11/2017
 */

// now we try to do dynamic mem allocation to facilitate the test of various N values
//#define N 502
#define MAX(a,b)  ( ( (a)>(b) ) ? (a) : (b) )

// To paralise the code, x[][], xnew[][] and solution[][] should not be global
//double x[N][N], xnew[N][N], solution[N][N];

double calcerror(double **x_matrix, int iter, double **solution_matrix, int N);
void matrix_initialise(double **x_matrix, double **xnew_matrix, double **solution_matrix, double h, int N);
void update_points(double **x_matrix, double **xnew_matrix, double omega, int N); 

int main(int argc, char *argv[]){
	int N;
	double tol=0.001, h, omega, error;
	double pi = (double)4.0*atan((double)1.0);
	int iter=0, i, j;
	// this "pointer to pointer" logic causes lots of troubles during the experiment, please see the
	// comments in update_points function.
	// Actually, x_ptr is unnecessary. Just use x. It would do the same work!
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
	matrix_initialise(x_ptr, xnew_ptr, solution_ptr, h, N);
	
	// start time
    calcerror_start = omp_get_wtime();

	// calculate the initial error
	error = calcerror(x_ptr, iter, solution_ptr, N);

	while(error >= tol){
		// update all the points in x if error is not acceptable
		// Increase N will affect the loading of xnew into cache (L1/L2/L3)
		update_points(x_ptr, xnew_ptr, omega, N);
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

void matrix_initialise(double **x_matrix, double **xnew_matrix, double **solution_matrix, double h, int N){
	int i, j;

	for(i=0; i<N; i++)
		x_matrix[i][N-1] = sin((double)i*h);
	for(i=0; i<N; i++)
		for(j=0; j<N-1; j++)
			x_matrix[i][j] = (double)j*h*x_matrix[i][N-1];
	for(i=0; i<N; i++)
		for(j=0; j<N; j++)
			solution_matrix[i][j] = sinh((double)j*h) * sin((double)i*h)/sinh(M_PI);

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

void update_points(double **x_matrix, double **xnew_matrix, double omega, int N){

	int i, j;

#if 0	
	/*
	 * The "pointer to pointer" logic here may hamper the speed lot. This was originally
	 * made to allow accesses to those 2-D matrices. However, experiments showed that 
	 * it will become very slow even trying to print out the value of x_matrix after
	 * calculating the xnew_matrix. 
	 * So let's switch back to the static array.
	 * Also: this program does not work with N=30. It will not terminate even though the error is less than tol.
	 */
	// a pointer to each diagonal strip, max length is N (the diagonal of an N*N matrix)
	// has to be a pointer to pointer(the starting point of that diagonal strip)
	double **dstrip_ptr;
	int n_task = 0; // number of non-boundary elements in the diagonal strip
	int counter; // counter to ignore boundary
	int dlength; // length of elements in the diagonal strip
	
	// An N*N square matrix has (2*n-1) diagnonal strips
	for(i=0; i<(2*N-1); i++) {
		counter = 1;
		printf("Slice %d:\n ", i);
		int shift = (i < N)? 0 : (i - N + 1);
		dlength = (i-shift) - shift + 1;
		printf("\t @@dlength: %d", dlength);
		for(j=shift; j<= (i-shift); j++){
			// ignore the boundary
			if ( (counter != 1) || (counter != dlength) ){
				xnew_matrix[i][j] = x_matrix[i][j]+0.25*omega*(xnew_matrix[i-1][j] + xnew_matrix[i][j-1] + x_matrix[i+1][j] + x_matrix[i][j+1] - (4*x_matrix[i][j]));
			}
			printf("%f ", x_matrix[j][i-j]);// This line is very time consuming!!!
			counter++;
		}
		printf("\n");
	}
#endif
	// update points in x, put it in the xnew and copy back to x 
	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++){
			xnew_matrix[i][j] = x_matrix[i][j]+0.25*omega*(xnew_matrix[i-1][j] + xnew_matrix[i][j-1] + x_matrix[i+1][j] + x_matrix[i][j+1] - (4*x_matrix[i][j]));
		}
	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			x_matrix[i][j] = xnew_matrix[i][j];
} 

double calcerror(double **x_matrix, int iter, double **solution_matrix, int N){
	int i,j;
	double error = 0.0;

	// TO-DO: implement parallel omp reduce
	for(i=1; i<N-1; i++)
		for(j=1; j<N-1; j++)
			error = MAX(error, fabs(solution_matrix[i][j] - x_matrix[i][j]));

	printf("On iteration %d error= %f\n",iter, error);
	return error;
}
