#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include<mpi.h>

void vector_update(double a[], double b[], double x[], int limit);
void matrix_update(double **l, double y[], int limit);
void vector_add(double a[], double b[], double x[], int limit);
double vec_sum(double x[], int limit); 
void trimv(double **l, double y[], double a[], int limit); 

/* SBLAS code */

void *main() {

  int n = 1024*60;
  int m = 512;
  int nsteps = 1000;
  int k, i;
  MPI_Status status;

  //double ll[m][m], **l;
  double **l;
  double a[n], b[n], x[n], y[m];

	// create a global sum
	double global_sum; 
  
  double res;
  double vector_start, sum_start, trimv_start, vec_update_start, matrix_update_start, total_start;

  double vec_update_time = 0.0;
  double vector_time = 0.0;
  double sum_time = 0.0;
  double matrix_update_time = 0.0;
  double trimv_time = 0.0;
  double total_time = 0.0;

  // MPI specific data
  int myid, numprocs, rc, ierr;
  int nx;

  // MPI initialisation

  ierr = MPI_Init(NULL, NULL); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  printf ("Process %d of %d is alive\n", myid, numprocs);
  
  // allocate (lower triangular) matrix as a contibuous block and set up pointers to allow l[][] access 
  l = (double **) malloc((unsigned) m*sizeof(double*));  
  for (i=0; i <= m-1; i++) {
    l[i] = (double *)malloc((unsigned) m * sizeof(double));
    //l[i] = ll[i];
  }

  // Loop over whole program to get reasonable execution time
  // (simulates a time-steping code)

  if (myid == 0) total_start = omp_get_wtime();

  for (k=1; k<=nsteps; k++ )  {
    //vec_update_start = omp_get_wtime();
    if (myid == 0) vector_update(a, b, x, n);
    //vec_update_time = vec_update_time + (omp_get_wtime() - vec_update_start);

    // Work out size of data on this processor (not a general solution!)

    nx = n/numprocs;

    vector_start = omp_get_wtime();
    // Distribute the vectors a and b to other processors

	// MASTER
	// send parts of the vectors a and b to other processors (x not reqd)
	ierr = MPI_Scatter(a, nx, MPI_DOUBLE_PRECISION, a, nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	ierr = MPI_Scatter(b, nx, MPI_DOUBLE_PRECISION, b, nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

    // vector addition operation: each processor adds its own parts 
    // Note: only master time is collected...

    // vector addition operation

    vector_add(a, b, x, nx);

    // Gather result vector on processor 0

	// MASTER
	// receive parts of result into the appropriate parts of my x

	ierr = MPI_Gather(x, nx, MPI_DOUBLE_PRECISION,	\
			x, nx, MPI_DOUBLE_PRECISION, \
			0, MPI_COMM_WORLD);

    // Rest of the code is performed on processor 0
    vector_time = vector_time + (omp_get_wtime() - vector_start);


	// sum operation
	
	// MASTER
	// send parts of the vectors x to other processes
	ierr = MPI_Scatter(x, nx, MPI_DOUBLE_PRECISION, x, nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

	sum_start = omp_get_wtime();
	res = vec_sum(x, nx);

	ierr = MPI_Reduce(&res, &global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

	sum_time = sum_time +(omp_get_wtime() - sum_start);

    if (myid == 0) {
      // 'Update' lower triangular matrix and y vector

      matrix_update_start = omp_get_wtime();
      matrix_update(l, y, m);
      matrix_update_time = matrix_update_time +(omp_get_wtime() - matrix_update_start);

      // triangular matrix times vector operation

      trimv_start = omp_get_wtime();
      trimv(l, y, a, m);
      trimv_time = trimv_time + (omp_get_wtime() - trimv_start);

    }

    // end time-step loop
  }
  
  if (myid == 0)  total_time = omp_get_wtime() - total_start;
  
  for (i=0;i<numprocs;i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == i) {
      if (myid == 0) { //MASTER
	printf("MASTER output\n");
	//printf( "vector update time = %f seconds \n \n", vec_update_time);
	printf( "x[0] = %f \n", x[0]);
	printf( "x[n-1] = %f \n", x[n-1]);
	printf( "vector add time = %f seconds \n \n", vector_time);
	printf( "Sum of x() = %f \n", global_sum);
	printf( "sum time = %f \n \n", sum_time);
	printf( "matrix update time = %f seconds \n \n", matrix_update_time);
	printf( "y[0] = %f \n", y[0]);
	printf( "y[m-1] = %f \n", y[m-1]);
	printf( "trimv time = %f seconds \n \n", trimv_time);
	printf( "Total time = %f seconds \n \n", total_time);
	printf("End MASTER output\n");
      }
      printf("myid=%d\n",myid);
      printf( "vector add time = %f seconds \n \n", vector_time);
    }
  }

  MPI_Finalize();
}

void vector_update(double a[], double b[], double x[] ,int limit) {
  int i;

  for (i=0; i<limit; i++ )  {
    x[i] = 0.0;
    a[i] = 3.142+i;
    b[i] = 3.142;
  }
}

void matrix_update(double **l, double y[], int limit) {
  int i, j;

  for (i=0; i<limit; i++ ) {
    y[i] = 0.0;
    for (j=0; j<=i; j++ ) { 
      l[i][j] = 2.0;
    }
  }
}

void vector_add(double a[], double b[], double x[], int limit) {
  int i;

  for (i=0; i<limit; i++ ) {
    x[i] = a[i] + b[i];
  }
}

double vec_sum(double x[], int limit) {
  int i;
  double sum = 0.0;

  for (i=0; i<limit; i++ ) {
    sum = sum + x[i];
  }
  return sum;
}

void trimv(double **l, double y[], double a[], int limit) {
  int i, j;

  for (i=0; i<limit; i++) {
    for (j=0; j<=i; j++) { 
      y[i] = y[i] + l[i][j]*a[j];
    }
  }
}

