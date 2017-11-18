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

  // Storage for this processors part of the data.
  // (set large enough for single processor to allow runs with 
  // different numbers of processors without recompilation. For
  // memory efficiency and scalability the space allocated should
  // be around 1/pth of the total memory requirement)

  int ns=n, ms=m;

  //double ll[m][m], **l;
  double **l;
  double a[ns], b[ns], x[ns], y[ms];
  double res, glob_res;
  
  double vector_start, sum_start, trimv_start, vec_update_start, matrix_update_start, total_start;

  double vec_update_time = 0.0;
  double vector_time = 0.0;
  double sum_time = 0.0;
  double matrix_update_time = 0.0;
  double trimv_time = 0.0;
  double total_time = 0.0;

  // MPI specific data
  int myid, numprocs, rc, ierr;

  // Local array sizes (on a particular processor)
  int nx, nnx, mx, nmx, modx, pdn, pdm;

  // MPI initialisation

  ierr = MPI_Init(NULL, NULL); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  printf ("Process %d of %d is alive\n", myid, numprocs);
  
  // allocate (lower triangular) matrix as a contibuous block and set up pointers to allow l[][] access 
  l = (double **) malloc((unsigned) ms*sizeof(double*));  
  for (i=0; i <= ms-1; i++) {
    l[i] = (double *)malloc((unsigned) ms * sizeof(double));
  }

  // Work out size of data on this processor and number of processors
  // with data

  nx = (n+numprocs-1)/numprocs;
  nnx = nx;
  pdn = (n+nx-1)/nx;
  if (myid == pdn-1) {
    modx = fmod(n,nx);
    if (modx != 0) nx = modx;
  } else if (myid >= pdn) {
    nx = 0;
  }

  mx = (m+numprocs-1)/numprocs;
  nmx = mx;
  pdm = (m+mx-1)/mx;
  if (myid == pdm-1) {
    modx = fmod(m,mx);
    if (modx != 0) mx = modx;
  } else if (myid >= pdm) {
    mx = 0;
  }

  // Loop over whole program to get reasonable execution time
  // (simulates a time-steping code)

  if (myid == 0) total_start = omp_get_wtime();

  for (k=1; k<=nsteps; k++ )  {
 
    vector_update(a, b, x, nx);
   
    // vector addition operation

    vector_start = omp_get_wtime();
    vector_add(a, b, x, nx);
    vector_time = vector_time + (omp_get_wtime() - vector_start);
    
    // sum operation
    
    sum_start = omp_get_wtime();
    res = vec_sum(x, nx);
    sum_time = sum_time +(omp_get_wtime() - sum_start);
    
    // Now calculate the sum of the partial results 
    
    glob_res = 0;
    ierr = MPI_Reduce(&res, &glob_res, 1, MPI_DOUBLE_PRECISION, \
		      MPI_SUM, 0, MPI_COMM_WORLD);

    // 'Update' lower triangular matrix and y vector
    
    matrix_update_start = omp_get_wtime();
    matrix_update(l, y, m); 
    matrix_update_time = matrix_update_time +(omp_get_wtime() - matrix_update_start);

    // triangular matrix times vector operation

    // perform the mv operation on the local data.
    // Note, array a was distributed block, but now all processors require
    // access to it. An all-to-all gather data exchange copies the data 
    // to everyone.

    // Get a copy of a for each processor.

    MPI_Gather(a, nnx, MPI_DOUBLE_PRECISION, a, nnx,	\
		  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    //MPI_Allgather(MPI_IN_PLACE, nnx, MPI_DOUBLE_PRECISION, a, nnx,	\
		  MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);

    // Broadcast a to all processes
    MPI_Bcast(a,ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);

    // Perform the matrix vector operation on the local data
    // Pass in the dimension of the local data and the global
    // size of the array (for interleaved scheduling)

    
    // triangular matrix times vector operation

    trimv_start = omp_get_wtime();
    trimv(l, y, a, m); 
    trimv_time = trimv_time + (omp_get_wtime() - trimv_start);

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
	printf( "Sum of x() = %f \n", glob_res);
	printf( "matrix update time = %f seconds \n \n", matrix_update_time);
	printf( "y[0] = %f \n", y[0]);
	printf( "Total time = %f seconds \n \n", total_time);
	printf("End MASTER output\n");
      }

      printf( "a[0] = %f myid=%d \n", a[0], myid);
      printf( "a[1] = %f myid=%d \n", a[1], myid);
      printf( "a[nx-1] = %f myid=%d \n", a[nx-1], myid);
      printf( "x[nx-1] = %f myid=%d \n", x[nx-1], myid);
      printf( "y[mx-1] = %f myid=%d \n", y[mx-1], myid);
      printf("myid=%d\n",myid);
      printf( "vector add time = %f seconds \n \n", vector_time);
      printf( "sum time = %f \n \n", sum_time);
      printf( "trimv time = %f seconds \n \n", trimv_time);
    }
  }

  MPI_Finalize();
}

void vector_update(double a[], double b[], double x[] ,int limit) {
  int i;
  int myid;
  MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  for (i=0; i<limit; i++ )  {
    x[i] = 0.0;
    a[i] = 3.142+i +myid*limit;
    b[i] = 3.142;
  }
}

void matrix_update(double **l, double y[], int m) {
  int i, j, ii;
  int myid, numprocs;

  MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

  // Use interleaved distribution (for load balance). Each processor
  // has m/p rows (max) of l and the corresponding elements of y
  
  ii=0;
  for (i=myid; i<m; i+=numprocs ) {
    y[ii] = 0.0;
    for (j=0; j<=i; j++ ) { 
      l[ii][j] = 2.0;
    }
    ii++;
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

void trimv(double **l, double y[], double a[], int m) {
  int i, j, ii;
  int myid, numprocs;
  
  MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

  ii=0;
  for (i=myid; i<m; i+=numprocs) {
    for (j=0; j<=i; j++) { 
      y[ii] = y[ii] + l[ii][j]*a[j];
    }
    ii++;
  }
}

