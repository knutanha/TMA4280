#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>

// Parallel toolboxes
#include <mpi.h>
#include <omp.h>

typedef double Real;
// GLOBAL VARS
#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286


// ALLOCATION FUNCTIONS
Real *createRealArray (int n);
Real **createReal2DArray (int n1, int n2);
void free2DArray(Real **input);

// TRANSPOSER
void transpose_paral(Real** bt, Real** b, int* sdispl,int* scount, int m, int dist, int* distribution, int size,  MPI_Comm communicator);

// FORTRAN FUNCTIONS
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

// I/O FUNCTIONS (Only to be used on very small systems)
void printRealArray(Real* input, int inputlength, int rank);
void print2DArray(Real** input, int rows, int cols, int rank);

// Function for filling right hand side
Real f(int j, int i, Real h)
{
    Real x,y;
    x = h + h*(double)i;
    y = h + h*(double)j;
    return 5.0*pi*pi * sin(pi*x) * sin(2.0*pi*y);
}
// Function for evaluating the solution
Real evaluate(int j, int i, Real u_ji, Real h)
{
    Real x,y,err;
    x = h + h*(double)i;
    y = h + h*(double)j;
    return fabs(sin(pi*x) * sin(2.0*pi*y) - u_ji);
}

int main(int argc, char **argv)
{
    //////////////////////////////////////////
    //                                      //
    //        INITIALIZING VARIABLES        //
    //                                      //
    //////////////////////////////////////////
    Real *diag, **b, **bt, **z;
    Real h, umax, globalMax, globalMin;
    int i, j, k, l, m, n, nn, omp_threads;
    
    // Parallel variables
    int size, rank, dist, leftovers, globalRank, distDispl;
    int* distribution,* scount, *sdispl;
    
    if( argc < 2 ) {
        printf("need a problem size\n");
        return 0;
    }
	
    // Problem size
    n = atoi( argv[argc - 1] );
    m  = n-1;
    nn = 4*n;
    h    = 1./(Real)n;
    
    //////////////////////////////////////////
    //                                      //
    //           INITIALIZING MPI           //
    //                                      //
    //////////////////////////////////////////
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Finding number of avaiable openMP threads 
    #pragma omp parallel
    {
        omp_threads = omp_get_num_threads();
    }
    
    //////////////////////////////////////////
    //                                      //
    //          ALLOCATING MEMORY           //
    //                                      //
    //////////////////////////////////////////
    distribution = (int*) malloc (size * sizeof(int));
    scount 		 = (int*) malloc (size * sizeof(int));
    sdispl 	     = (int*) malloc (size * sizeof(int));
	
    // Distribution in between MPI threads
    dist =  m / size;
    leftovers = m - dist * size;
    for (i = 0; i < size; i++)
    {
	    distribution[i] = (dist + (i < leftovers ? 1:0) );
    }
    
    for (i = 0, distDispl = 0; i < rank; i++)
    {
        distDispl += distribution[i];
    }
    dist = distribution[rank];
	
	// Allocating memory for the numerical system
    diag = createRealArray(m);
    b = createReal2DArray(dist,m);
    bt = createReal2DArray(dist,m);
    
    // Allocating memory for each omp_thread to avoid
    // memory conflicts when parallilizing fst with 
    // openMP
    z = createReal2DArray(omp_threads,nn);
    
    
    //////////////////////////////////////////
    //                                      //
    //        FILL RIGHT HAND SIDE          //
    //                                      //
    //////////////////////////////////////////
    #pragma omp parallel for schedule(dynamic) private(i)
    for (j = 0; j < dist; j++)
    {
        for (i = 0; i < m; i++)
        {
            b[j][i] = h*h * f(j+distDispl,i,h);
        }
    }
    
    //////////////////////////////////////////
    //                                      //
    //      START OF NUMERICAL METHOD       //
    //                                      //
    //////////////////////////////////////////
    
    // Set eigenvalues
    for (i=0; i < m; i++) {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }

    // Perform fst on all columns
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {
        fst_(b[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // Transpose entire matrix
    transpose_paral(bt,b,sdispl, scount, m, dist, distribution, size, MPI_COMM_WORLD);
    
    // Perform fstinv on all columns of transformed matrix
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {	
        fstinv_(bt[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // Divide by eigenvalues
    #pragma omp parallel for schedule(dynamic) private(i)
    for (j=0; j < dist; j++) 
    {
        for (i=0; i < m; i++) 
        {
            bt[j][i] = bt[j][i]/(diag[j+distDispl]+diag[i]);
        }
    }
    
    // Perform fst on on all columns
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {
        fst_(bt[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // Transpose entire matrix
    transpose_paral(b,bt,sdispl, scount, m, dist, distribution, size,  MPI_COMM_WORLD);
    
    // Perform fstinv on all columns
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    { 
        fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    //////////////////////////////////////////
    //                                      //
    //          FIND MAXIMUM ERROR          //
    //                                      //
    //////////////////////////////////////////
    
    umax = 0.0;
    #pragma omp parallel for schedule(dynamic) private(i) shared(umax) 
    for (j=0; j < dist; j++) 
    {
        Real temp;
        for (i=0; i < m; i++) 
        {
            temp = evaluate(j+distDispl,i,b[j][i],h);
            if (temp > umax) umax = temp;
        }
    }
    MPI_Reduce(&umax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // if (rank == 0) printf("Maximum error = %e\n", globalMax);
    
    //////////////////////////////////////////
    //                                      //
    //        FIND MAXIMUM U VALUE          //
    //                                      //
    //////////////////////////////////////////
    
    umax = 0.0;
    #pragma omp parallel for schedule(dynamic) private(i) shared(umax)
    for (j=0; j < dist; j++) 
    {
        for (i=0; i < m; i++) 
        {
            if (b[j][i] > umax) umax = b[j][i];
        }
    }
    MPI_Reduce(&umax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("globalMax = %e\n",globalMax);
    
    
    //////////////////////////////////////////
    //                                      //
    //          CLEANING UP MEMORY          //
    //                                      //
    //////////////////////////////////////////
    
    free2DArray(b);
    free2DArray(bt);
    free(z);
    free(diag);
    free(sdispl);
    free(scount);
    free(distribution);
    
    MPI_Finalize();
}



// ------ PARALLELL TRANSPOSE FUNCTION ------ //
void transpose_paral(Real** bt, Real** b, int* sdispl,int* scount, int m, int dist, int* distribution, int size,  MPI_Comm communicator)
{
    int i,j;
    // Generating required meta data for MPI_Alltoallv()
	sdispl[0] = 0;
    scount[0] = distribution[0]*dist;
    for (i=1; i < size; i++)
    {
        scount[i] = distribution[i] * dist;
	    sdispl[i] = scount[i-1] + sdispl[i-1];
    }
    // Copying chunks of memory to save some time
    #pragma omp parallel for schedule(dynamic) private(i)
    for (j = 0; j < size; j++)
	{
		for (i = 0; i < dist; i++)
		{
			memcpy(bt[0] + sdispl[j] + distribution[j]*i, b[i] + sdispl[j]/dist, distribution[j]*sizeof(double));
		}
	}
	// Sending data in between MPI threads
	MPI_Alltoallv(bt[0], scount, sdispl, MPI_DOUBLE, b[0], scount, sdispl, MPI_DOUBLE, communicator);
	
	// Inserting 
	#pragma omp parallel for schedule(dynamic) private(j)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < dist; j++)
		{
			bt[j][i] = *(b[0] + i*dist + j);
		}
	}
}

//-------------- UTILITY -------------------//
        
Real *createRealArray (int n){
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++) 
    {
        a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0 ,n*sizeof(Real));
    return (a);
}

void free2DArray(Real **input) {
    free(input[0]);
    free(input);
}


void printRealArray(Real* input, int inputlength, int rank)
{
    int i;
	printf("Rank %i = {", rank);
	for (i=0; i<inputlength; i++)
	{
		printf("%e",input[i]);
		if (i < inputlength-1)
			printf(", ");
		else
			printf("}\n");
	}
}      


void print2DArray(Real** input, int rows, int cols, int rank)
{
    int i,j;
	printf("Rank%i = [\n", rank);
	for (i = 0; i < rows; i++)
	{
	    printf("[");
	    for (j = 0; j < cols; j++)
	    {
            printf( "%f",input[j][i] );
            if(j < cols-1) printf(" ");
        }
        printf("]\n");
    }
    printf("]\n\n");
}

