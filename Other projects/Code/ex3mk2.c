//EXERCISE 4, TASK 3
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

// pi/6 from Wolfram-Alpha (http://www.wolframalpha.com/)
#define S_TRUE 1.6449340668482264364724151666460251892189499012067984

double pow( double base, double exp );

void printHeader(int N, int P, int size)
{
    printf("\nROOT:\tN = 2^k = %i, P = 2^p = %i\n\tFull cluster size: %i\n",N,P,size);
}

void printResult(double * sum, double time)
{
    printf("\nS - S_n = %1.16E \n\nTime Elapsed: %f\n\n", S_TRUE-*sum,time);
}

double * createVector(int iStart, int size)
{
    double * vector   = (double*)malloc(size*sizeof(double));
    int i;
    for(i=0; i<size; i++)
    {
        vector[i] = (((double)1.0)/pow(iStart+i+1,2));
    }
    return vector;
}

double sumVector(double * vector, int size)
{
    double sum = 0.0;
    int i;
    for(i=0;i<size;i++)
    {
        sum += vector[i];
    }
    return sum;
}

int main(int argc, char** argv)
{
    if (argc < 3) 
    {
        printf("Error: Need k and p.\n");
        return 1;
    }
    
    int size, rank, N, P, n, iStart, i;
    double s_n, s_n_sum,time;
    
    time = MPI_Wtime();
    
    // --------- MPI --------- //
    MPI_Status status;
    MPI_Comm cluster;
    MPI_Init(&argc, &argv);
    
    // Figure out the number of processes and our rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    N = pow(2,atoi(argv[1]));
    P = pow(2,atoi(argv[2]));
    
    if ( N < P )
    {
        printf("Error: N must be larger or equal to P!\n");
        MPI_Finalize();
        return 1;
    }
                   
    // Sort out superfluous processors.
    // In the case where size > P;
    if (size >= P) MPI_Comm_split(MPI_COMM_WORLD, rank < P, rank, &cluster);
    else cluster = MPI_COMM_WORLD;
    
    // Stop processes with rank >= P
    if (rank >= P){MPI_Finalize();return 0;}
    
    // Header print
    if (rank == 0) printHeader(N,P,size);
    
    // ------------ COMPUTING ------------- //
    n = N/P;
    iStart = rank*n;
    double * v = createVector(iStart,n);
    s_n = sumVector(v,n);
    MPI_Reduce(&s_n, &s_n_sum, 1, MPI_DOUBLE, MPI_SUM, 0, cluster);
    
    // PRINT OUT
    time = MPI_Wtime()-time;
    if (rank == 0) printResult(&s_n_sum,time);
    free(v);
    
    MPI_Finalize();
    
    return 0;
}
