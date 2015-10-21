//EXERCISE 4, TASK 3
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define S 1.6449340668482264364724151666460251892189499012067984

double pow( double base, double exp );

int main(int argc, char** argv)
{
    if (argc < 3) 
    {
        printf("Error: Need k and p.\n");
        return 1;
    }
    // --------- MPI --------- //
    MPI_Status status;
    MPI_Init(&argc, &argv);
    
    int size, rank;
    // Figure out the number of processes and our rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int N = pow(2,atoi(argv[1]));
    int P = pow(2,atoi(argv[2]));
    
    // If running on a larger system than P
    if (rank >= P)
    {
        MPI_Finalize();
        return 0;
    }
    // Less elements in system than processors.
    else if (N<P)
    {
        printf("Error: N must be larger or equal to P!\n");
        MPI_Finalize();
        return 1;
    }
        
    // ----------- MPI SENDING ---------- //
    int i;
    int n = N/P;
    int tag = 100;
    int iStart = rank*n;
    
    if (rank == 0)
    {
        printf("CPU0: \nK = %i, N = 2^k = %i\nNo. of CPUs: P = 2^p = %i\n", atoi(argv[1]),N,P);
    }
    
    // ------------ COMPUTING ------------- //
    // Generate vector v dynamic size
    double * v   = (double*)malloc(n*sizeof(double));
    double * S_n = (double*)malloc(sizeof(double));
    // Fill vector
    for(i=0; i<n; i++)
    {
        v[i] = (((double)1.0)/pow(iStart+i+1,2));
    }
    //Calculate sum
    *S_n = 0.0;
    for(i=0;i<n;i++)
    {
        *S_n += v[i];
    }
    
    
    // ----------- MPI COLLECTING ---------- //
    if (rank == 0)
    {
        double S_n_sum = *S_n;
        for (i=1;i<P;i++)
        {
            MPI_Recv(S_n,1,MPI_DOUBLE,i,tag,MPI_COMM_WORLD, &status);
            S_n_sum += *S_n;
        }
        // PRINT OUT
        printf("\nS_n = %1.16E\nS - S_n = %1.16E \n\n", S_n_sum, S-S_n_sum);
    }
    else
    {
        MPI_Send(S_n,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
    }
    
    free(v);
    free(S_n);
    
    MPI_Finalize();
    
    return 0;
}
