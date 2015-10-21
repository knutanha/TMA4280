//EXERCISE 4, TASK 1
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define S 1.6449340668482264364724151666460251892189499012067984
double pow( double base, double exp );

int main(int argc, char** argv)
{

    if (argc < 2) 
    {
        printf("Error: Need one input parameter k, such that n = 2^k.\n");
        return 1;
    }
    int N=pow(2,atoi(argv[1]));
    printf("\nK = %i, N = 2^k = %i\n",atoi(argv[1]),N);
    
    // Generate vector v dynamic size
    double *v = (double*)malloc(N*sizeof(double));
    // Fill vector
#pragma omp parallel for
    for(int i=0;i<N;i++)
    {
        v[i] = (((double)1.0)/pow(i+1,2));
        if(i==0 || i==N-1)
            printf("(%p): v[%i] = %1.16f\n",v+i,i,v[i]);
    }
    //Calculate sum
    double S_n = 0.0;
#pragma omp parallel for schedule(static) reduction(+:S_n)
    for(int i=0;i<N;i++)
    {
        S_n += v[i];
        int th_id = omp_get_thread_num();
        printf("Hello World from thread %d\n", th_id);
    }
    
    printf("\nS_n = %1.16E\nS - S_n = %1.16E \n\n",S_n, S-S_n);
    return 0;
} 


