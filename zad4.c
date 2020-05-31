#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>


double gamma_Euler_sum(int n, int m)
{
    double sum = 0;
    for(int k = n;k < m; k++)
        sum +=  1./ k;
    return sum;
}

int main(int argc, char *argv[])
{

	int npes;
	int myrank;
	double start, end;


	MPI_Init(&argc, &argv);
	
	
    MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int N, n, r; //N ilość iteracji, n ilość iteracji dla pojedynczego procesu, r reszta
    if(!myrank)
    {
        N = atoi(argv[1]);
        n = N/npes;
        r = N - n * npes;
    }
    
    start = MPI_Wtime();
    
    int Bcast_buff[] = {n,r};
	MPI_Bcast( Bcast_buff, 2, MPI_INT,0, MPI_COMM_WORLD );
    n = Bcast_buff[0];
    r = Bcast_buff[1];

	double localsum = 0, globalsum = 0;
  
    if(myrank < (npes - 1))
        localsum = gamma_Euler_sum(1+ n*myrank,1+ n*myrank+n);
	else
        localsum = gamma_Euler_sum(1+ n*myrank,1+ n*myrank+n+r);

	//printf("Localsum %lf\n", localsum);

	MPI_Reduce( &localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

	end = MPI_Wtime();
    double wynik = globalsum - log(N);
	if(!myrank){
        printf("Suma: %lf\n",globalsum);
        printf("Wynik: %lf\n", wynik);
		printf("obliczono w: %lfs\n", end-start);
	}
	
	MPI_Finalize();
	return 0;
}
