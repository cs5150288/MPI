#include <stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<math.h>


int IsEqual(double *A, double *B, int m, int n){
	int i, j = 0;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			if(A[i*n + j] != B[i*n + j])
				return 0;

	return 1;
}

void Matrix_Multiply(double *A, double *B, double *C, int m, int n, int p){
	int i, j, k = 0;

	for(i = 0; i < m; i++){
		for(j = 0; j < p; j++){
			C[i*p + j] = 0;

			for(k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k]*B[k*p + j];

		}
	}
}

int main(int argc, char ** argv){
	double *A, *B, *C, *C_serial;

	int rank;
	int size;
	int N = 100;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf( "Hello world from process %d of %d\n", rank, size);

	A = (double*)malloc(sizeof(double)*N*32);
	B = (double*)malloc(sizeof(double)*32*N);
	C = (double*)malloc(sizeof(double)*N*N);
	C_serial = (double*)malloc(sizeof(double)*N*N);

	if(rank == 0){
		for(int i = 0; i < N; i++)
			for(int j = 0; j < 32; j++)
				A[i*32+j] = drand48();

		for(int i = 0; i < 32; i++)
			for(int j = 0; j < N; j++)
				B[i*N+j] = drand48();
	}
	// Barrier, wait for thread 0 to initialize A[N][32] and B[32][N]
	MPI_Barrier(MPI_COMM_WORLD);

	//Perform matrix multiplication for AXB using communication and store it in C
	//Communication between process 0, and the rest



	MPI_Barrier(MPI_COMM_WORLD);

	//Check the correctness
	if(rank == 0){
		Matrix_Multiply(A, B, C_serial, N, 32, N);
		IsEqual(C, C_serial, N, N);
	}


	MPI_Finalize();
	return 0;
}