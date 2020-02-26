#include <stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<math.h>


int IsEqual(double *A, double *B, int m, int n){
	int i, j = 0;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			if(A[i*n + j] != B[i*n + j]){
				printf("Failure at index i : %d, j = %d\n", i, j);
				return 0;
			}

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

void dump(double *C_serial, double *C, int M, int N, FILE *f1, FILE *f2){

	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			fprintf(f1, "%lf ", C_serial[i*N+j]);
			fprintf(f2, "%lf ", C[i*N+j]);
		}

		fprintf(f1, "\n");
		fprintf(f2, "\n");
	}

}

int main(int argc, char ** argv){
	double *A, *B, *C, *C_serial;

	int rank;
	int size;
	int N = 100;
	int tag = 0;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int width = N / (size-1);

	if(N % (size-1) != 0)
		width++;

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

	MPI_Bcast(A, N*32, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B, 32*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	if(rank != 0){
		//Compute partial rows, of the matrix
		int start_ind = width*(rank-1);
		int end_ind = start_ind + width;
		
		if(end_ind > N)
			end_ind = N;

		for(int i = start_ind; i < end_ind; i++){
			for(int j = 0; j < N; j++){
				C[i*N + j] = 0;
				for(int k = 0; k < 32; k++)
					C[i*N + j] += A[i*32 + k]*B[k*N + j];
			}
		}

		printf("Work done by proc %d, is %d\n", rank, end_ind-start_ind);
		//Send Message to rank 0
		MPI_Send(C+start_ind*N, (end_ind-start_ind)*N, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}

	if(rank == 0){
		int start_ind;
		int end_ind;

		for(int i = 1; i < size; i++){
			start_ind = width*(i-1);
			end_ind = start_ind + width;

			if(end_ind > N)
				end_ind = N;

			MPI_Recv(C+start_ind*N, (end_ind-start_ind)*N, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Check the correctness
	if(rank == 0){
		Matrix_Multiply(A, B, C_serial, N, 32, N);
		int result = IsEqual(C, C_serial, N, N);

		FILE *f1, *f2;

		f1 = fopen("c_serial.txt", "w");
		f2 = fopen("c.txt", "w");

		dump(C_serial, C, N, N, f1, f2);

		if(result == 0){
			printf("Matrix Multiplication in MPI FAILURE\n");
		}
		else{
			printf("Matrix Multiplication in MPI SUCCESS\n");
		}
	}


	MPI_Finalize();
	return 0;
}