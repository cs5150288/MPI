#include<stdio.h>
#include<pthread.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<tgmath.h>
#include<string.h>


#define MATRIX_SIZE 4000

pthread_barrier_t barrier;

// double a[MATRIX_SIZE][MATRIX_SIZE];
// double L[MATRIX_SIZE][MATRIX_SIZE];
// double U[MATRIX_SIZE][MATRIX_SIZE];
// int pi[MATRIX_SIZE];

double *a[MATRIX_SIZE];
double *b[MATRIX_SIZE];
double *L[MATRIX_SIZE];
double *U[MATRIX_SIZE];
double *result[MATRIX_SIZE];
int pi[MATRIX_SIZE];


int NUM_THREADS;
int *k_;
double *max_array;
double overall_max;
int overall_ind;
int ZERO_FLAG = 0;


void swap(int *a, int *b){
	int t = *a;
	*a = *b;
	*b = t;
}

void fswap(double *a, double *b){
	double t = *a;
	*a = *b;
	*b = t;
}

void print_array(double * a, int N){

	printf("\n**************************\n");
	for(int i = 0; i < N; i++){
		printf("%lf ", a[i]);
	}
	printf("\n---------------------------\n");
}


void* divide(void* arg){
	int tid = (int) arg;
	printf("Thread number : %d\n", tid);
	int n = MATRIX_SIZE;
	int k = 0;
	int i = 0, j = 0;

	for(k = 0; k < n; k++){

		k_[tid] = 0; max_array[tid] = 0;
	
		for(i = k+tid; i < n; i += NUM_THREADS){
			if(fabs(a[i][k]) > max_array[tid]){
				max_array[tid] = fabs(a[i][k]);
				k_[tid] = i;
			}
		}

		pthread_barrier_wait(&barrier);	

		if(tid == 0){
			
			overall_max = max_array[0];
			overall_ind = k_[0];
			for(int t = 1; t < NUM_THREADS; t++){
				if(overall_max < max_array[t]){
					overall_max = max_array[t];
					overall_ind = k_[t];
				}
			}

			swap(&pi[k], &pi[overall_ind]);
			U[k][k] = a[overall_ind][k];
		}

		pthread_barrier_wait(&barrier);

		if(overall_max == 0){
			ZERO_FLAG = 1;
			pthread_exit(NULL);
		}


		for(i = k+tid; i < n; i += NUM_THREADS)
			fswap(&a[k][i], &a[overall_ind][i]);

		for(i = tid; i < k; i += NUM_THREADS)
			fswap(&L[k][i], &L[overall_ind][i]);

		pthread_barrier_wait(&barrier);

		for(i = k+1+tid; i < n; i += NUM_THREADS){
			L[i][k] = a[i][k] / U[k][k];
			U[k][i] = a[k][i];
		}

		pthread_barrier_wait(&barrier);

		for(i = k + 1; i < n; i++){
			for(j = k + 1 + tid; j < n; j += NUM_THREADS){
				a[i][j] = a[i][j] - L[i][k]*U[k][j];
			}
		}

		pthread_barrier_wait(&barrier);

	}
		
	pthread_exit(NULL);
}


void LU_factorization(int n, int num_threads){
	int k = 0;
	int i = 0, j = 0;
	
	pthread_t tids[num_threads];
	// pthread_attr_t attr;	// pthread_attr_init(&attr);
	pthread_barrier_init(&barrier, NULL, num_threads);


	k_ = (int*)malloc(sizeof(int)*num_threads);
	max_array = (double*)malloc(sizeof(double)*num_threads);

	for(i = 0; i < num_threads; i++){
		int er = pthread_create(&tids[i], NULL, divide, (void*) i);
		if(er != 0){
			printf("Thread %d was not created\n", i);
			return;
		}
	}
	for(i = 0; i < num_threads; i++)
		pthread_join(tids[i], NULL);

	if(ZERO_FLAG == 1){
		printf("Error, SINGULAR MATRIX\n");
		return;
	}

}

void print_matrix(int n, double a[][n]){
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("%lf ", a[i][j]);
		}
		printf("\n");
	}

}

void check_error(int n, double **a){

	for(int i = 0; i < MATRIX_SIZE; i++){
		for(int j = 0; j < MATRIX_SIZE; j++){
			double res = 0;
			for(int k = 0; k < MATRIX_SIZE; k++){
				res += L[i][k]*U[k][j];
			}
			result[i][j] = res;
		}
	}


/*
	print_matrix(MATRIX_SIZE, result);
	printf("------------------------------------\n");
	
	for(int i = 0; i < MATRIX_SIZE; i++)
		printf("%d ", pi[i]);
	printf("\n------------------------------------\n");
*/


	long double error[MATRIX_SIZE];
	
	for(int i = 0; i < MATRIX_SIZE; i++)
		error[i] = 0;


	// print_matrix(MATRIX_SIZE, a);
	// printf("----------------------------------\n");

	for(int i = 0; i < MATRIX_SIZE; i++){
		int row_A = pi[i];
		for(int c = 0; c < MATRIX_SIZE; c++){
			error[c] += powl((long double)a[row_A][c] - (long double)result[i][c], 2);
			// printf("%lf ", a[row_A][c] - result[i][c]);
		}
		// printf("\n");
	}

	long double err = 0;
	
	for(int i = 0; i < MATRIX_SIZE; i++){
		err += sqrt(error[i]);
	}

	printf("Error of the LU decomp is : %Lf\n",err);
}

void load_matrix(char* filename){
	FILE *fp;
	fp = fopen(filename, "r");

	if(fp == NULL){
		perror("Error opening the file.\n");
		return;
	}
	else{
		for(int i = 0; i < MATRIX_SIZE; i++){
			for(int j = 0; j < MATRIX_SIZE; j++){
				fscanf(fp, "%lf", &a[i][j]);
			}
		}
	}

	fclose(fp);
}


void dump_PLU(char* p, char* l, char* u, int num_threads, int size){
// x_size_nthreads.txt
	char pf[16], lf[16], uf[16], buf[16];
	memset(pf, 0, 16);
	memset(lf, 0, 16);
	memset(uf, 0, 16);

	sprintf(buf, "%d_%d.txt", size, num_threads);
	strcat(pf, p); strcat(pf, buf);
	strcat(lf, l); strcat(lf, buf);
	strcat(uf, u); strcat(uf, buf);

	// printf("%s\n", pf);
	// printf("%s\n", lf);
	// printf("%s\n", uf);

	FILE *f1, *f2, *f3;

	f1 = fopen(pf, "w"); f2 = fopen(lf, "w"); f3 = fopen(uf, "w");

	for(int i = 0; i < MATRIX_SIZE; i++){
		for(int j = 0; j < MATRIX_SIZE; j++){
						
			fprintf(f2, "%lf", L[i][j]);
			fprintf(f3, "%lf", U[i][j]);
			
			if(j != MATRIX_SIZE-1){
				fprintf(f2, " ");
				fprintf(f3, " ");
			}
		}
		fprintf(f2, "\n");
		fprintf(f3, "\n");
	}


	for(int i = 0; i < MATRIX_SIZE; i++){
		
		for(int j = 0; j < MATRIX_SIZE; j++){
			if(pi[i] == j)
				fprintf(f1, "%d", 1);
			// (j != MATRIX_SIZE-1)
			else
				fprintf(f1, "%d", 0);

			if(j != MATRIX_SIZE - 1)
				fprintf(f1, " ");

		}
		
		fprintf(f1, "\n");
	}

	fclose(f1); fclose(f2); fclose(f3);
}

int main(int argc, char** argv){

	// printf("Came here");

	int load = 1;

	for(int it = 0; it < MATRIX_SIZE; it++){
		a[it] = (double*)malloc(MATRIX_SIZE*sizeof(double));
		L[it] = (double*)malloc(MATRIX_SIZE*sizeof(double));
		U[it] = (double*)malloc(MATRIX_SIZE*sizeof(double));
		b[it] = (double*)malloc(MATRIX_SIZE*sizeof(double));
		result[it] = (double*)malloc(MATRIX_SIZE*sizeof(double));
	}


	int n = MATRIX_SIZE;
	NUM_THREADS = atoi(argv[1]);
	load = atoi(argv[2]);
	time_t seconds;

	printf("Number of threads : %d\n", NUM_THREADS);

	for(int i = 0; i < n; i++)
		pi[i] = i;

	// double b[MATRIX_SIZE][MATRIX_SIZE];
	// double b[][MATRIX_SIZE] = {{1, 1, 1}, {4, 3, -1}, {3, 5, 3}};


	if(load == 1)
		load_matrix("B_4000.txt");


	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(i == j)
				L[i][j] = 1.0;
			else if(i < j)
				L[i][j] = 0;
			else
				U[i][j] = 0;

			if(load == 0)
				a[i][j] = drand48();
			b[i][j] = a[i][j];
			// a[i][j] = b[i][j];
			
		}
	}

	
	// print_matrix(n, a);
	// printf("------------------------------------\n");
	

	time(&seconds);
	printf("Time start : %ld\n", seconds);
	LU_factorization(n, NUM_THREADS);
	time(&seconds);
	printf("Time end : %ld\n", seconds);

	dump_PLU("P_", "L_", "U_", NUM_THREADS, MATRIX_SIZE);

/*
	print_matrix(n, L);
	printf("------------------------------------\n");
	
	print_matrix(n, U);
	printf("------------------------------------\n");
*/

	// check_error(MATRIX_SIZE, b);
}
