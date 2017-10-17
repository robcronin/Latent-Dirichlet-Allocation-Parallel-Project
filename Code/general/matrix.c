#include "matrix.h"

// allocates memory for a matrix
void allocMat(double ***A, int n, int m){
	double *temp;
	temp = malloc(n*m*sizeof(double));
	*A = malloc(n*sizeof(double *));

	int i;
	for(i=0; i<n; i++){
		(*A)[i] = &temp[i*m];
	}

	return;
}

// frees memory
void freeMat(double **A){
	free(A[0]);
	free(A);
	
	return;
}

// fills a matrix either from a file or else randomly
void fillMat(double ***A, int *n, int *m, char *filename){
	int i, j;

	// if from a file
	if(filename != NULL){

		// open up file
		FILE *file;
		if( (file = fopen(filename, "r")) == NULL){
			fprintf(stderr, "*** ERROR ***\nCouldn't open %s\n\n",filename);
			exit(1);
		}
		else{
			char lineBuffer[1024];
			char *p, *end;

			// get dimensions
			if( fgets(lineBuffer, 1024, file) == NULL ){
				fprintf(stderr, "*** ERROR ***\nToo few lines in %s\n\n", filename);
				exit(1);
			}
			p = lineBuffer;
			*n = strtod(p, &end);
			if(p == end){
				fprintf(stderr, "*** ERROR ***\nDimensions in %s not entered correctly at top of file\n\n", filename);
				exit(1);
			}	
			p = end;
			*m = strtod(p, &end);

			// allocate memory
			allocMat(A, *n, *m);

			// fill from file while error checking
			for(i = 0; i < *n; i++){
				if( fgets(lineBuffer, 1024, file) == NULL ){
					fprintf(stderr, "*** ERROR ***\nToo few lines in %s\n\n", filename);
					exit(1);
				}
				p = lineBuffer;
				for(j = 0; j < *m; j++){
					(*A)[i][j] = strtod(p, &end);
					if(p == end){
						fprintf(stderr, "*** ERROR ***\nLine %d of %s does not contain enough values\n\n", i+1, filename);
						exit(1);
					}	
					p = end;
				}
			}
		}		
		fclose(file);
	}

	// if no file specified then fill randomly
	else{
		// allocate memory
		allocMat(A, *n, *m);

		// randomly fill from (-2, 2)
		for(i = 0; i < *n; i++){
			for(j = 0; j < *m; j++){
				(*A)[i][j] = 4*(drand48() - 0.5);
			}
		}
	}
	
	return;
}

// trivial matrix multiplication 
void multMat(double **A, double **B, double **C, int n, int m, int q){
	int i, j, k;
	double sum;
	for(i = 0; i < n; i++){
		for(j = 0; j < q; j++){
			sum = 0;
			for(k = 0; k < m; k++){
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}
	}
	
	return;
}

// prints matrix
void printMat(double **A, int n, int m){
	int i, j;
	for(i = 0; i<n; i++){
		for(j=0; j<m; j++){
			printf("%lf\t", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	return;
}


// swaps two columns of a matrix	
void swapCol(double **A, int c1, int c2, int n, int m){
	if(c1 > m || c2 > m){
		fprintf(stderr, "\tERROR\nInvalid column entry for swapCol\n");
		return;
	}
	
	int i;
	double temp;
	// swaps value from each row
	for(i = 0; i < n; i++){
		temp = A[i][c1];
		A[i][c1] = A[i][c2];
		A[i][c2] = temp;
	}

	return;
} 

// swaps the diagonal entries of a matrix
void swapDiag(double **A, int d1, int d2){
	double temp = A[d1][d1];
	A[d1][d1] = A[d2][d2];
	A[d2][d2] = temp;
	
	return;
}


// transposes a matrix and stores the result in At
void transpose(double **A, double **At, int n, int m){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			At[j][i] = A[i][j];
		}
	}

	return;
}


// flips the sign of a given column (used for eigenvectors mostly)
void flipSignColMat(double **A, int c, int n, int m){
	if(c >= m){
		fprintf(stderr, "*** ERROR ***\nInvalid column entry for flipSignCol\n\n");
		return;
	}

	int i;
	for(i = 0; i < n; i++){
		A[i][c] *= -1;
	}
	
	return;
}


// normalises the entries in a matrix
void normMat(double **A, int n, int m){
	int i, j;
	double sum;
	for(i = 0; i < n; i++){
		sum = 0;
		for(j = 0; j < m; j++){
			sum += A[i][j] * A[i][j];
		}
		for(j = 0; j < m; j++){
			A[i][j] /= sum;
		}
	}

	return;
}


// calculates A * At without storing At (from Jacobi method)
void calcAAt(double **A, double **AAt, int n, int m){
	int i, j, k;
	double sum;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			sum = 0;
			for(k=0; k<m; k++){
				sum += A[i][k] * A[j][k];
			}
			AAt[i][j] = sum;
		}
	}
	
	return;
}


// calculates At * A without storing At (from Jacobi method)
void calcAtA(double **A, double **AtA, int n, int m){
	int i, j, k;
	double sum;
	for(i=0; i<m; i++){
		for(j=0; j<m; j++){
			sum = 0;
			for(k=0; k<n; k++){
				sum += A[k][i] * A[k][j];
			}
			AtA[i][j] = sum;
		}
	}
	
	return;
}

// ensures the row of a matrix sums to 1
void sumRowOneMat(double **A, int r, int n, int m){
	if(r >= n){
		fprintf(stderr, "*** ERROR ***\nInvalid row passed to normRowMat\n\n");
		return;
	}

	double sum = 0;
	int i;
	for(i = 0; i < m; i++){
		sum += A[r][i];
	}
	if(sum != 1){
		for(i = 0; i < m; i++){
			A[r][i] /= sum;
		}
	}

	return;
}




// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

//					INTEGERS


// allocates memory for a matrix
void  allocMatInt(int ***A, int n, int m){
	int *temp;
	if((temp = malloc(n*m*sizeof(int))) == NULL){
		printf("FAIL FAIL FAIL\n");	
	}
	*A = malloc(n*sizeof(int *));

	int i;
	for(i=0; i<n; i++){
		(*A)[i] = &temp[i*m];
	}

	return;
}

// allocates memory for a matrix and initialises to 0
void  callocMatInt(int ***A, int n, int m){
	int *temp;
	temp = malloc(n*m*sizeof(int));
	*A = malloc(n*sizeof(int *));

	int i;
	for(i=0; i<n; i++){
		(*A)[i] = &temp[i*m];
	}

	zeroMatInt(*A, n, m);

	return;
}

// frees memory
void  freeMatInt(int **A){
	free(A[0]);
	free(A);
	
	return;
}

// allocates and zeroes an integer array
void	callocArrInt(int **A, int n){
	*A = malloc(n * sizeof(int));
	zeroArrInt(*A, n);

	return;
}

// zeroes all entries in matrix
void	zeroMatInt(int **A, int n, int m){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			A[i][j] = 0;
		}
	}

	return;
}

// zeroes all entries in an array
void	zeroArrInt(int *A, int n){
	int i;
	for(i = 0; i < n; i++){
		A[i] = 0;
	}

	return;
}
