#include "svd.h"

// reorders a matrix of eigenvectors and of corresponding eigenvalues, in order of eigenvalues
void reorderEigen(double **A, double **D, int n){
	int i, j;
	double max;
	int maxInd;
	for(i = 0; i < n-1; i++){
		max = D[i][i];
		maxInd = i;
		for(j = i+1; j < n; j++){
			if(D[j][j] > max){
				max = D[j][j];
				maxInd = j;
			}
		}
		if(maxInd != i){
			swapCol(A, i, maxInd, n, n);
			swapDiag(D, i, maxInd);
		}
	} 
}

// prints the various matrices from a SVD decomposition
void printSvdResults(double **A, double **U, double **V, double **S, int n, int m, int verbose){
	
	if(!verbose){
		return;
	}
	printf("A\n");
	printMat(A, n, m);
	
	printf("U\n");
	printMat(U, n, n);
	
	printf("V\n");
	printMat(V, m, m);

	printf("S\n");
	printMat(S, n, m);

	return;
}


// runs SVD using Jacobi method
void svdJacobi(double **A, double **U, double **V, double **S, int n, int m, int verbose, int errorCheck){
	int i, j;
	int errors;

	double **AAt, **AtA, **Eu, **Ev;

	// allocates matrices and does eigenvalue decomp for AAt
	allocMat(&AAt, n, n);
	calcAAt(A, AAt, n, m);
	allocMat(&Eu, n, n);
	jacobiMet(AAt, &U, Eu, n);
	freeMat(AAt);
	reorderEigen(U, Eu, n);

	// frees if not needed anymore
	if(!errorCheck){freeMat(Eu);}


	// does same for AtA
	allocMat(&AtA, m, m);
	calcAtA(A, AtA, n, m);
	allocMat(&Ev, m, m);
	jacobiMet(AtA, &V, Ev, m);
	freeMat(AtA);
	reorderEigen(V, Ev, m);

	int minDim = minInt(m, n);
	// error checks that the eigenvalues match up for both decompisitions
	if(errorCheck){
		errors = 0;
		// TEST EIGENVALUES
		for(i = 0; i < minDim; i++){
			if( fabs(Eu[i][i] - Ev[i][i]) > 1.E-6){
				errors ++;
			}
		}
		fprintf(stderr, "Eigenvalues Compared\t\t%d Errors\n", errors);
		freeMat(Eu);
	}

	// FILLS S with sqrt of found eigenvalues
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			S[i][j] = 0;
		}
	}	
	for(i = 0; i < minDim; i++){
		S[i][i] = sqrt(Ev[i][i]);
	}

	freeMat(Ev);

	


	// FIX SIGNS so that U * S * Vt = A
	// checks that first entry of each col of US = first entry of each col of AV	(inverse of V is Vt)
	double sum;

	// checks up to min of n and m
	for(i = 0; i < minDim; i++){
		sum = 0;
		// calcs first entry of AV column i
		for(j = 0; j < m; j++){
			sum += A[0][j] * V[j][i];
		}
		// if they don't equal then flip the sign of the relevant eigenvector of V
		if(fabs(sum - U[0][i] * S[i][i]) > 1.E-6){
			flipSignColMat(V, i, m, m);
		}
	}

	// prints results if verbose turned on
	printSvdResults(A, U, V, S, n, m, verbose);


	// error checks decomposition
	if(errorCheck){
		// TRANSPOSE V
		double **Vt;
		allocMat(&Vt, m, m);
		transpose(V, Vt, m, m);

		// TEST DECOMPOSITION
		errors = 0;
		double **newA1, **newA2;
		allocMat(&newA1, n, m);
		allocMat(&newA2, n, m);
		multMat(U, S, newA1, n, n, m);
		multMat(newA1, Vt, newA2, n, m, m);
		for(i = 0; i < n; i++){
			for(j = 0; j < m; j++){
				if(fabs(A[i][j] - newA2[i][j]) > 1.E-6){
					errors++;
				}
			}
		}
		fprintf(stderr, "Decomposition Compared\t\t%d Errors\n", errors);

		// prints U * S * Vt
		if(verbose){
			printf("newA\n");
			printMat(newA2, n, m);
		}

		freeMat(Vt);freeMat(newA1);freeMat(newA2);

	}
	
	return;

}
