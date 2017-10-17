// Functions for running the Jacobi method to find eigenvalues and eigenvectors
// Method mostly adopted from the following paper
// 	http://www.cmi.ac.in/~ksutar/NLA2013/iterativemethods.pdf




#include "jacobiEigen.h"

// finds the index for the max value off the diagonal of a matrix
void findMaxIndex(double **D, int *maxi, int *maxj, int n){
	int i, j;
	double max = fabs(D[0][1]);
	*maxi = 0; *maxj = 1;

	// loops through all non-diagonal entries
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i != j){
				if(fabs(D[i][j]) > max){
					max = fabs(D[i][j]);
					*maxi = i;
					*maxj = j;
				}
			}
		}
	}
	
	return;
}

// calculates the theta value
double findTheta(double **D, int i, int j){
	if(D[i][i] == D[j][j]){
		if(D[i][j]>0){
			return M_PI/4;
		}
		else{
			return -M_PI/4;
		}
	}
	
	return (0.5 * atan( (2*D[i][j]) / (D[j][j] - D[i][i]) ));
}

// fills in the Givens rotation matrix
void fillS1(double **S1, int maxi, int maxj, double theta, int n){
	int i, j;

	// start with Identity matrix
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			S1[i][j] = 0;
		}
		S1[i][i] = 1;
	}

	// fill in cos and sin values
	S1[maxi][maxi] = cos(theta);
	S1[maxj][maxj] = cos(theta);
	S1[maxi][maxj] = sin(theta);
	S1[maxj][maxi] = -sin(theta);

	// fix if j > i
	if(maxj>maxi){
		S1[maxi][maxj] *= -1;
		S1[maxj][maxi] *= -1;
	}

	return;
}

// updates the Givens rotation matrix based on the old i and j values
void updateS1(double **S1, int maxi, int maxj, int oldMaxi, int oldMaxj, double theta){
	// returns old part to identity
	S1[oldMaxi][oldMaxi] = 1; 
	S1[oldMaxj][oldMaxj] = 1; 
	S1[oldMaxi][oldMaxj] = 0; 
	S1[oldMaxj][oldMaxi] = 0; 

	// fills in new part
	S1[maxi][maxi] = cos(theta);
	S1[maxj][maxj] = cos(theta);
	S1[maxi][maxj] = sin(theta);
	S1[maxj][maxi] = -sin(theta);

	if(maxj>maxi){
		S1[maxi][maxj] *= -1;
		S1[maxj][maxi] *= -1;
	}

	return;
}

// runs the jacobi Method
void jacobiMet(double **A, double ***S, double **D, int n){
	int i, j;

	// Set up matrices
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			D[i][j] = A[i][j];	// initialise D to A
			(*S)[i][j] = 0;		// S starts as identity
		}
		(*S)[i][i] = 1;
	}

	int maxi, maxj, oldMaxi, oldMaxj;
	double theta;
	double **S1, **temp; // **tp;
	
	// allocates S1 and a temp matrix
	allocMat(&S1, n, n);
	allocMat(&temp, n, n);

	// finds first set of max indices and theta values and initialses S1
	findMaxIndex(D, &maxi, &maxj, n);
	oldMaxi = maxi; oldMaxj = maxj;
	theta = findTheta(D, maxi, maxj);
	fillS1(S1, maxi, maxj, theta, n);
		
	int counter = 0;	

	// loops until D is practically diagonal
	while(fabs(D[maxi][maxj]) > 1.E-6){


		// calculates S = S * S1
		multMat(*S, S1, temp, n, n, n);
		/*
		tp = *S;
		*S = temp;
		temp = tp;
		*/
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				(*S)[i][j] = temp[i][j];
			}
		}
	
		// calculates D = S1t * D * S1
		multMat(D, S1, temp, n, n, n);
		S1[maxi][maxj] *= -1;
		S1[maxj][maxi] *= -1;
		multMat(S1, temp, D, n, n, n);


		counter ++;

		// updates indices and theta etc
		findMaxIndex(D, &maxi, &maxj, n);
		theta = findTheta(D, maxi, maxj);
		fillS1(S1, maxi, maxj, theta, n);

		
	/*
 		A LOOOOTTT SLOWER FOR SOME REASON

		if(oldMaxi != maxi || oldMaxj != maxj){
			theta = findTheta(D, maxi, maxj);
			updateS1(S1, maxi, maxj, oldMaxi, oldMaxj, theta);
			oldMaxi = maxi;
			oldMaxj = maxj;
		}
		else{
			S1[maxi][maxj] *= -1;
			S1[maxj][maxi] *= -1;

		}
	*/
		
	}

	// frees matrices used
	freeMat(S1);
	freeMat(temp);


	return;


}
