// Functions for running the Jacobi method to find eigenvalues and eigenvectors
//Method mostly adopted from the following paper
//	http://www.cmi.ac.in/~ksutar/NLA2013/iterativemethods.pdf


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../general/matrix.h"

void 	jacobiMet(double **A, double ***S, double **D, int n);



// void 	findMaxIndex(double **D, int *maxi, int *maxj, int n);
// double 	findTheta(double **D, int i, int j);
// void 	fillS1(double **S1, int maxi, int maxj, double theta, int n);
// void		updateS1(double **S1, int maxi, int maxj, int oldMaxi, int oldMaxj, double theta);
