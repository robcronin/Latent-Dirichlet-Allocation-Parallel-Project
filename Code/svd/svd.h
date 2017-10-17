// Runs SVD decomposition for various methods


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../general/matrix.h"
#include "jacobiEigen.h"
#include "../general/extra.h"

void 	printSvdResults(double **A, double **U, double **V, double **S, int n, int m, int verbose);
void 	svdJacobi(double **A, double **U, double **V, double **S, int n, int m, int verbose, int errorCheck);


// void 	reorderEigen(double **A, double **D, int n);
