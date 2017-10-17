/*	Functions relating to matrices and some arrays for ints and doubles
 *		Assigns memory contiguously for MPI_Send/MPI_Recv
 *		Frees memory
 *		Resets to zero
 *		Multiply, transpose, swap cols etc
 *		More efficient methods of specific multiplications
 */



#ifndef MATRIXH
#define MATRIXH

// Some functions dealing with matrix operations


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

void 	allocMat(double ***A, int n, int m);
void	freeMat(double **A);
void 	fillMat(double ***A, int *n, int *m, char *filename);
void 	multMat(double **A, double **B, double **C, int n, int m, int q);
void 	printMat(double **A, int n, int m);
void 	swapCol(double **A, int c1, int c2, int n, int m);
void 	swapDiag(double **A, int d1, int d2);
void 	transpose(double **A, double **At, int n, int m);
void 	flipSignColMat(double **A, int c, int n, int m);
void 	normMat(double **A, int n, int m);
void 	calcAAt(double **A, double **AAt, int n, int m);
void 	calcAtA(double **A, double **AtA, int n, int m);
void 	sumRowOneMat(double **A, int r, int n, int m);




void 	allocMatInt(int ***A, int n, int m);
void 	callocMatInt(int ***A, int n, int m);
void	freeMatInt(int **A);
void	callocArrInt(int **A, int n);
void	zeroMatInt(int **A, int n, int m);
void	zeroArrInt(int *A, int n);


#endif
