/*	Quicksort functions for array and matrix
 *		Matrices can be sorted by any column
 *	Adapted from Marina Marinkovic's MA5613 2016/2017 TCD Dublin
 */



#ifndef QUICKSORT
#define QUICKSORT


#include <stdio.h>
#include "matrix.h"

void 	quicksortArr(int *arr, int n);
void	quicksortMat(int **A, int n, int m, int sortCol);


#endif
