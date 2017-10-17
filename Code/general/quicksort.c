#include "quicksort.h"

// Quicksort function adapted from Marina Marinkovic's MA5613 2016/2017 TCD Dublin
void 	recursiveSortArray(int *arr, int low, int high){
	int pivot, i, j, temp;
	if(low < high){
		pivot = low;
		i = low;
		j = high;
		while(i < j){
			while(i <= high && arr[i] <= arr[pivot]){
				i ++;
			}
			while(j >= low && arr[j] > arr[pivot]){
				j --;
			}
			if(i < j){
				temp = arr[i];
				arr[i] = arr[j];
				arr[j] = temp;
			}
		}
		temp = arr[j];
		arr[j] = arr[pivot];
		arr[pivot] = temp;
		recursiveSortArray(arr, low, j-1);
		recursiveSortArray(arr, j+1, high);
	}

	return;
}



void 	quicksortArray(int *arr, int n){
	recursiveSortArray(arr, 0, n-1);

	return;
}



// Quicksort function adapted from Marina Marinkovic's MA5613 2016/2017 TCD Dublin
void 	recursiveSortMatrix(int **A, int cols, int low, int high, int sortCol){
	int pivot, i, j, temp, loop;
	if(low < high){
		pivot = low;
		i = low;
		j = high;
		while(i < j){
			while(i <= high && A[i][sortCol] <= A[pivot][sortCol]){
				i ++;
			}
			while(j >= low && A[j][sortCol] > A[pivot][sortCol]){
				j --;
			}
			if(i < j){
				for(loop = 0; loop < cols; loop ++){
					temp = A[i][loop];
					A[i][loop] = A[j][loop];
					A[j][loop] = temp;
				}
			}
		}
		for(loop = 0; loop < cols; loop ++){
			temp = A[j][loop];
			A[j][loop] = A[pivot][loop];
			A[pivot][loop] = temp;
		}
		recursiveSortMatrix(A, cols, low, j-1, sortCol);
		recursiveSortMatrix(A, cols, j+1, high, sortCol);
	}

	return;
}



void 	quicksortMat(int **A, int n, int m, int sortCol){
	recursiveSortMatrix(A, m, 0, n-1, sortCol);

	return;
}



#ifdef MAIN
int main(){
	int i;
	int seed = 1501762216;//time(NULL);
	seed = time(NULL);
	srand(seed);
	printf("seed is %d\n", seed);

	int array[9]={1,12,5,26,7,14,3,7,2};
	for(i = 0; i < 9; i++){
		array[i] = (int) (rand() % 20);
	}
	
	quicksortArray(array, 9);
	for(i=0;i<9;i++){
		printf("%d ",array[i]);
	}
	printf("\n\n");


	int **mat;
	int n = 10;	
	int m = 3;
	allocMatInt(&mat, n, m);
	int j;
	for(i = 0; i < n; i++){
		mat[i][0] = (int) (rand() % 20);
		mat[i][1] = (int) (rand() % 20);
		for(j = 2; j < m; j++){
			mat[i][j] = (int) (rand() % 2 + 1);
		}
	} 

	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			printf("%2d  ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n"); 

	quicksortMat(mat, n, m, 2);

	int loc = 0;
	while(mat[loc][2] == 1){
		loc ++;
	}
	printf("loc is %d\n", loc);













	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			printf("%2d  ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
#endif
