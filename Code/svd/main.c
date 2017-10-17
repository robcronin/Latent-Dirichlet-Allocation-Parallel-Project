#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../general/matrix.h"
#include "jacobiEigen.h"
#include "../general/extra.h"
#include "svd.h"

void parseArgs(int argc, char *argv[], int *n, int *m, int *verbose, char **filename, int *errorCheck, int *seed);
void printHelp();


int main(int argc, char *argv[]){

	//set default arguments
	int n = 4;
	int m = 4;
	int verbose = 0;
	int errorCheck = 0;
	char *filename = NULL;
	int seed = time(NULL);

	parseArgs(argc, argv, &n, &m, &verbose, &filename, &errorCheck, &seed);
	printf("Seed is %d\n", seed);
	srand48(seed);


	double **A, **U, **V, **S;
	fillMat(&A, &n, &m, filename);
	allocMat(&S, n, m);
	allocMat(&U, n, n);
	allocMat(&V, m, m);

	svdJacobi(A, U, V, S, n, m, verbose, errorCheck);


	freeMat(A);freeMat(S);freeMat(U);freeMat(V);

	return 0;
}


// parses command line arguments
void parseArgs(int argc, char *argv[], int *n, int *m, int *verbose, char **filename, int *errorCheck, int *seed){
	int opt;
	while((opt=getopt(argc,argv,"ef:hn:m:s:v"))!=-1){
		switch(opt){
			case 'e':
				*errorCheck = 1;
				break;
			case 'f':
				*filename =  optarg;
				break;
			case 'h':
				printHelp();
				exit(1);
			case 'n':
				*n = atoi(optarg);
				break;
			case 'm':
				*m = atoi(optarg);
				break;
			case 's':
				*seed = atoi(optarg);
				break;
			case 'v':
				*verbose = 1;
				break;
			default:
				printHelp();
				exit(EXIT_FAILURE);
		}
	}
	return;
}

void printHelp(){
	printf("-e\t\tRuns error checking\n");
	printf("-f [FILE]\tUses matrix in given file\n");
	printf("-h\t\tPrints help\n");
	printf("-n [INT]\tSets rows of matrix (default: 4)\n");
	printf("-m [INT]\tSets cols of matrix (default: 4)\n");
	printf("-s [INT]\tSets seed (default: random)\n");
	printf("-v\t\tTurns on verbose\n");

	return;
}
