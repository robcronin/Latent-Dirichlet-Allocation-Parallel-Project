/*	All the various implementations of Gibbs Sampling for PLDA+
 *		Standard method with BST
 *		Condensed matrix for Cword
 *		Condensed matrix for Cword and Inverted Index
 *		RMA
 *		Various "beta" versions for pthreads
 */



#ifndef GIBBSPLDAPLUS
#define GIBBSPLDAPLUS


#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "../binaryTree/binaryTreeCword.h"
#include "../binaryTree/binaryTreeInv.h"
#include "pthreadFunctions.h"
#include "../general/extra.h"



void    runPLDAplusGibbs(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm);

void    runPLDAplusGibbsMatrix(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm mainComm);

void    runPLDAplusGibbsMatrixInv(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, int **invertedIndexMatrix, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **wordLocationMat, MPI_Comm mainComm);

void    runPLDAplusGibbsMatrixRMA(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm);

void    runPLDAplusGibbsMatrixInvRMA(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, int **invertedIndexMatrix, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, int **wordLocationMat, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm);


// ISSUES STILL EXIST WITH THIS FUNCTION, IF RUN FOR LONG ENOUGH IT WILL HANG
void    runPLDAplusGibbsPthreadsWorkload(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm);




// THE FOLLOWING FUNCTIONS STILL CONTAINS SOME BUGS AND WILL OCCASSIONALLY SEGFAULT (see runPLDAplusGibbsPthreadsWorkload for best pthreads version)
void    runPLDAplusGibbsPthreads(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm);

void    runPLDAplusGibbsPthreadsBoth(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm);

void    runPLDAplusGibbsPthreadsDeadline(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *gammaDir, double gammaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm);




#endif
