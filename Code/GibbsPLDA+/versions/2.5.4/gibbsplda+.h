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

#include "../Gibbs-LDA-MPI/gibbsLDA-MPI.h"
#include "../binaryTree/binaryTreeCword.h"
#include "../binaryTree/binaryTreeInv.h"
#include "../general/quicksort.h"
#include "../linkedList/linkedList.h"
#include "pthreadFunctions.h"



void    runPLDAplusGibbs(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm);
void    runPLDAplusGibbsMatrix(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm mainComm);
void    runPLDAplusGibbsMatrixInv(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, int **invertedIndexMatrix, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **wordLocationMat, MPI_Comm mainComm);
void    runPLDAplusGibbsMatrixRMA(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, MPI_Win **winArr, MPI_Comm mainComm);
void    runPLDAplusGibbsMatrixRMAFixed(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm);
void    runPLDAplusGibbsPthreads(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm);
void    runPLDAplusGibbsPthreadsBoth(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm);
void    runPLDAplusGibbsPthreadsDeadline(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm);
void    runPLDAplusGibbsPthreadsDeadlineNew(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm);



#endif
