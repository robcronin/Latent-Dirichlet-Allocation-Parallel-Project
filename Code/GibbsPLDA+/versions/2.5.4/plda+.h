#ifndef PLDAPLUS
#define PLDAPLUS


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

#include "gibbsplda+.h"
#include "pthreadFunctions.h"



void	procAssignment(int rank, int size, int noFiles, MPI_Comm *pdComm, MPI_Comm *pwComm, double pwToPdRatio, int *pdSize, int *pwSize, int *prank, char *procType, MPI_Comm mainComm);
void	fillM(char procType, int vocabSize, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int *wordSpread, MPI_Comm pdComm);
void	distWordsToPw(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int *wordLocation, int noTopics, struct treeRootCword **wordCountPerTopicTree, MPI_Comm pdComm, MPI_Comm mainComm);
void	distWordsToPwMatrix(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int **wordLocationMat, int noTopics, int ***wordCountPerTopic, int *CwordSize, int wordLocCols, MPI_Comm pdComm, MPI_Comm mainComm);
void	createComms(MPI_Comm **commArr, int pwSize, int pdSize, MPI_Comm mainComm);
void	createWindowsNew(MPI_Win *wordCountWin, MPI_Win *topicCountWin, char procType, MPI_Comm mainComm, int CwordSize, int noTopics, int **wordCountPerTopic, int *localTopicCountTotal);
void	buildInvertedIndex(char procType, struct treeRootInv **invertedIndexTree, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int prank, MPI_Comm pdComm);
void	buildInvertedIndexMatrix(char procType, int ***invertedIndexMatrix, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int prank, int vocabSize, int **wordLocationMat, int *invertedIndexSize, MPI_Comm pdComm);
void	readInDocsAll(char procType, int prank, int pdSize, int noFiles, int *locNoFiles, int parallelReadIn, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults,  char **pathNames, int *totalWords, struct timeval *t1Sec, struct timeval *t2Sec, int gutenFiles, MPI_Comm pdComm);
void	readInDocsAllNew(char procType, int prank, int pdSize, int noFiles, int *locNoFiles, int parallelReadIn, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults,  char **pathNames, int *totalWords, struct timeval *t1Sec, struct timeval *t2Sec, int gutenFiles, char *dictfile,  MPI_Comm pdComm);
void	initTopicsPLDAplus(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int *wordLocation, int **topicCountPerDoc, int pdSize, int size, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm);
void	initTopicsPLDAplusMatrix(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int **wordCountPerTopic, int *localTopicCountTotal, int rank, MPI_Comm mainComm);
void	initTopicsPLDAplusMatrixRMA(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win **winArr, MPI_Comm mainComm);
void	initTopicsPLDAplusMatrixRMAFixed(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm);
void	createWordBundles(char procType, struct treeRootInv *invertedIndexTree, int *wordLocation, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int *totNoBundles, int prank, MPI_Comm pdComm);
void	createWordBundlesMatrix(char procType, struct treeRootInv *invertedIndexTree, int **wordLocationMat, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int ***bundles2dConv, int *totNoBundles, int prank, MPI_Comm pdComm);
void	createWordBundlesMatrixInv(char procType, int **invertedIndexMatrix, int **wordLocationMat, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int ***bundles2dConv, int *totNoBundles, int prank, int invertedIndexSize, int locNoFiles, MPI_Comm pdComm);
void	circleQueue(char procType, int vocabSize, int pdSize, int prank, int *circleStartBundle, int **bundles2d);
void	collectAndPrint(int rank, int vocabSize, int noTopics, int pwSize, char procType, struct treeRootCword *wordCountPerTopicTree, int printResults, double *betaDir, char **vocabulary, MPI_Comm mainComm);
void	collectAndPrintMatrix(int rank, int vocabSize, int noTopics, char procType, int **wordCountPerTopic, int printResults, double *betaDir, char **vocabulary, int **wordLocationMat, int pdSize, int size, MPI_Comm mainComm);



#endif
