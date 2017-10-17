/*	Functions associated with AD-LDA method of LDA
 *	(Re-uses many serial LDA functions where applicable)
 *		Parallel or Serial read in of documents
 *		Distribution of docMatrix, word counts
 *		Creating counters
 *		Altering communicators
 *		Gibbs Sampling via AD-LDA
 */


#ifndef GIBBSMPI
#define GIBBSMPI



#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../Gibbs-LDA2/gibbsLDA2.h"
#include "../general/extra.h"
#include "../general/matrix.h"




void	filesPerProc(int **localFiles, int rank, int size, int noFiles, int *locNoFiles, MPI_Comm mainComm);

void	alterComm(int rank, int *size, int noFiles, MPI_Comm *mainComm);

void	readInDocsSerial(int rank, int size, char **pathNames, int noFiles, char *stopfile, int *vocabSize, char ***vocabulary, int printResults, int **locUniqueCount, int **locTotalCount, int locNoFiles, int *localFiles, int ****localDocMatrix, int *totalWords, int gutenFiles, MPI_Comm mainComm);

void	distWordCounts(int **locUniqueCount, int **locTotalCount, int *uniqueWordCounts, int *totalWordCounts, int locNoFiles, int rank, int size, int *localFiles, int **totalUniques, MPI_Comm mainComm); 

void	distDocMatrix(int ****localDocMatrix, int locNoFiles, int *localFiles, int *locUniqueCount, int rank, int size, int ***wholeDocMatrix, int *totalUniques, MPI_Comm mainComm);

void	readInDocsParallel(int rank, int size, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults, char **pathNames, int *localFiles, int locNoFiles, int *totalWords, int gutenFiles, MPI_Comm mainComm);

void	readInDocsParallelDict(int rank, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults, char **pathNames, int *localFiles, int locNoFiles, int *totalWords, int gutenFiles, char *dictfile, MPI_Comm mainComm);

void	sendPhi(struct alpha *phi, int src, int dest, int rank, MPI_Comm mainComm);

void	bcastPhi(struct alpha *phi, int src, int rank, MPI_Comm mainComm);

void	initialiseCountersMPI(int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **topicCountPerDoc, int *topicCountTotal, int *localTopicCountTotal, int **wordCountPerTopic, int **localWordCountPerTopic);

void	gibbsSamplerMPI(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int ***topicIndex, int noTopics, int noFiles, int ***docMatrix, int *uniqueWordCounts, double *alphaDir, double *gammaDir, double gammaSum, int iterations, int **localWordCountPerTopic, int *localTopicCountTotal, int vocabSize, MPI_Comm mainComm);


#endif
