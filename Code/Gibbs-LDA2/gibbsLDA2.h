/* Functions associated with running Gibbs Sampling version of Latent Dirichlet Allocation(LDA)
 * 	Create vocab, docMatrix, topicIndex, Counters, Dirichlet priors
 * 	Freeing of these
 * 	Gibbs sampling function
 * 	Print results
 * 	Information function that prints to screen
 */


#ifndef GIBBSLDA
#define GIBBSLDA


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../tfidf/tfidf.h"
#include "../general/matrix.h"
#include "../general/extra.h"
#include "../general/quicksort.h"

void 	createVocab(struct alpha *phi, char ***vocabulary);
void	freeVocab(char ***vocabulary, int vocabSize);
void	allocDocMatrix(int ****docMatrix, int noFiles, int *uniqueWordCounts);
void	freeDocMatrix(int ****docMatrix, int noFiles);
void	fillDocMatrix(int ***docMatrix, char **vocabulary, struct alpha **all, int noFiles);
void	fillDocMatrixNew(int ***docMatrix, char **vocabulary, struct alpha **all, int noFiles, struct alpha *phi);
void	gibbsSampler2(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int ***topicIndex, int noTopics, int noFiles, int ***docMatrix, int *uniqueWordCounts, double *alphaDir, double *gammaDir, double gammaSum, int iterations);


void 	printTopXLDA2(double **phiMatrix, int noTopics, int vocabSize, char **vocabulary, int **wordCountPerTopic, int x);

void	createPhi(double ***phiMatrix, int vocabSize, int noTopics, int **wordCountPerTopic, double *gammaDir);

void 	createTopicIndex(int ****topicIndex, int *totalWordCounts, int noFiles, int *uniqueWordCounts, int ***docMatrix);
void	freeTopicIndex(int ***topicIndex, int noFiles);
void 	createCounters2(int ***topicCountPerDoc, int **topicCountTotal, int ***wordCountPerTopic, int noFiles, int noTopics, int vocabSize);
void	freeCounters2(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic);
void	createWordCounts(int **uniqueWordCounts, int **totalWordCounts, int noFiles, struct alpha **all);
void 	initialiseCounters(int ***topicIndex, int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int noFiles, int noTopics, int *uniqueWordCounts, int ***docMatrix);
void	initialiseAlpha(double **alpha, int noTopics);
void	initialiseGamma(double **gamma, int noWords, double *gammaSum);
void	printInfoLDA(int rank, int index, int LDAmethod, int seed, char *stopfile, int noTopics, int noFiles, int vocabSize, int iterations, int size, int parallelReadIn, int totalWords, int pdSize, int pwSize, int CwordMethod, int rma, int threads, int bundleSize, int invMethod, int sortRead);


#endif
