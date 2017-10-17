#ifndef PTHREADFCNS
#define PTHREADFCNS



#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "../general/matrix.h"
#include "../binaryTree/binaryTreeInv.h"

// 	structure to hold arguments for pthreads
struct	argStruct{
	int 	*stop;
	int 	stop2;
	int	*topicCountTotal;
	int 	noTopics;
	int	pdSize;
	int	size;
	int	**queueMat;
	int	queueSize;
	int	**bundles2d;
	int	bundleSize;
	int	**bundles2dConv;
	int	**CwordBuffer;
	int	**CtopicBuffer;
	int	locNoFiles;
	int	***localDocMatrix;
	int	***topicIndex;
	int	**topicCountPerDoc;
	double	*alphaDir;
	double	*betaDir;
	double	betaSum;
	int 	iterations;
	int 	totNoBundles;
	int	rank;
	int	*localTopicCountTotalAct;
	int	**wordCountPerTopic;
	int	sleepTime;
	int	*deadlines;
	int	*workload;
	int	overload;
	int	**invertedIndexMatrix;
	int	**wordLocationMat;
	int	loopsPerSample;
	int	*topicSignal;
	pthread_mutex_t *locks;
	struct treeRootInv *invertedIndexTree;
	MPI_Comm mainComm;
};



void	*fetchCtopic(void *arg);
void	*fetchCword(void *arg);
void	*runGibbs(void *arg);
void	*updateCword(void *arg);
void	*recvFetchCtopic(void *arg);
void	*recvFetchCword(void *arg);
void	*recvUpdateCword(void *arg);
void	*stopProg(void *arg);
void	*fetchCtopicDeadline(void *arg);
void	*fetchCwordDeadline(void *arg);
void	*runGibbsDeadline(void *arg);
void	*updateCwordDeadline(void *arg);
void	*fetchCtopicDeadlineNew(void *arg);
void	*fetchCwordDeadlineNew(void *arg);
void	*runGibbsDeadlineNew(void *arg);
void	*updateCwordDeadlineNew(void *arg);



#endif
