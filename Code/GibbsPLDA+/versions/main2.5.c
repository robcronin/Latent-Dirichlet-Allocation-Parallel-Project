/*
 *
 *
 * 		STILL NEED TO CHECK AND FREE ALL MALLOCS
 * 		STILL GETTING ERRORS WITHOUT ACCTIME DUE TO WRONG TAGS RECEIVING
 *
 */



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

#define N 10000
#define	M 10000
#define Z -1023

// 	structure to hold arguments for pthreads
struct	argStruct{
	int 	*stop;
	int 	stop2;
	int		*topicCountTotal;
	int 	noTopics;
	int		pdSize;
	int		size;
	int		**queueMat;
	int		queueSize;
	int		**bundles2d;
	int		bundleSize;
	int		**bundles2dConv;
	int		**CwordBuffer;
	int		**CtopicBuffer;
	int		locNoFiles;
	int		***localDocMatrix;
	int		***topicIndex;
	int		**topicCountPerDoc;
	double	*alphaDir;
	double	*betaDir;
	double	betaSum;
	int 	iterations;
	int 	totNoBundles;
	int		rank;
	int		*localTopicCountTotalAct;
	int		**wordCountPerTopic;
	int		sleepTime;
	int		*deadlines;
	int		*workload;
	int		overload;
	int		**invertedIndexMatrix;
	int		**wordLocationMat;
	int		loopsPerSample;
	int		*topicSignal;
	pthread_mutex_t *locks;
	struct treeRootInv *invertedIndexTree;
	MPI_Comm mainComm;
};



void 	printHelp();
void 	parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma, int *bundleSize, double *pwToPdRatio, int *invMethod, int *sortRead);



//	function to assign processors as Pd or Pw
void	procAssignment(int rank, int size, int noFiles, MPI_Comm *pdComm, MPI_Comm *pwComm, double pwToPdRatio, int *pdSize, int *pwSize, int *prank, char *procType, MPI_Comm mainComm){

	// returns error if only one processor as need at least one of each
	if(size < 2){
		printf("*** ERROR ***\nPLDA+ requires at least 2 processes\n\n");
		MPI_Finalize();
		exit(1);
	}

	// calculates number of each based on ratio
	*pwSize = (int) rint(pwToPdRatio * size);
	*pdSize = size - *pwSize;

	// reassigns if less files than Pd processes
	if(*pdSize > noFiles){
		*pdSize = noFiles;
		*pwSize = size - *pdSize;
	}

	// ensures neither are put at 0
	if(*pdSize == 0){
		*pdSize = 1;
		*pwSize = size - 1;
	}
	else if(*pwSize == 0){
		*pwSize = 1;
		*pdSize = size - 1;
	}

	// creates group of all processes
	MPI_Group origGroup, pdGroup, pwGroup;
	MPI_Comm_group(mainComm, &origGroup);

	// creates rank matrix to assign first pdSize processes as Pd and rest as Pw
	int *ranks = malloc((*pdSize) * sizeof(int));		// FREED
	int i;
	for(i = 0; i < *pdSize; i++){
			ranks[i] = i;
	}

	// create pd and pw Groups
	MPI_Group_incl(origGroup, *pdSize, ranks, &pdGroup);
	MPI_Group_excl(origGroup, *pdSize, ranks, &pwGroup);
	free(ranks);

	// create pd and pw Comms
	MPI_Comm_create(mainComm, pdGroup, pdComm);
	MPI_Comm_create(mainComm, pwGroup, pwComm);

	// assign letter identifier to each process
	if(rank < *pdSize){
		MPI_Comm_rank(*pdComm, prank);
		*procType = 'd';
	}
	else{
		MPI_Comm_rank(*pwComm, prank);
		*procType = 'w';
	}	
		
	return;
}



// fills the wordSpread array (known as m in PLDA+ paper)
void	fillM(char procType, int vocabSize, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int *wordSpread, MPI_Comm pdComm){

	// only run on Pd processes
	if(procType == 'd'){
		int doc, word;
		int *locWordSpread;
		callocArrInt(&locWordSpread, vocabSize);		// FREED

		// loops through every word in local documents and marks them as present on local m
		for(doc = 0; doc < locNoFiles; doc++){
			for(word = 0; word < locUniqueCount[doc]; word++){
				locWordSpread[ localDocMatrix[doc][word][0] ] = 1;
			}
		}

		// reduces across all Pds to create a global m
		MPI_Allreduce(locWordSpread, wordSpread, vocabSize, MPI_INT, MPI_SUM, pdComm);
		free(locWordSpread);
	}
	
	return;
}


// Distributes words from Pds to be stored on various Pws (BST METHOD)
void	distWordsToPw(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int *wordLocation, int noTopics, struct treeRootCword **wordCountPerTopicTree, MPI_Comm pdComm, MPI_Comm mainComm){ 

	struct timeval t1Sec;
	startTimeMPI(&t1Sec, mainComm);

	// All words need to be assigned to a specific Pw
	// 		- each Pd distributes a portion of the total words not its own words

	// DISTIBUTE WORDS
	if(procType == 'd'){

		// splits the words evenly across Pds
		int gap = (int)vocabSize/pdSize;
		int startPoint = rank*gap;
		int endPoint = (rank+1)*gap;
		if(prank == pdSize - 1){endPoint = vocabSize;}
			
		int len, word;
		int *locWordLocation;
		callocArrInt(&locWordLocation, vocabSize);		// FREED

		// starts at a different Pw for each Pd
		int dest = prank % pwSize + pdSize;

		// distribute from largest possible spread down to ensure each Pw has roughly even workload
		for(len = pdSize; len >= 0; len --){
			// loops through the words and finds the one with the current len
			// Not the most efficient method but quick function in overall program
			for(word = startPoint; word < endPoint; word ++){
				if(wordSpread[word] == len){

					// sends that word index to the current dest
					MPI_Send(&word, 1, MPI_INT, dest, 2, mainComm);
					locWordLocation[word] = dest;		// marks the location off
						
					// loops through the Pws
					dest ++;
					if(dest == size){dest = pdSize;}		
				}
				// catches any errors
				else if(wordSpread[word] < 0 || wordSpread[word] > pdSize){
					printf("%d is %d\n", word, wordSpread[word]);
				}
			}
		}
		// send finish signal to all Pw's
		int i;
		for(i = pdSize; i < size; i++){
			MPI_Send(&word, 1, MPI_INT, i, 0, mainComm);
		}
	
		// consolidate wordLocation across all Pds
		MPI_Allreduce(locWordLocation, wordLocation, vocabSize, MPI_INT, MPI_SUM, pdComm);
		free(locWordLocation);
	}
	// RECEIVE WORDS
	else if(procType == 'w'){
		int index;
		MPI_Status stat;

		// creates a BST to store the word topic counts for its designated words
		createTreeCword(wordCountPerTopicTree, noTopics);
		int noFinished = 0;
		int total = 0;

		// receives from any Pd until all have sent finish signals
		while(noFinished < pdSize){
			// receives word and adds it to the tree
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 2){
				addValueCword(*wordCountPerTopicTree, index);
				total ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// checks for receiving wrong tag
			else{
				printf("\t\t\tERROR DIST WORDS RECEIVING WRONG THING\n");
			}
		}
		// rebalances the BST for better recall performance
		rebalanceTreeCword(*wordCountPerTopicTree);

		#ifdef DEBUG
		int totalTreeSize;
		MPI_Reduce(&((*wordCountPerTopicTree) -> size), &totalTreeSize, 1, MPI_INT, MPI_SUM, 0, pwComm);
		if(prank == 0){printf("Totoal tree is %d\n", totalTreeSize);}
		#endif
	}

	endTimeMPI(t1Sec, "Cw TREE (Pw)", rank, mainComm);

	return;
}




// Distributes words from Pds to be stored on various Pws (CONDENSED MATRIX METHOD)
void	distWordsToPwMatrix(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int **wordLocationMat, int noTopics, int ***wordCountPerTopic, int *CwordSize, int wordLocCols, MPI_Comm pdComm, MPI_Comm mainComm){ 

	struct timeval t1Sec;
	startTimeMPI(&t1Sec, mainComm);
	MPI_Status stat;

	
	// All words need to be assigned to a specific Pw
	// 		- each Pd distributes a portion of the total words not its own words

	// DISTIBUTE WORDS
	if(procType == 'd'){

		// splits the words evenly across Pds
		int gap = (int)vocabSize/pdSize;
		int startPoint = rank*gap;
		int endPoint = (rank+1)*gap;
		if(prank == pdSize - 1){endPoint = vocabSize;}
			
		int len, word;
		int **locWordLocation;
		callocMatInt(&locWordLocation, vocabSize, wordLocCols);		// FREED

		// starts at a different Pw for each Pd
		int dest = prank % pwSize + pdSize;

		// distribute from largest spread down
		for(len = pdSize; len >= 0; len --){

			// loops through the words and finds the one with the current len
			// Not the most efficient method but quick function in overall program
			for(word = startPoint; word < endPoint; word ++){
				if(wordSpread[word] == len){

					// sends that word index to the current dest
					MPI_Send(&word, 1, MPI_INT, dest, 2, mainComm);
		
					// marks the rank dest in the first column
					locWordLocation[word][0] = dest;

					// Receives the condensed matrix location in the 2nd column
					MPI_Recv(&(locWordLocation[word][1]), 1, MPI_INT, dest, 2, mainComm, &stat);

					// loops through Pws to send to
					dest ++;
					if(dest == size){dest = pdSize;}		
				}
				// catches any errors
				else if(wordSpread[word] < 0 || wordSpread[word] > pdSize){
					printf("%d is %d\n", word, wordSpread[word]);
				}
			}
		}
		// send finish signal to all Pw's
		int i;
		for(i = pdSize; i < size; i++){
			MPI_Send(&word, 1, MPI_INT, i, 0, mainComm);
		}
	
		// consolidate wordLocation across all Pds
		MPI_Allreduce(locWordLocation[0], wordLocationMat[0], vocabSize*wordLocCols, MPI_INT, MPI_SUM, pdComm);
		freeMatInt(locWordLocation);
	}
	// RECEIVE WORDS
	else if(procType == 'w'){
		int index;
		int noFinished = 0;
		int matIndex = 0;

		// receives from any Pd until all have sent finish signals
		while(noFinished < pdSize){

			// receives a request to store a word
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			
			// sends back a location index for that word
			if(stat.MPI_TAG == 2){
				MPI_Send(&matIndex, 1, MPI_INT, stat.MPI_SOURCE, 2, mainComm);
				matIndex ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			else{
				printf("\t\t\tERROR DIST WORDS RECEIVING WRONG THING\n");
			}
		}

		// Create Cword matrix for the number of words it received requests for
		callocMatInt(wordCountPerTopic, matIndex, noTopics);
		*CwordSize = matIndex;

	}

	endTimeMPI(t1Sec, "Cw MATRIX (Pw)", rank, mainComm);

	return;
}



// NOT USED ANYMORE FIXME
void	createComms(MPI_Comm **commArr, int pwSize, int pdSize, MPI_Comm mainComm){
	int i;

	// orig group
	MPI_Group origGroup;
	MPI_Comm_group(mainComm, &origGroup);

	// create comm array
	*commArr = malloc(pwSize * sizeof(MPI_Comm));

	// fill in all Pd's into rank
	int *ranks = malloc((pdSize + 1)*sizeof(int));
	for(i = 0; i < pdSize; i++){	
		ranks[i] = i;
	}
	for(i = 0; i < pwSize; i++){

		// create new group
		ranks[pdSize] = i + pdSize;
		MPI_Group	newGroup;
		MPI_Group_incl(origGroup, pdSize + 1, ranks, &newGroup);

		// create comm
		MPI_Comm_create(mainComm, newGroup, &((*commArr)[i]));

		MPI_Group_free(&newGroup);
	}

	free(ranks);
	MPI_Group_free(&origGroup);
	return;
}
	


// NOT USED ANYMORE FIXME
void	createWindows(MPI_Win ***winArr, int pwSize, char procType, MPI_Comm *commArr, int prank, int CwordSize, int noTopics, int **wordCountPerTopic, int *localTopicCountTotal){
	
	int i;
	// create window matrix
	MPI_Win *temp = malloc(pwSize * 2 * sizeof(MPI_Win));
	*winArr = malloc(pwSize * sizeof(MPI_Win*));
	for(i = 0; i < pwSize; i++){
		(*winArr)[i] = &temp[i*2];
	}

	for(i = 0; i < pwSize; i++){
		if(procType == 'd'){
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, commArr[i], &((*winArr)[i][0])); 
			MPI_Win_fence(0, (*winArr)[i][0]);
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, commArr[i], &((*winArr)[i][1])); 
			MPI_Win_fence(0, (*winArr)[i][1]);
		}
		else if(procType == 'w' && prank == i){
			MPI_Win_create(&(wordCountPerTopic[0][0]), CwordSize*noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, commArr[i], &((*winArr)[i][0]));	
			MPI_Win_fence(0, (*winArr)[i][0]);
			MPI_Win_create(localTopicCountTotal, noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, commArr[i], &((*winArr)[i][1]));	
			MPI_Win_fence(0, (*winArr)[i][1]);
		}		
	}

	return;
}


	

// FIXME NAME
// creates two RMA windows for the data stored on Pw
// 		- one for the Cword matrix and one for the local Ctopic
void	createWindowsNew(MPI_Win *wordCountWin, MPI_Win *topicCountWin, char procType, MPI_Comm mainComm, int CwordSize, int noTopics, int **wordCountPerTopic, int *localTopicCountTotal){
	
	// Pd creates the two windows but doesn't expose any data (will be using MPI_Get)
	if(procType == 'd'){
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, mainComm, wordCountWin); 
			MPI_Win_fence(0, *wordCountWin);
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, mainComm, topicCountWin); 
			MPI_Win_fence(0, *topicCountWin);
	}
	// Pw also creates the two windows but exposes the relevant memory
	else if(procType == 'w'){
			MPI_Win_create(&(wordCountPerTopic[0][0]), CwordSize*noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, mainComm, wordCountWin);	
			MPI_Win_fence(0, *wordCountWin);
			MPI_Win_create(localTopicCountTotal, noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, mainComm, topicCountWin);	
			MPI_Win_fence(0, *topicCountWin);
	}		

	return;
}


	
// builds the inverted index BST on Pd
void	buildInvertedIndex(char procType, struct treeRootInv **invertedIndexTree, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int prank, MPI_Comm pdComm){
	struct timeval t1Sec;


	// Pd word organisation
	if(procType == 'd'){
		startTimeMPI(&t1Sec, pdComm);

		// set up tree
		createTreeInv(invertedIndexTree, locNoFiles);

		// fill in tree from DocMatrix
		int doc, word;
		for(doc = 0; doc < locNoFiles; doc ++){
			for(word = 0; word < locUniqueCount[doc]; word ++){
				addValueInv(*invertedIndexTree, localDocMatrix[doc][word][0], doc, word, localDocMatrix[doc][word][1]);
			}
			// rebalances after each document
			rebalanceTreeInv(*invertedIndexTree);
		}

		#ifdef DEBUG
		if(prank == 0){
			int testIndex = localDocMatrix[2][435][0];
			struct nodeInv *search;
			search = searchTreeInv(*invertedIndexTree, testIndex);
			for(doc = 0; doc < locNoFiles; doc++){
				printf("%d\n", search->invIndex[doc]);
			}
		}
		#endif
		endTimeMPI(t1Sec, "INV INDEX (Pd)", prank, pdComm);
		
	}

	return;
}



void	buildInvertedIndexMatrix(char procType, int ***invertedIndexMatrix, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int prank, int vocabSize, int **wordLocationMat, int *invertedIndexSize, MPI_Comm pdComm){
	struct timeval t1Sec;

	// Pd word organisation
	if(procType == 'd'){
		startTimeMPI(&t1Sec, pdComm);

		// find local number words
		int *vocabMarker;
		callocArrInt(&vocabMarker, vocabSize);
		int totNoWords = 0;
		int doc, word;
		for(doc = 0; doc < locNoFiles; doc ++){
			for(word = 0; word < locUniqueCount[doc]; word ++){
				if( vocabMarker[ localDocMatrix[doc][word][0] ] == 0){
					totNoWords ++;
					// mark the first document used (+1 so as not to confuse with 0)
					vocabMarker[ localDocMatrix[doc][word][0] ] = doc + 1;
				}
			}
		}
	
		*invertedIndexSize = totNoWords;

		// set up matrix
		// first locNoFiles indicate location in that file
		// next indicates total occurences
		// last indicates actual wordIndex
		allocMatInt(invertedIndexMatrix, totNoWords, locNoFiles+2);
		int i, j;
		for(i = 0; i < totNoWords; i ++){
			for(j = 0; j < locNoFiles; j++){
				(*invertedIndexMatrix)[i][j] = -1;
			}
			(*invertedIndexMatrix)[i][locNoFiles] = 0;
		}

		// fill in tree from DocMatrix
		int index = 0;
		int wordIndex, tempIndex;
		for(doc = 0; doc < locNoFiles; doc ++){
			for(word = 0; word < locUniqueCount[doc]; word ++){
				wordIndex = localDocMatrix[doc][word][0];
				
				// if first occurrence
				if(vocabMarker[wordIndex] == doc + 1){
					(*invertedIndexMatrix)[index][doc] = word;
					(*invertedIndexMatrix)[index][locNoFiles] += localDocMatrix[doc][word][1];
					(*invertedIndexMatrix)[index][locNoFiles + 1] = wordIndex;
					wordLocationMat[wordIndex][2] = index;
					index ++;
				}
				else{
					tempIndex = wordLocationMat[wordIndex][2];	
					(*invertedIndexMatrix)[tempIndex][doc] = word;
					(*invertedIndexMatrix)[tempIndex][locNoFiles] += localDocMatrix[doc][word][1];
				}
			}
		}
		#ifdef DEBUG
		if(prank == 0){	
			int testIndex = localDocMatrix[2][435][0];
			tempIndex = wordLocationMat[testIndex][2];
			for(doc = 0; doc < locNoFiles; doc++){
				printf("%d\n", (*invertedIndexMatrix)[tempIndex][doc]);
			}
		}
		#endif
		free(vocabMarker);
		endTimeMPI(t1Sec, "INV INDEX MATRIX(Pd)", prank, pdComm);
		
	}

	return;
}



// does relevant set up to read in documents and then uses either serial or parallel readin method
void	readInDocsAll(char procType, int prank, int pdSize, int noFiles, int *locNoFiles, int parallelReadIn, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults,  char **pathNames, int *totalWords, struct timeval *t1Sec, struct timeval *t2Sec, int gutenFiles, MPI_Comm pdComm){
	if(procType == 'd'){	
		startTimeMPI(t1Sec, pdComm);

		int *localFiles;
		
		// DETERMINE NUMBER OF FILES PER RANK
		filesPerProc(&localFiles, prank, pdSize, noFiles, locNoFiles, pdComm);

		// EACH PROC READS IN A PORTION AND THEN COLLABORATE
		if(parallelReadIn){
			readInDocsParallel(prank, pdSize, stopfile, vocabSize, vocabulary, locUniqueCount, locTotalCount, localDocMatrix, printResults, pathNames, localFiles, *locNoFiles, totalWords, gutenFiles, pdComm);
		}
		// READ IN DOCS ON RANK 0 AND DISTRIBUTE INFO
		else{
			readInDocsSerial(prank, pdSize, pathNames, noFiles, stopfile, vocabSize, vocabulary, printResults, locUniqueCount, locTotalCount, *locNoFiles, localFiles, localDocMatrix, totalWords, gutenFiles, pdComm);
		}
		free(localFiles);
		startTimeMPI(t2Sec, pdComm);
	}
	
	return;
}


// initialises and distributes topics to relevant locations (BST METHOD)
void	initTopicsPLDAplus(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int *wordLocation, int **topicCountPerDoc, int pdSize, int size, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm){ 

	int 	sendBuff2[2];
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	// Pd initialises and sends
	if(procType == 'd'){
		int d, w, occ;
		int topic, word;

		// loops through all occurrences of all words of all local docs
		for(d = 0; d < locNoFiles; d ++){
			for(w = 0; w < locUniqueCount[d]; w ++){
				word = localDocMatrix[d][w][0];
				for(occ = 0; occ < localDocMatrix[d][w][1]; occ ++){

					// randomises topic and updates Cdoc and topicIndex
					topic = (int) (drand48() * noTopics);
					topicIndex[d][w][occ] = topic;		// z
					topicCountPerDoc[d][topic] ++;		// Cdoc

					// then sends the word index and the new topic to the relevant location (from wordLocation)
					sendBuff2[0] = word;
					sendBuff2[1] = topic;
					MPI_Send(&sendBuff2, 2, MPI_INT, wordLocation[word], 3, mainComm);

				}
			}
		}
		int i;
		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&sendBuff2, 2, MPI_INT, i, 0, mainComm);
		}
	}
	// Receive Cw
	else if(procType == 'w'){
		int noFinished = 0;
		struct nodeCword *current;
		MPI_Status stat;

		// loops through until all Pd's have sent finish signal
		while(noFinished < pdSize){

			// receives word and topic
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 3){

				// finds word in Cword tree and updates
				current = searchTreeCword(wordCountPerTopicTree, sendBuff2[0]);
				current -> topicArray[sendBuff2[1]] ++;

				// also updates the local Ctopic
				wordCountPerTopicTree -> topicCountTotal[sendBuff2[1]] ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// catches wrong tag being sent
			else{
				printf("\t\tERROR receive Cw recving wrong tag\n");
			}
		}

		#ifdef DEBUG
		int total = countTopicsCword(wordCountPerTopicTree);
		printf("total is %d\n", total);
		#endif
	}

	endTimeMPI(t1, "INIT TOPICS", rank, mainComm);

	return;
}



// initialises and distributes topics to relevant locations (MATRIX METHOD)
void	initTopicsPLDAplusMatrix(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int **wordCountPerTopic, int *localTopicCountTotal, int rank, MPI_Comm mainComm){ 

	int 	sendBuff2[2];
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	// Pd initialises and sends
	if(procType == 'd'){
		int d, w, occ;
		int topic, word;

		// loops through all occurrences of all words of all local docs
		for(d = 0; d < locNoFiles; d ++){
			for(w = 0; w < locUniqueCount[d]; w ++){
				word = localDocMatrix[d][w][0];
				for(occ = 0; occ < localDocMatrix[d][w][1]; occ ++){
		
					// randomises topic and updates Cdoc and topicIndex
					topic = (int) (drand48() * noTopics);
					topicIndex[d][w][occ] = topic;		// z
					topicCountPerDoc[d][topic] ++;		// Cdoc

					// send CONDENSED matrix location and topic
					// Recall the matrix on Pw has random words and the index for that matrix is stored in the 2nd column of wordLocationMat
					sendBuff2[0] = wordLocationMat[word][1];
					sendBuff2[1] = topic;
					MPI_Send(&sendBuff2, 2, MPI_INT, wordLocationMat[word][0], 3, mainComm);

				}
			}
		}
		int i;
		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&sendBuff2, 2, MPI_INT, i, 0, mainComm);
		}
	}
	// Receive Cw
	else if(procType == 'w'){
		int noFinished = 0;
		MPI_Status stat;

		// loops through until all Pd's have sent finish signal
		while(noFinished < pdSize){

			// receives word location in wordCountPerTopic and topic
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// updates in matrix and in local Ctopic
			if(stat.MPI_TAG == 3){
				wordCountPerTopic[sendBuff2[0]][sendBuff2[1]] ++;
				localTopicCountTotal[sendBuff2[1]] ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// catches any errors
			else{
				printf("\t\tERROR receive Cw recving wrong tag\n");
			}
		}

		#ifdef DEBUG
		int wordTot = 0;
		int topicTot = 0;
		int i, j;
		for(i = 0; i < 26865; i++){
			for(j = 0; j < noTopics; j++){
				wordTot += wordCountPerTopic[i][j];
			}
		}
		for(i = 0; i < noTopics; i++){
			topicTot += localTopicCountTotal[i];
		}
		printf("rank %d: wordTot %d, topicTot %d\n", rank, wordTot, topicTot);
		#endif

	}
	endTimeMPI(t1, "INIT TOPICS MATRIX", rank, mainComm);

	return;
}



// NOT USED ANYMORE FIXME
void	initTopicsPLDAplusMatrixRMA(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win **winArr, MPI_Comm mainComm){ 
	int 	sendBuff2[2];
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		int d, w, occ;
		int topic, word;
		int one = 1;
		int **Ctopictemp;
		callocMatInt(&Ctopictemp, pwSize, noTopics);
		for(d = 0; d < locNoFiles; d ++){
			for(w = 0; w < locUniqueCount[d]; w ++){
				word = localDocMatrix[d][w][0];
				for(occ = 0; occ < localDocMatrix[d][w][1]; occ ++){
					topic = (int) (drand48() * noTopics);
					topicIndex[d][w][occ] = topic;		// z
					topicCountPerDoc[d][topic] ++;		// Cdoc
					//localTopicCountTotal[topic] ++;		// Ctopic

					// send matrix location and topic
					sendBuff2[0] = wordLocationMat[word][1];
					sendBuff2[1] = topic;
					//MPI_Send(&sendBuff2, 2, MPI_INT, wordLocationMat[word][0], 3, mainComm);
					int win = wordLocationMat[word][0] - pdSize;
					int disp = wordLocationMat[word][1] * noTopics + topic;

					// update Cword
					MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize, 0, winArr[win][0]);
					MPI_Accumulate(&one, 1, MPI_INT, pdSize, disp, 1, MPI_INT, MPI_SUM, winArr[win][0]);
					MPI_Win_unlock(pdSize, winArr[win][0]);

					/*
					// update Ctopic 
					MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize, 0, winArr[win][1]);
					MPI_Accumulate(&one, 1, MPI_INT, pdSize, topic, 1, MPI_INT, MPI_SUM, winArr[win][1]);
					MPI_Win_unlock(pdSize, winArr[win][1]);
					*/
					Ctopictemp[win][topic]++;
				}
			}
		}
		int i;
		// update complete Ctopic
		for(i = 0; i < pwSize; i++){
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize, 0, winArr[i][1]);
				MPI_Accumulate(Ctopictemp[i], noTopics, MPI_INT, pdSize, 0, noTopics, MPI_INT, MPI_SUM, winArr[i][1]);
				MPI_Win_unlock(pdSize, winArr[i][1]);
		}
		freeMatInt(Ctopictemp);

		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&sendBuff2, 2, MPI_INT, i, 0, mainComm);
		}
	}
	// Receive Cw
	else if(procType == 'w'){
		int noFinished = 0;
		MPI_Status stat;
		while(noFinished < pdSize){
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, 0, mainComm, &stat);
			noFinished ++;
		}

		#ifdef DEBUG
		int wordTot = 0;
		int topicTot = 0;
		int i, j;
		for(i = 0; i < 26865; i++){
			for(j = 0; j < noTopics; j++){
				wordTot += wordCountPerTopic[i][j];
			}
		}
		for(i = 0; i < noTopics; i++){
			topicTot += localTopicCountTotal[i];
		}
		printf("rank %d: wordTot %d, topicTot %d\n", rank, wordTot, topicTot);
		#endif

	}
	endTimeMPI(t1, "INIT TOPICS RMA", rank, mainComm);

	return;
}



// FIXME NAME
// initialises and distributes topics to relevant locations (MATRIX METHOD with RMA)
void	initTopicsPLDAplusMatrixRMAFixed(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm){ 

	struct	timeval t1;
	startTimeMPI(&t1, mainComm);
	int sendBuff2[2];

	// Pd initialises and sends
	if(procType == 'd'){
		int d, w, occ;
		int topic, word;
		int dest, disp;
		int one = 1;
		int **Ctopictemp;
		callocMatInt(&Ctopictemp, pwSize, noTopics);		// FREED

		// loops through all occurrences of all words of all local docs
		for(d = 0; d < locNoFiles; d ++){
			for(w = 0; w < locUniqueCount[d]; w ++){
				word = localDocMatrix[d][w][0];
				for(occ = 0; occ < localDocMatrix[d][w][1]; occ ++){

					// randomises topic and updates Cdoc and topicIndex
					topic = (int) (drand48() * noTopics);
					topicIndex[d][w][occ] = topic;		// z
					topicCountPerDoc[d][topic] ++;		// Cdoc

					// finds destination and displacement into window		
					dest = wordLocationMat[word][0];
					disp = wordLocationMat[word][1] * noTopics + topic;

					// update Cword by locking window on dest and adding one by accumulating
					MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dest, 0, wordCountWin);
					MPI_Accumulate(&one, 1, MPI_INT, dest, disp, 1, MPI_INT, MPI_SUM, wordCountWin);
					MPI_Win_unlock(dest, wordCountWin);

					// updates the local Ctopic for the relevant Pw
					Ctopictemp[dest - pdSize][topic]++;
				}
			}
		}
		int i;
		// update complete Ctopic count to each of the Pws
		for(i = 0; i < pwSize; i++){
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize + i, 0, topicCountWin);
				MPI_Accumulate(Ctopictemp[i], noTopics, MPI_INT, pdSize + i, 0, noTopics, MPI_INT, MPI_SUM, topicCountWin);
				MPI_Win_unlock(pdSize + i, topicCountWin);
		}
		freeMatInt(Ctopictemp);

		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&sendBuff2, 2, MPI_INT, i, 0, mainComm);
		}
	}
	// Receive Cw Pw doesn't have to do anything with RMA, just wait for finish signal
	else if(procType == 'w'){
		int noFinished = 0;
		MPI_Status stat;
		while(noFinished < pdSize){
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			// Counts number of Pd's who send finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// catches any errors
			else{
				printf("\t\tERROR receive Cw recving wrong tag\n");
			}
		}

		#ifdef DEBUG
		int wordTot = 0;
		int topicTot = 0;
		int i, j;
		for(i = 0; i < 26865; i++){
			for(j = 0; j < noTopics; j++){
				wordTot += wordCountPerTopic[i][j];
			}
		}
		for(i = 0; i < noTopics; i++){
			topicTot += localTopicCountTotal[i];
		}
		printf("rank %d: wordTot %d, topicTot %d\n", rank, wordTot, topicTot);
		#endif

	}
	endTimeMPI(t1, "INIT TOPICS RMA", rank, mainComm);

	return;
}


// creates the word bundles on Pd(BST METHOD)
void	createWordBundles(char procType, struct treeRootInv *invertedIndexTree, int *wordLocation, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int *totNoBundles, int prank, MPI_Comm pdComm){


	// CREATE WORD BUNDLES
	if(procType == 'd'){
		int i, j;
		struct timeval t1;
		startTimeMPI(&t1, pdComm);

		// Sort by occurences
		int 	**occMat;
		int	totalWords = invertedIndexTree -> size;
		allocMatInt(&occMat, totalWords, 3);		// FREED
		fillOccMat(occMat, invertedIndexTree, wordLocation);
		

		// sort by index (for ability to do circular queue)
		quicksortMat(occMat, totalWords, 3, 0);

		// sort by location
		quicksortMat(occMat, totalWords, 3, 2);

		// find break points between location change
		int *breakers = malloc((pwSize+1) * sizeof(int));	// FREED
		breakers[0] = 0;
		breakers[pwSize] = totalWords;
		int loc = 0;
		for(i = 0; i < pwSize - 1; i++){
			while(occMat[loc][2] == i + pdSize){
				loc ++;
			}
			breakers[i+1] = loc;
		}

		// sort by occ within loc 
		for(i = 0; i < pwSize; i++){
			quicksortMat( &(occMat[breakers[i]]), breakers[i+1] - breakers[i], 3, 1);
		}

		
		// Find how many bundles are located on each Pw (and total bundles at last index)
		int *noBundles = malloc((pwSize + 1) * sizeof(int));		// FREED
		noBundles[pwSize] = 0;
		for(i = 0; i < pwSize; i++){
			noBundles[i] = (breakers[i+1] - breakers[i]) / bundleSize;
			if( (breakers[i+1] - breakers[i]) % bundleSize != 0){
				noBundles[i] ++;
			}
			noBundles[pwSize] += noBundles[i];
		}

		// Allocate bundle array (2d version for use and 3d version for sorting)
		// 			3d version has:
		// 				 location index 
		// 				 then bundle number at that index 
		// 				 and then the words in that bundle
		// 			2d version has:
		// 				bundles number
		// 				then the words in that bundle
		// 			All bundles have their location stored in the last index
		// 				eg if bundleSize = 4, a bundle stored at location 5 would look like:
		// 						812  23  293  901  5
		*bundles2d = malloc(noBundles[pwSize] * sizeof(int *));;
		int *temp = malloc( pwSize * noBundles[pwSize] * (bundleSize+1) * sizeof(int));
		int totalIndex = 0;
		int bunNo = 0;
		int ***bundles3d = malloc(pwSize * sizeof(int *));
		for(i = 0; i < pwSize; i++){
			bundles3d[i] = malloc(noBundles[i] * sizeof(int *));
			for(j = 0; j < noBundles[i]; j++){
				bundles3d[i][j] = &temp[totalIndex];
				(*bundles2d)[bunNo] = &temp[totalIndex];

				totalIndex += bundleSize + 1;
				bunNo ++;
			}
		}


		// fill bundle array
		int bun;
		int high, low;
	
		// fill for each set of bundles (each set referring to the Pw they are stored
		for(loc = 0; loc < pwSize; loc++){
			// low is lowest occurring word for current loc, high is highest occ
			low = breakers[loc];
			high = breakers[loc + 1] - 1;

			// fills for the number of bundles located on loc Pw
			for(bun = 0; bun < noBundles[loc] - 1; bun++){
				
				// fills highest word
				bundles3d[loc][bun][0] = occMat[high][0];
				high --;

				// and then remainder are from lowest occuring words
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = occMat[low][0];
					low ++;
				}
				
				// stores location of bundle at last index
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// LEFTOVER BUNDLE
			// if one word left over
			if(high == low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = -1;
				}
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// >1 but <bundleSize left
			else if(high > low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					if(low < high){
						bundles3d[loc][bun][i] = occMat[low][0];
						low ++;
					}
					else{
						bundles3d[loc][bun][i] = -1;
					}
				}
				bundles3d[loc][bun][bundleSize] = loc;
					
			}
		}

		// FREE 3D MEMORY AND START USING 2D
		for(i = 0; i < pwSize; i++){
			free(bundles3d[i]);
		}
		free(breakers);
		free(bundles3d);
		freeMatInt(occMat); 
		*totNoBundles = noBundles[pwSize];
		free(noBundles);

		
		// Sort by index of first word in bundle to randomise location of bundles
		quicksortMat(*bundles2d, *totNoBundles, bundleSize + 1, 0); 

		endTimeMPI(t1, "WORD BUNDLES (BST)", prank, pdComm);
	}

	return;
}



// creates the word bundles on Pd(MATRIX METHOD)
void	createWordBundlesMatrix(char procType, struct treeRootInv *invertedIndexTree, int **wordLocationMat, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int ***bundles2dConv, int *totNoBundles, int prank, MPI_Comm pdComm){


	// CREATE WORD BUNDLES
	if(procType == 'd'){
		struct timeval t1;
		startTimeMPI(&t1, pdComm);
		int i, j;

		// Sort by occurences
		int 	**occMat;
		int	totalWords = invertedIndexTree -> size;
		allocMatInt(&occMat, totalWords, 3);		// FREED
		fillOccMatMatrix(occMat, invertedIndexTree, wordLocationMat);
		

		// sort by index (for ability to do circular queue)
		quicksortMat(occMat, totalWords, 3, 0);

		// sort by location
		quicksortMat(occMat, totalWords, 3, 2);

		// find break points between location change
		int *breakers = malloc((pwSize+1) * sizeof(int));
		breakers[0] = 0;
		breakers[pwSize] = totalWords;
		int loc = 0;
		for(i = 0; i < pwSize - 1; i++){
			while(occMat[loc][2] == i + pdSize){
				loc ++;
			}
			breakers[i+1] = loc;
		}

		// sort by occ within loc
		for(i = 0; i < pwSize; i++){
			quicksortMat( &(occMat[breakers[i]]), breakers[i+1] - breakers[i], 3, 1);
		}

		
		// Find how many bundles are located on each Pw (and total bundles at last index)
		int *noBundles = malloc((pwSize + 1) * sizeof(int));		// FREED
		noBundles[pwSize] = 0;
		for(i = 0; i < pwSize; i++){
			noBundles[i] = (breakers[i+1] - breakers[i]) / bundleSize;
			if( (breakers[i+1] - breakers[i]) % bundleSize != 0){
				noBundles[i] ++;
			}
			noBundles[pwSize] += noBundles[i];
		}

		// Allocate bundle array (2d version for use and 3d version for sorting)
		// 			3d version has:
		// 				 location index 
		// 				 then bundle number at that index 
		// 				 and then the words in that bundle
		// 			2d version has:
		// 				bundles number
		// 				then the words in that bundle
		// 			All bundles have their location stored in the last index
		// 				eg if bundleSize = 4, a bundle stored at location 5 would look like:
		// 						812  23  293  901  5
		*bundles2d = malloc(noBundles[pwSize] * sizeof(int *));;
		int *temp = malloc( pwSize * noBundles[pwSize] * (bundleSize+1) * sizeof(int));

		int totalIndex = 0;
		int bunNo = 0;
		int ***bundles3d = malloc(pwSize * sizeof(int *));
		for(i = 0; i < pwSize; i++){
			bundles3d[i] = malloc(noBundles[i] * sizeof(int *));
			for(j = 0; j < noBundles[i]; j++){
				bundles3d[i][j] = &temp[totalIndex];
				(*bundles2d)[bunNo] = &temp[totalIndex];

				totalIndex += bundleSize + 1;
				bunNo ++;
			}
		}


		// fill bundle array
		int bun;
		int high, low;

		// fill for each set of bundles (each set referring to the Pw they are stored
		for(loc = 0; loc < pwSize; loc++){

			// low is lowest occ for current loc, high is highest occ
			low = breakers[loc];
			high = breakers[loc + 1] - 1;

			// fills for the number of bundles located on loc Pw
			for(bun = 0; bun < noBundles[loc] - 1; bun++){

				// fills highest word
				bundles3d[loc][bun][0] = occMat[high][0];
				high --;

				// and then remainder are from the lowest occuring words
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = occMat[low][0];
					low ++;
				}
		
				//  stores location of bundle at last index
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// LEFTOVER BUNDLE
			// if one word left over
			if(high == low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = -1;
				}
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// >1 but <bundleSize left
			else if(high > low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					if(low < high){
						bundles3d[loc][bun][i] = occMat[low][0];
						low ++;
					}
					else{
						bundles3d[loc][bun][i] = -1;
					}
				}
				bundles3d[loc][bun][bundleSize] = loc;
					
			}
		}

		// FREE 3D MEMORY AND START USING 2D
		for(i = 0; i < pwSize; i++){
			free(bundles3d[i]);
		}
		free(breakers);
		free(bundles3d);
		freeMatInt(occMat);
		*totNoBundles = noBundles[pwSize];

		
		// Sort by index to randomise location
		quicksortMat(*bundles2d, *totNoBundles, bundleSize + 1, 0); 

		// Fill in converted bundle indexes instead of actual word indexes
		allocMatInt(bundles2dConv, noBundles[loc], bundleSize);
		for(bun = 0; bun < noBundles[loc]; bun ++){
			for(i = 0; i < bundleSize; i++){
				if((*bundles2d)[bun][i] != -1){
					(*bundles2dConv)[bun][i] = wordLocationMat[(*bundles2d)[bun][i]][1];
				}
				else{
					(*bundles2dConv)[bun][i] = -1;
				}
			}
		}
		free(noBundles);
		endTimeMPI(t1, "WORD BUNDLES (MATRIX)", prank, pdComm);
	}

	return;
}



void	createWordBundlesMatrixInv(char procType, int **invertedIndexMatrix, int **wordLocationMat, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int ***bundles2dConv, int *totNoBundles, int prank, int invertedIndexSize, int locNoFiles, MPI_Comm pdComm){


	// CREATE WORD BUNDLES
	if(procType == 'd'){
		struct timeval t1;
		startTimeMPI(&t1, pdComm);
		int i, j;

		// Sort by occurences
		int 	**occMat;
		int	totalWords = invertedIndexSize;
		allocMatInt(&occMat, totalWords, 3);		// FREED

		// fill OccMatrix
		for(i = 0; i < totalWords; i++){
			occMat[i][0] = invertedIndexMatrix[i][locNoFiles+1];
			occMat[i][1] = invertedIndexMatrix[i][locNoFiles];
			occMat[i][2] = wordLocationMat[ occMat[i][0] ][0];
		}
		

		// sort by index (for ability to do circular queue)
		quicksortMat(occMat, totalWords, 3, 0);

		// sort by location
		quicksortMat(occMat, totalWords, 3, 2);

		// find break points between location change
		int *breakers = malloc((pwSize+1) * sizeof(int));
		breakers[0] = 0;
		breakers[pwSize] = totalWords;
		int loc = 0;
		for(i = 0; i < pwSize - 1; i++){
			while(occMat[loc][2] == i + pdSize){
				loc ++;
			}
			breakers[i+1] = loc;
		}

		// sort by occ within loc
		for(i = 0; i < pwSize; i++){
			quicksortMat( &(occMat[breakers[i]]), breakers[i+1] - breakers[i], 3, 1);
		}

		
		// Find how many bundles are located on each Pw (and total bundles at last index)
		int *noBundles = malloc((pwSize + 1) * sizeof(int));		// FREED
		noBundles[pwSize] = 0;
		for(i = 0; i < pwSize; i++){
			noBundles[i] = (breakers[i+1] - breakers[i]) / bundleSize;
			if( (breakers[i+1] - breakers[i]) % bundleSize != 0){
				noBundles[i] ++;
			}
			noBundles[pwSize] += noBundles[i];
		}

		// Allocate bundle array (2d version for use and 3d version for sorting)
		// 			3d version has:
		// 				 location index 
		// 				 then bundle number at that index 
		// 				 and then the words in that bundle
		// 			2d version has:
		// 				bundles number
		// 				then the words in that bundle
		// 			All bundles have their location stored in the last index
		// 				eg if bundleSize = 4, a bundle stored at location 5 would look like:
		// 						812  23  293  901  5
		*bundles2d = malloc(noBundles[pwSize] * sizeof(int *));;
		int *temp = malloc( pwSize * noBundles[pwSize] * (bundleSize+1) * sizeof(int));

		int totalIndex = 0;
		int bunNo = 0;
		int ***bundles3d = malloc(pwSize * sizeof(int *));
		for(i = 0; i < pwSize; i++){
			bundles3d[i] = malloc(noBundles[i] * sizeof(int *));
			for(j = 0; j < noBundles[i]; j++){
				bundles3d[i][j] = &temp[totalIndex];
				(*bundles2d)[bunNo] = &temp[totalIndex];

				totalIndex += bundleSize + 1;
				bunNo ++;
			}
		}


		// fill bundle array
		int bun;
		int high, low;

		// fill for each set of bundles (each set referring to the Pw they are stored
		for(loc = 0; loc < pwSize; loc++){

			// low is lowest occ for current loc, high is highest occ
			low = breakers[loc];
			high = breakers[loc + 1] - 1;

			// fills for the number of bundles located on loc Pw
			for(bun = 0; bun < noBundles[loc] - 1; bun++){

				// fills highest word
				bundles3d[loc][bun][0] = occMat[high][0];
				high --;

				// and then remainder are from the lowest occuring words
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = occMat[low][0];
					low ++;
				}
		
				//  stores location of bundle at last index
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// LEFTOVER BUNDLE
			// if one word left over
			if(high == low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = -1;
				}
				bundles3d[loc][bun][bundleSize] = loc;
			}
			// >1 but <bundleSize left
			else if(high > low){
				bundles3d[loc][bun][0] = occMat[high][0];
				for(i = 1; i < bundleSize; i++){
					if(low < high){
						bundles3d[loc][bun][i] = occMat[low][0];
						low ++;
					}
					else{
						bundles3d[loc][bun][i] = -1;
					}
				}
				bundles3d[loc][bun][bundleSize] = loc;
					
			}
		}

		// FREE 3D MEMORY AND START USING 2D
		for(i = 0; i < pwSize; i++){
			free(bundles3d[i]);
		}
		free(breakers);
		free(bundles3d);
		freeMatInt(occMat);
		*totNoBundles = noBundles[pwSize];

		
		// Sort by index to randomise location
		quicksortMat(*bundles2d, *totNoBundles, bundleSize + 1, 0); 

		// Fill in converted bundle indexes instead of actual word indexes
		allocMatInt(bundles2dConv, noBundles[loc], bundleSize);
		for(bun = 0; bun < noBundles[loc]; bun ++){
			for(i = 0; i < bundleSize; i++){
				if((*bundles2d)[bun][i] != -1){
					(*bundles2dConv)[bun][i] = wordLocationMat[(*bundles2d)[bun][i]][1];
				}
				else{
					(*bundles2dConv)[bun][i] = -1;
				}
			}
		}
		free(noBundles);
		endTimeMPI(t1, "WORD BUNDLES (INV MATRIX)", prank, pdComm);
	}

	return;
}



// assigns different starting point to each Pd on the circular queue
void	circleQueue(char procType, int vocabSize, int pdSize, int prank, int *circleStartBundle, int **bundles2d){

	if(procType == 'd'){
		int gap = vocabSize / pdSize;
		int wordIndex = prank * gap;

		// finds the first bundle with starting word near the break point	
		*circleStartBundle = 0;
		while(bundles2d[*circleStartBundle][0] < wordIndex){
			(*circleStartBundle) ++;
		} 
	}

	return;
}



// base model of PLDA+ Gibbs sampler(BST METHOD)
void    runPLDAplusGibbs(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm){
	
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	// Pd runs program
	if(procType == 'd'){
		int i, j, k;
		int iter;
		int bunNo, bun, dest;
		int *topicCountTotal, *localTopicCountTotal;
		MPI_Status stat;

		// SET UP C TOPIC
		callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
		callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed

		// Fetch Ctopic from all Pw processes
		zeroArrInt(topicCountTotal, noTopics);
		for(i = pdSize; i < size; i++){
		
			// send tag 4 to request and then receive that Pw's local Ctopic
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);

			// add each Pw's count into overall topicCountTotal
			for(j = 0; j < noTopics; j++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

		int d, w, occ;
		int word, topic;
		int matLoc;
		struct nodeInv *current;

		// prob dist used for sampling
		double *probDist = malloc(noTopics * sizeof(double));	// FREED
		double probTotal, sampler;

		// CwordBuffer used to receive current values from Pw
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
	
		// CtopicBuffer used to tally changes made to Ctopic (just changes not overall count)
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);					// FREED

		// Run Gibbs
		for(iter = 0; iter < iterations; iter++){
			// loops through all bundles
			for(bunNo = 0; bunNo < totNoBundles; bunNo ++){
				
				// adapts the bundle based on circle start point
				bun = bunNo + circleStartBundle;
				if(bun >= totNoBundles){bun -= totNoBundles;}
				dest = bundles2d[bun][bundleSize] + pdSize;
				
				// fetch Cword using tag 5
				MPI_Send(&iter, 1, MPI_INT, dest, 5, mainComm);
				MPI_Send(bundles2d[bun], bundleSize, MPI_INT, dest, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);
				zeroArrInt(CtopicBuffer, noTopics);
				
				// run Gibbs on the bundle
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];

					// checks it is not an empty word
					if(word != -1){			

						// finds in the inverted index and loops through all local docs
						current = searchTreeInv(invertedIndexTree, word);
						for(d = 0; d < locNoFiles; d++){
							matLoc = current -> invIndex[d];

							// runs Gibbs on all docs containing that word
							if(matLoc != -1){
								for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
									topic = topicIndex[d][matLoc][occ];
									
									// get -i values
									topicCountPerDoc[d][topic] --;
									CwordBuffer[w][topic] --;
									topicCountTotal[topic] --;
									CtopicBuffer[topic] --;
	
									// update probs
									probTotal = 0;
									for(k = 0; k < noTopics; k++){
										probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
										probTotal += probDist[k];
									}
			
									// sample new topic
									sampler = drand48() * probTotal;
									probTotal = probDist[0];
									topic = 0;
									while(probTotal < sampler && topic < noTopics-1){
										topic ++;
										probTotal += probDist[topic];
									}


									// update z
									topicIndex[d][matLoc][occ] = topic;

									// update counters
									topicCountPerDoc[d][topic] ++;
									CwordBuffer[w][topic] ++;
									topicCountTotal[topic] ++;
									CtopicBuffer[topic] ++;
								}	
							}
						}
					}
				}
				
		
				// update Cword with 6 tag
				MPI_Send(&iter, 1, MPI_INT, dest, 6, mainComm);
				MPI_Send(bundles2d[bun], bundleSize, MPI_INT, dest, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest, 6, mainComm);

			}
			// Fetch Ctopic after each full iteration with 4 tag
			zeroArrInt(topicCountTotal, noTopics);
			for(i = pdSize; i < size; i++){
				MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
				MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
				for(j = 0; j < noTopics; j++){
					topicCountTotal[j] += localTopicCountTotal[j];
				}
			}
		}

		// Send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}

		free(probDist);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;
		struct nodeCword *current;

		// wordIndex used to receive locations of words
		int *wordIndexes = malloc(bundleSize * sizeof(int));	// FREED

		// CwordBuffer used to store results before sending
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED

		// CtopicBuffer used to store Ctopic results coming in
		int *CtopicBuffer = malloc(noTopics * sizeof(int));		// FREED
	
		// loops through until finish signal processing requests
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(wordCountPerTopicTree -> topicCountTotal, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);

				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					current = searchTreeCword(wordCountPerTopicTree, wordIndexes[i]);
					if(current != NULL){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = current -> topicArray[j];
						}
					}
				}	

				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);

			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					current = searchTreeCword(wordCountPerTopicTree, wordIndexes[i]);
					if(current != NULL){
						for(j = 0; j < noTopics; j++){
							current -> topicArray[j] = CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					wordCountPerTopicTree -> topicCountTotal[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}
	
		free(wordIndexes);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
		
	}
	
	endTimeMPI(t1, "GIBBS SAMPLING", rank, mainComm);
	return;
}



// base model of PLDA+ Gibbs sampler (MATRIX METHOD)
void    runPLDAplusGibbsMatrix(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm mainComm){
	
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	// Pd runs program
	if(procType == 'd'){
		int i, j, k;
		int iter;
		int bunNo, bun, dest;
		int *topicCountTotal, *localTopicCountTotal;
		MPI_Status stat;

		// SET UP C TOPIC
		callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
		callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed

		// Fetch Ctopic from all Pw processes
		zeroArrInt(topicCountTotal, noTopics);
		for(i = pdSize; i < size; i++){
			
			// send tag 4 to request and then receive that Pw's local Ctopic
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);

			// add each Pw's count into overall topicCountTotal
			for(j = 0; j < noTopics; j++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

		int d, w, occ;
		int word, topic;
		int matLoc;
		struct nodeInv *current;

		// prob dist used for sampling
		double *probDist = malloc(noTopics * sizeof(double));		// FREED
		double probTotal, sampler;

		// CwordBuffer used to receive current values from Pw
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
	
		// CtopicBuffer used to tally changes made to Ctopic (just changes not overall count)
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);					// FREED

		// Run Gibbs
		for(iter = 0; iter < iterations; iter++){
			// loops through all bundles
			for(bunNo = 0; bunNo < totNoBundles; bunNo ++){

				// adapts the bundle based on circle start point
				bun = bunNo + circleStartBundle;
				if(bun >= totNoBundles){bun -= totNoBundles;}
				dest = bundles2d[bun][bundleSize] + pdSize;
				
				// fetch Cword using tag 5
				MPI_Send(&iter, 1, MPI_INT, dest, 5, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);
				zeroArrInt(CtopicBuffer, noTopics);
				
				// run Gibbs on the bundle
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];

					// checks it is not an empty word
					if(word != -1){			

						// finds in the inverted index and loops through all local docs
						current = searchTreeInv(invertedIndexTree, word);
						for(d = 0; d < locNoFiles; d++){
							matLoc = current -> invIndex[d];

							// runs Gibbs on all docs containing that word
							if(matLoc != -1){
								for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
									topic = topicIndex[d][matLoc][occ];
									
									// get -i values
									topicCountPerDoc[d][topic] --;
									CwordBuffer[w][topic] --;
									topicCountTotal[topic] --;
									CtopicBuffer[topic] --;
	
									// update probs
									probTotal = 0;
									for(k = 0; k < noTopics; k++){
										probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
										probTotal += probDist[k];
									}
			
									// sample new topic
									sampler = drand48() * probTotal;
									probTotal = probDist[0];
									topic = 0;
									while(probTotal < sampler && topic < noTopics-1){
										topic ++;
										probTotal += probDist[topic];
									}


									// update z
									topicIndex[d][matLoc][occ] = topic;

									// update counters
									topicCountPerDoc[d][topic] ++;
									CwordBuffer[w][topic] ++;
									topicCountTotal[topic] ++;
									CtopicBuffer[topic] ++;
								}	
							}
						}
					}
				}
				
		
				// update Cword with 6 tag
				MPI_Send(&iter, 1, MPI_INT, dest, 6, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest, 6, mainComm);

			}
			// Fetch Ctopic after a full iteration with 4 tag
			zeroArrInt(topicCountTotal, noTopics);
			for(i = pdSize; i < size; i++){
				MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
				MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
				for(j = 0; j < noTopics; j++){
					topicCountTotal[j] += localTopicCountTotal[j];
				}
			}
		}

		// Send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	
		free(probDist);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		// wordIndex used to receive locations of words
		int *wordIndexes = malloc(bundleSize * sizeof(int));	//FREED

		// CwordBuffer used to store results before sending
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED

		// CtopicBuffer used to store Ctopic results coming in
		int *CtopicBuffer = malloc(noTopics * sizeof(int));		// FREED
	
		// loops through until finish signal processing requests
		// uses wordCountPerTopic matrix here instead of BST tree
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							wordCountPerTopic[wordIndexes[i]][j] = CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}
		free(wordIndexes);
		free(CwordBuffer);	
		free(CtopicBuffer);

	}


	
	endTimeMPI(t1, "GIBBS SAMPLING MATRIX", rank, mainComm);
	return;
}



// base model of PLDA+ Gibbs sampler (MATRIX INV METHOD)
void    runPLDAplusGibbsMatrixInv(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, int **invertedIndexMatrix, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **wordLocationMat, MPI_Comm mainComm){
	
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	// Pd runs program
	if(procType == 'd'){
		int i, j, k;
		int iter;
		int bunNo, bun, dest;
		int *topicCountTotal, *localTopicCountTotal;
		MPI_Status stat;

		// SET UP C TOPIC
		callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
		callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed

		// Fetch Ctopic from all Pw processes
		zeroArrInt(topicCountTotal, noTopics);
		for(i = pdSize; i < size; i++){
			
			// send tag 4 to request and then receive that Pw's local Ctopic
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);

			// add each Pw's count into overall topicCountTotal
			for(j = 0; j < noTopics; j++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

		int d, w, occ;
		int word, topic;
		int matLoc;
		int invIndex;

		// prob dist used for sampling
		double *probDist = malloc(noTopics * sizeof(double));		// FREED
		double probTotal, sampler;

		// CwordBuffer used to receive current values from Pw
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
	
		// CtopicBuffer used to tally changes made to Ctopic (just changes not overall count)
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);					// FREED

/*
		int *invertedConversion = malloc(vocabSize * sizeof(int));
		for(w = 0; w < vocabSize; w++){
			invertedConversion[w] = wordLocationMat[w][2];
		}
*/
		// Run Gibbs
		for(iter = 0; iter < iterations; iter++){
			// loops through all bundles
			for(bunNo = 0; bunNo < totNoBundles; bunNo ++){

				// adapts the bundle based on circle start point
				bun = bunNo + circleStartBundle;
				if(bun >= totNoBundles){bun -= totNoBundles;}
				dest = bundles2d[bun][bundleSize] + pdSize;
				
				// fetch Cword using tag 5
				MPI_Send(&iter, 1, MPI_INT, dest, 5, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);
				zeroArrInt(CtopicBuffer, noTopics);
				
				// run Gibbs on the bundle
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];

					// checks it is not an empty word
					if(word != -1){			

						// finds in the inverted index and loops through all local docs
						invIndex = wordLocationMat[word][2];
						//invIndex = invertedConversion[word];
						for(d = 0; d < locNoFiles; d++){
							matLoc = invertedIndexMatrix[invIndex][d];

							// runs Gibbs on all docs containing that word
							if(matLoc != -1){
								for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
									topic = topicIndex[d][matLoc][occ];
									
									// get -i values
									topicCountPerDoc[d][topic] --;
									CwordBuffer[w][topic] --;
									topicCountTotal[topic] --;
									CtopicBuffer[topic] --;
	
									// update probs
									probTotal = 0;
									for(k = 0; k < noTopics; k++){
										probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
										probTotal += probDist[k];
									}
			
									// sample new topic
									sampler = drand48() * probTotal;
									probTotal = probDist[0];
									topic = 0;
									while(probTotal < sampler && topic < noTopics-1){
										topic ++;
										probTotal += probDist[topic];
									}


									// update z
									topicIndex[d][matLoc][occ] = topic;

									// update counters
									topicCountPerDoc[d][topic] ++;
									CwordBuffer[w][topic] ++;
									topicCountTotal[topic] ++;
									CtopicBuffer[topic] ++;
								}	
							}
						}
					}
				}
				
		
				// update Cword with 6 tag
				MPI_Send(&iter, 1, MPI_INT, dest, 6, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest, 6, mainComm);

			}
			// Fetch Ctopic after a full iteration with 4 tag
			zeroArrInt(topicCountTotal, noTopics);
			for(i = pdSize; i < size; i++){
				MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
				MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
				for(j = 0; j < noTopics; j++){
					topicCountTotal[j] += localTopicCountTotal[j];
				}
			}
		}

		// Send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	
		free(probDist);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		// wordIndex used to receive locations of words
		int *wordIndexes = malloc(bundleSize * sizeof(int));	//FREED

		// CwordBuffer used to store results before sending
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED

		// CtopicBuffer used to store Ctopic results coming in
		int *CtopicBuffer = malloc(noTopics * sizeof(int));		// FREED
	
		// loops through until finish signal processing requests
		// uses wordCountPerTopic matrix here instead of BST tree
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							wordCountPerTopic[wordIndexes[i]][j] = CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}
		free(wordIndexes);
		free(CwordBuffer);	
		free(CtopicBuffer);

	}


	
	endTimeMPI(t1, "GIBBS SAMPLING MATRIX INV", rank, mainComm);
	return;
}



// NOT USED ANYMORE FIXME
void    runPLDAplusGibbsMatrixRMA(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, MPI_Win **winArr, MPI_Comm mainComm){
	
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	// Pd runs program
	if(procType == 'd'){
		int i, j, k;
		int iter;
		int bunNo, bun, dest;
		int *topicCountTotal, *localTopicCountTotal;

		// SET UP C TOPIC
		callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
		callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed


		// Fetch Ctopic
		zeroArrInt(topicCountTotal, noTopics);
		for(i = 0; i < pwSize; i++){
			MPI_Win_lock(MPI_LOCK_SHARED, pdSize, 0, winArr[i][1]);
			MPI_Get(localTopicCountTotal, noTopics, MPI_INT, pdSize, 0, noTopics, MPI_INT, winArr[i][1]);
			MPI_Win_unlock(pdSize, winArr[i][1]);
			for(j = 0; j < noTopics; j++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

		int d, w, occ;
		int word, topic;
		int matLoc;
		struct nodeInv *current;

		// prob dist
		double *probDist = malloc(noTopics * sizeof(double));		// NOT freed
		double probTotal, sampler;

		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed
	
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);

		// Run Gibbs
		for(iter = 0; iter < iterations; iter++){
			for(bunNo = 0; bunNo < totNoBundles; bunNo ++){
				bun = bunNo + circleStartBundle;
				if(bun >= totNoBundles){bun -= totNoBundles;}

				dest = bundles2d[bun][bundleSize];
				
				// fetch Cword
				/*
				MPI_Send(&iter, 1, MPI_INT, dest+pdSize, 5, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest+pdSize, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest+pdSize, 5, mainComm, &stat);
				*/
				//*
				MPI_Win_lock(MPI_LOCK_SHARED, pdSize, 0, winArr[dest][0]);
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						MPI_Get(CwordBuffer[w], noTopics, MPI_INT, pdSize, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, winArr[dest][0]);
					}
				}
				MPI_Win_unlock(pdSize, winArr[dest][0]);
				//*/	
					
				
				zeroArrInt(CtopicBuffer, noTopics);
				// run Gibbs
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						current = searchTreeInv(invertedIndexTree, word);
						for(d = 0; d < locNoFiles; d++){
							matLoc = current -> invIndex[d];
							if(matLoc != -1){
								for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
									topic = topicIndex[d][matLoc][occ];
									
									// get -i values
									topicCountPerDoc[d][topic] --;
									CwordBuffer[w][topic] --;
									topicCountTotal[topic] --;
									CtopicBuffer[topic] --;
	
									// update probs
									probTotal = 0;
									for(k = 0; k < noTopics; k++){
										probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
										probTotal += probDist[k];
									}
			
									// sample new topic
									sampler = drand48() * probTotal;
									probTotal = probDist[0];


									// PUT BACK
									topic = 0;
									while(probTotal < sampler && topic < noTopics-1){
										topic ++;
										probTotal += probDist[topic];
									}
									if(probTotal < sampler && topic == noTopics - 1){
										FILE *errorlog = fopen("errorlog.txt", "a");
										fprintf(errorlog, "probTotal is %lf, sampler is %lf\n", probTotal, sampler);
										fclose(errorlog);
									}


									// update z
									topicIndex[d][matLoc][occ] = topic;

									// update counters
									topicCountPerDoc[d][topic] ++;
									CwordBuffer[w][topic] ++;
									topicCountTotal[topic] ++;
									CtopicBuffer[topic] ++;
								}	
							}
						}
					}
				}
				
		
				// update Cword
				/*
				MPI_Send(&iter, 1, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest+pdSize, 6, mainComm);
				*/

				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize, 0, winArr[dest][0]);
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						MPI_Put(CwordBuffer[w], noTopics, MPI_INT, pdSize, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, winArr[dest][0]);
					}
				}
				MPI_Win_unlock(pdSize, winArr[dest][0]);
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize, 0, winArr[dest][1]);
				MPI_Accumulate(CtopicBuffer, noTopics, MPI_INT, pdSize, 0, noTopics, MPI_INT, MPI_SUM, winArr[dest][1]);
				MPI_Win_unlock(pdSize, winArr[dest][1]);

			}
			// Fetch Ctopic
			zeroArrInt(topicCountTotal, noTopics);
			for(i = 0; i < pwSize; i++){
				MPI_Win_lock(MPI_LOCK_SHARED, pdSize, 0, winArr[i][1]);
				MPI_Get(localTopicCountTotal, noTopics, MPI_INT, pdSize, 0, noTopics, MPI_INT, winArr[i][1]);
				MPI_Win_unlock(pdSize, winArr[i][1]);
				for(j = 0; j < noTopics; j++){
					topicCountTotal[j] += localTopicCountTotal[j];
				}
			}
		}

		// Send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							wordCountPerTopic[wordIndexes[i]][j] = CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}

	}
	
	endTimeMPI(t1, "GIBBS SAMPLING RMA", rank, mainComm);
	return;
}


// COMMENTS DONE TO HERE FIXME
// FIXME NAME
void    runPLDAplusGibbsMatrixRMAFixed(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int circleStartBundle, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm){
	
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	// Pd runs program
	if(procType == 'd'){
		int i, j, k;
		int iter;
		int bunNo, bun, dest;
		int *topicCountTotal, *localTopicCountTotal;

		// SET UP C TOPIC
		callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
		callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed


		// Fetch Ctopic
		zeroArrInt(topicCountTotal, noTopics);
		for(i = 0; i < pwSize; i++){
			MPI_Win_lock(MPI_LOCK_SHARED, pdSize+i, 0, topicCountWin);
			MPI_Get(localTopicCountTotal, noTopics, MPI_INT, pdSize+i, 0, noTopics, MPI_INT, topicCountWin);
			MPI_Win_unlock(pdSize+i, topicCountWin);
			for(j = 0; j < noTopics; j++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

		int d, w, occ;
		int word, topic;
		int matLoc;
		struct nodeInv *current;

		// prob dist
		double *probDist = malloc(noTopics * sizeof(double));		// NOT freed
		double probTotal, sampler;

		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed
	
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);

		// Run Gibbs
		for(iter = 0; iter < iterations; iter++){
			for(bunNo = 0; bunNo < totNoBundles; bunNo ++){
				bun = bunNo + circleStartBundle;
				if(bun >= totNoBundles){bun -= totNoBundles;}

				dest = bundles2d[bun][bundleSize];
				
				// fetch Cword
				/*
				MPI_Send(&iter, 1, MPI_INT, dest+pdSize, 5, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest+pdSize, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest+pdSize, 5, mainComm, &stat);
				*/
				//*
				MPI_Win_lock(MPI_LOCK_SHARED, pdSize + dest, 0, wordCountWin);
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						MPI_Get(CwordBuffer[w], noTopics, MPI_INT, pdSize + dest, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, wordCountWin);
					}
				}
				MPI_Win_unlock(pdSize + dest, wordCountWin);
				//*/	
					
				
				zeroArrInt(CtopicBuffer, noTopics);
				// run Gibbs
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						current = searchTreeInv(invertedIndexTree, word);
						for(d = 0; d < locNoFiles; d++){
							matLoc = current -> invIndex[d];
							if(matLoc != -1){
								for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
									topic = topicIndex[d][matLoc][occ];
									
									// get -i values
									topicCountPerDoc[d][topic] --;
									CwordBuffer[w][topic] --;
									topicCountTotal[topic] --;
									CtopicBuffer[topic] --;
	
									// update probs
									probTotal = 0;
									for(k = 0; k < noTopics; k++){
										probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
										probTotal += probDist[k];
									}
			
									// sample new topic
									sampler = drand48() * probTotal;
									probTotal = probDist[0];


									// PUT BACK
									topic = 0;
									while(probTotal < sampler && topic < noTopics-1){
										topic ++;
										probTotal += probDist[topic];
									}
									if(probTotal < sampler && topic == noTopics - 1){
										FILE *errorlog = fopen("errorlog.txt", "a");
										fprintf(errorlog, "probTotal is %lf, sampler is %lf\n", probTotal, sampler);
										fclose(errorlog);
									}


									// update z
									topicIndex[d][matLoc][occ] = topic;

									// update counters
									topicCountPerDoc[d][topic] ++;
									CwordBuffer[w][topic] ++;
									topicCountTotal[topic] ++;
									CtopicBuffer[topic] ++;
								}	
							}
						}
					}
				}
				
		
				// update Cword
				/*
				MPI_Send(&iter, 1, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest+pdSize, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest+pdSize, 6, mainComm);
				*/

				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize + dest, 0, wordCountWin);
				for(w = 0; w < bundleSize; w++){
					word = bundles2d[bun][w];
					if(word != -1){			
						MPI_Put(CwordBuffer[w], noTopics, MPI_INT, pdSize + dest, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, wordCountWin);
					}
				}
				MPI_Win_unlock(pdSize + dest, wordCountWin);
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pdSize + dest, 0, topicCountWin);
				MPI_Accumulate(CtopicBuffer, noTopics, MPI_INT, pdSize + dest, 0, noTopics, MPI_INT, MPI_SUM, topicCountWin);
				MPI_Win_unlock(pdSize + dest, topicCountWin);

			}
			// Fetch Ctopic
			zeroArrInt(topicCountTotal, noTopics);
			for(i = 0; i < pwSize; i++){
				MPI_Win_lock(MPI_LOCK_SHARED, pdSize + i, 0, topicCountWin);
				MPI_Get(localTopicCountTotal, noTopics, MPI_INT, pdSize + i, 0, noTopics, MPI_INT, topicCountWin);
				MPI_Win_unlock(pdSize + i, topicCountWin);
				for(j = 0; j < noTopics; j++){
					topicCountTotal[j] += localTopicCountTotal[j];
				}
			}
		}

		// Send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							wordCountPerTopic[wordIndexes[i]][j] = CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}

	}
	
	endTimeMPI(t1, "GIBBS SAMPLING RMA", rank, mainComm);
	return;
}



void	collectAndPrint(int rank, int vocabSize, int noTopics, int pwSize, char procType, struct treeRootCword *wordCountPerTopicTree, int printResults, double *betaDir, char **vocabulary, MPI_Comm mainComm){

	int **wordCountPerTopic;
	struct timeval t1;
	startTimeMPI(&t1, mainComm);

	// SET UP OVERALL wordCountPerTopic MATRIX on Pd rank 0
	if(rank == 0){
		int i;
		callocMatInt(&wordCountPerTopic, vocabSize, noTopics);		// NOT FREED

		int noFinished = 0;	
		int index;
		int *topicArray = malloc(noTopics * sizeof(int));			// NOT FREED
		MPI_Status stat;
		while(noFinished < pwSize){
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, 7, mainComm, &stat);
			if(index != -1){
				MPI_Recv(topicArray, noTopics, MPI_INT, stat.MPI_SOURCE, 7, mainComm, &stat);
				for(i = 0; i < noTopics; i++){
					wordCountPerTopic[index][i] = topicArray[i];
				}
			}
			else{
				noFinished++;
			}
		}
		
	}

	if(procType == 'w'){
		sendForPrintCword(wordCountPerTopicTree, mainComm);	

	}



	// PRINT RESULTS
	if(printResults && rank == 0){
		double **phiMatrix;
		createPhi(&phiMatrix, vocabSize, noTopics, wordCountPerTopic, betaDir);		// freed
		printTopXLDA2(phiMatrix, noTopics, vocabSize, vocabulary, wordCountPerTopic, 10);
		freeVocab(&vocabulary, vocabSize);
		freeMat(phiMatrix);
	}

	endTimeMPI(t1, "PRINTING", rank, mainComm);
	return;
}



void	collectAndPrintMatrix(int rank, int vocabSize, int noTopics, char procType, int **wordCountPerTopic, int printResults, double *betaDir, char **vocabulary, int **wordLocationMat, int pdSize, int size, MPI_Comm mainComm){

	int **wordCountPerTopicFull;
	struct timeval t1;
	startTimeMPI(&t1, mainComm);

	// SET UP OVERALL wordCountPerTopic MATRIX on Pd rank 0
	if(rank == 0){
		int 	i, k;
		MPI_Status stat;
		callocMatInt(&wordCountPerTopicFull, vocabSize, noTopics);		// NOT FREED
		int	*topicArray = malloc(noTopics * sizeof(int));
		int 	dest, matIndex;
		// request Cword for every word from respective Pw
		for(i = 0; i < vocabSize; i++){
			dest = wordLocationMat[i][0];
			matIndex = wordLocationMat[i][1];
			if(matIndex == -1){printf("\t\tBINGO\n");}
			
			MPI_Send(&matIndex, 1, MPI_INT, dest, 7, mainComm);
			MPI_Recv(topicArray, noTopics, MPI_INT, dest, 7, mainComm, &stat);
			for(k = 0; k < noTopics; k++){
				wordCountPerTopicFull[i][k] = topicArray[k];
			}
		}
		int end = -1;
		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&end, 1, MPI_INT, i, 7, mainComm);
		}
		
	}

	if(procType == 'w'){
		int noFinished = 0;	
		int matIndex;
		MPI_Status stat;
		while(noFinished < 1){
			MPI_Recv(&matIndex, 1, MPI_INT, 0, 7, mainComm, &stat);
			if(matIndex != -1){
				MPI_Send(wordCountPerTopic[matIndex], noTopics, MPI_INT, 0, 7, mainComm);
			}
			else{
				noFinished ++;
			}
		}
	}



	// PRINT RESULTS
	if(printResults && rank == 0){
		double **phiMatrix;
		createPhi(&phiMatrix, vocabSize, noTopics, wordCountPerTopicFull, betaDir);		// freed
		printTopXLDA2(phiMatrix, noTopics, vocabSize, vocabulary, wordCountPerTopicFull, 10);
		freeVocab(&vocabulary, vocabSize);
		freeMat(phiMatrix);
	}

	endTimeMPI(t1, "PRINTING", rank, mainComm);
	return;
}

void	*fetchCtopic(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int		size = (*pargs) -> size;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;


	int 	*localTopicCountTotal, *holder;
	callocArrInt(&localTopicCountTotal, noTopics);
	callocArrInt(&holder, noTopics);

	int i, j;
	while((*stop) == 0){
		usleep(50000);
		zeroArrInt(holder, noTopics);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
			for(j = 0; j < noTopics; j ++){
				holder[j] += localTopicCountTotal[j];
			}
		}

		// update overall topicCountTotal (no need to lock as completely accurate counts aren't vital
		for(i = 0; i < noTopics; i++){
			topicCountTotal[i] = holder[i];
		}
	}

	free(localTopicCountTotal);
	free(holder);

	return NULL;
}


void	*fetchCword(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;

	int 	queueIndex = 0;
	int		bun, dest;
	while((*stop) == 0){
		// fetch
		if(queueMat[queueIndex][0] == 0){
			bun = queueMat[queueIndex][1];
			dest = bundles2d[bun][bundleSize] + pdSize;

			// request and read in the values
			MPI_Send(&bun, 1, MPI_INT, dest, 5, mainComm);
			MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
			MPI_Recv(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);

			// change operation to Gibbs
			queueMat[queueIndex][0] = 1;

			#ifdef DEBUG
			int i, j;
			for(i = queueIndex * bundleSize; i < (queueIndex + 1)*bundleSize; i++){
				for(j = 0; j < noTopics; j++){
					printf("%d  ", CwordBuffer[i][j]);
				}
				printf("\n");
			}
			#endif

		}
		
		// update queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}


	return NULL;
}


void	*runGibbs(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		locNoFiles = (*pargs) -> locNoFiles;
	int		***localDocMatrix = (*pargs) -> localDocMatrix;
	int		***topicIndex = (*pargs) -> topicIndex;
	int		**topicCountPerDoc = (*pargs) -> topicCountPerDoc;
	double	*alphaDir = (*pargs) -> alphaDir;
	double	*betaDir = (*pargs) -> betaDir;
	double	betaSum = (*pargs) -> betaSum;
	struct treeRootInv *invertedIndexTree = (*pargs) -> invertedIndexTree;

	int w, d, i, j, k;
	int word, bun, matLoc, occ, topic;
	struct nodeInv *current;

	double *probDist = malloc(noTopics * sizeof(double));
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);

	int 	queueIndex = 0;
	while((*stop) == 0){
		// run Gibbs
		if(queueMat[queueIndex][0] == 1){
			// zero the changes
			zeroMatInt(CwordChanges, bundleSize, noTopics);
			for(i = 0; i < noTopics; i++){
				CtopicBuffer[queueIndex][topic] = 0;
			}
	
			// find the bun and sample new topics
			bun = queueMat[queueIndex][1];
			for(w = 0; w < bundleSize; w++){
				word = bundles2d[bun][w];
				if(word != -1){
					current = searchTreeInv(invertedIndexTree, word);
					for(d = 0; d < locNoFiles; d++){
						matLoc = current -> invIndex[d];
						if(matLoc != -1){
							for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
								topic = topicIndex[d][matLoc][occ];

								// get -i values
								topicCountPerDoc[d][topic] --;
								CwordBuffer[queueIndex * bundleSize + w][topic] --;
								topicCountTotal[topic] --;
								CtopicBuffer[queueIndex][topic] --;
								CwordChanges[w][topic] --;

								// update probs
								probTotal = 0;
								for(k = 0; k < noTopics; k++){
									probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[queueIndex * bundleSize + w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
									probTotal += probDist[k];
								}

								// sample new topic
								sampler = drand48() * probTotal;
								probTotal = probDist[0];
								topic = 0;
								while(probTotal < sampler && topic < noTopics -1){
									topic ++;
									probTotal += probDist[topic];
								}

								// update z
								topicIndex[d][matLoc][occ] = topic;
				
								// update counters
								topicCountPerDoc[d][topic] ++;
								CwordBuffer[queueIndex * bundleSize + w][topic] ++;
								topicCountTotal[topic] ++;
								CtopicBuffer[queueIndex][topic] ++;
								CwordChanges[w][topic] ++;
							}
						}
					}
				}
			}
		
			// fill in CwordBuffer with CwordChanges (for sending)
			for(i = 0; i < bundleSize; i++){
				for(j = 0; j < noTopics; j++){
					CwordBuffer[queueIndex * bundleSize + i][j] = CwordChanges[i][j];
				}
			}

			// change operation to update Cword
			queueMat[queueIndex][0] = 2;
		}
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}
	
	return NULL;
}


void	*updateCword(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		iterations = (*pargs) -> iterations;
	int		totNoBundles = (*pargs) -> totNoBundles;
	MPI_Comm mainComm = (*pargs) -> mainComm;


	int 	stopCount = 0;
	int 	queueIndex = 0;
	int		bun, dest;
	int		currentIter = 0;
	int		currentBun = queueSize - 1;
	while((*stop) == 0){
		// update Cword
		if(queueMat[queueIndex][0] == 2){
			bun = queueMat[queueIndex][1];
			dest = bundles2d[bun][bundleSize] + pdSize;

			// request and read in the values
			MPI_Send(&bun, 1, MPI_INT, dest, 6, mainComm);
			MPI_Send(&(bundles2dConv[bun][0]), bundleSize, MPI_INT, dest, 6, mainComm);
			MPI_Send(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
			MPI_Send(&(CtopicBuffer[queueIndex][0]), noTopics, MPI_INT, dest, 6, mainComm);

			// change operation to Cword fetch for new bundle
			currentBun ++;
			if(currentBun < totNoBundles){
				queueMat[queueIndex][1] = currentBun;
				queueMat[queueIndex][0] = 0;
			}
			else{
				currentIter ++;
				if(currentIter < iterations){
					currentBun = 0;
					queueMat[queueIndex][1] = currentBun;
					queueMat[queueIndex][0] = 0;
					
				}
				else{
					stopCount++;
					queueMat[queueIndex][0] = -1;
					if(stopCount == queueSize){
						(*stop) = 1;
						return NULL;
					}
				}
			}
				

		}
			
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}

	
	return NULL;
}


void    runPLDAplusGibbsPthreads(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 2);
		int i;
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);
		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> topicCountTotal = topicCountTotal;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueMat = queueMat;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> CwordBuffer = CwordBuffer;
		pargs -> CtopicBuffer = CtopicBuffer;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		


		// Initial fetch Ctopic
		int 	*localTopicCountTotal;
		callocArrInt(&localTopicCountTotal, noTopics);
		MPI_Status stat;
		int j;
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
			for(j = 0; j < noTopics; j ++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}


	
		MPI_Barrier(pdComm);
		pthread_create(&(handles[0]), NULL, fetchCtopic, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, fetchCword, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, runGibbs, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, updateCword, (void *)&pargs);
		MPI_Barrier(pdComm);

		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		MPI_Barrier(pdComm);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	}
	else{
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							//FIXME NOTICE THE += CHANGE
							wordCountPerTopic[wordIndexes[i]][j] += CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS", rank, mainComm);
	return;
}


void	*recvFetchCtopic(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*localTopicCountTotalAct = (*pargs) -> localTopicCountTotalAct;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;



	int i = 0;
	int noFinished = 0;
	int counter = 0;
	while(noFinished < pdSize){
		counter ++;
		if(counter % 10 == 0){printf("\trecvFetchCtopic looping\n");}
		//printf("fetchtop recv\n");
		MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, 4, mainComm, &stat);
		if(i != -1){
			MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
		}
		else{
		//	printf("BREAKING FROM recvFetchCtopic\n");
			printf("fetchCtopic from %d\n", stat.MPI_SOURCE);
			noFinished ++;
		}
	}


	return NULL;
}



void	*recvFetchCword(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**wordCountPerTopic = (*pargs) -> wordCountPerTopic;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;

	int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
	int **CwordBuffer;
	allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed
	
	int i, j;
	int src;
	int noFinished = 0;

	int counter = 0;

	while(noFinished < pdSize){
		counter ++;
		if(counter % M == 0){printf("\trecvFetchCword looping\n");}
		MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, 5, mainComm, &stat);
		if(i != -1){
			src = stat.MPI_SOURCE;
			MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
			// fill in CwordBuffer
			for(i = 0; i < bundleSize; i++){
				if(wordIndexes[i] != -1){
					for(j = 0; j < noTopics; j++){
						CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
					}
				}
			}		
			// Send CwordBuffer
			MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
		}
		else{
			noFinished ++;
			printf("fetchCword from %d\n", stat.MPI_SOURCE);
			//printf("BREAKING FROM recvFetchCword\n");
		}

	}


	return NULL;
}



void	*recvUpdateCword(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**wordCountPerTopic = (*pargs) -> wordCountPerTopic;
	int 	*localTopicCountTotalAct = (*pargs) -> localTopicCountTotalAct;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;

	int i, j;
	int src;
	MPI_Status stat;


	int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
	int **CwordBuffer;
	allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed
	int *CtopicBuffer = malloc(noTopics * sizeof(int));
	int noFinished = 0;

	int counter = 0;
	while(noFinished < pdSize){
		counter ++;
		if(counter %M == 0){printf("\t recvUpdateCword looping\n");}
		MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, 6, mainComm, &stat);
		if(i != -1){
			//printf("AbundleSize is %d, noTopics uis %d\n", bundleSize, noTopics);
			src = stat.MPI_SOURCE;
			MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
			MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
			MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);

			// MAYBE PUT IN THE IRecv type things here??


			// Update topicArray
			for(i = 0; i < bundleSize; i++){
				if(wordIndexes[i] != -1){
					for(j = 0; j < noTopics; j++){
						//FIXME NOTICE THE += CHANGE
						wordCountPerTopic[wordIndexes[i]][j] += CwordBuffer[i][j];
					}
				}
			}	
			// update topicCountTotal
			for(i = 0; i < noTopics; i++){
				localTopicCountTotalAct[i] += CtopicBuffer[i];
			}


		}
		else{
			noFinished ++;
			printf("updateCword from %d\n", stat.MPI_SOURCE);
			//printf("recvUpdateCword BROKE\n");
		}
	}

	
	return NULL;
}



void	*stopProg(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	
	int noFinished = 0;
	int i;
	MPI_Status stat;
	
	while(noFinished < pdSize){
		MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, 0, mainComm, &stat);
		noFinished ++;
	}
	printf("Setting stop\n");
	(*stop) = 1;

	return NULL;
}




void    runPLDAplusGibbsPthreadsBoth(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 2);
		int i;
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);
		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> topicCountTotal = topicCountTotal;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueMat = queueMat;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> CwordBuffer = CwordBuffer;
		pargs -> CtopicBuffer = CtopicBuffer;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		


		// Initial fetch Ctopic
		int 	*localTopicCountTotal;
		callocArrInt(&localTopicCountTotal, noTopics);
		MPI_Status stat;
		int j;
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
			for(j = 0; j < noTopics; j ++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}


	
		MPI_Barrier(pdComm);
		pthread_create(&(handles[0]), NULL, fetchCtopic, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, fetchCword, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, runGibbs, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, updateCword, (void *)&pargs);
		MPI_Barrier(pdComm);

		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		MPI_Barrier(pdComm);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	}
	else{
		// set up threads
		pthread_t *handles;
		int numThreads = 3;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));

		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		pargs -> localTopicCountTotalAct = localTopicCountTotalAct;
		pargs -> wordCountPerTopic = wordCountPerTopic;



	//	MPI_Barrier(pwComm);
		pthread_create(&(handles[0]), NULL, recvFetchCtopic, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, recvFetchCword, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, recvUpdateCword, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, stopProg, (void *)&pargs);
	//	MPI_Barrier(pwComm);


		int i;
		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}



		/*
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							//FIXME NOTICE THE += CHANGE
							wordCountPerTopic[wordIndexes[i]][j] += CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}


		}
		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS BOTH", rank, mainComm);
	return;
}




void	*fetchCtopicDeadline(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int		size = (*pargs) -> size;
	int		queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		rank = (*pargs) -> rank;
	int		loopsPerSample = (*pargs) -> loopsPerSample;
	int		*topicSignal = (*pargs) -> topicSignal;
	pthread_mutex_t *locks = (*pargs) -> locks;	
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;


	int 	*localTopicCountTotal, *holder;
	callocArrInt(&localTopicCountTotal, noTopics);
	callocArrInt(&holder, noTopics);




	int i, j;
	int topicCounter = 0;
	while((*stop) == 0){
		// fetch Ctopic
		//usleep(sleepTime);
		//

		if((*topicSignal) == 1){
			zeroArrInt(holder, noTopics);
			for(i = pdSize; i < size; i++){
				MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
				MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
				for(j = 0; j < noTopics; j ++){
					holder[j] += localTopicCountTotal[j];
				}
			}

			// update overall topicCountTotal (no need to lock as completely accurate counts aren't vital
			for(i = 0; i < noTopics; i++){
				topicCountTotal[i] = holder[i];
			}
			topicCounter ++;
			(*topicSignal) = 0;
		}




		// XXX
		// RUN ADMIN
		for(i = 0; i < queueSize; i++){
			
			// if deadline up
			if(queueMat[i][2] == 0){
				pthread_mutex_lock(&(locks[i]));
				if(queueMat[i][2] == 0){
				//	printf("Deadline Killed\n");
					queueMat[i][2] = -1;
					pthread_mutex_unlock(&(locks[i]));
				
					// reset  it
					queueMat[i][0] = 3;
				}
				else{
					pthread_mutex_unlock(&(locks[i]));
				}
			}
			// else if still waiting
			else if(queueMat[i][2] > 0){
				pthread_mutex_lock(&(locks[i]));
				if(queueMat[i][2] > 0){
					queueMat[i][2] --;
					pthread_mutex_unlock(&(locks[i]));
				}
				else{
					pthread_mutex_unlock(&(locks[i]));
				}
			}
			else if(queueMat[i][2] < -1){
				printf("\t\t\tFATAL ERROR %d\n", queueMat[i][2]);
			}
		}


	}

	free(localTopicCountTotal);
	free(holder);
	printf("%d had %d topic samples\n", rank, topicCounter);

	return NULL;
}




void	*fetchCwordDeadline(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		*deadlines = (*pargs) -> deadlines;
	pthread_mutex_t *locks = (*pargs) -> locks;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;

	int 	queueIndex = 0;
	int		bun, dest;
	while((*stop) == 0){
		// fetch
		if(queueMat[queueIndex][0] == 0){
			if(queueMat[queueIndex][2] > 0){
				pthread_mutex_lock(&(locks[queueIndex]));
				if(queueMat[queueIndex][2] > 0 && queueMat[queueIndex][0] == 0){
					queueMat[queueIndex][2] = -1;
					pthread_mutex_unlock(&(locks[queueIndex]));

					bun = queueMat[queueIndex][1];
					dest = bundles2d[bun][bundleSize] + pdSize;

					// request and read in the values
					MPI_Send(&bun, 1, MPI_INT, dest, 5, mainComm);
					MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
					MPI_Recv(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);

					// change operation to Gibbs
					queueMat[queueIndex][2] = deadlines[1];
					queueMat[queueIndex][0] = 1;

					#ifdef DEBUG
					int i, j;
					for(i = queueIndex * bundleSize; i < (queueIndex + 1)*bundleSize; i++){
						for(j = 0; j < noTopics; j++){
							printf("%d  ", CwordBuffer[i][j]);
						}
						printf("\n");
					}
					#endif
				}
				else{
					pthread_mutex_unlock(&(locks[queueIndex]));
				}
			}

		}
		
		// update queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}


	return NULL;
}


void	*runGibbsDeadline(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		locNoFiles = (*pargs) -> locNoFiles;
	int		***localDocMatrix = (*pargs) -> localDocMatrix;
	int		***topicIndex = (*pargs) -> topicIndex;
	int		**topicCountPerDoc = (*pargs) -> topicCountPerDoc;
	int		*deadlines = (*pargs) -> deadlines;
	int		**invertedIndexMatrix = (*pargs) -> invertedIndexMatrix;
	int		**wordLocationMat = (*pargs) -> wordLocationMat;
	double	*alphaDir = (*pargs) -> alphaDir;
	double	*betaDir = (*pargs) -> betaDir;
	double	betaSum = (*pargs) -> betaSum;
	pthread_mutex_t *locks = (*pargs) -> locks;

	int w, d, i, j, k;
	int word, bun, matLoc, occ, topic;
	int	invIndex;

	double *probDist = malloc(noTopics * sizeof(double));
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);

	int 	queueIndex = 0;
	while((*stop) == 0){
		// run Gibbs
		if(queueMat[queueIndex][0] == 1){
			if(queueMat[queueIndex][2] > 0){
				pthread_mutex_lock(&(locks[queueIndex]));
				if(queueMat[queueIndex][2] > 0 && queueMat[queueIndex][0] == 1){
					queueMat[queueIndex][2] = -1;
					pthread_mutex_unlock(&(locks[queueIndex]));

					// zero the changes
					zeroMatInt(CwordChanges, bundleSize, noTopics);
					for(i = 0; i < noTopics; i++){
						CtopicBuffer[queueIndex][i] = 0;
					}
			
					// find the bun and sample new topics
					bun = queueMat[queueIndex][1];
					for(w = 0; w < bundleSize; w++){
						word = bundles2d[bun][w];
						if(word != -1){
							invIndex = wordLocationMat[word][2];
							for(d = 0; d < locNoFiles; d++){
								matLoc = invertedIndexMatrix[invIndex][d];
								if(matLoc != -1){
									for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
										topic = topicIndex[d][matLoc][occ];

										// get -i values
										topicCountPerDoc[d][topic] --;
										CwordBuffer[queueIndex * bundleSize + w][topic] --;
										topicCountTotal[topic] --;
										CtopicBuffer[queueIndex][topic] --;
										CwordChanges[w][topic] --;

										// update probs
										probTotal = 0;
										for(k = 0; k < noTopics; k++){
											probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[queueIndex * bundleSize + w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
											probTotal += probDist[k];
										}

										// sample new topic
										sampler = drand48() * probTotal;
										probTotal = probDist[0];
										topic = 0;
										while(probTotal < sampler && topic < noTopics -1){
											topic ++;
											probTotal += probDist[topic];
										}

										// update z
										topicIndex[d][matLoc][occ] = topic;
						
										// update counters
										topicCountPerDoc[d][topic] ++;
										CwordBuffer[queueIndex * bundleSize + w][topic] ++;
										topicCountTotal[topic] ++;
										CtopicBuffer[queueIndex][topic] ++;
										CwordChanges[w][topic] ++;
									}
								}
							}
						}
					}
				
					// fill in CwordBuffer with CwordChanges (for sending)
					for(i = 0; i < bundleSize; i++){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[queueIndex * bundleSize + i][j] = CwordChanges[i][j];
						}
					}

					// change operation to update Cword
					queueMat[queueIndex][2] = deadlines[2];
					queueMat[queueIndex][0] = 2;
				}
				else{
					pthread_mutex_unlock(&(locks[queueIndex]));
				}
			}
		}
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}
	
	return NULL;
}


void	*updateCwordDeadline(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		totNoBundles = (*pargs) -> totNoBundles;
	int		iterations = (*pargs) -> iterations;
	int		*deadlines = (*pargs) -> deadlines;
	int		*workload = (*pargs) -> workload;
	int		overload = (*pargs) -> overload;
	int		rank = (*pargs) -> rank;
	int		*topicSignal = (*pargs) -> topicSignal;
	pthread_mutex_t *locks = (*pargs) -> locks;
	MPI_Comm mainComm = (*pargs) -> mainComm;

	printf("overload is %d\n", overload);


	/*
		bun = queueMat[queueIndex][1];
		dest = bundles2d[bun][bundleSize] + pdSize;
	*/

	int 	queueIndex = 0;
	int		bun, dest;
	int		stopCount = 0;
	int		currentIter = 0;
	int		currentBun = queueSize - 1;
	int		skippedCount = 0;
	pthread_mutex_t stopLock;
	pthread_mutex_init(&stopLock, NULL);
	while((*stop) == 0){
	//	printf("WORK %d %d \n", workload[0], workload[1]);
		// update Cword
		
	//	if(currentIter >= iterations -1){
		//	printf("stopCount is %d, iterations is %d\n", stopCount, iterations);
	//	}
		if(queueMat[queueIndex][0] == 2){
			if(queueMat[queueIndex][2] > 0){
				pthread_mutex_lock(&(locks[queueIndex]));
				if(queueMat[queueIndex][2] > 0 && queueMat[queueIndex][0] == 2){
					queueMat[queueIndex][2] = -1;
					pthread_mutex_unlock(&(locks[queueIndex]));

					bun = queueMat[queueIndex][1];
					dest = bundles2d[bun][bundleSize] + pdSize;

					// request and read in the values
					MPI_Send(&bun, 1, MPI_INT, dest, 6, mainComm);
					MPI_Send(&(bundles2dConv[bun][0]), bundleSize, MPI_INT, dest, 6, mainComm);
					MPI_Send(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
					MPI_Send(&(CtopicBuffer[queueIndex][0]), noTopics, MPI_INT, dest, 6, mainComm);

					// change operation to Cword fetch for new bundle
					workload[ bundles2d[ queueMat[queueIndex][1] ][bundleSize] ] --;
					currentBun ++;
					while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
						currentBun++;
						skippedCount ++;
					}
					if(currentBun < totNoBundles){
						queueMat[queueIndex][1] = currentBun;
						queueMat[queueIndex][2] = deadlines[0];
						queueMat[queueIndex][0] = 0;
						workload[ bundles2d[currentBun][bundleSize] ] ++;
					}
					else{
						currentIter ++;
						if(currentIter < iterations){
							*topicSignal = 1;
			//				printf("%d iteration no %d\n", rank, currentIter);
							currentBun = 0;
							while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
								currentBun++;
								skippedCount ++;
							}
							queueMat[queueIndex][1] = currentBun;
							queueMat[queueIndex][2] = deadlines[0];
							queueMat[queueIndex][0] = 0;
							workload[ bundles2d[currentBun][bundleSize] ] ++;
							
						}
						else{
							pthread_mutex_lock(&stopLock);
							stopCount++;
							pthread_mutex_unlock(&stopLock);
							queueMat[queueIndex][0] = -1;
							if(stopCount == queueSize){
								(*stop) = 1;
								printf("%d skipped %d of %d\n", rank, skippedCount, totNoBundles*iterations);
								return NULL;
							}
						}
					}
				}
				else{
					pthread_mutex_unlock(&(locks[queueIndex]));
				}
			}

		}
		else if(queueMat[queueIndex][0] == 3){
			// change operation to Cword fetch for new bundle
			workload[ bundles2d[ queueMat[queueIndex][1] ][bundleSize] ] --;
			currentBun ++;
			skippedCount ++;	// extra count since bundle was killed
			while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
				currentBun++;
				skippedCount ++;
			//	printf("SKIPPING\n");
			}
			if(currentBun < totNoBundles){
				queueMat[queueIndex][1] = currentBun;
				queueMat[queueIndex][2] = deadlines[0];
				queueMat[queueIndex][0] = 0;
				workload[ bundles2d[currentBun][bundleSize] ] ++;
			}
			else{
				currentIter ++;
				if(currentIter < iterations){
					*topicSignal = 1;
			//		printf("%d iteration no %d\n", rank, currentIter);
					currentBun = 0;
					skippedCount ++;	// extra count since bundle was killed
					while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
						currentBun++;
						skippedCount ++;
					//	printf("SKIPPING\n");
					}
					queueMat[queueIndex][1] = currentBun;
					queueMat[queueIndex][2] = deadlines[0];
					queueMat[queueIndex][0] = 0;
					workload[ bundles2d[currentBun][bundleSize] ] ++;
					
				}
				else{
					pthread_mutex_lock(&stopLock);
					stopCount++;
					pthread_mutex_unlock(&stopLock);
					queueMat[queueIndex][0] = -1;
					if(stopCount == queueSize){
						(*stop) = 1;
						printf("%d skipped %d of %d\n", rank, skippedCount, totNoBundles*iterations);
						return NULL;
					}
				}
			}
		}
			
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}

	
	return NULL;
}




void    runPLDAplusGibbsPthreadsDeadline(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		int i;

		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int queueSize = 100;
		int loopsPerSample = 100;
		handles = malloc(numThreads * sizeof(pthread_t));

		// XXX

		// deadlines
		int		sleepTime = 5000;
		int		*deadlines = malloc(3 * sizeof(int));
		deadlines[0] = 15000*loopsPerSample;//25000 / sleepTime;
		deadlines[1] = 100*loopsPerSample;//15000 / sleepTime;
		deadlines[2] = 100*loopsPerSample;//10000 / sleepTime;
		

		int pwSize = size - pdSize; 		// FIXME
		int *workload;
		callocArrInt(&workload, pwSize);
		int overload = 1.5*(queueSize/pwSize);


		// locks
		pthread_mutex_t *locks = malloc(queueSize * sizeof(pthread_mutex_t));
		for(i = 0; i < queueSize; i++){
			pthread_mutex_init(&(locks[i]), NULL);
		}

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 3);
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
			queueMat[i][2] = deadlines[0];
			workload[ bundles2d[i][bundleSize] ] ++;	
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);
		int stop = 0;
		int topicSignal = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> topicCountTotal = topicCountTotal;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueMat = queueMat;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> CwordBuffer = CwordBuffer;
		pargs -> CtopicBuffer = CtopicBuffer;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		pargs -> deadlines = deadlines;
		pargs -> sleepTime = sleepTime;
		pargs -> locks = locks;
		pargs -> workload = workload;
		pargs -> overload = overload;
		pargs -> invertedIndexMatrix = invertedIndexMatrix;
		pargs -> wordLocationMat = wordLocationMat;
		pargs -> loopsPerSample = loopsPerSample;
		pargs -> topicSignal = &topicSignal;
		
	// XXX


		// Initial fetch Ctopic
		int 	*localTopicCountTotal;
		callocArrInt(&localTopicCountTotal, noTopics);
		MPI_Status stat;
		int j;
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
			for(j = 0; j < noTopics; j ++){
				topicCountTotal[j] += localTopicCountTotal[j];
			}
		}


	
		MPI_Barrier(pdComm);
		pthread_create(&(handles[0]), NULL, fetchCtopicDeadline, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, fetchCwordDeadline, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, runGibbsDeadline, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, updateCwordDeadline, (void *)&pargs);
		MPI_Barrier(pdComm);

		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		MPI_Barrier(pdComm);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	}
	else{
/*
		// set up threads
		pthread_t *handles;
		int numThreads = 3;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));

		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		pargs -> localTopicCountTotalAct = localTopicCountTotalAct;
		pargs -> wordCountPerTopic = wordCountPerTopic;



	//	MPI_Barrier(pwComm);
		pthread_create(&(handles[0]), NULL, recvFetchCtopic, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, recvFetchCword, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, recvUpdateCword, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, stopProg, (void *)&pargs);
	//	MPI_Barrier(pwComm);


		int i;
		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		*/
//		/*

		MPI_Status stat;
		int noFinished = 0;	// FIXME
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
				//printf("%d recvd finish from %d\n", rank, stat.MPI_SOURCE);
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							//FIXME NOTICE THE += CHANGE
							wordCountPerTopic[wordIndexes[i]][j] += CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}

		}
//		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS w/ deadline", rank, mainComm);
	return;
}



void	*fetchCtopicDeadlineNew(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int		size = (*pargs) -> size;
	int		rank = (*pargs) -> rank;
	int		loopsPerSample = (*pargs) -> loopsPerSample;
	int		*topicSignal = (*pargs) -> topicSignal;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;


	int 	*localTopicCountTotal, *holder;
	callocArrInt(&localTopicCountTotal, noTopics);
	callocArrInt(&holder, noTopics);




	int i, j;
	int topicCounter = 0;
	int counter = 0;
	while((*stop) == 0){
		counter ++;
		//if(counter % N == Z){printf("fetchCtopicDeadlineNew looping\n");}
		// fetch Ctopic
		//usleep(sleepTime);
		//

		if((*topicSignal) == 1){
			zeroArrInt(holder, noTopics);
			for(i = pdSize; i < size; i++){
				MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
				MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
				for(j = 0; j < noTopics; j ++){
					holder[j] += localTopicCountTotal[j];
				}
			}

			// update overall topicCountTotal (no need to lock as completely accurate counts aren't vital
			for(i = 0; i < noTopics; i++){
				topicCountTotal[i] = holder[i];
			}
			topicCounter ++;
			(*topicSignal) = 0;
		}




		// XXX

	}

	free(localTopicCountTotal);
	free(holder);
	printf("%d had %d topic samples\n", rank, topicCounter);

	return NULL;
}




void	*fetchCwordDeadlineNew(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;

	int 	queueIndex = 0;
	int		bun, dest;
	int 	counter = 0;
	while((*stop) == 0){
		counter ++;
	//	if(counter % N == Z){printf("fetchCwordDeadlineNew looping\n");}
		// fetch
		if(queueMat[queueIndex][0] == 0){

					bun = queueMat[queueIndex][1];
					dest = bundles2d[bun][bundleSize] + pdSize;
					//printf("bun is %d, dest is %d\n", bun, dest);

					// request and read in the values
				//	printf("bundleSize is %d, noTopics is %d\n", bundleSize, noTopics);
					MPI_Send(&bun, 1, MPI_INT, dest, 5, mainComm);
					MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
					MPI_Recv(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);

					// change operation to Gibbs
					queueMat[queueIndex][0] = 1;

					#ifdef DEBUG
					int i, j;
					for(i = queueIndex * bundleSize; i < (queueIndex + 1)*bundleSize; i++){
						for(j = 0; j < noTopics; j++){
							printf("%d  ", CwordBuffer[i][j]);
						}
						printf("\n");
					}
					#endif

		}
		
		// update queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}


	return NULL;
}


void	*runGibbsDeadlineNew(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		locNoFiles = (*pargs) -> locNoFiles;
	int		***localDocMatrix = (*pargs) -> localDocMatrix;
	int		***topicIndex = (*pargs) -> topicIndex;
	int		**topicCountPerDoc = (*pargs) -> topicCountPerDoc;
	int		**invertedIndexMatrix = (*pargs) -> invertedIndexMatrix;
	int		**wordLocationMat = (*pargs) -> wordLocationMat;
	double	*alphaDir = (*pargs) -> alphaDir;
	double	*betaDir = (*pargs) -> betaDir;
	double	betaSum = (*pargs) -> betaSum;

	int w, d, i, j, k;
	int word, bun, matLoc, occ, topic;
	int	invIndex;

	double *probDist = malloc(noTopics * sizeof(double));
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);

	int counter = 0;
	int 	queueIndex = 0;
	while((*stop) == 0){
		counter ++;
		//if(counter % N == Z){printf("runGibbsDeadlineNew looping\n");}
		// run Gibbs
		if(queueMat[queueIndex][0] == 1){

					// zero the changes
					zeroMatInt(CwordChanges, bundleSize, noTopics);
					for(i = 0; i < noTopics; i++){
						CtopicBuffer[queueIndex][i] = 0;
					}
			
					// find the bun and sample new topics
					bun = queueMat[queueIndex][1];
					for(w = 0; w < bundleSize; w++){
						word = bundles2d[bun][w];
						if(word != -1){
							invIndex = wordLocationMat[word][2];
							for(d = 0; d < locNoFiles; d++){
								matLoc = invertedIndexMatrix[invIndex][d];
								if(matLoc != -1){
									for(occ = 0; occ < localDocMatrix[d][matLoc][1]; occ++){
										topic = topicIndex[d][matLoc][occ];

										// get -i values
										topicCountPerDoc[d][topic] --;
										CwordBuffer[queueIndex * bundleSize + w][topic] --;
										topicCountTotal[topic] --;
										CtopicBuffer[queueIndex][topic] --;
										CwordChanges[w][topic] --;

										// update probs
										probTotal = 0;
										for(k = 0; k < noTopics; k++){
											probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (CwordBuffer[queueIndex * bundleSize + w][k] + betaDir[word]) / (topicCountTotal[k] + betaSum);
											probTotal += probDist[k];
										}

										// sample new topic
										sampler = drand48() * probTotal;
										probTotal = probDist[0];
										topic = 0;
										while(probTotal < sampler && topic < noTopics -1){
											topic ++;
											probTotal += probDist[topic];
										}

										// update z
										topicIndex[d][matLoc][occ] = topic;
						
										// update counters
										topicCountPerDoc[d][topic] ++;
										CwordBuffer[queueIndex * bundleSize + w][topic] ++;
										topicCountTotal[topic] ++;
										CtopicBuffer[queueIndex][topic] ++;
										CwordChanges[w][topic] ++;
									}
								}
							}
						}
					}
				
					// fill in CwordBuffer with CwordChanges (for sending)
					for(i = 0; i < bundleSize; i++){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[queueIndex * bundleSize + i][j] = CwordChanges[i][j];
						}
					}

					// change operation to update Cword
					queueMat[queueIndex][0] = 2;
		}
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}
	
	return NULL;
}


void	*updateCwordDeadlineNew(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int 	queueSize = (*pargs) -> queueSize;
	int		**queueMat = (*pargs) -> queueMat;
	int		**bundles2d = (*pargs) -> bundles2d;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**bundles2dConv = (*pargs) -> bundles2dConv;
	int		**CwordBuffer = (*pargs) -> CwordBuffer;
	int		**CtopicBuffer = (*pargs) -> CtopicBuffer;
	int		totNoBundles = (*pargs) -> totNoBundles;
	int		iterations = (*pargs) -> iterations;
	int		*workload = (*pargs) -> workload;
	int		overload = (*pargs) -> overload;
	int		rank = (*pargs) -> rank;
	int		*topicSignal = (*pargs) -> topicSignal;
	MPI_Comm mainComm = (*pargs) -> mainComm;

	printf("overload is %d\n", overload);


	/*
		bun = queueMat[queueIndex][1];
		dest = bundles2d[bun][bundleSize] + pdSize;
	*/

	int 	queueIndex = 0;
	int		bun, dest;
	int		stopCount = 0;
	int		currentIter = 0;
	int		currentBun = queueSize - 1;
	int		skippedCount = 0;
	pthread_mutex_t stopLock;
	pthread_mutex_init(&stopLock, NULL);
	int counter = 0;
	while((*stop) == 0){
		counter ++;
		//if(counter % N == Z){printf("updateCwordDeadlineNew looping\n");}
	//	printf("WORK %d %d \n", workload[0], workload[1]);
		// update Cword
		
	//	if(currentIter >= iterations -1){
		//	printf("stopCount is %d, iterations is %d\n", stopCount, iterations);
	//	}
		if(queueMat[queueIndex][0] == 2){
					bun = queueMat[queueIndex][1];
					dest = bundles2d[bun][bundleSize] + pdSize;

					// request and read in the values
					//printf("sending %d\n", bun);
					if(bun == -1){printf("\tBINGO\n");}
					MPI_Send(&bun, 1, MPI_INT, dest, 6, mainComm);
					//printf("sent %d\n", bun);
					MPI_Send(&(bundles2dConv[bun][0]), bundleSize, MPI_INT, dest, 6, mainComm);
					MPI_Send(&(CwordBuffer[bundleSize*queueIndex][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
					MPI_Send(&(CtopicBuffer[queueIndex][0]), noTopics, MPI_INT, dest, 6, mainComm);

					// change operation to Cword fetch for new bundle
					workload[ bundles2d[ queueMat[queueIndex][1] ][bundleSize] ] --;
					currentBun ++;
					while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
						currentBun++;
						skippedCount ++;
					}
					if(currentBun < totNoBundles){
						queueMat[queueIndex][1] = currentBun;
						queueMat[queueIndex][0] = 0;
						workload[ bundles2d[currentBun][bundleSize] ] ++;
					}
					else{
						currentIter ++;
						if(currentIter < iterations){
							*topicSignal = 1;
							printf("%d iteration no %d\n", rank, currentIter);
							currentBun = 0;
							while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
								currentBun++;
								skippedCount ++;
							}
							queueMat[queueIndex][1] = currentBun;
							queueMat[queueIndex][0] = 0;
							workload[ bundles2d[currentBun][bundleSize] ] ++;
							
						}
						else{
							pthread_mutex_lock(&stopLock);
							stopCount++;
							pthread_mutex_unlock(&stopLock);
							queueMat[queueIndex][0] = -1;
							if(stopCount == queueSize){
								(*stop) = 1;
								printf("%d skipped %d of %d\n", rank, skippedCount, totNoBundles*iterations);
								return NULL;
							}
						}
					}
		}
		else if(queueMat[queueIndex][0] == 3){
			// change operation to Cword fetch for new bundle
			workload[ bundles2d[ queueMat[queueIndex][1] ][bundleSize] ] --;
			currentBun ++;
			skippedCount ++;	// extra count since bundle was killed
			while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
				currentBun++;
				skippedCount ++;
			//	printf("SKIPPING\n");
			}
			if(currentBun < totNoBundles){
				queueMat[queueIndex][1] = currentBun;
				queueMat[queueIndex][0] = 0;
				workload[ bundles2d[currentBun][bundleSize] ] ++;
			}
			else{
				currentIter ++;
				if(currentIter < iterations){
					*topicSignal = 1;
					printf("%d iteration no %d\n", rank, currentIter);
					currentBun = 0;
					skippedCount ++;	// extra count since bundle was killed
					while(currentBun < totNoBundles - 1 && workload[ bundles2d[currentBun][bundleSize]] + 1 > overload){
						currentBun++;
						skippedCount ++;
					//	printf("SKIPPING\n");
					}
					queueMat[queueIndex][1] = currentBun;
					queueMat[queueIndex][0] = 0;
					workload[ bundles2d[currentBun][bundleSize] ] ++;
					
				}
				else{
					pthread_mutex_lock(&stopLock);
					stopCount++;
					pthread_mutex_unlock(&stopLock);
					queueMat[queueIndex][0] = -1;
					if(stopCount == queueSize){
						(*stop) = 1;
						printf("%d skipped %d of %d\n", rank, skippedCount, totNoBundles*iterations);
						return NULL;
					}
				}
			}
		}
			
		
		// next queueIndex
		queueIndex ++;
		if(queueIndex == queueSize){
			queueIndex = 0;
		}
	}

	
	return NULL;
}




void    runPLDAplusGibbsPthreadsDeadlineNew(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);
//	printf("RUNNING NEW\n");


	int queueSize = 200;

	if(procType == 'd'){
		int i;

		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int loopsPerSample = 100;
		handles = malloc(numThreads * sizeof(pthread_t));

		// XXX

		// deadlines
		int		sleepTime = 5000;
		int		*deadlines = malloc(3 * sizeof(int));
		deadlines[0] = 150000*loopsPerSample;//25000 / sleepTime;
		deadlines[1] = 10000*loopsPerSample;//15000 / sleepTime;
		deadlines[2] = 10000*loopsPerSample;//10000 / sleepTime;
		

		int pwSize = size - pdSize; 		// FIXME
		int *workload;
		callocArrInt(&workload, pwSize);
		int overload = 1.3*(queueSize/pwSize);


		// locks
		pthread_mutex_t *locks = malloc(queueSize * sizeof(pthread_mutex_t));
		for(i = 0; i < queueSize; i++){
			pthread_mutex_init(&(locks[i]), NULL);
		}

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 3);
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
			queueMat[i][2] = deadlines[0];
			workload[ bundles2d[i][bundleSize] ] ++;	
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);
		int stop = 0;
		int topicSignal = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> topicCountTotal = topicCountTotal;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueMat = queueMat;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> CwordBuffer = CwordBuffer;
		pargs -> CtopicBuffer = CtopicBuffer;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		pargs -> sleepTime = sleepTime;
		pargs -> locks = locks;
		pargs -> workload = workload;
		pargs -> overload = overload;
		pargs -> invertedIndexMatrix = invertedIndexMatrix;
		pargs -> wordLocationMat = wordLocationMat;
		pargs -> loopsPerSample = loopsPerSample;
		pargs -> topicSignal = &topicSignal;
		
	// XXX


		// Initial fetch Ctopic
		int 	*localTopicCountTotal;
		callocArrInt(&localTopicCountTotal, noTopics);
		MPI_Status stat;
		int j;
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);		
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
			for(j = 0; j < noTopics; j ++){
					topicCountTotal[j] += localTopicCountTotal[j];
			}
		}

	
		MPI_Barrier(pdComm);
		pthread_create(&(handles[0]), NULL, fetchCtopicDeadlineNew, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, fetchCwordDeadlineNew, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, runGibbsDeadlineNew, (void *)&pargs);
		pthread_create(&(handles[3]), NULL, updateCwordDeadlineNew, (void *)&pargs);
		MPI_Barrier(pdComm);

		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		MPI_Barrier(pdComm);

	//	/*
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	//	*/
/*
		// Send Finish signal
			for(i = pdSize; i < size; i++){
				//MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
				j = -1;
				MPI_Send(&j, 1, MPI_INT, i, 4, mainComm);
				MPI_Send(&j, 1, MPI_INT, i, 5, mainComm);
				MPI_Send(&j, 1, MPI_INT, i, 6, mainComm);
			
			}
*/
	}
	else{
/*
		// set up threads
		pthread_t *handles;
		int numThreads = 3;
		handles = malloc(numThreads * sizeof(pthread_t));

		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));
		pargs -> stop = &stop;
		pargs -> noTopics = noTopics;
		pargs -> pdSize = pdSize;
		pargs -> size = size;
		pargs -> mainComm = mainComm;
		pargs -> queueSize = queueSize;
		pargs -> bundles2d = bundles2d;
		pargs -> bundleSize = bundleSize;
		pargs -> bundles2dConv = bundles2dConv;
		pargs -> invertedIndexTree = invertedIndexTree;
		pargs -> locNoFiles = locNoFiles;
		pargs -> localDocMatrix = localDocMatrix;
		pargs -> topicIndex = topicIndex;
		pargs -> topicCountPerDoc = topicCountPerDoc;
		pargs -> alphaDir = alphaDir;
		pargs -> betaDir = betaDir;
		pargs -> betaSum = betaSum;
		pargs -> iterations = iterations;
		pargs -> totNoBundles = totNoBundles;
		pargs -> rank = rank;
		pargs -> localTopicCountTotalAct = localTopicCountTotalAct;
		pargs -> wordCountPerTopic = wordCountPerTopic;



	//	MPI_Barrier(pwComm);
		pthread_create(&(handles[0]), NULL, recvFetchCtopic, (void *)&pargs);
		pthread_create(&(handles[1]), NULL, recvFetchCword, (void *)&pargs);
		pthread_create(&(handles[2]), NULL, recvUpdateCword, (void *)&pargs);
	//	pthread_create(&(handles[3]), NULL, stopProg, (void *)&pargs);
	//	MPI_Barrier(pwComm);


		int i;
		for(i = 0; i < numThreads; i++){
			pthread_join(handles[i], NULL);
		}
		*/
//		/*

		MPI_Status stat;
		int noFinished = 0;	// FIXME
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// NOT FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// NOT freed

		int *CtopicBuffer = malloc(noTopics * sizeof(int));
	
		while(noFinished < pdSize){
			MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);

			// Finish signal
			if(stat.MPI_TAG == 0){
				noFinished ++;
				//printf("%d recvd finish from %d\n", rank, stat.MPI_SOURCE);
			}
			// Ctopic being requested
			else if(stat.MPI_TAG == 4){
				MPI_Send(localTopicCountTotalAct, noTopics, MPI_INT, stat.MPI_SOURCE, 4, mainComm);
			}
			// Cword being requested
			else if(stat.MPI_TAG == 5){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 5, mainComm, &stat);
				// fill in CwordBuffer
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							CwordBuffer[i][j] = wordCountPerTopic[wordIndexes[i]][j];
						}
					}
				}	
				// Send CwordBuffer
				MPI_Send(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 5, mainComm);
			}
			// Cword being updated
			else if(stat.MPI_TAG == 6){
				src = stat.MPI_SOURCE;
				MPI_Recv(wordIndexes, bundleSize, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(&(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, src, 6, mainComm, &stat);
				MPI_Recv(CtopicBuffer, noTopics, MPI_INT, src, 6, mainComm, &stat);
		
				// MAYBE PUT IN THE IRecv type things here??


				// Update topicArray
				for(i = 0; i < bundleSize; i++){
					if(wordIndexes[i] != -1){
						for(j = 0; j < noTopics; j++){
							//FIXME NOTICE THE += CHANGE
							wordCountPerTopic[wordIndexes[i]][j] += CwordBuffer[i][j];
						}
					}
				}	
				// update topicCountTotal
				for(i = 0; i < noTopics; i++){
					localTopicCountTotalAct[i] += CtopicBuffer[i];
				}
			}
			else{
				printf("\t\tERROR MAIN Recving wrong tag\n");
			}

		}
//		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS w/ deadline NEW", rank, mainComm);
	return;
}






int 	main(int argc, char *argv[]){


	// ************************************************
	// ************************************************
	// ***********      VARIABLES          ************
	// ************************************************
	// ************************************************



	
	// CREATE VARIABLES AND INITIAL VALUES
	int 	seed = time(NULL);
	int 	noTopics = 10;
	char 	*directory = malloc(1024);		// freed
	sprintf(directory, "../../Inputs/Gutenberg");
	char	*stopfile = malloc(1024);		// freed
	sprintf(stopfile, "../../Inputs/stopwords.txt");
	int		printResults = 0;
	int 	maxFiles = 0;
	int		iterations = 100;
	int 	parallelReadIn = 0;
	int 	LDAmethod = 3;
	double	pwToPdRatio = 0.6;
	int		bundleSize = 4;
	int		gutenFiles = 0;

	struct 	timeval t1, t1Sec, t2Sec;
	int 	noFiles, vocabSize, totalWords;
	int 	locNoFiles;
	char	**vocabulary;
	int	**topicCountPerDoc, ***topicIndex;
	double 	*alphaDir, *betaDir, betaSum;
	int 	*locUniqueCount, *locTotalCount;
	int 	***localDocMatrix;
	char 	**pathNames;
	int 	rank, size;


	// PLDA+
	MPI_Comm pdComm, pwComm;
	int		pdSize, pwSize, prank;
	char 	procType;
	int		*wordSpread, *wordLocation;
	int		**wordLocationMat;
	int		**wordCountPerTopic;
	struct 	treeRootCword *wordCountPerTopicTree;
	struct	treeRootInv *invertedIndexTree;
	int		**bundles2d, **bundles2dConv, totNoBundles, circleStartBundle;
	int		CwordMethod = 2;	// 1 for bst, 2 for matrix
	int		invMethod = 2;		// 1 for bst, 2 for matrix
	int		sortRead = 1;
	int		*localTopicCountTotal;
	int		rma = 0;
	int		wordLocCols;
	//MPI_Win	**winArr;
	//MPI_Comm	*commArr;
	int		CwordSize = 0;
	MPI_Win	wordCountWin, topicCountWin;
	int 	threads = 0;

	// need extra column in wordLocationMat if invMethod is 2
	if(invMethod == 1){
		wordLocCols = 2;
	}
	else{
		wordLocCols = 3;
	}
	
	#ifdef THREADS
	threads = 1;
	#endif


	// ************************************************
	// ************************************************
	// ***********       SETUP         ****************
	// ************************************************
	// ************************************************
	
	// MPI SETUP
	
	if(threads == 1){
		int requested = MPI_THREAD_MULTIPLE, provided;
		MPI_Init_thread(&argc, &argv, requested, &provided);
		if(provided != MPI_THREAD_MULTIPLE){
			printf("*** ERROR ***\nMPI_THREAD_MULTIPLE not supported\n\n");
			MPI_Finalize();
			return 0;
		}
	}		
	else{
		MPI_Init(&argc, &argv);
	}
	MPI_Comm mainComm = MPI_COMM_WORLD;
	MPI_Comm_rank(mainComm, &rank);
	MPI_Comm_size(mainComm, &size);

	// START TIMING
	startTimeMPI(&t1, mainComm);

	// PARSE ARGS
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn, &CwordMethod, &rma, &bundleSize, &pwToPdRatio, &invMethod, &sortRead);


	// FIND NO FILES
	findNoFiles(directory, &noFiles, maxFiles);

	// ASSIGN PW and PD
	procAssignment(rank, size, noFiles, &pdComm, &pwComm, pwToPdRatio, &pdSize, &pwSize, &prank, &procType, mainComm);

	// FIND PATH NAMES
	if(sortRead == 1){
		findPathNamesSorted(directory, &pathNames, &noFiles, maxFiles, pdSize);
	}
	else{
		findPathNames(directory, &pathNames, &noFiles, maxFiles);
	}
	free(directory);


	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma, threads, bundleSize, invMethod, sortRead);

	// READ IN FILES
	readInDocsAll(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, gutenFiles, pdComm);

	// FREE
	freePathNames(pathNames);
	free(stopfile);

	// PRINT INFO
	printInfoLDA(rank, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma, threads, bundleSize, invMethod, sortRead);
	if(rank == 0){printTime(t1Sec, t2Sec, "READING IN DOCS");}








	// ************************************************
	// ************************************************
	// ************      m, l & others     ************
	// ************************************************
	// ************************************************
	if(rank == 0){gettimeofday(&t1Sec, NULL);}
	
	// Pd
	if(procType == 'd'){
		// m
		callocArrInt(&wordSpread, vocabSize);

		if(CwordMethod == 1){
			callocArrInt(&wordLocation, vocabSize);
		}
		else{
			callocMatInt(&wordLocationMat, vocabSize, wordLocCols);
		}

		// ALPHA & BETA
		initialiseAlpha(&alphaDir, noTopics);		// freed
		initialiseBeta(&betaDir, vocabSize, &betaSum);	// freed

		// SET UP COUNTERS
		callocMatInt(&topicCountPerDoc, locNoFiles, noTopics);		// local C_doc	-- freed
		createTopicIndex(&topicIndex, locTotalCount, locNoFiles, locUniqueCount, localDocMatrix); // local z  -- freed
		free(locTotalCount);

		// FILL M
		fillM(procType, vocabSize, locNoFiles, locUniqueCount, localDocMatrix, wordSpread, pdComm);
	}
	else if(procType == 'w'){
		callocArrInt(&localTopicCountTotal, noTopics);
	}










	// ************************************************
	// ************************************************
	// *********     DISTRIBUTE WORDS       ***********
	// ************************************************
	// ************************************************
	
	if(CwordMethod == 1){
		distWordsToPw(procType, vocabSize, pdSize, pwSize, rank, size, prank, wordSpread, wordLocation, noTopics, &wordCountPerTopicTree, pdComm, mainComm);
	}	
	else if(CwordMethod == 2){
		distWordsToPwMatrix(procType, vocabSize, pdSize, pwSize, rank, size, prank, wordSpread, wordLocationMat, noTopics, &wordCountPerTopic, &CwordSize, wordLocCols, pdComm, mainComm);
	}

	if(rma > 0){
		//createComms(&commArr, pwSize, pdSize, mainComm);
		//createWindows(&winArr, pwSize, procType, commArr, prank, CwordSize, noTopics, wordCountPerTopic, localTopicCountTotal);	
		createWindowsNew(&wordCountWin, &topicCountWin, procType, mainComm, CwordSize, noTopics, wordCountPerTopic, localTopicCountTotal);	
	}


	// ************************************************
	// ************************************************
	// ***********     INVERTED INDEX      ************
	// ************************************************
	// ************************************************

	int **invertedIndexMatrix;
	int	invertedIndexSize;
	if(invMethod == 1 || CwordMethod == 1){
		buildInvertedIndex(procType, &invertedIndexTree, locNoFiles, locUniqueCount, localDocMatrix, prank, pdComm);
	}
	else{
		buildInvertedIndexMatrix(procType, &invertedIndexMatrix, locNoFiles, locUniqueCount, localDocMatrix, prank, vocabSize, wordLocationMat, &invertedIndexSize, pdComm);
	}
	


	// ************************************************
	// ************************************************
	// ************     INIT TOPICS       *************
	// ************************************************
	// ************************************************
	
	if(CwordMethod == 1){
		initTopicsPLDAplus(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocation, topicCountPerDoc, pdSize, size, wordCountPerTopicTree, rank, mainComm);
	}
	else if(CwordMethod == 2){
		if(rma == 0 || rma == 2){
			initTopicsPLDAplusMatrix(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocationMat, topicCountPerDoc, pdSize, size, wordCountPerTopic, localTopicCountTotal, rank, mainComm);
		}
		else{
		//	initTopicsPLDAplusMatrixRMA(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocationMat, topicCountPerDoc, pdSize, size, rank, pwSize, winArr, mainComm);
			initTopicsPLDAplusMatrixRMAFixed(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocationMat, topicCountPerDoc, pdSize, size, rank, pwSize, wordCountWin, topicCountWin, mainComm);
		}
		#ifdef DEBUG
		int r = 1120, c;
		if(rank == 2){
			int i, j;
			for(i = 0; i < CwordSize; i++){
				for(j = 0; j < noTopics; j++){
				//	wordCountPerTopic[i][j] = 100*i + j;
				}
			}
			for(c = 0; c < noTopics; c++){
			//	printf("On Pw entry %d is %d\n", c, wordCountPerTopic[r][c]);
				printf("On Pw entry %d is %d\n", c, localTopicCountTotal[c]);
			}
		}
		if(rank == 0){
			int *test;
			callocArrInt(&test, noTopics);
			MPI_Win_lock(MPI_LOCK_SHARED, pdSize, 0, winArr[0][1]); 
			//MPI_Get(test, noTopics, MPI_INT, pdSize, (r*noTopics), noTopics, MPI_INT, winArr[0][0]);
			MPI_Get(test, noTopics, MPI_INT, pdSize, 0, noTopics, MPI_INT, winArr[0][1]);
			MPI_Win_unlock(pdSize, winArr[0][1]);
			for(c = 0; c < noTopics; c++){
				printf("On Pd entry %d is %d\n", c, test[c]);
			}
		}
		if(rank == 1){
		}
		#endif
	}

	// ************************************************
	// ************************************************
	// ***********     WORD BUNDLES       *************
	// ************************************************
	// ************************************************
	
	if(CwordMethod == 1){
		createWordBundles(procType, invertedIndexTree, wordLocation, pwSize, pdSize, bundleSize, &bundles2d, &totNoBundles, prank, pdComm);
	}
	else if(CwordMethod == 2){
		if(invMethod == 1){
			createWordBundlesMatrix(procType, invertedIndexTree, wordLocationMat, pwSize, pdSize, bundleSize, &bundles2d, &bundles2dConv, &totNoBundles, prank, pdComm);
		}
		else{
			createWordBundlesMatrixInv(procType, invertedIndexMatrix, wordLocationMat, pwSize, pdSize, bundleSize, &bundles2d, &bundles2dConv, &totNoBundles, prank, invertedIndexSize, locNoFiles, pdComm);
		}
	}









	// ************************************************
	// ************************************************
	// ***********     CIRCULAR QUEUE       ***********
	// ************************************************
	// ************************************************
	circleQueue(procType, vocabSize, pdSize, prank, &circleStartBundle, bundles2d);




	// WORD BUNDLES MAY NEED TO BE ALTERED TO ACCOUNT FOR NEW matIndex


MPI_Barrier(mainComm);

	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************
	if(CwordMethod == 1){
		runPLDAplusGibbs(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopicTree, rank, mainComm);
	}
	else if(CwordMethod == 2){
		if(threads == 1){
			//runPLDAplusGibbsPthreads(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pdComm, mainComm);
			//runPLDAplusGibbsPthreadsBoth(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pdComm, mainComm);
			//runPLDAplusGibbsPthreadsDeadline(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
			runPLDAplusGibbsPthreadsDeadlineNew(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
		}
		else if(invMethod == 2){
			runPLDAplusGibbsMatrixInv(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexMatrix, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, wordLocationMat, mainComm);
		}
		else if(rma == 0 || rma == 1){
			runPLDAplusGibbsMatrix(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, mainComm);
		}
		else{
			//runPLDAplusGibbsMatrixRMA(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pwSize, winArr, mainComm);
			runPLDAplusGibbsMatrixRMAFixed(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pwSize, wordCountWin, topicCountWin, mainComm);
		}
	}

	if(procType == 'd'){
		free(alphaDir);
		freeTopicIndex(topicIndex, locNoFiles);
		freeDocMatrix(&localDocMatrix, locNoFiles, locUniqueCount);
		free(locUniqueCount);
	}




	// ************************************************
	// ************************************************
	// *******    COLLECT & PRINT RESULTS     *********     
	// ************************************************
	// ************************************************
	
	MPI_Barrier(mainComm);
	// Barrier to avoid interfering with Gibbs sampling requests etc
	if(CwordMethod == 1){
		collectAndPrint(rank, vocabSize, noTopics, pwSize, procType, wordCountPerTopicTree, printResults, betaDir, vocabulary, mainComm);
	}
	else if(CwordMethod == 2){
		collectAndPrintMatrix(rank, vocabSize, noTopics, procType, wordCountPerTopic, printResults, betaDir, vocabulary, wordLocationMat, pdSize, size, mainComm);
	}



	// ************************************************
	// ************************************************
	// *************     CLEAN UP      ****************     
	// ************************************************
	// ************************************************


	// FINISH TIMING
	endTimeMPI(t1, "WHOLE PROGRAM", rank, mainComm);


	// FREE REMAINING
	if(mainComm != MPI_COMM_WORLD){MPI_Comm_free(&mainComm);}
	if(procType == 'd'){
		free(betaDir);
	}


	MPI_Finalize();
	return 0;
}




void printHelp(){
	printf("-a [DOUBLE]\tSet Pw:Pd ratio\t\t\t\t\tDefault: 0.6\n");
	printf("-b [INT]\tSet Bundle Size\t\t\t\t\tDefault: 4\n");
	printf("-c [INT]\tSet Cword Method: 1=bst, 2=matrix\t\tDefault: 2\n");
	printf("-d [DIR]\tSet Directory\t\t\t\t\tDefault: ../../Gutenberg\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs\t\tDefault: 100\n");
	printf("-m [INT]\tSet max no files\t\t\t\tDefault: All Files\n");
	printf("-o\t\tSet random file read in\t\t\t\tDefault: Balances Load\n");
	printf("-p\t\tPrint Results\n");
	printf("-r\t\tSet Parallel Read In\t\t\t\tDefault: Serial & Bcast\n");
	printf("-s [INT]\tSet seed\t\t\t\t\tDefault: Random\n");
	printf("-t [INT]\tSet number of Topics\t\t\t\tDefault: 10\n");
	printf("-v [INT]\tSet Inv Method: 1=bst, 2=matrix\t\t\tDefault: 2 (Must use Cword Method 2)\n");
	printf("-w [INT]\tSet RMA Method: 0=off, 1=init, 2=gibbs, 3=both\tDefault: 0 (Must use Cword Method 2)\n");

	return;
}



void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma, int *bundleSize, double *pwToPdRatio, int *invMethod, int *sortRead){
	int pArgs[11];
	double pArgsD[1];
	int endProg = 0;

	//parse command line arguments on rank 0
	if(rank == 0){
		int opt;
		while((opt=getopt(argc,argv,"a:b:c:d:hi:m:oprs:t:v:w:"))!=-1){
			switch(opt){
				case 'a':
					*pwToPdRatio = strtod(optarg, NULL);
					if((*pwToPdRatio) <= 0 || (*pwToPdRatio) >= 1){
						printf("*** ERROR ***\nPw:Pd ratio must be in (0, 1)\n\n");
						endProg = 1;
					}
					break;
				case 'b':
					*bundleSize = atoi(optarg);
					if((*bundleSize) <= 0){
						printf("*** ERROR ***\nBundle Size must be non-zero positive integer\n\n");
						endProg = 1;
					}
					break;
				case 'c':
					*CwordMethod = atoi(optarg);
					if(*CwordMethod != 1 && *CwordMethod != 2){
						printf("*** ERROR ***\nCword Method must be in {1, 2}\n\n");
						endProg = 1;
					}
					break;					
				case 'd':
					//*directory = optarg;
					sprintf(*directory, optarg);
					break;
				case 'h':
					endProg = 1;
					break;
				case 'i':
					*iterations = atoi(optarg);
					break;
				case 'm':
					*maxFiles = atoi(optarg);
					break;
				case 'o':
					*sortRead = 0;
					break;
				case 'p':
					*printResults = 1;
					break;
				case 'r':
					*parallelReadIn = 1;
					break;
				case 's':
					*seed = atoi(optarg);
					break;
				case 't':
					*noTopics = atoi(optarg);
					break;
				case 'v':
					*invMethod = atoi(optarg);
					if(*invMethod != 1 && *invMethod != 2){
						printf("*** ERROR ***\nInv Method must be in {1, 2}\n\n");
						endProg = 1;
					}
					break;					
				case 'w':
					*rma = atoi(optarg);
					if(*rma < 0 || *rma > 3){
						printf("*** ERROR ***\nRMA Method must be in {0, 1, 2, 3}\n\n");
						endProg = 1;
					}
					break;					
				default:
					endProg = 1;
					break;
			}
		}
		if(*CwordMethod != 2 && *rma != 0){
			printf("*** ERROR ***\nRMA only avilable with matrix version of Cword\n\n");
			endProg = 1;
		}
		if(*CwordMethod != 2 && *invMethod != 1){
			printf("*** ERROR ***\nMatrix method for Inverted Index only avilable with matrix version of Cword\n\n");
			endProg = 1;
		}
		

		pArgs[0] = *seed;
		pArgs[1] = *noTopics;
		pArgs[2] = *iterations;
		pArgs[3] = *maxFiles;
		pArgs[4] = endProg;
		pArgs[5] = *parallelReadIn;
		pArgs[6] = *CwordMethod;
		pArgs[7] = *rma;
		pArgs[8] = *bundleSize;
		pArgs[9] = *invMethod;
		pArgs[10] = *sortRead;
		pArgsD[0] = *pwToPdRatio;
	}

	// Send relevant args to rest of processes
	MPI_Bcast(pArgs, 11, MPI_INT, 0, mainComm);
	MPI_Bcast(pArgsD, 1, MPI_DOUBLE, 0, mainComm);

	*seed = pArgs[0];
	*noTopics = pArgs[1];
	*iterations = pArgs[2];
	*maxFiles = pArgs[3];
	endProg = pArgs[4];
	*parallelReadIn = pArgs[5];
	*CwordMethod = pArgs[6];
	*rma = pArgs[7];
	*bundleSize = pArgs[8];
	*invMethod = pArgs[9];
	*sortRead = pArgs[10];
	*pwToPdRatio = pArgsD[0];

	srand48(78*rank + (*seed) );
	MPI_Bcast(*directory, 1024, MPI_CHAR, 0, mainComm);


	if(endProg){
		if(rank == 0){
			printHelp();
		}
		MPI_Finalize();
		exit(0);
	}
	
	return;
}




//# vim: ts=4
