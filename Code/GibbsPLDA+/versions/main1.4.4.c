/*
 *
 *
 * 		STILL NEED TO CHECK AND FREE ALL MALLOCS
 * 		STILL GETTING ERRORS WITHOUT ACCTIME DUE TO WRONG TAGS RECEIVING
 *
 */



#include <mpi.h>
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

void 	printHelp();
void 	parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma);



struct	wordBundle{
	int	*indexes;
	int	location;
};


void	procAssignment(int rank, int size, int noFiles, MPI_Comm *pdComm, MPI_Comm *pwComm, double pwToPdRatio, int *pdSize, int *pwSize, int *prank, char *procType, MPI_Comm mainComm){
	if(size < 2){
		printf("*** ERROR ***\nPLDA+ requires at least 2 processes\n\n");
		MPI_Finalize();
		exit(1);
	}
	*pwSize = (int) rint(pwToPdRatio * size);
	*pdSize = size - *pwSize;

	if(*pdSize > noFiles){
		*pdSize = noFiles;
		*pwSize = size - *pdSize;
	}

	MPI_Group origGroup, pdGroup, pwGroup;
	MPI_Comm_group(mainComm, &origGroup);
	int *ranks = malloc((*pdSize) * sizeof(int));		// freed2
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



void	fillM(char procType, int vocabSize, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int *wordSpread, MPI_Comm pdComm){
	// INITIALISE M
	if(procType == 'd'){
		int doc, word;
		int *locWordSpread;
		callocArrInt(&locWordSpread, vocabSize);
		for(doc = 0; doc < locNoFiles; doc++){
			for(word = 0; word < locUniqueCount[doc]; word++){
				locWordSpread[ localDocMatrix[doc][word][0] ] = 1;
			}
		}
		MPI_Allreduce(locWordSpread, wordSpread, vocabSize, MPI_INT, MPI_SUM, pdComm);
		free(locWordSpread);
	}
	
	return;
}



void	distWordsToPw(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int *wordLocation, int noTopics, struct treeRootCword **wordCountPerTopicTree, MPI_Comm pdComm, MPI_Comm mainComm){ 

	struct timeval t1Sec;
	startTimeMPI(&t1Sec, mainComm);

	

	// DISTIBUTE WORDS
	if(procType == 'd'){
		int gap = (int)vocabSize/pdSize;
		int startPoint = rank*gap;
		int endPoint = (rank+1)*gap;
		if(prank == pdSize - 1){endPoint = vocabSize;}
			
		int len, word;
		int *locWordLocation;
		callocArrInt(&locWordLocation, vocabSize);
		int dest = prank % pwSize + pdSize;
		// distribute from largest spread down
		for(len = pdSize; len >= 0; len --){
			for(word = startPoint; word < endPoint; word ++){
				if(wordSpread[word] == len){
					MPI_Send(&word, 1, MPI_INT, dest, 2, mainComm);
					locWordLocation[word] = dest;
					// loops through Pws to send to
					dest ++;
					if(dest == size){dest = pdSize;}		
				}
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
	
		// consolidate wordLocation
		MPI_Allreduce(locWordLocation, wordLocation, vocabSize, MPI_INT, MPI_SUM, pdComm);
		free(locWordLocation);
	}
	// RECEIVE WORDS
	else if(procType == 'w'){
		int index;
		MPI_Status stat;
		createTreeCword(wordCountPerTopicTree, noTopics);
		int noFinished = 0;
		int total = 0;
		while(noFinished < pdSize){
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 2){
				addValueCword(*wordCountPerTopicTree, index);
				total ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			else{
				printf("\t\t\tERROR DIST WORDS RECEIVING WRONG THING\n");
			}
		}
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




void	distWordsToPwMatrix(char procType, int vocabSize, int pdSize, int pwSize, int rank, int size, int prank, int *wordSpread, int **wordLocationMat, int noTopics, int ***wordCountPerTopic, int *CwordSize, MPI_Comm pdComm, MPI_Comm mainComm){ 

	struct timeval t1Sec;
	startTimeMPI(&t1Sec, mainComm);
	MPI_Status stat;

	

	// DISTIBUTE WORDS
	if(procType == 'd'){
		int gap = (int)vocabSize/pdSize;
		int startPoint = rank*gap;
		int endPoint = (rank+1)*gap;
		if(prank == pdSize - 1){endPoint = vocabSize;}
			
		int len, word;
		int **locWordLocation;
		callocMatInt(&locWordLocation, vocabSize, 2);
		int dest = prank % pwSize + pdSize;
		// distribute from largest spread down
		for(len = pdSize; len >= 0; len --){
			for(word = startPoint; word < endPoint; word ++){
				if(wordSpread[word] == len){
					MPI_Send(&word, 1, MPI_INT, dest, 2, mainComm);
					locWordLocation[word][0] = dest;
					MPI_Recv(&(locWordLocation[word][1]), 1, MPI_INT, dest, 2, mainComm, &stat);
				//	printf("%d: mat loc is %d\n", rank, locWordLocation[word][1]);
					// loops through Pws to send to
					dest ++;
					if(dest == size){dest = pdSize;}		
				}
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
	
		// consolidate wordLocation
		MPI_Allreduce(locWordLocation[0], wordLocationMat[0], vocabSize*2, MPI_INT, MPI_SUM, pdComm);
		freeMatInt(locWordLocation);
	}
	// RECEIVE WORDS
	else if(procType == 'w'){
		int index;
		int noFinished = 0;
		int matIndex = 0;
		while(noFinished < pdSize){
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
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

		// Create Cword matrix
		callocMatInt(wordCountPerTopic, matIndex, noTopics);
		*CwordSize = matIndex;

	}

	endTimeMPI(t1Sec, "Cw MATRIX (Pw)", rank, mainComm);

	return;
}




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
	

//XXX
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


	
void	createWindowsNew(MPI_Win *wordCountWin, MPI_Win *topicCountWin, char procType, MPI_Comm mainComm, int CwordSize, int noTopics, int **wordCountPerTopic, int *localTopicCountTotal){
	
	if(procType == 'd'){
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, mainComm, wordCountWin); 
			MPI_Win_fence(0, *wordCountWin);
			MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, mainComm, topicCountWin); 
			MPI_Win_fence(0, *topicCountWin);
	}
	else if(procType == 'w'){
			MPI_Win_create(&(wordCountPerTopic[0][0]), CwordSize*noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, mainComm, wordCountWin);	
			MPI_Win_fence(0, *wordCountWin);
			MPI_Win_create(localTopicCountTotal, noTopics*sizeof(int), sizeof(int), MPI_INFO_NULL, mainComm, topicCountWin);	
			MPI_Win_fence(0, *topicCountWin);
	}		

	return;
}


	
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
			rebalanceTreeInv(*invertedIndexTree);
		}
	
		endTimeMPI(t1Sec, "INV INDEX (Pd)", prank, pdComm);
		
	}

	return;
}



void	readInDocsAll(char procType, int prank, int pdSize, int noFiles, int *locNoFiles, int parallelReadIn, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults,  char **pathNames, int *totalWords, struct timeval *t1Sec, struct timeval *t2Sec, MPI_Comm pdComm){
	if(procType == 'd'){	
		startTimeMPI(t1Sec, pdComm);

		int *localFiles;
		
		// DETERMINE NUMBER OF FILES PER RANK
		filesPerProc(&localFiles, prank, pdSize, noFiles, locNoFiles, pdComm);


		if(parallelReadIn){
			readInDocsParallel(prank, pdSize, stopfile, vocabSize, vocabulary, locUniqueCount, locTotalCount, localDocMatrix, printResults, pathNames, localFiles, *locNoFiles, totalWords, pdComm);
		}
		else{
			// READ IN DOCS ON RANK 0 AND DISTRIBUTE INFO
			readInDocsSerial(prank, pdSize, pathNames, noFiles, stopfile, vocabSize, vocabulary, printResults, locUniqueCount, locTotalCount, *locNoFiles, localFiles, localDocMatrix, totalWords, pdComm);
		}
		free(localFiles);
		startTimeMPI(t2Sec, pdComm);
	}
	
	return;
}



void	initTopicsPLDAplus(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int *wordLocation, int **topicCountPerDoc, int pdSize, int size, struct treeRootCword *wordCountPerTopicTree, int rank, MPI_Comm mainComm){ 
	int 	sendBuff2[2];
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		int d, w, occ;
		int topic, word;
		for(d = 0; d < locNoFiles; d ++){
			for(w = 0; w < locUniqueCount[d]; w ++){
				word = localDocMatrix[d][w][0];
				for(occ = 0; occ < localDocMatrix[d][w][1]; occ ++){
					topic = (int) (drand48() * noTopics);
					topicIndex[d][w][occ] = topic;		// z
					topicCountPerDoc[d][topic] ++;		// Cdoc
					//localTopicCountTotal[topic] ++;		// Ctopic

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
		while(noFinished < pdSize){
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 3){
				current = searchTreeCword(wordCountPerTopicTree, sendBuff2[0]);
				current -> topicArray[sendBuff2[1]] ++;
				wordCountPerTopicTree -> topicCountTotal[sendBuff2[1]] ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
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



void	initTopicsPLDAplusMatrix(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int **wordCountPerTopic, int *localTopicCountTotal, int rank, MPI_Comm mainComm){ 
	int 	sendBuff2[2];
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		int d, w, occ;
		int topic, word;
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
		while(noFinished < pdSize){
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 3){
				wordCountPerTopic[sendBuff2[0]][sendBuff2[1]] ++;
				localTopicCountTotal[sendBuff2[1]] ++;
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
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
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 3){
			//	wordCountPerTopic[sendBuff2[0]][sendBuff2[1]] ++;
			//	localTopicCountTotal[sendBuff2[1]] ++;
			//	NOT USED FOR RMA
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
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



void	initTopicsPLDAplusMatrixRMAFixed(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm){ 
	struct	timeval t1;
	startTimeMPI(&t1, mainComm);
	int sendBuff2[2];
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

					int win = wordLocationMat[word][0] - pdSize;
					int dest = wordLocationMat[word][0];
					int disp = wordLocationMat[word][1] * noTopics + topic;

					// update Cword
					MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dest, 0, wordCountWin);
					MPI_Accumulate(&one, 1, MPI_INT, dest, disp, 1, MPI_INT, MPI_SUM, wordCountWin);
					MPI_Win_unlock(dest, wordCountWin);

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
	// Receive Cw
	else if(procType == 'w'){
		int noFinished = 0;
		MPI_Status stat;
		while(noFinished < pdSize){
			MPI_Recv(&sendBuff2, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 3){
			//	wordCountPerTopic[sendBuff2[0]][sendBuff2[1]] ++;
			//	localTopicCountTotal[sendBuff2[1]] ++;
			//	NOT USED FOR RMA
			}
			// Counts number of Pd's who send finish signal
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
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



void	createWordBundles(char procType, struct treeRootInv *invertedIndexTree, int *wordLocation, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int *totNoBundles, int prank, MPI_Comm pdComm){


	// CREATE WORD BUNDLES
	if(procType == 'd'){
		struct timeval t1;
		startTimeMPI(&t1, pdComm);

		int i, j;
		// Sort by occurences
		int 	**occMat;
		int	totalWords = invertedIndexTree -> size;
		allocMatInt(&occMat, totalWords, 3);		// NOT freed
		fillOccMat(occMat, invertedIndexTree, wordLocation);
		

		// sort by index (for ability to do circular queue)
		quicksortMat(occMat, totalWords, 3, 0);

		// sort by location
		quicksortMat(occMat, totalWords, 3, 2);

		// find break points between location change
		int *breakers = malloc((pwSize+1) * sizeof(int));	// NOT freed
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

		
		// Find how many bundles
		int *noBundles = malloc((pwSize + 1) * sizeof(int));		// NOT freed
		noBundles[pwSize] = 0;
		for(i = 0; i < pwSize; i++){
			noBundles[i] = (breakers[i+1] - breakers[i]) / bundleSize;
			if( (breakers[i+1] - breakers[i]) % bundleSize != 0){
				noBundles[i] ++;
			}
			noBundles[pwSize] += noBundles[i];
		}

		// Create bundle array
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
		for(loc = 0; loc < pwSize; loc++){
			// low is lowest occ for current loc, high is highest occ
			low = breakers[loc];
			high = breakers[loc + 1] - 1;
			for(bun = 0; bun < noBundles[loc] - 1; bun++){
				bundles3d[loc][bun][0] = occMat[high][0];
				high --;
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = occMat[low][0];
					low ++;
				}
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
		free(bundles3d);
		*totNoBundles = noBundles[pwSize];
		free(noBundles);

		
		// Sort by index to randomise location
		quicksortMat(*bundles2d, *totNoBundles, bundleSize + 1, 0); 

		endTimeMPI(t1, "WORD BUNDLES", prank, pdComm);
	}

	return;
}



void	createWordBundlesMatrix(char procType, struct treeRootInv *invertedIndexTree, int **wordLocationMat, int pwSize, int pdSize, int bundleSize, int ***bundles2d, int ***bundles2dConv, int *totNoBundles, int prank, MPI_Comm pdComm){


	// CREATE WORD BUNDLES
	if(procType == 'd'){
		struct timeval t1;
		startTimeMPI(&t1, pdComm);

		int i, j;
		// Sort by occurences
		int 	**occMat;
		int	totalWords = invertedIndexTree -> size;
		allocMatInt(&occMat, totalWords, 3);		// NOT freed
		fillOccMatMatrix(occMat, invertedIndexTree, wordLocationMat);
		

		// sort by index (for ability to do circular queue)
		quicksortMat(occMat, totalWords, 3, 0);

		// sort by location
		quicksortMat(occMat, totalWords, 3, 2);

		// find break points between location change
		int *breakers = malloc((pwSize+1) * sizeof(int));	// NOT freed
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

		
		// Find how many bundles
		int *noBundles = malloc((pwSize + 1) * sizeof(int));		// NOT freed
		noBundles[pwSize] = 0;
		for(i = 0; i < pwSize; i++){
			noBundles[i] = (breakers[i+1] - breakers[i]) / bundleSize;
			if( (breakers[i+1] - breakers[i]) % bundleSize != 0){
				noBundles[i] ++;
			}
			noBundles[pwSize] += noBundles[i];
		}

		// Create bundle array
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
		for(loc = 0; loc < pwSize; loc++){
			// low is lowest occ for current loc, high is highest occ
			low = breakers[loc];
			high = breakers[loc + 1] - 1;
			for(bun = 0; bun < noBundles[loc] - 1; bun++){
				bundles3d[loc][bun][0] = occMat[high][0];
				high --;
				for(i = 1; i < bundleSize; i++){
					bundles3d[loc][bun][i] = occMat[low][0];
					low ++;
				}
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
		free(bundles3d);
		*totNoBundles = noBundles[pwSize];

		
		// Sort by index to randomise location
		quicksortMat(*bundles2d, *totNoBundles, bundleSize + 1, 0); 

		// Fill in converted bundle indexes
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
		endTimeMPI(t1, "WORD BUNDLES", prank, pdComm);
	}

	return;
}



void	circleQueue(char procType, int vocabSize, int pdSize, int prank, int *circleStartBundle, int **bundles2d){
	if(procType == 'd'){
		int gap = vocabSize / pdSize;
		int wordIndex = prank * gap;
	
		*circleStartBundle = 0;
		while(bundles2d[*circleStartBundle][0] < wordIndex){
			(*circleStartBundle) ++;
		} 
	}

	return;
}



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

		// Fetch Ctopic
		zeroArrInt(topicCountTotal, noTopics);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
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
				dest = bundles2d[bun][bundleSize] + pdSize;
				
				// fetch Cword
				MPI_Send(&iter, 1, MPI_INT, dest, 5, mainComm);
				MPI_Send(bundles2d[bun], bundleSize, MPI_INT, dest, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);
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
				MPI_Send(&iter, 1, MPI_INT, dest, 6, mainComm);
				MPI_Send(bundles2d[bun], bundleSize, MPI_INT, dest, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest, 6, mainComm);

				/*
				// Fetch Ctopic
				zeroArrInt(topicCountTotal, noTopics);
				for(i = pdSize; i < size; i++){
					MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
					MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
					for(j = 0; j < noTopics; j++){
						topicCountTotal[j] += localTopicCountTotal[j];
					}
				}*/
			}
			// Fetch Ctopic
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

	}

	// Pw processes requests
	if(procType == 'w'){
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;
		struct nodeCword *current;

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

	}
	
	endTimeMPI(t1, "GIBBS SAMPLING", rank, mainComm);
	return;
}



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

		// Fetch Ctopic
		zeroArrInt(topicCountTotal, noTopics);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
			MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
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
				dest = bundles2d[bun][bundleSize] + pdSize;
				
				// fetch Cword
				MPI_Send(&iter, 1, MPI_INT, dest, 5, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 5, mainComm);
				MPI_Recv( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 5, mainComm, &stat);
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
				MPI_Send(&iter, 1, MPI_INT, dest, 6, mainComm);
				MPI_Send(bundles2dConv[bun], bundleSize, MPI_INT, dest, 6, mainComm);
				MPI_Send( &(CwordBuffer[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
				MPI_Send(CtopicBuffer, noTopics, MPI_INT, dest, 6, mainComm);

				/*
				// Fetch Ctopic
				zeroArrInt(topicCountTotal, noTopics);
				for(i = pdSize; i < size; i++){
					MPI_Send(&i, 1, MPI_INT, i, 4, mainComm);
					MPI_Recv(localTopicCountTotal, noTopics, MPI_INT, i, 4, mainComm, &stat);
					for(j = 0; j < noTopics; j++){
						topicCountTotal[j] += localTopicCountTotal[j];
					}
				}*/
			}
			// Fetch Ctopic
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
	
	endTimeMPI(t1, "GIBBS SAMPLING MATRIX", rank, mainComm);
	return;
}



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
			MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 7){
				MPI_Recv(topicArray, noTopics, MPI_INT, stat.MPI_SOURCE, 7, mainComm, &stat);
				for(i = 0; i < noTopics; i++){
					wordCountPerTopic[index][i] = topicArray[i];
				}
			}
			else if(stat.MPI_TAG == 0){
				noFinished++;
			}
			else{
				printf("\t\tERROR PRINT RECVING OTHER FEATURE\n");
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
			
			MPI_Send(&matIndex, 1, MPI_INT, dest, 7, mainComm);
			MPI_Recv(topicArray, noTopics, MPI_INT, dest, 7, mainComm, &stat);
			for(k = 0; k < noTopics; k++){
				wordCountPerTopicFull[i][k] = topicArray[k];
			}
		}
		// send finish signal
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
		
	}

	if(procType == 'w'){
		int noFinished = 0;	
		int matIndex;
		MPI_Status stat;
		while(noFinished < 1){
			MPI_Recv(&matIndex, 1, MPI_INT, 0, MPI_ANY_TAG, mainComm, &stat);
			if(stat.MPI_TAG == 7){
				MPI_Send(wordCountPerTopic[matIndex], noTopics, MPI_INT, 0, 7, mainComm);
			}
			else if(stat.MPI_TAG == 0){
				noFinished ++;
			}
			else{
				printf("\t\tERROR PRINT RECVING OTHER FEATURE\n");
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
	int	printResults = 0;
	int 	maxFiles = 0;
	int	iterations = 100;
	int 	parallelReadIn = 0;
	int 	LDAmethod = 3;
	double	pwToPdRatio = 0.6;
	int	bundleSize = 4;

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
	int		*localTopicCountTotal;
	int		rma = 1;
	//MPI_Win	**winArr;
	//MPI_Comm	*commArr;
	int		CwordSize = 0;
	MPI_Win	wordCountWin, topicCountWin;




	// ************************************************
	// ************************************************
	// ***********       SETUP         ****************
	// ************************************************
	// ************************************************
	
	// MPI SETUP
	MPI_Init(&argc, &argv);
	MPI_Comm mainComm = MPI_COMM_WORLD;
	MPI_Comm_rank(mainComm, &rank);
	MPI_Comm_size(mainComm, &size);

	// START TIMING
	startTimeMPI(&t1, mainComm);

	// PARSE ARGS
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn, &CwordMethod, &rma);

	// FIND NO FILES & PATH NAMES
	findPathNames(directory, &pathNames, &noFiles, maxFiles);
	free(directory);

	// ASSIGN PW and PD
	procAssignment(rank, size, noFiles, &pdComm, &pwComm, pwToPdRatio, &pdSize, &pwSize, &prank, &procType, mainComm);

	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma);

	// READ IN FILES
	readInDocsAll(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, pdComm);

	// FREE
	freePathNames(pathNames);
	free(stopfile);

	// PRINT INFO
	printInfoLDA(rank, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma);
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
			callocMatInt(&wordLocationMat, vocabSize, 2);
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
		distWordsToPwMatrix(procType, vocabSize, pdSize, pwSize, rank, size, prank, wordSpread, wordLocationMat, noTopics, &wordCountPerTopic, &CwordSize, pdComm, mainComm);
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

	buildInvertedIndex(procType, &invertedIndexTree, locNoFiles, locUniqueCount, localDocMatrix, prank, pdComm);



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
		createWordBundlesMatrix(procType, invertedIndexTree, wordLocationMat, pwSize, pdSize, bundleSize, &bundles2d, &bundles2dConv, &totNoBundles, prank, pdComm);
	}









	// ************************************************
	// ************************************************
	// ***********     CIRCULAR QUEUE       ***********
	// ************************************************
	// ************************************************
	circleQueue(procType, vocabSize, pdSize, prank, &circleStartBundle, bundles2d);




	// WORD BUNDLES MAY NEED TO BE ALTERED TO ACCOUNT FOR NEW matIndex




	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************
	if(CwordMethod == 1){
		runPLDAplusGibbs(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopicTree, rank, mainComm);
	}
	else if(CwordMethod == 2){
		if(rma == 0 || rma == 1){
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
	printf("-c [INT]\tSet Cword Method: 1=bst, 2=matrix (default: 2)\n");
	printf("-d [DIR]\tSet Directory (default: ../../Gutenberg)\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs (default: 1000)\n");
	printf("-m [INT]\tSet max no files (default: all files)\n");
	printf("-p\t\tPrint Results\n");
	printf("-s [INT]\tSet seed (default: random)\n");
	printf("-t [INT]\tSet number of Topics (default: 10)\n");
	printf("-w [INT]\tSet RMA Method: 0=none, 1=init, 2=gibbs, 3=both (default: 1) Must use Cword Method 2\n");

	return;
}



void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma){
	int pArgs[8];
	int endProg = 0;

	//parse command line arguments on rank 0
	if(rank == 0){
		int opt;
		while((opt=getopt(argc,argv,"c:d:hi:m:prs:t:w:"))!=-1){
			switch(opt){
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

		pArgs[0] = *seed;
		pArgs[1] = *noTopics;
		pArgs[2] = *iterations;
		pArgs[3] = *maxFiles;
		pArgs[4] = endProg;
		pArgs[5] = *parallelReadIn;
		pArgs[6] = *CwordMethod;
		pArgs[7] = *rma;
	}

	// Send relevant args to rest of processes
	MPI_Bcast(pArgs, 8, MPI_INT, 0, mainComm);
	*seed = pArgs[0];
	*noTopics = pArgs[1];
	*iterations = pArgs[2];
	*maxFiles = pArgs[3];
	endProg = pArgs[4];
	*parallelReadIn = pArgs[5];
	*CwordMethod = pArgs[6];
	*rma = pArgs[7];
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
