#include "initTopics.h"



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




// initialises and distributes topics to relevant locations (MATRIX METHOD with RMA)
void	initTopicsPLDAplusMatrixRMA(char procType, int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **wordLocationMat, int **topicCountPerDoc, int pdSize, int size, int rank, int pwSize, MPI_Win wordCountWin, MPI_Win topicCountWin, MPI_Comm mainComm){ 

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


//# vim: ts=4
