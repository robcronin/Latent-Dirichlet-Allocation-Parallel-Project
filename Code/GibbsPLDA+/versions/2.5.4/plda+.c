#include "plda+.h"




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


// does relevant set up to read in documents and then uses either serial or parallel readin method
void	readInDocsAllNew(char procType, int prank, int pdSize, int noFiles, int *locNoFiles, int parallelReadIn, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults,  char **pathNames, int *totalWords, struct timeval *t1Sec, struct timeval *t2Sec, int gutenFiles, char *dictfile,  MPI_Comm pdComm){
	if(procType == 'd'){	
		startTimeMPI(t1Sec, pdComm);

		int *localFiles;
		
		// DETERMINE NUMBER OF FILES PER RANK
		filesPerProc(&localFiles, prank, pdSize, noFiles, locNoFiles, pdComm);

		// EACH PROC READS IN A PORTION AND THEN COLLABORATE
		if(parallelReadIn && dictfile == NULL){
			readInDocsParallel(prank, pdSize, stopfile, vocabSize, vocabulary, locUniqueCount, locTotalCount, localDocMatrix, printResults, pathNames, localFiles, *locNoFiles, totalWords, gutenFiles, pdComm);
		}
		else if(parallelReadIn){
			readInDocsParallelNew(prank, stopfile, vocabSize, vocabulary, locUniqueCount, locTotalCount, localDocMatrix, printResults, pathNames, localFiles, *locNoFiles, totalWords, gutenFiles, dictfile, pdComm);
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



