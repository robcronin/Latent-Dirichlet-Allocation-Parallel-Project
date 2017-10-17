#include "gibbsLDA-MPI.h"


// Determines the number of files to be stored on each rank
void	filesPerProc(int **localFiles, int rank, int size, int noFiles, int *locNoFiles, MPI_Comm mainComm){
	int i;

	// stores amount for each process
	(*localFiles) = malloc(size * sizeof(int));	// freed

	// rank 0 calculates and then distributes
	if(rank == 0){
		// starts by getting the minimum per proc
		*locNoFiles = noFiles/size;
		for(i = 0; i < size; i++){
			(*localFiles)[i] = *locNoFiles;
		}

		// then adds on any remainders
		int total = (*locNoFiles)*size;
		i = 0;
		while(total < noFiles){
			(*localFiles)[i] ++;
			total ++;
			i++;
		}
	}
	MPI_Bcast(*localFiles, size, MPI_INT, 0, mainComm);
	*locNoFiles = (*localFiles)[rank];

	return;

}



// In the case of more processes than files (only for testing small scale)
void	alterComm(int rank, int *size, int noFiles, MPI_Comm *mainComm){
	MPI_Comm newComm;
	int i;

	if(*size > noFiles){
		// set new size
		*size = noFiles;

		// create new group
		MPI_Group originalGroup, mainGroup;
		MPI_Comm_group(*mainComm, &originalGroup);
		int *ranks = malloc((*size) * sizeof(int));	// freed
		for(i = 0; i < (*size); i++){
			ranks[i] = i;
		}
		MPI_Group_incl(originalGroup, *size, ranks, &mainGroup);

		// create new comm
		MPI_Comm_create(*mainComm, mainGroup, &newComm);	// freed
		*mainComm = newComm;

		free(ranks);
		
	}
	// kill any excess processes
	if(rank >= (*size)){
		printf("Killing Process %d\n", rank);
		MPI_Finalize();
		exit(0);
	}

	return;
}



// Reads in all documents on rank 0 and then distributes chunks to all other processes
void	readInDocsSerial(int rank, int size, char **pathNames, int noFiles, char *stopfile, int *vocabSize, char ***vocabulary, int printResults, int **locUniqueCount, int **locTotalCount, int locNoFiles, int *localFiles, int ****localDocMatrix, int *totalWords, int gutenFiles, MPI_Comm mainComm){
	struct 	alpha *phi, **all;
	int 	***wholeDocMatrix;
	int	*uniqueWordCounts, *totalWordCounts, *totalUniques;

	if(rank == 0){
	
		// SETUP CORPUS VOCABULARY
		// Reads in all documents and creates an overall vocabulary of all used words
		initialiseLocalTFIDF(pathNames, &all, noFiles, stopfile, 0, gutenFiles);	// freed
		create_dict(&phi);		// freed
		createBeta(all, noFiles, phi);	// phi in literature not beta	-- freed
		*vocabSize = phi -> uniqueCount;
		createVocab(phi, vocabulary);	// freed
		freeBeta(&phi);

		// SET UP OVERALL DOC MATRIX
		// Counts unique and total amounts per doc to allocate docMatrix correctly and then fills
		createWordCounts(&uniqueWordCounts, &totalWordCounts, noFiles, all);	// freed
		*totalWords = totalWordCounts[noFiles];
		allocDocMatrix(&wholeDocMatrix, noFiles, uniqueWordCounts);		// freed
		fillDocMatrix(wholeDocMatrix, *vocabulary, all, noFiles);
		freeLocalTFIDF(&all, noFiles);
	
		// All words now referred to as indexes (only need vocabulary if need to get words back to print)
		if(!printResults){
			freeVocab(vocabulary, *vocabSize);
		}
	}
	MPI_Bcast(vocabSize, 1, MPI_INT, 0, mainComm);

	// DISTRIBUTE WORD COUNTS
	distWordCounts(locUniqueCount, locTotalCount, uniqueWordCounts, totalWordCounts, locNoFiles, rank, size, localFiles, &totalUniques, mainComm);
	if(rank == 0){free(totalWordCounts);}
	
	// DISTRIBUTE DOC MATRIX
	distDocMatrix(localDocMatrix, locNoFiles, localFiles, *locUniqueCount, rank, size, wholeDocMatrix, totalUniques, mainComm);
	free(totalUniques);

	// free main doc matrix
	if(rank == 0){
		freeDocMatrix(&wholeDocMatrix, noFiles);
		free(uniqueWordCounts);
	}

	return;
}



// Distributes word counts to processes (after having read in on rank 0)
void	distWordCounts(int **locUniqueCount, int **locTotalCount, int *uniqueWordCounts, int *totalWordCounts, int locNoFiles, int rank, int size, int *localFiles, int **totalUniques, MPI_Comm mainComm){ 
	int i, j;
	int index;
	MPI_Status stat;

	(*locUniqueCount) = malloc( (locNoFiles+1) * sizeof(int) );	// freed
	(*locTotalCount) = malloc( (locNoFiles+1) * sizeof(int) );	// freed
	if(rank == 0){
		// send word counts to each proc
		index = localFiles[0];
		for(i = 1; i < size; i++){
			MPI_Send(&(uniqueWordCounts[index]), localFiles[i], MPI_INT, i, 0, mainComm);
			MPI_Send(&(totalWordCounts[index]), localFiles[i], MPI_INT, i, 0, mainComm);
			index += localFiles[i];
		}

		// fill in rank 0 local versions
		(*locUniqueCount)[locNoFiles] = 0;
		(*locTotalCount)[locNoFiles] = 0;
		for(i = 0; i < locNoFiles; i++){
			(*locUniqueCount)[i] = uniqueWordCounts[i];
			(*locUniqueCount)[locNoFiles] += (*locUniqueCount)[i];
			(*locTotalCount)[i] = totalWordCounts[i];
			(*locTotalCount)[locNoFiles] += (*locTotalCount)[i];
		}
	}
	else{
		// receive into every other rank
		MPI_Recv(*locUniqueCount, localFiles[rank], MPI_INT, 0, 0, mainComm, &stat);
		MPI_Recv(*locTotalCount, localFiles[rank], MPI_INT, 0, 0, mainComm, &stat);

		// find local totals
		(*locUniqueCount)[locNoFiles] = 0;
		(*locTotalCount)[locNoFiles] = 0;
		for(i = 0; i < locNoFiles; i++){
			(*locUniqueCount)[locNoFiles] += (*locUniqueCount)[i];
			(*locTotalCount)[locNoFiles] += (*locTotalCount)[i];
		}
	}

	// FIND TOTAL UNIQUES PER PROC
	*totalUniques = malloc(size * sizeof(int));	// freed
	if(rank == 0){
		index = 0;
		for(i = 0; i < size; i++){
			(*totalUniques)[i] = 0;
			for(j = 0; j < localFiles[i]; j++){
				(*totalUniques)[i] += uniqueWordCounts[index];
				index ++;
			}
		}
	}

	return;

}



// Distributes doc matrix in chunks
void	distDocMatrix(int ****localDocMatrix, int locNoFiles, int *localFiles, int *locUniqueCount, int rank, int size, int ***wholeDocMatrix, int *totalUniques, MPI_Comm mainComm){
	int i, j, index;
	MPI_Status stat;

	// all procs can allocate their local docMatrix now having received total/unique counts
	allocDocMatrix(localDocMatrix, locNoFiles, locUniqueCount);	// freed
	if(rank == 0){
		// rank 0 distributes docMatrix based on localFiles (number of docs per rank)
		index = localFiles[0];
		for(i = 1; i < size; i++){
			MPI_Send( &(wholeDocMatrix[index][0][0]), 2 * totalUniques[i], MPI_INT, i, 0, mainComm);
			index += localFiles[i];
		}
		
		// fills in its own local version
		for(i = 0; i < locNoFiles; i++){
			for(j = 0; j < locUniqueCount[i]; j++){
				(*localDocMatrix)[i][j][0] = wholeDocMatrix[i][j][0];
				(*localDocMatrix)[i][j][1] = wholeDocMatrix[i][j][1];
			}
		}


	}
	// Rest of ranks wait to receive their docMatrix
	else{
		MPI_Recv(&((*localDocMatrix)[0][0][0]), 2 * locUniqueCount[locNoFiles], MPI_INT, 0, 0, mainComm, &stat);
	}

	return;

}



// Reads in local files on each rank and sends around an overall vocab so everyone can index their words
void	readInDocsParallel(int rank, int size, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults, char **pathNames, int *localFiles, int locNoFiles, int *totalWords, int gutenFiles, MPI_Comm mainComm){
	struct 	alpha *phi, **localAll;
	int i;

	
	// find start index in path names for given rank
	int startIndex = 0;
	for(i = 0; i < rank; i++){
		startIndex += localFiles[i];
	}

	// fill in local files
	initialiseLocalTFIDF(pathNames, &localAll, locNoFiles, stopfile, startIndex, gutenFiles);	// freed

	create_dict(&phi);	// freed

	createBeta(localAll, locNoFiles, phi);	// phi not beta -- freed

	// 0 fills sends to 1, who fills and sends to 2 and so on
	for(i = 0; i < size-1; i++){
		MPI_Barrier(mainComm);
		sendPhi(phi, i, (i+1)%size, rank, mainComm);
	}

	// then size-1 has full vocab and bcasts to everyone
	bcastPhi(phi, size-1, rank, mainComm);
	*vocabSize = phi -> uniqueCount;
	createVocab(phi, vocabulary);	// freed
	freeBeta(&phi);

	// find local/total counts and fill in docMatrix
	createWordCounts(locUniqueCount, locTotalCount, locNoFiles, localAll);	// freed
	allocDocMatrix(localDocMatrix, locNoFiles, *locUniqueCount);		// freed
	fillDocMatrix(*localDocMatrix, *vocabulary, localAll, locNoFiles);
	freeLocalTFIDF(&localAll, locNoFiles);

	// reduces to find the total overall number of words
	*totalWords = 0;
	MPI_Allreduce(&((*locTotalCount)[locNoFiles]), totalWords, 1, MPI_INT, MPI_SUM, mainComm);

	if(!printResults){
		freeVocab(vocabulary, *vocabSize);
	}

	return;
}



// reads in in parallel but uses a common dictionary to index words to avoid costly communication
// 	** Prooves very costly memory wise as counters need to be created for all words, many of which never occur
void	readInDocsParallelDict(int rank, char *stopfile, int *vocabSize, char ***vocabulary, int **locUniqueCount, int **locTotalCount, int ****localDocMatrix, int printResults, char **pathNames, int *localFiles, int locNoFiles, int *totalWords, int gutenFiles, char *dictfile, MPI_Comm mainComm){
	struct 	alpha *phi, **localAll;
	int i;

	
	// find start index in path names for given rank
	int startIndex = 0;
	for(i = 0; i < rank; i++){
		startIndex += localFiles[i];
	}


	create_dict(&phi);	// freed

	
	// fill in dictionary from file
	struct alpha *holder;
	create_dict(&holder);
	fill_dict(dictfile, &phi, holder, 0);
	free_dict(&holder);


	// create common vocab to index against
	*vocabSize = phi -> uniqueCount;
	createVocab(phi, vocabulary);	// freed


	// fill in local files but ignore any words not in given dictionary
	initialiseLocalTFIDFDict(pathNames, &localAll, locNoFiles, stopfile, startIndex, gutenFiles, phi);	// freed
	freeBeta(&phi);

	// find local/total counts and fill in docMatrix
	createWordCounts(locUniqueCount, locTotalCount, locNoFiles, localAll);	// freed
	allocDocMatrix(localDocMatrix, locNoFiles, *locUniqueCount);		// freed
	fillDocMatrix(*localDocMatrix, *vocabulary, localAll, locNoFiles);
	freeLocalTFIDF(&localAll, locNoFiles);

	// reduces to find the total overall number of words
	*totalWords = 0;
	MPI_Allreduce(&((*locTotalCount)[locNoFiles]), totalWords, 1, MPI_INT, MPI_SUM, mainComm);

	if(!printResults){
		freeVocab(vocabulary, *vocabSize);
	}

	return;
}



// Sends the vocab phi from src to dest
void	sendPhi(struct alpha *phi, int src, int dest, int rank, MPI_Comm mainComm){
	int sendCount;
	int allocLength = 30;		// if you change this also change 28 in add_word (tfidf.c) and below
	char **phiBuffer;
	MPI_Status stat;

	// On the source end
	if(rank == src){
		// allocates and fills a contiguous phi buffer
		sendCount = phi -> uniqueCount;
		allocPhiBuffer(&phiBuffer, sendCount, allocLength);	// freed
		fillPhiBuffer(phi, phiBuffer, allocLength);
		
		// sends the amount and then the buffer
		MPI_Send(&sendCount, 1, MPI_INT, dest, 0, mainComm);
		MPI_Send(&(phiBuffer[0][0]), sendCount*allocLength, MPI_CHAR, dest, 0, mainComm);
		freePhiBuffer(phiBuffer);
	}
	// At destination
	else if(rank == dest){
		// receives the count and allocates accordingly before receiving
		MPI_Recv(&sendCount, 1, MPI_INT, src, 0, mainComm, &stat);
		allocPhiBuffer(&phiBuffer, sendCount, allocLength);	// freed
		MPI_Recv(&(phiBuffer[0][0]), sendCount*allocLength, MPI_CHAR, src, 0, mainComm, &stat);

		// fills its own phi from this buffer
		fillBetaFromBuffer(phi, phiBuffer, sendCount);
		freePhiBuffer(phiBuffer);
	}

	return;
}



// Broadcasts a vocab phi 
void	bcastPhi(struct alpha *phi, int src, int rank, MPI_Comm mainComm){
	int sendCount;
	int allocLength = 30;		// if you change this also change 28 in add_word (tfidf.c) and above
	char **phiBuffer;

	// Source allocates and fills contiguous buffer
	if(rank == src){
		sendCount = phi -> uniqueCount;
		allocPhiBuffer(&phiBuffer, sendCount, allocLength);	// freed
		fillPhiBuffer(phi, phiBuffer, allocLength);
	}

	// broadcasts size
	MPI_Bcast(&sendCount, 1, MPI_INT, src, mainComm);

	// everyone else allocates appropriate memory
	if(rank != src){
		allocPhiBuffer(&phiBuffer, sendCount, allocLength);	// freed
	}
	
	// then the buffer is broadcast
	MPI_Bcast(&(phiBuffer[0][0]), sendCount*allocLength, MPI_CHAR, src, mainComm);

	// everyone else fills in from this buffer
	if(rank != src){
		fillBetaFromBuffer(phi, phiBuffer, sendCount);
	}
	freePhiBuffer(phiBuffer);

	return;
}



// initialises the counters by randomising topics to every word
void	initialiseCountersMPI(int locNoFiles, int *locUniqueCount, int ***localDocMatrix, int noTopics, int ***topicIndex, int **topicCountPerDoc, int *topicCountTotal, int *localTopicCountTotal, int **wordCountPerTopic, int **localWordCountPerTopic){
	int d, w, occ;
	int word, topic;

	// loops through all occurrences of every word in every document
	for(d = 0; d < locNoFiles; d++){
		for(w = 0; w < locUniqueCount[d]; w++){
			word = localDocMatrix[d][w][0];
			for(occ = 0; occ < localDocMatrix[d][w][1]; occ++){
					
				// randomly assigns a topic and fills in all the counters accordingly 
				topic = (int) (drand48() * noTopics);
				topicIndex[d][w][occ] = topic;
				topicCountPerDoc[d][topic] ++;

				topicCountTotal[topic] ++;
				localTopicCountTotal[topic] ++;

				wordCountPerTopic[word][topic] ++;
				localWordCountPerTopic[word][topic] ++;
			}
		}
	}

	return;

}


// runs gibbs sampling for AD-LDA
void	gibbsSamplerMPI(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int ***topicIndex, int noTopics, int noFiles, int ***docMatrix, int *uniqueWordCounts, double *alphaDir, double *betaDir, double betaSum, int iterations, int **localWordCountPerTopic, int *localTopicCountTotal, int vocabSize, MPI_Comm mainComm){
	int rank;
	MPI_Comm_rank(mainComm, &rank);

	int i, d, w, occ, k;
	int wordIndex, topic;


	// prob dist
	double *probDist = malloc(noTopics * sizeof(double));		// freed
	double probTotal, sampler;
	

	// Each iteration every occurrence of every word in every document is looped
	for(i = 0; i < iterations; i++){	
		for(d = 0; d < noFiles; d++){
			for(w = 0; w < uniqueWordCounts[d]; w++){
				wordIndex = docMatrix[d][w][0];
				for(occ = 0; occ < docMatrix[d][w][1]; occ++){

					// finds current topic
					topic = topicIndex[d][w][occ];	
					
					// get -i values
					topicCountPerDoc[d][topic] --;
					wordCountPerTopic[wordIndex][topic] --;
					topicCountTotal[topic] --;
					localWordCountPerTopic[wordIndex][topic] --;
					localTopicCountTotal[topic] --;

					// update probs
					probTotal = 0;
					for(k = 0; k < noTopics; k++){
						probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (wordCountPerTopic[wordIndex][k] + betaDir[wordIndex]) / (topicCountTotal[k] + betaSum);
						probTotal += probDist[k];
					}
		
					// sample new topic
					sampler = drand48() * probTotal;
					probTotal = probDist[0];
					topic = 0;
					while(probTotal < sampler && topic < noTopics - 1){
						topic ++;
						probTotal += probDist[topic];
					}

					// update z
					topicIndex[d][w][occ] = topic;

					// update counters
					topicCountPerDoc[d][topic] ++;
					wordCountPerTopic[wordIndex][topic] ++;
					topicCountTotal[topic] ++;
					localWordCountPerTopic[wordIndex][topic] ++;
					localTopicCountTotal[topic] ++;
				}
			}
		}
		// After iteration of all words Cword and Ctopic are reduced to everyone to update their global counters
		MPI_Allreduce(localWordCountPerTopic[0], wordCountPerTopic[0], vocabSize*noTopics, MPI_INT, MPI_SUM, mainComm);
		MPI_Allreduce(localTopicCountTotal, topicCountTotal, noTopics, MPI_INT, MPI_SUM, mainComm);
	}



	free(probDist);
					
	return;
}
