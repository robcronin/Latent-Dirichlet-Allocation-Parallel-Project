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

void 	printHelp();
void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn);



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
	int	pdSize, pwSize, prank;
	char 	procType;
	int	*wordSpread, *wordLocation;
	struct 	treeRootCword *wordCountPerTopicTree;
	struct	treeRootInv *invertedIndexTree;
	int	**bundles2d, totNoBundles, circleStartBundle;




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
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn);

	// FIND NO FILES & PATH NAMES
	findPathNames(directory, &pathNames, &noFiles, maxFiles);
	free(directory);

	// ASSIGN PW and PD
	procAssignment(rank, size, noFiles, &pdComm, &pwComm, pwToPdRatio, &pdSize, &pwSize, &prank, &procType, mainComm);

	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize);

	// READ IN FILES
	readInDocsAll(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, pdComm);

	// FREE
	freePathNames(pathNames);
	free(stopfile);

	// PRINT INFO
	printInfoLDA(rank, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize);
	if(rank == 0){printTime(t1Sec, t2Sec, "READING IN DOCS");}











	// ************************************************
	// ************************************************
	// ************      m, l & others     ************
	// ************************************************
	// ************************************************
	if(rank == 0){gettimeofday(&t1Sec, NULL);}
	
	// Pd
	if(procType == 'd'){
		// m and l
		callocArrInt(&wordSpread, vocabSize);
		callocArrInt(&wordLocation, vocabSize);

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










	// ************************************************
	// ************************************************
	// *********     DISTRIBUTE WORDS       ***********
	// ************************************************
	// ************************************************
	
	distWordsToPw(procType, vocabSize, pdSize, pwSize, rank, size, prank, wordSpread, wordLocation, noTopics, &wordCountPerTopicTree, pdComm, mainComm);










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
	
	initTopicsPLDAplus(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocation, topicCountPerDoc, pdSize, size, wordCountPerTopicTree, rank, mainComm);





	// ************************************************
	// ************************************************
	// ***********     WORD BUNDLES       *************
	// ************************************************
	// ************************************************
	
	createWordBundles(procType, invertedIndexTree, wordLocation, pwSize, pdSize, bundleSize, &bundles2d, &totNoBundles, prank, pdComm);






	// ************************************************
	// ************************************************
	// ***********     CIRCULAR QUEUE       ***********
	// ************************************************
	// ************************************************
	circleQueue(procType, vocabSize, pdSize, prank, &circleStartBundle, bundles2d);







	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************
	runPLDAplusGibbs(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopicTree, rank, mainComm);

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
	collectAndPrint(rank, vocabSize, noTopics, pwSize, procType, wordCountPerTopicTree, printResults, betaDir, vocabulary, mainComm);





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
	printf("-d [DIR]\tChange Directory (default: ../../Gutenberg)\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs (default: 1000)\n");
	printf("-m [INT]\tSet max no files (default: all files)\n");
	printf("-p\t\tPrint Results\n");
	printf("-s [INT]\tSet seed (default: random)\n");
	printf("-t [INT]\tSet number of Topics (default: 10)\n");

	return;
}



void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn){
	int pArgs[6];
	int endProg = 0;

	//parse command line arguments on rank 0
	if(rank == 0){
		int opt;
		while((opt=getopt(argc,argv,"d:hi:m:prs:t:"))!=-1){
			switch(opt){
				case 'd':
					//*directory = optarg;
					sprintf(*directory, optarg);
					break;
				case 'h':
					printHelp();
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
				default:
					printHelp();
					endProg = 1;
					break;
			}
		}

		pArgs[0] = *seed;
		pArgs[1] = *noTopics;
		pArgs[2] = *iterations;
		pArgs[3] = *maxFiles;
		pArgs[4] = endProg;
		pArgs[5] = *parallelReadIn;
	}

	// Send relevant args to rest of processes
	MPI_Bcast(pArgs, 6, MPI_INT, 0, mainComm);
	*seed = pArgs[0];
	*noTopics = pArgs[1];
	*iterations = pArgs[2];
	*maxFiles = pArgs[3];
	endProg = pArgs[4];
	*parallelReadIn = pArgs[5];
	srand48(78*rank + (*seed) );
	MPI_Bcast(*directory, 1024, MPI_CHAR, 0, mainComm);


	if(endProg){
		MPI_Finalize();
		exit(0);
	}
	
	return;
}

//# vim: ts=4
