#include "gibbsplda+.h"




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
		int **CwordBuffer, **CwordChanges;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
		allocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED
	
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
				zeroMatInt(CwordChanges, bundleSize, noTopics);
				
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
									CwordChanges[w][topic] --;
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
									CwordChanges[w][topic] ++;
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
				MPI_Send( &(CwordChanges[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
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
		freeMatInt(CwordChanges);
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
							current -> topicArray[j] += CwordBuffer[i][j];
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
		int **CwordBuffer, **CwordChanges;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
		allocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED
	
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
				zeroMatInt(CwordChanges, bundleSize, noTopics);
				
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
									CwordChanges[w][topic] --;
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
									CwordChanges[w][topic] ++;
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
				MPI_Send( &(CwordChanges[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
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
		freeMatInt(CwordChanges);
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
		int **CwordBuffer, **CwordChanges;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
		callocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED
	
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
				zeroMatInt(CwordChanges, bundleSize, noTopics);
				
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
									CwordChanges[w][topic] --;
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
									CwordChanges[w][topic] ++;
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
				MPI_Send( &(CwordChanges[0][0]), bundleSize * noTopics, MPI_INT, dest, 6, mainComm);
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
		freeMatInt(CwordChanges);
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
		double *probDist = malloc(noTopics * sizeof(double));		// FREED
		double probTotal, sampler;

		int **CwordBuffer, **CwordChanges;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
		allocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED		
	
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);					// FREED

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
				zeroMatInt(CwordChanges, bundleSize, noTopics);
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
									CwordChanges[w][topic] --;
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
									CwordChanges[w][topic] ++;
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
						MPI_Accumulate(CwordChanges[w], noTopics, MPI_INT, pdSize, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, MPI_SUM, winArr[dest][0]);
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
		
		
		freeMatInt(CwordBuffer);
		freeMatInt(CwordChanges);
		free(CtopicBuffer);
		free(probDist);

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
		double *probDist = malloc(noTopics * sizeof(double));		// FREED
		double probTotal, sampler;

		int **CwordBuffer, **CwordChanges;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);		// FREED
		allocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED
	
		int *CtopicBuffer;
		callocArrInt(&CtopicBuffer, noTopics);					// FREED

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
				zeroMatInt(CwordChanges, bundleSize, noTopics);
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
									CwordChanges[w][topic] --;
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
									CwordChanges[w][topic] ++;
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
						//MPI_Put(CwordBuffer[w], noTopics, MPI_INT, pdSize + dest, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, wordCountWin);
						MPI_Accumulate(CwordChanges[w], noTopics, MPI_INT, pdSize + dest, bundles2dConv[bun][w]*noTopics, noTopics, MPI_INT, MPI_SUM, wordCountWin);
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

		freeMatInt(CwordBuffer);
		freeMatInt(CwordChanges);
		free(CtopicBuffer);
		free(probDist);	

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
	
	endTimeMPI(t1, "GIBBS SAMPLING RMA", rank, mainComm);
	return;
}



void    runPLDAplusGibbsPthreads(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));		// FREED

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);				// FREED
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 2);					// FREED
		int i;
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);	// FREED
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);				// FREED
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


		free(handles);
		free(topicCountTotal);
		freeMatInt(queueMat);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);

		MPI_Barrier(pdComm);
		for(i = pdSize; i < size; i++){
			MPI_Send(&i, 1, MPI_INT, i, 0, mainComm);
		}
	}
	else{
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED

		int *CtopicBuffer = malloc(noTopics * sizeof(int));			// FREED
	
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

		free(wordIndexes);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS", rank, mainComm);
	return;
}



void    runPLDAplusGibbsPthreadsBoth(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);

	if(procType == 'd'){
		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));					// FREED

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
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);		// FREED
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);					// FREED
		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));			// FREED
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
		callocArrInt(&localTopicCountTotal, noTopics);						// FREED
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
		
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
		free(localTopicCountTotal);
		free(handles);
		free(pargs);
	}
	else{
		// set up threads
		pthread_t *handles;
		int numThreads = 3;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));			// FREED

		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));	// FREED
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
		
		free(handles);
		free(pargs);


		/*
		MPI_Status stat;
		int noFinished = 0;
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED

		int *CtopicBuffer = malloc(noTopics * sizeof(int));			// FREED
	
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
	
		free(wordIndexes);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS BOTH", rank, mainComm);
	return;
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
		handles = malloc(numThreads * sizeof(pthread_t));		// FREED

		// XXX

		// deadlines
		int		sleepTime = 5000;
		int		*deadlines = malloc(3 * sizeof(int));			// FREED
		deadlines[0] = 1500000;
		deadlines[1] = 10000;
		deadlines[2] = 10000;
		

		int pwSize = size - pdSize; 
		int *workload;
		callocArrInt(&workload, pwSize);						// FREED
		int overload = 1.5*(queueSize/pwSize);


		// locks
		pthread_mutex_t *locks = malloc(queueSize * sizeof(pthread_mutex_t));		// FREED
		for(i = 0; i < queueSize; i++){
			pthread_mutex_init(&(locks[i]), NULL);									// FREED
		}

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);			// FREED
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 3);				// FREED
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
			queueMat[i][2] = deadlines[0];
			workload[ bundles2d[i][bundleSize] ] ++;	
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);	// FREED
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);				// FREED
		int stop = 0;
		int topicSignal = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));		// FREED
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

		for(i = 0; i < queueSize; i++){
			pthread_mutex_destroy(&(locks[i]));
		}
		free(locks);
		free(handles);
		free(deadlines);
		free(workload);
		free(topicCountTotal);
		freeMatInt(queueMat);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
		free(pargs);
	
	}
	else{
/*
		// set up threads
		pthread_t *handles;
		int numThreads = 3;
		int queueSize = 10;
		handles = malloc(numThreads * sizeof(pthread_t));			// FREED

		int stop = 0;
		
		// fill argStruct
		struct argStruct *pargs = malloc(sizeof(struct argStruct));	// FREED
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

		free(handles);
		free(pargs);
		*/
//		/*

		MPI_Status stat;
		int noFinished = 0;	// FIXME
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED

		int *CtopicBuffer = malloc(noTopics * sizeof(int));			// FREED
	
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
		free(wordIndexes);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
//		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS w/ deadline", rank, mainComm);
	return;
}



void    runPLDAplusGibbsPthreadsDeadlineNew(char procType, int noTopics, int pdSize, int size, int bundleSize, int iterations, int totNoBundles, int **bundles2d, int **bundles2dConv, struct treeRootInv *invertedIndexTree, int locNoFiles, int ***localDocMatrix, int ***topicIndex, int **topicCountPerDoc, double *alphaDir, double *betaDir, double betaSum, int **wordCountPerTopic, int *localTopicCountTotalAct, int rank, int **invertedIndexMatrix, int **wordLocationMat, MPI_Comm pdComm, MPI_Comm mainComm){


	struct	timeval t1;
	startTimeMPI(&t1, mainComm);


	int queueSize = 100;

	if(procType == 'd'){
		int i;

		// set up threads
		pthread_t *handles;
		int numThreads = 4;
		int loopsPerSample = 100;
		handles = malloc(numThreads * sizeof(pthread_t));		// FREED

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
		pthread_mutex_t *locks = malloc(queueSize * sizeof(pthread_mutex_t));	// FREED
		for(i = 0; i < queueSize; i++){
			pthread_mutex_init(&(locks[i]), NULL);								// FREED
		}

		// variables
		int 	*topicCountTotal;
		callocArrInt(&topicCountTotal, noTopics);		// FREED
		int 	**queueMat;
		callocMatInt(&queueMat, queueSize, 3);			// FREED
		for(i = 0; i < queueSize; i++){
			queueMat[i][0] = 0;
			queueMat[i][1] = i;
			queueMat[i][2] = deadlines[0];
			workload[ bundles2d[i][bundleSize] ] ++;	
		}
		int 	**CwordBuffer;
		allocMatInt(&CwordBuffer, queueSize * bundleSize, noTopics);		// FREED
		int		**CtopicBuffer;
		callocMatInt(&CtopicBuffer, queueSize, noTopics);					// FREED
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
		callocArrInt(&localTopicCountTotal, noTopics);		// FREED
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

		for(i = 0; i < queueSize; i++){
			pthread_mutex_destroy(&(locks[i]));
		}
		free(locks);
		free(handles);
		free(topicCountTotal);
		freeMatInt(queueMat);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
		free(pargs);
		free(localTopicCountTotal);
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
		free(handles);
		free(pargs);
		*/
//		/*

		MPI_Status stat;
		int noFinished = 0;	// FIXME
		int i, j, src;

		int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
		int **CwordBuffer;
		allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED

		int *CtopicBuffer = malloc(noTopics * sizeof(int));			// FREED
	
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

		free(wordIndexes);
		freeMatInt(CwordBuffer);
		free(CtopicBuffer);
//		*/
	}

	endTimeMPI(t1, "GIBBS SAMPLING PTHREADS w/ deadline NEW", rank, mainComm);
	return;
}
