#include "pthreadFunctions.h"



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

	double *probDist = malloc(noTopics * sizeof(double));		// FREED
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);			// FREED

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

	freeMatInt(CwordChanges);
	free(CtopicBuffer);
	
	
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



void	*recvFetchCtopic(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
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
	int		noTopics = (*pargs) -> noTopics;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**wordCountPerTopic = (*pargs) -> wordCountPerTopic;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;

	int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
	int **CwordBuffer;
	allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED
	
	int i, j;
	int src;
	int noFinished = 0;

	int counter = 0;

	while(noFinished < pdSize){
		counter ++;
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

	free(wordIndexes);
	freeMatInt(CwordBuffer);


	return NULL;
}



void	*recvUpdateCword(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int		noTopics = (*pargs) -> noTopics;
	int		bundleSize = (*pargs) -> bundleSize;
	int		**wordCountPerTopic = (*pargs) -> wordCountPerTopic;
	int 	*localTopicCountTotalAct = (*pargs) -> localTopicCountTotalAct;
	int		pdSize = (*pargs) -> pdSize;
	MPI_Comm mainComm = (*pargs) -> mainComm;

	int i, j;
	int src;
	MPI_Status stat;


	int *wordIndexes = malloc(bundleSize * sizeof(int));		// FREED
	int **CwordBuffer;
	allocMatInt(&CwordBuffer, bundleSize, noTopics);			// FREED
	int *CtopicBuffer = malloc(noTopics * sizeof(int));			// FREED
	int noFinished = 0;

	int counter = 0;
	while(noFinished < pdSize){
		counter ++;
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

	free(wordIndexes);
	freeMatInt(CwordBuffer);
	free(CtopicBuffer);

	
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
	int		*topicSignal = (*pargs) -> topicSignal;
	pthread_mutex_t *locks = (*pargs) -> locks;	
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;


	int 	*localTopicCountTotal, *holder;
	callocArrInt(&localTopicCountTotal, noTopics);				// FREED
	callocArrInt(&holder, noTopics);							// FREED




	int i, j;
	int topicCounter = 0;
	while((*stop) == 0){
		// fetch Ctopic
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

	double *probDist = malloc(noTopics * sizeof(double));				// FREED
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);					// FREED

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

	freeMatInt(CwordChanges);
	free(probDist);
	
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

	int 	queueIndex = 0;
	int		bun, dest;
	int		stopCount = 0;
	int		currentIter = 0;
	int		currentBun = queueSize - 1;
	int		skippedCount = 0;
	pthread_mutex_t stopLock;
	pthread_mutex_init(&stopLock, NULL);			// FREED
	while((*stop) == 0){
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
								pthread_mutex_destroy(&stopLock);
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
						pthread_mutex_destroy(&stopLock);
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

	
	pthread_mutex_destroy(&stopLock);
	return NULL;
}



void	*fetchCtopicDeadlineNew(void *arg){
	struct argStruct **pargs = (struct argStruct **)arg;
	int 	*stop = (*pargs) -> stop;
	int 	*topicCountTotal = (*pargs) -> topicCountTotal;
	int		noTopics = (*pargs) -> noTopics;
	int		pdSize = (*pargs) -> pdSize;
	int		size = (*pargs) -> size;
	int		rank = (*pargs) -> rank;
	int		*topicSignal = (*pargs) -> topicSignal;
	MPI_Comm mainComm = (*pargs) -> mainComm;
	MPI_Status stat;


	int 	*localTopicCountTotal, *holder;
	callocArrInt(&localTopicCountTotal, noTopics);		// FREED
	callocArrInt(&holder, noTopics);					// FREED




	int i, j;
	int topicCounter = 0;
	int counter = 0;
	while((*stop) == 0){
		counter ++;
		// fetch Ctopic
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

	double *probDist = malloc(noTopics * sizeof(double));	// FREED
	double probTotal, sampler;

	int **CwordChanges;
	callocMatInt(&CwordChanges, bundleSize, noTopics);		// FREED

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

	free(probDist);
	freeMatInt(CwordChanges);
	
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
		// update Cword
		if(queueMat[queueIndex][0] == 2){
					bun = queueMat[queueIndex][1];
					dest = bundles2d[bun][bundleSize] + pdSize;

					// request and read in the values
					//printf("sending %d\n", bun);
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
								pthread_mutex_destroy(&stopLock);
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
						pthread_mutex_destroy(&stopLock);
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

	
	pthread_mutex_destroy(&stopLock);
	return NULL;
}


//# vim: ts=4
