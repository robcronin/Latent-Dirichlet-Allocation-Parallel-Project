#include "gibbsLDA2.h"


// creates an array of words from an alpha struct
void 	createVocab(struct alpha *phi, char ***vocabulary){
	int 	noWords = phi -> uniqueCount;

	// allocates for number of words
	*vocabulary = malloc(noWords * sizeof(char *));
	int i;
	int index = 0;
	struct wordholder *search;

	// loops through the struct and adds each word
	for(i = 0; i < 27; i++){
		search = phi -> letters[i] -> head;
		while(search != NULL){
			(*vocabulary)[index] = malloc(strlen(search -> word) + 2);
			strcpy( (*vocabulary)[index], search -> word);
			index ++;
			search = search -> next;
		}
	}

	return;
}



// frees a vocab array
void	freeVocab(char ***vocabulary, int vocabSize){
	int i;
	for(i = 0; i < vocabSize; i++){
		free((*vocabulary)[i]);
	}
	free(*vocabulary);

	return;
}



// allocates a docMatrix based on unique number of words per files
void	allocDocMatrix(int ****docMatrix, int noFiles, int *uniqueWordCounts){
	int i, j;
	int totalIndex = 0;
	int docIndex;

	// total size is 2 times the sum of unique counts from all docs (stored at noFiles index of uniqueWordCounts)
	int *temp = malloc(uniqueWordCounts[noFiles] * 2 * sizeof(int));

	// allocates pointers for each file
	(*docMatrix) = malloc(noFiles * sizeof(int *));

	// loops though files
	for(i = 0; i < noFiles; i++){
		docIndex = 0;
	
		// allocates pointers to each word in a document	
		(*docMatrix)[i] = malloc(uniqueWordCounts[i] * sizeof(int *));

		// loops through unique number of words to point at correct part of temp
		for(j = 0; j < uniqueWordCounts[i]; j++){
			(*docMatrix)[i][j] = &temp[totalIndex + docIndex];
			docIndex += 2;
		}
		totalIndex += 2*uniqueWordCounts[i];
	}

	return;
}



// frees a doc matrix
void	freeDocMatrix(int ****docMatrix, int noFiles){
	int i;
	
	free( (*docMatrix)[0][0] );
	for(i = 0; i < noFiles; i++){
		free( (*docMatrix)[i] );
	}
	free(*docMatrix);	

	return;
}




void	fillDocMatrix(int ***docMatrix, char **vocabulary, struct alpha **all, int noFiles){
	struct wordholder *search;
	int i, j;

	int dictIndex;
	int docIndex;
	for(i = 0; i < noFiles; i++){
		dictIndex = 0;
		docIndex = 0;
		for(j = 0; j < 27; j++){
			search = all[i] -> letters[j] -> head;
			while(search != NULL){
				while(strcmp(search -> word, vocabulary[dictIndex]) != 0){
					dictIndex ++;
				}
				docMatrix[i][docIndex][0] = dictIndex;
				docMatrix[i][docIndex][1] = search -> count;

				search = search -> next;			
				dictIndex ++;
				docIndex ++;
			}
		}
	}
	return;
}

void	fillDocMatrixNew(int ***docMatrix, char **vocabulary, struct alpha **all, int noFiles, struct alpha *phi){
	struct wordholder *search;
	int i, j;

	int dictIndex;
	int docIndex;
	for(i = 0; i < noFiles; i++){
		dictIndex = 0;
		docIndex = 0;
		for(j = 0; j < 27; j++){
			search = all[i] -> letters[j] -> head;
			while(search != NULL){
				if(is_present(phi, search -> word)){
					while(strcmp(search -> word, vocabulary[dictIndex]) != 0){
						dictIndex ++;
					}
					docMatrix[i][docIndex][0] = dictIndex;
					docMatrix[i][docIndex][1] = search -> count;

					dictIndex ++;
					docIndex ++;
				}
				search = search -> next;			
			}
		}
	}
	return;
}

void	gibbsSampler2(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int ***topicIndex, int noTopics, int noFiles, int ***docMatrix, int *uniqueWordCounts, double *alphaDir, double *gammaDir, double gammaSum, int iterations){

	int i, d, w, occ, k;
	int wordIndex, topic;

	// prob dist
	double *probDist = malloc(noTopics * sizeof(double));		// freed
	double probTotal, sampler;

	for(i = 0; i < iterations; i++){	
		for(d = 0; d < noFiles; d++){
			for(w = 0; w < uniqueWordCounts[d]; w++){
				wordIndex = docMatrix[d][w][0];
				for(occ = 0; occ < docMatrix[d][w][1]; occ++){
					topic = topicIndex[d][w][occ];	
					
					// get -i values
					topicCountPerDoc[d][topic] --;
					wordCountPerTopic[wordIndex][topic] --;
					topicCountTotal[topic] --;

					// update probs
					probTotal = 0;
					for(k = 0; k < noTopics; k++){
						probDist[k] = (topicCountPerDoc[d][k] + alphaDir[k]) * (wordCountPerTopic[wordIndex][k] + gammaDir[wordIndex]) / (topicCountTotal[k] + gammaSum);
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
				//	printf("topic is %d\n", topic);

					// update z
					topicIndex[d][w][occ] = topic;

					// update counters
					topicCountPerDoc[d][topic] ++;
					wordCountPerTopic[wordIndex][topic] ++;
					topicCountTotal[topic] ++;
				}
			}
		}
	}

	free(probDist);
					
	return;
}



void 	printTopXLDA2(double **phiMatrix, int noTopics, int vocabSize, char **vocabulary, int **wordCountPerTopic, int x){
	int i, j;

	// set up matrices	
	double **max;
	int **topWords;
	allocMat(&max, noTopics, x);		// freed
	allocMatInt(&topWords, noTopics, x);	// freed
	for(i = 0; i < noTopics; i++){
		for(j = 0; j < x; j++){
			max[i][j] = 0;
			topWords[i][j] = -1;
		}
	}

	// loop through and find max
	int k, w;
	for(w = 0; w < vocabSize; w++){
		for(k = 0; k < noTopics; k++){
			if(phiMatrix[w][k] > max[k][x-1]){
				max[k][x-1] = phiMatrix[w][k];
				topWords[k][x-1] = w;
				for(i = x-2; i >= 0; i--){
					if(phiMatrix[w][k] > max[k][i]){
						max[k][i+1] = max[k][i];
						topWords[k][i+1] = topWords[k][i];
						max[k][i] = phiMatrix[w][k];
						topWords[k][i] = w;
					}
					else{
						break;
					}
				}
			}
		}
	}

	// print results
	for(k = 0; k < noTopics; k++){
		if(wordCountPerTopic[ topWords[k][0] ][k] > 3){
			printf("Topic %d\n", k+1);
			for(i = 0; i < x; i++){
				printf("%s with index %d with count %d\n",vocabulary[topWords[k][i]], topWords[k][i], wordCountPerTopic[topWords[k][i]][k] );
			}
			printf("\n");
		}
	}

	// free memory
	freeMat(max);
	freeMatInt(topWords);

	return;
}



void	createPhi(double ***phiMatrix, int vocabSize, int noTopics, int **wordCountPerTopic, double *gammaDir){
	int i, k;

	// set up phi
	allocMat(&(*phiMatrix), vocabSize, noTopics);
	double *sum = malloc(noTopics * sizeof(double));

	// reset sum
	for(i = 0; i < noTopics; i++){
		sum[i] = 0;
	}

	// fill initial values
	for(i = 0; i < vocabSize; i++){
		for(k = 0; k < noTopics; k++){
			(*phiMatrix)[i][k] = wordCountPerTopic[i][k] + gammaDir[i];
			sum[k] += (*phiMatrix)[i][k];
		}
	}

	// sum to 1
	for(i = 0; i < vocabSize; i++){
		for(k = 0; k < noTopics; k++){
			(*phiMatrix)[i][k] /= sum[k];
		}
	} 

	// free
	free(sum);

	return;
}

void 	createTopicIndex(int ****topicIndex, int *totalWordCounts, int noFiles, int *uniqueWordCounts, int ***docMatrix){
	int i, j;
	
	// set up matrix
	int *temp = malloc(totalWordCounts[noFiles] * sizeof(int));
	zeroArrInt(temp, totalWordCounts[noFiles]);
	*topicIndex = malloc(noFiles * sizeof(int *));
	
	// set up pointers correctly
	int totalIndex = 0;
	for(i = 0; i < noFiles; i++){
		(*topicIndex)[i] = malloc(totalWordCounts[i] * sizeof(int *));
		for(j = 0; j < uniqueWordCounts[i]; j++){
			(*topicIndex)[i][j] = &temp[totalIndex];
			//if(&(docMatrix[i][j][1]) == NULL){printf("\t\t\tCATASTROPIC\n");}
			totalIndex += docMatrix[i][j][1];
		}
	}

	return;

}

void	freeTopicIndex(int ***topicIndex, int noFiles){
	int i;
	free(topicIndex[0][0]);
	for(i = 0; i < noFiles; i++){
		free(topicIndex[i]);
	}
	free(topicIndex);

	return;
}


void 	createCounters2(int ***topicCountPerDoc, int **topicCountTotal, int ***wordCountPerTopic, int noFiles, int noTopics, int vocabSize){
	
	// CREATE COUNTERS
	allocMatInt(&(*topicCountPerDoc), noFiles, noTopics);		// n_d,k
	*topicCountTotal = malloc(noTopics * sizeof(int));		// n_k
	allocMatInt(&(*wordCountPerTopic), vocabSize, noTopics);		// n_w,k 

	// ZERO COUNTERS
	zeroMatInt(*topicCountPerDoc, noFiles, noTopics);
	zeroArrInt(*topicCountTotal, noTopics);
	zeroMatInt(*wordCountPerTopic, vocabSize, noTopics);

}

void	freeCounters2(int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic){
	freeMatInt(topicCountPerDoc);
	free(topicCountTotal);
	freeMatInt(wordCountPerTopic);

	return;

}


void	createWordCounts(int **uniqueWordCounts, int **totalWordCounts, int noFiles, struct alpha **all){
	int i;

	*uniqueWordCounts = malloc((noFiles+1) * sizeof(int));
	*totalWordCounts = malloc((noFiles+1) * sizeof(int));
	(*uniqueWordCounts)[noFiles] = 0;
	(*totalWordCounts)[noFiles] = 0;
	for(i = 0; i < noFiles; i++){
		(*uniqueWordCounts)[i] = all[i] -> uniqueCount;
		(*uniqueWordCounts)[noFiles] += (*uniqueWordCounts)[i];

		(*totalWordCounts)[i] = all[i] -> total_count;
		(*totalWordCounts)[noFiles] += (*totalWordCounts)[i];
	}

	return;

}


void 	initialiseCounters(int ***topicIndex, int **topicCountPerDoc, int *topicCountTotal, int **wordCountPerTopic, int noFiles, int noTopics, int *uniqueWordCounts, int ***docMatrix){
	int i, j, k;
	int topic, word;


	for(i = 0; i < noFiles; i++){
		for(j = 0; j < uniqueWordCounts[i]; j++){
			word = docMatrix[i][j][0];
			for(k = 0; k < docMatrix[i][j][1]; k++){
				topic = (int) (drand48() * noTopics);
				topicIndex[i][j][k] = topic;
				topicCountPerDoc[i][topic] ++;
				topicCountTotal[topic] ++;
				wordCountPerTopic[word][topic] ++;
			}
		}
	}

	return;

}

void	initialiseAlpha(double **alpha, int noTopics){
	(*alpha) = malloc(noTopics * sizeof(double));
	
	int i;
	for(i = 0; i < noTopics; i++){
		(*alpha)[i] = 1.0;
	}

	return;
}


void	initialiseGamma(double **gamma, int noWords, double *gammaSum){
	(*gamma) = malloc(noWords * sizeof(double));
	*gammaSum = 0;

	int i;
	for(i = 0; i < noWords; i++){
		(*gamma)[i] = 1.0;
		*gammaSum += (*gamma)[i];
	}

	return;
}

void	printInfoLDA(int rank, int index, int LDAmethod, int seed, char *stopfile, int noTopics, int noFiles, int vocabSize, int iterations, int size, int parallelReadIn, int totalWords, int pdSize, int pwSize, int CwordMethod, int rma, int threads, int bundleSize, int invMethod, int sortRead){

	char 	*datestr, *timestr;

	if(rank == 0 && (index == 0 || index == 1)){
		time_t	rawtime;
		struct	tm *timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		char	day[5], num[5], mon[5], year[5], time[20];
		strftime(day, sizeof(day), "%a", timeinfo);
		strftime(num, sizeof(num), "%d", timeinfo);
		strftime(mon, sizeof(mon), "%b", timeinfo);
		strftime(year, sizeof(year), "%Y", timeinfo);
		strftime(time, sizeof(time), "%T", timeinfo);
		datestr = malloc(20 * sizeof(char));
		timestr = malloc(20 * sizeof(char));
		sprintf(datestr, "%s %s %s", num, mon, year);
		sprintf(timestr, "%s", time);
	}

	if(rank == 0 && (index == 0 || index == 1)){
		printf("\n");
		printf("***********************************\n");
		printf("            INFORMATION\n");
		printf("***********************************\n");


		printf("Date:\t\t\t%s\n", datestr);
		printf("Time:\t\t\t%s\n", timestr);
		free(datestr);
		free(timestr);

		if(LDAmethod == 1){printf("LDA Method:\t\tSERIAL\n");}
		else if(LDAmethod == 2){printf("LDA Method:\t\tADLDA/PLDA\n");}
		else if(LDAmethod == 3){printf("LDA Method:\t\tPLDA+\n");}
		else{printf("LDA Method:\t\tUNKNOWN METHOD\n");}

		if(parallelReadIn == 1){printf("Read In Method:\t\tPARALLEL\n");}
		else if (parallelReadIn == 0){printf("Read In Method:\t\tSERIAL\n");}
		if(sortRead == 1){printf("Read in Order:\t\tSORTED\n");}
		else if(sortRead == 0){printf("Read in Order:\t\tRANDOM\n");}
		

		if(LDAmethod == 3){
			if(CwordMethod == 1){
				printf("Cword Method:\t\tBINARY SEARCH TREE\n");
			}
			else{
				printf("Cword Method:\t\tCONDENSED MATRIX\n");
			}
			if(invMethod == 1){
				printf("Inv Method:\t\tBINARY SEARCH TREE\n");
			}
			else{
				printf("Inv Method:\t\tCONDENSED MATRIX\n");
			}
			if(rma == 0){
				printf("RMA Method:\t\tNONE\n");
			}
			else if(rma == 1){
				printf("RMA Method:\t\tINIT TOPICS\n");
			}
			else if(rma == 2){
				printf("RMA Method:\t\tGIBBS SAMPLING\n");
			}
			else if(rma == 3){
				printf("RMA Method:\t\tBOTH\n");
			}
		}

		if(LDAmethod != 1){printf("Processes:\t\t%d\n", size);}
		if(LDAmethod == 3){
			printf("Pd Processes:\t\t%d\n", pdSize);
			printf("Pw Processes:\t\t%d\n", pwSize);
			if(threads == 1){
				printf("Threads Per Process:\t4\n");
			}
			printf("Bundle Size:\t\t%d\n", bundleSize);
		}

		if(stopfile != NULL){printf("Stop Words:\t\tYES\n");}
		else{printf("Stop Words\t\tNO\n");}

		printf("Seed:\t\t\t%d\n", seed);
		printf("Topics:\t\t\t%d\n", noTopics);
		printf("Gibbs Iterations:\t%d\n", iterations);
		printf("Files:\t\t\t%d\n", noFiles);
	}
	if(rank == 0 && (index == 0 || index == 2)){
		printf("Vocab Size:\t\t%d\n", vocabSize);
		printf("Total Words:\t\t%d\n", totalWords);

		printf("***********************************\n");
		printf("\n");
	}

	return;

}

