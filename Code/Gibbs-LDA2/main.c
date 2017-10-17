#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "gibbsLDA2.h"
#include "../general/extra.h"

void 	printHelp();
void 	parseArgs(int argc, char *argv[], char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations);

int 	main(int argc, char *argv[]){
	
	// CREATE VARIABLES AND INITIAL VALUES
	int 	seed = time(NULL);
	int 	noTopics = 10;
	char 	*directory = malloc(1024);
	sprintf(directory, "../../Inputs/Gutenberg");		// freed
	char	*stopfile = malloc(1024);			// freed
	sprintf(stopfile, "../../Inputs/stopwords.txt");
	int	printResults = 1;
	int 	maxFiles = 0;
	int	iterations = 100;
	int 	parallelReadIn = -1;
	int 	LDAmethod = 1;
	int	gutenFiles = 1;


	int 	size = 1;
	int 	noFiles, vocabSize = 0, totalWords = 0;
	struct 	alpha *phi, **all;
	int	**topicCountPerDoc, *topicCountTotal, **wordCountPerTopic, ***topicIndex;
	double 	*alphaDir, *gammaDir;
	double	gammaSum;
	char	**vocabulary;
	int 	*uniqueWordCounts, *totalWordCounts;
	int 	***docMatrix;
	double 	**phiMatrix;
	char	**pathNames;
	struct 	timeval t1, t2, t1Sec, t2Sec;



	gettimeofday(&t1, NULL);
	parseArgs(argc, argv, &directory, &maxFiles, &seed, &noTopics, &iterations);
	srand48(seed);




	gettimeofday(&t1Sec, NULL);
	// SETUP CORPUS VOCABULARY
	findPathNames(directory, &pathNames, &noFiles, maxFiles);
	printInfoLDA(0, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, 0, 0, 0, 0, 0, 0, 0, 0);
	free(directory);	
	initialiseLocalTFIDF(pathNames, &all, noFiles, stopfile, 0, gutenFiles);
	freePathNames(pathNames);
	create_dict(&phi);		// freed
	createBeta(all, noFiles, phi);	// phi not beta	-- freed
	vocabSize = phi -> uniqueCount;
	createVocab(phi, &vocabulary);	// freed
	freeBeta(&phi);



	// SET UP DOC MATRIX
	createWordCounts(&uniqueWordCounts, &totalWordCounts, noFiles, all);
	totalWords = totalWordCounts[noFiles];
	allocDocMatrix(&docMatrix, noFiles, uniqueWordCounts);
	fillDocMatrix(docMatrix, vocabulary, all, noFiles);

	if(!printResults){freeVocab(&vocabulary, vocabSize);}
	freeLocalTFIDF(&all, noFiles);
	
	gettimeofday(&t2Sec, NULL);
	printInfoLDA(0, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, 0, 0, 0, 0, 0, 0, 0, 0);
	free(stopfile);
	printTime(t1Sec, t2Sec, "Reading In Docs");


	gettimeofday(&t1Sec, NULL);
	// SET UP COUNTERS
	createCounters2(&topicCountPerDoc, &topicCountTotal, &wordCountPerTopic, noFiles, noTopics, vocabSize);		// freed
	createTopicIndex(&topicIndex, totalWordCounts, noFiles, uniqueWordCounts, docMatrix);		// freed
	initialiseCounters(topicIndex, topicCountPerDoc, topicCountTotal, wordCountPerTopic, noFiles, noTopics, uniqueWordCounts, docMatrix);

	free(totalWordCounts);


	// ALPHA & BETA
	initialiseAlpha(&alphaDir, noTopics);
	initialiseGamma(&gammaDir, uniqueWordCounts[noFiles], &gammaSum);

	gettimeofday(&t2Sec, NULL);
	printTime(t1Sec, t2Sec, "Setting Up");



	gettimeofday(&t1Sec, NULL);
	// RUN PROGRAM
	gibbsSampler2(topicCountPerDoc, topicCountTotal, wordCountPerTopic, topicIndex, noTopics, noFiles, docMatrix, uniqueWordCounts, alphaDir, gammaDir, gammaSum, iterations);

	gettimeofday(&t2Sec, NULL);
	printTime(t1Sec, t2Sec, "Gibbs Sampling");

	// FREE MEMORY
	freeDocMatrix(&docMatrix, noFiles);
	free(uniqueWordCounts);
	freeMatInt(topicCountPerDoc);
	free(topicCountTotal);
	freeTopicIndex(topicIndex, noFiles);
	free(alphaDir);
	

	if(printResults){
		gettimeofday(&t1Sec, NULL);
		// UPDATE PHI
		createPhi(&phiMatrix, vocabSize, noTopics, wordCountPerTopic, gammaDir);

		// PRINT RESULTS
		printTopXLDA2(phiMatrix, noTopics, vocabSize, vocabulary, wordCountPerTopic, 10);

		gettimeofday(&t2Sec, NULL);
		printTime(t1Sec, t2Sec, "Printing");

		freeMat(phiMatrix);
		freeVocab(&vocabulary, vocabSize);

	}


	gettimeofday(&t2, NULL);
	printTime(t1, t2, "Whole Program");


	// FREE MEMORY
	freeMatInt(wordCountPerTopic);
	free(gammaDir);

	return 0;
}


void printHelp(){
	printf("-d [DIR]\tChange Directory (default: ../testfiles/smallTestDocs)\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs (default: 1000)\n");
	printf("-m [INT]\tSet max no files (default: all files)\n");
	printf("-s [INT]\tSet seed (default: random)\n");
	printf("-t [INT]\tSet number of Topics (default: 10)\n");

	return;
}


void parseArgs(int argc, char *argv[], char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations){

	//parse command line arguments
	int opt;
	while((opt=getopt(argc,argv,"d:hi:m:s:t:"))!=-1){
		switch(opt){
			case 'd':
				sprintf(*directory, optarg);
				break;
			case 'h':
				printHelp();
				exit(1);
				break;
			case 'i':
				*iterations = atoi(optarg);
				break;
			case 'm':
				*maxFiles = atoi(optarg);
				break;
			case 's':
				*seed = atoi(optarg);
				break;
			case 't':
				*noTopics = atoi(optarg);
				break;
			default:
				printHelp();
				exit(EXIT_FAILURE);
		}
	}
	return;
}
