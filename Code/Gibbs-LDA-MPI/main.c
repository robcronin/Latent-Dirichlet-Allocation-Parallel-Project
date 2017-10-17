#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "gibbsLDA-MPI.h"

void 	printHelp();
void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn);




int 	main(int argc, char *argv[]){


	// ************************************************
	// ************************************************
	// ***********         SETUP           ************
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
	int 	LDAmethod = 2;
	int	sortRead = 1;
	int	gutenFiles = 1;

	struct 	timeval t1, t2, t1Sec, t2Sec;
	int 	noFiles, vocabSize, totalWords;
	int 	locNoFiles, *localFiles;
	char	**vocabulary;
	int	**topicCountPerDoc, *topicCountTotal, **wordCountPerTopic, ***topicIndex;
	int 	**localWordCountPerTopic, *localTopicCountTotal;
	double 	*alphaDir, *gammaDir, gammaSum;
	int 	*locUniqueCount, *locTotalCount;
	int 	***localDocMatrix;
	double 	**phiMatrix;
	char 	**pathNames;


	// MPI SETUP
	int 	rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm mainComm = MPI_COMM_WORLD;
	MPI_Comm_rank(mainComm, &rank);
	MPI_Comm_size(mainComm, &size);


	// START TIMING
	if(rank == 0)
		gettimeofday(&t1, NULL);

	// PARSE ARGS
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn);



	// ************************************************
	// ************************************************
	// ***********     READ IN FILES       ************
	// ************************************************
	// ************************************************



	if(rank == 0){gettimeofday(&t1Sec, NULL);}

	// FIND NO FILES
	findNoFiles(directory, &noFiles, maxFiles);

	// FIND PATH NAMES
	if(sortRead == 1){
		findPathNamesSorted(directory, &pathNames, &noFiles, maxFiles, size);
	}
	else{
		findPathNames(directory, &pathNames, &noFiles, maxFiles);
	}
	free(directory);
	
	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, 0, 0, 0, 0, 0, 0, 0, sortRead);

	// ALTER COMMUNICATOR IF NECESSARY
	alterComm(rank, &size, noFiles, &mainComm);

	// DETERMINE NUMBER OF FILES PER RANK
	filesPerProc(&localFiles, rank, size, noFiles, &locNoFiles, mainComm);


	if(parallelReadIn){
		readInDocsParallel(rank, size, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, localFiles, locNoFiles, &totalWords, gutenFiles, mainComm);

	}
	else{
		// READ IN DOCS ON RANK 0 AND DISTRIBUTE INFO
		readInDocsSerial(rank, size, pathNames, noFiles, stopfile, &vocabSize, &vocabulary, printResults, &locUniqueCount, &locTotalCount, locNoFiles, localFiles, &localDocMatrix, &totalWords, gutenFiles, mainComm);

	}

	freePathNames(pathNames);
	free(localFiles);

	if(rank == 0){gettimeofday(&t2Sec, NULL);}

	// PRINT INFO
	printInfoLDA(rank, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, 0, 0, 0, 0, 0, 0, 0, sortRead);
	if(rank == 0){printTime(t1Sec, t2Sec, "Reading In Docs");}
	free(stopfile);









	// ************************************************
	// ************************************************
	// *********    SET UP COUNTERS ETC      **********
	// ************************************************
	// ************************************************
	if(rank == 0){gettimeofday(&t1Sec, NULL);}


	

	// ALPHA & BETA
	initialiseAlpha(&alphaDir, noTopics);		// freed
	initialiseGamma(&gammaDir, vocabSize, &gammaSum);	// freed

	// SET UP COUNTERS
	callocMatInt(&topicCountPerDoc, locNoFiles, noTopics);		// local C_doc	-- freed
	createTopicIndex(&topicIndex, locTotalCount, locNoFiles, locUniqueCount, localDocMatrix); // local z  -- freed
	free(locTotalCount);
	callocMatInt(&wordCountPerTopic, vocabSize, noTopics);		// global C_word  -- freed
	callocMatInt(&localWordCountPerTopic, vocabSize, noTopics);	// local C_word  -- freed
	callocArrInt(&topicCountTotal, noTopics);			// global C_topic  -- freed
	callocArrInt(&localTopicCountTotal, noTopics);			// local C_topic  -- freed

	// INITIALISE COUNTERS
	initialiseCountersMPI(locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, topicCountPerDoc, topicCountTotal, localTopicCountTotal, wordCountPerTopic, localWordCountPerTopic);	// freed


	// INITAL ALL REDUCE OP
	MPI_Allreduce(localWordCountPerTopic[0], wordCountPerTopic[0], vocabSize*noTopics, MPI_INT, MPI_SUM, mainComm);
	MPI_Allreduce(localTopicCountTotal, topicCountTotal, noTopics, MPI_INT, MPI_SUM, mainComm);

	
	if(rank == 0){gettimeofday(&t2Sec, NULL);}
	if(rank == 0){printTime(t1Sec, t2Sec, "Setting Up");}






	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************




	if(rank == 0){gettimeofday(&t1Sec, NULL);}
	// RUN PROGRAM
	gibbsSamplerMPI(topicCountPerDoc, topicCountTotal, wordCountPerTopic, topicIndex, noTopics, locNoFiles, localDocMatrix, locUniqueCount, alphaDir, gammaDir, gammaSum, iterations, localWordCountPerTopic, localTopicCountTotal, vocabSize, mainComm);
	if(rank == 0){gettimeofday(&t2Sec, NULL);}
	if(rank == 0){printTime(t1Sec, t2Sec, "Gibbs Sampling");}

	free(alphaDir);
	freeTopicIndex(topicIndex, locNoFiles);
	freeDocMatrix(&localDocMatrix, locNoFiles);
	free(locUniqueCount);
	if(mainComm != MPI_COMM_WORLD){MPI_Comm_free(&mainComm);}

	
	// CREATE PHI
	if(rank == 0)

	// PRINT RESULTS
	if(rank == 0){gettimeofday(&t1Sec, NULL);}
	if(printResults && rank == 0){
		createPhi(&phiMatrix, vocabSize, noTopics, wordCountPerTopic, gammaDir);		// freed
		printTopXLDA2(phiMatrix, noTopics, vocabSize, vocabulary, wordCountPerTopic, 10);
		freeVocab(&vocabulary, vocabSize);
		freeMat(phiMatrix);
		if(rank == 0){gettimeofday(&t2Sec, NULL);}
		if(rank == 0){printTime(t1Sec, t2Sec, "PRINTING");}
	}

	// FREE REMAINING
	free(gammaDir);
	freeCounters2(topicCountPerDoc, topicCountTotal, wordCountPerTopic);
	freeMatInt(localWordCountPerTopic);
	free(localTopicCountTotal);


	// FINISH TIMING
	if(rank == 0){gettimeofday(&t2, NULL);}
	if(rank == 0){printTime(t1, t2, "Whole Program");}






	MPI_Finalize();
	return 0;

}


void printHelp(){
	printf("-d [DIR]\tChange Directory (default: ../../Gutenberg)\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs (default: 1000)\n");
	printf("-m [INT]\tSet max no files (default: all files)\n");
	printf("-p\t\tPrint Results\n");
	printf("-r\t\tSet Parallel Read In\n");
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
