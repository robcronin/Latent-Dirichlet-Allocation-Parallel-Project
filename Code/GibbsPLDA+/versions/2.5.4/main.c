#include "plda+.h"
void	printHelp();
void 	parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma, int *bundleSize, double *pwToPdRatio, int *invMethod, int *sortRead);




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
	char	*dictfile = malloc(1024);
	//sprintf(dictfile, "../../Inputs/unixDict.txt");
	dictfile = NULL;
	int		printResults = 0;
	int 	maxFiles = 0;
	int		iterations = 100;
	int 	parallelReadIn = 0;
	int 	LDAmethod = 3;
	double	pwToPdRatio = 0.6;
	int		bundleSize = 4;
	int		gutenFiles = 1;

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
	int		invMethod = 2;		// 1 for bst, 2 for matrix
	int		sortRead = 1;
	int		*localTopicCountTotal;
	int		rma = 0;
	int		wordLocCols;
	//MPI_Win	**winArr;
	//MPI_Comm	*commArr;
	int		CwordSize = 0;
	MPI_Win	wordCountWin, topicCountWin;
	int 	threads = 0;

	// need extra column in wordLocationMat if invMethod is 2
	if(invMethod == 1){
		wordLocCols = 2;
	}
	else{
		wordLocCols = 3;
	}
	
	#ifdef THREADS
	threads = 1;
	#endif


	// ************************************************
	// ************************************************
	// ***********       SETUP         ****************
	// ************************************************
	// ************************************************
	
	// MPI SETUP
	
	if(threads == 1){
		int requested = MPI_THREAD_MULTIPLE, provided;
		MPI_Init_thread(&argc, &argv, requested, &provided);
		if(provided != MPI_THREAD_MULTIPLE){
			printf("*** ERROR ***\nMPI_THREAD_MULTIPLE not supported\n\n");
			MPI_Finalize();
			return 0;
		}
	}		
	else{
		MPI_Init(&argc, &argv);
	}
	MPI_Comm mainComm = MPI_COMM_WORLD;
	MPI_Comm_rank(mainComm, &rank);
	MPI_Comm_size(mainComm, &size);

	// START TIMING
	startTimeMPI(&t1, mainComm);

	// PARSE ARGS
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn, &CwordMethod, &rma, &bundleSize, &pwToPdRatio, &invMethod, &sortRead);


	// FIND NO FILES
	findNoFiles(directory, &noFiles, maxFiles);

	// ASSIGN PW and PD
	procAssignment(rank, size, noFiles, &pdComm, &pwComm, pwToPdRatio, &pdSize, &pwSize, &prank, &procType, mainComm);

	// FIND PATH NAMES
	if(sortRead == 1){
		findPathNamesSorted(directory, &pathNames, &noFiles, maxFiles, pdSize);
	}
	else{
		findPathNames(directory, &pathNames, &noFiles, maxFiles);
	}
	free(directory);


	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma, threads, bundleSize, invMethod, sortRead);

	// READ IN FILES
	//readInDocsAll(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, gutenFiles, pdComm);

	// SLOWER
	readInDocsAllNew(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, gutenFiles, dictfile, pdComm);
	// FREE
	freePathNames(pathNames);
	free(stopfile);
	free(dictfile);

	// PRINT INFO
	printInfoLDA(rank, 2, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma, threads, bundleSize, invMethod, sortRead);
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
		callocArrInt(&wordSpread, vocabSize);							// FREED

		if(CwordMethod == 1){
			callocArrInt(&wordLocation, vocabSize);						// FREED
		}
		else{
			callocMatInt(&wordLocationMat, vocabSize, wordLocCols);		// FREED
		}

		// ALPHA & BETA
		initialiseAlpha(&alphaDir, noTopics);							// FREED
		initialiseBeta(&betaDir, vocabSize, &betaSum);					// FREED

		// SET UP COUNTERS
		callocMatInt(&topicCountPerDoc, locNoFiles, noTopics);			// local C_doc	-- FREED
		createTopicIndex(&topicIndex, locTotalCount, locNoFiles, locUniqueCount, localDocMatrix); // local z  -- FREED
		free(locTotalCount);

		// FILL M
		fillM(procType, vocabSize, locNoFiles, locUniqueCount, localDocMatrix, wordSpread, pdComm);
	}
	else if(procType == 'w'){
		callocArrInt(&localTopicCountTotal, noTopics);					// FREED
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
		distWordsToPwMatrix(procType, vocabSize, pdSize, pwSize, rank, size, prank, wordSpread, wordLocationMat, noTopics, &wordCountPerTopic, &CwordSize, wordLocCols, pdComm, mainComm);
	}

	if(procType == 'd'){free(wordSpread);}

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

	int **invertedIndexMatrix;
	int	invertedIndexSize;
	if(invMethod == 1 || CwordMethod == 1){
		buildInvertedIndex(procType, &invertedIndexTree, locNoFiles, locUniqueCount, localDocMatrix, prank, pdComm);
	}
	else{
		buildInvertedIndexMatrix(procType, &invertedIndexMatrix, locNoFiles, locUniqueCount, localDocMatrix, prank, vocabSize, wordLocationMat, &invertedIndexSize, pdComm);
	}
	


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
		if(procType == 'd'){free(wordLocation);}
	}
	else if(CwordMethod == 2){
		if(invMethod == 1){
			createWordBundlesMatrix(procType, invertedIndexTree, wordLocationMat, pwSize, pdSize, bundleSize, &bundles2d, &bundles2dConv, &totNoBundles, prank, pdComm);
		}
		else{
			createWordBundlesMatrixInv(procType, invertedIndexMatrix, wordLocationMat, pwSize, pdSize, bundleSize, &bundles2d, &bundles2dConv, &totNoBundles, prank, invertedIndexSize, locNoFiles, pdComm);
		}
	}









	// ************************************************
	// ************************************************
	// ***********     CIRCULAR QUEUE       ***********
	// ************************************************
	// ************************************************
	circleQueue(procType, vocabSize, pdSize, prank, &circleStartBundle, bundles2d);




	// WORD BUNDLES MAY NEED TO BE ALTERED TO ACCOUNT FOR NEW matIndex


	// To avoid cross contamination of sends a barrier is thrown before sampling
	MPI_Barrier(mainComm);

	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************
	if(CwordMethod == 1){
		runPLDAplusGibbs(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopicTree, rank, mainComm);
	}
	else if(CwordMethod == 2){
		if(threads == 1){
			//runPLDAplusGibbsPthreads(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pdComm, mainComm);
			//runPLDAplusGibbsPthreadsBoth(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, pdComm, mainComm);
			//runPLDAplusGibbsPthreadsDeadline(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
			runPLDAplusGibbsPthreadsDeadlineNew(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
		}
		else if(invMethod == 2){
			runPLDAplusGibbsMatrixInv(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexMatrix, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, betaDir, betaSum, wordCountPerTopic, localTopicCountTotal, rank, wordLocationMat, mainComm);
		}
		else if(rma == 0 || rma == 1){
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
		freeMatInt(topicCountPerDoc);
	}
	else if(procType == 'w'){
		free(localTopicCountTotal);
	}




	// ************************************************
	// ************************************************
	// *******    COLLECT & PRINT RESULTS     *********     
	// ************************************************
	// ************************************************
	
	MPI_Barrier(mainComm);
	// Barrier to avoid interfering with Gibbs sampling requests etc
	if(CwordMethod == 1){
		collectAndPrint(rank, vocabSize, noTopics, pwSize, procType, wordCountPerTopicTree, printResults, betaDir, vocabulary, mainComm);
	}
	else if(CwordMethod == 2){
		collectAndPrintMatrix(rank, vocabSize, noTopics, procType, wordCountPerTopic, printResults, betaDir, vocabulary, wordLocationMat, pdSize, size, mainComm);
		if(procType == 'd'){freeMatInt(wordLocationMat);}
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
	printf("-a [DOUBLE]\tSet Pw:Pd ratio\t\t\t\t\tDefault: 0.6\n");
	printf("-b [INT]\tSet Bundle Size\t\t\t\t\tDefault: 4\n");
	printf("-c [INT]\tSet Cword Method: 1=bst, 2=matrix\t\tDefault: 2\n");
	printf("-d [DIR]\tSet Directory\t\t\t\t\tDefault: ../../Gutenberg\n");
	printf("-h\t\tPrint this help\n");
	printf("-i [INT]\tSet number of Iterations for Gibbs\t\tDefault: 100\n");
	printf("-m [INT]\tSet max no files\t\t\t\tDefault: All Files\n");
	printf("-o\t\tSet random file read in\t\t\t\tDefault: Balances Load\n");
	printf("-p\t\tPrint Results\n");
	printf("-r\t\tSet Parallel Read In\t\t\t\tDefault: Serial & Bcast\n");
	printf("-s [INT]\tSet seed\t\t\t\t\tDefault: Random\n");
	printf("-t [INT]\tSet number of Topics\t\t\t\tDefault: 10\n");
	printf("-v [INT]\tSet Inv Method: 1=bst, 2=matrix\t\t\tDefault: 2 (Must use Cword Method 2)\n");
	printf("-w [INT]\tSet RMA Method: 0=off, 1=init, 2=gibbs, 3=both\tDefault: 0 (Must use Cword Method 2)\n");

	return;
}



void parseArgs(int argc, char *argv[], int rank, MPI_Comm mainComm, char **directory, int *maxFiles, int *seed, int *noTopics, int *iterations, int *printResults, int *parallelReadIn, int *CwordMethod, int *rma, int *bundleSize, double *pwToPdRatio, int *invMethod, int *sortRead){
	int pArgs[11];
	double pArgsD[1];
	int endProg = 0;

	//parse command line arguments on rank 0
	if(rank == 0){
		int opt;
		while((opt=getopt(argc,argv,"a:b:c:d:hi:m:oprs:t:v:w:"))!=-1){
			switch(opt){
				case 'a':
					*pwToPdRatio = strtod(optarg, NULL);
					if((*pwToPdRatio) <= 0 || (*pwToPdRatio) >= 1){
						printf("*** ERROR ***\nPw:Pd ratio must be in (0, 1)\n\n");
						endProg = 1;
					}
					break;
				case 'b':
					*bundleSize = atoi(optarg);
					if((*bundleSize) <= 0){
						printf("*** ERROR ***\nBundle Size must be non-zero positive integer\n\n");
						endProg = 1;
					}
					break;
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
				case 'o':
					*sortRead = 0;
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
				case 'v':
					*invMethod = atoi(optarg);
					if(*invMethod != 1 && *invMethod != 2){
						printf("*** ERROR ***\nInv Method must be in {1, 2}\n\n");
						endProg = 1;
					}
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
		if(*CwordMethod != 2 && *invMethod != 1){
			printf("*** ERROR ***\nMatrix method for Inverted Index only avilable with matrix version of Cword\n\n");
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
		pArgs[8] = *bundleSize;
		pArgs[9] = *invMethod;
		pArgs[10] = *sortRead;
		pArgsD[0] = *pwToPdRatio;
	}

	// Send relevant args to rest of processes
	MPI_Bcast(pArgs, 11, MPI_INT, 0, mainComm);
	MPI_Bcast(pArgsD, 1, MPI_DOUBLE, 0, mainComm);

	*seed = pArgs[0];
	*noTopics = pArgs[1];
	*iterations = pArgs[2];
	*maxFiles = pArgs[3];
	endProg = pArgs[4];
	*parallelReadIn = pArgs[5];
	*CwordMethod = pArgs[6];
	*rma = pArgs[7];
	*bundleSize = pArgs[8];
	*invMethod = pArgs[9];
	*sortRead = pArgs[10];
	*pwToPdRatio = pArgsD[0];

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
