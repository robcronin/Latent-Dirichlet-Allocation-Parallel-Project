#include "plda+.h"



int 	main(int argc, char *argv[]){

	// PROGRAM SETTINGS
	// 	values which affect performance have been set to current known optimal values(for large scale cases)
	int 	seed = time(NULL);
	int		iterations = 100;
	int 	noTopics = 10;
	int 	parallelReadIn = 0;
	int		sortRead = 1;
	int		gutenFiles = 1;
	int 	maxFiles = 0;
	int 	LDAmethod = 3;
	double	pwToPdRatio = 0.6;
	int		bundleSize = 8;
	int		CwordMethod = 2;	// 1 for bst, 2 for matrix
	int		invMethod = 2;		// 1 for bst, 2 for matrix
	int		rma = 0;
	int 	threads = 0;
	int		printResults = 0;
	char 	*directory = malloc(1024);		// FREED
	char	*stopfile = malloc(1024);		// FREED
	char	*dictfile = malloc(1024);		// FREED
	sprintf(directory, "../../Inputs/Gutenberg");
	sprintf(stopfile, "../../Inputs/stopwords.txt");
	// If you want to use an overall dictionary (SLOWER - requires parallel read in also)
	//sprintf(dictfile, "../../Inputs/unixDict.txt");
	dictfile = NULL;

	#ifdef THREADS
	threads = 1;
	#endif






	// VARIABLES AND COUNTERS etc
	struct 	timeval t1, t1Sec, t2Sec;
	int 	noFiles, vocabSize, totalWords;
	int 	locNoFiles;
	char	**vocabulary;
	int		**topicCountPerDoc, ***topicIndex;
	double 	*alphaDir, *gammaDir, gammaSum;
	int 	*locUniqueCount, *locTotalCount;
	int 	***localDocMatrix;
	char 	**pathNames;
	int 	rank, size;
	int		pdSize, pwSize, prank;
	char 	procType;
	int		*wordSpread, *wordLocation;
	int		**wordLocationMat;
	int		**wordCountPerTopic;
	struct 	treeRootCword *wordCountPerTopicTree;
	struct	treeRootInv *invertedIndexTree;
	int		**bundles2d, **bundles2dConv, totNoBundles, circleStartBundle;
	int		*localTopicCountTotal;
	int		wordLocCols;
	int		CwordSize = 0;
	MPI_Comm pdComm, pwComm;
	MPI_Win	wordCountWin, topicCountWin;

	


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
	parseArgs(argc, argv, rank, mainComm, &directory, &maxFiles, &seed, &noTopics, &iterations, &printResults, &parallelReadIn, &CwordMethod, &rma, &bundleSize, &pwToPdRatio, &invMethod, &sortRead, &gutenFiles);

	if(threads == 1 && invMethod == 1){
		if(rank == 0){
			printf("*** ERROR ***\nMatrix version of Inverted method required for threads\n\n");
			printHelp();
		}
		MPI_Finalize();
		exit(0);
	}

	// need extra column in wordLocationMat if invMethod is 2
	if(invMethod == 1){wordLocCols = 2;}
	else{wordLocCols = 3;}

	// FIND NO FILES
	findNoFiles(directory, &noFiles, maxFiles);

	// ASSIGN PW and PD
	procAssignment(rank, size, noFiles, &pdComm, &pwComm, pwToPdRatio, &pdSize, &pwSize, &prank, &procType, mainComm);

	// PRINT INFO
	printInfoLDA(rank, 1, LDAmethod, seed, stopfile, noTopics, noFiles, vocabSize, iterations, size, parallelReadIn, totalWords, pdSize, pwSize, CwordMethod, rma, threads, bundleSize, invMethod, sortRead);

	// FIND PATH NAMES
	if(sortRead == 1){findPathNamesSorted(directory, &pathNames, &noFiles, maxFiles, pdSize);}
	else{findPathNames(directory, &pathNames, &noFiles, maxFiles);}
	free(directory);


	// READ IN FILES
	readInDocsAll(procType, prank, pdSize, noFiles, &locNoFiles, parallelReadIn, stopfile, &vocabSize, &vocabulary, &locUniqueCount, &locTotalCount, &localDocMatrix, printResults, pathNames, &totalWords, &t1Sec, &t2Sec, gutenFiles, dictfile, pdComm);

	// FREE
	freePathNames(pathNames);
	free(stopfile);
	free(dictfile);

	// PRINT REST OF INFO
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
		initialiseGamma(&gammaDir, vocabSize, &gammaSum);				// FREED

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
		createWindows(&wordCountWin, &topicCountWin, procType, mainComm, CwordSize, noTopics, wordCountPerTopic, localTopicCountTotal);	
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
			initTopicsPLDAplusMatrixRMA(procType, locNoFiles, locUniqueCount, localDocMatrix, noTopics, topicIndex, wordLocationMat, topicCountPerDoc, pdSize, size, rank, pwSize, wordCountWin, topicCountWin, mainComm);
		}
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
	circleQueue(procType, vocabSize, pdSize, prank, &circleStartBundle, bundles2d, totNoBundles);




	// WORD BUNDLES MAY NEED TO BE ALTERED TO ACCOUNT FOR NEW matIndex


	// To avoid cross contamination of sends a barrier is thrown before sampling
	MPI_Barrier(mainComm);

	// ************************************************
	// ************************************************
	// ************    RUN PROGRAM        *************
	// ************************************************
	// ************************************************
	if(CwordMethod == 1){
		runPLDAplusGibbs(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopicTree, rank, mainComm);
	}
	else if(CwordMethod == 2){
		if(threads == 1){
			//runPLDAplusGibbsPthreadsDeadline(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
			runPLDAplusGibbsPthreadsWorkload(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, invertedIndexMatrix,wordLocationMat, pdComm, mainComm);
		}
		else if(invMethod == 2){
			if(rma == 0 || rma == 1){
				runPLDAplusGibbsMatrixInv(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexMatrix, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, wordLocationMat, mainComm);
			}
			else{
				runPLDAplusGibbsMatrixInvRMA(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexMatrix, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, pwSize, wordLocationMat, wordCountWin, topicCountWin, mainComm);
			}
		}
		else if(rma == 0 || rma == 1){
			runPLDAplusGibbsMatrix(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, mainComm);
		}
		else{
			runPLDAplusGibbsMatrixRMA(procType, noTopics, pdSize, size, bundleSize, iterations, totNoBundles, circleStartBundle, bundles2d, bundles2dConv, invertedIndexTree, locNoFiles, localDocMatrix, topicIndex, topicCountPerDoc, alphaDir, gammaDir, gammaSum, wordCountPerTopic, localTopicCountTotal, rank, pwSize, wordCountWin, topicCountWin, mainComm);
		}
	}



	// FREE
	if(procType == 'd'){
		free(alphaDir);
		freeTopicIndex(topicIndex, locNoFiles);
		freeDocMatrix(&localDocMatrix, locNoFiles);
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
		collectAndPrint(rank, vocabSize, noTopics, pwSize, procType, wordCountPerTopicTree, printResults, gammaDir, vocabulary, mainComm);
	}
	else if(CwordMethod == 2){
		collectAndPrintMatrix(rank, vocabSize, noTopics, procType, wordCountPerTopic, printResults, gammaDir, vocabulary, wordLocationMat, pdSize, size, mainComm);
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
		free(gammaDir);
	}


	MPI_Finalize();
	return 0;
}








//# vim: ts=4
