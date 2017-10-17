#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include "tfidf.h"


void 	parse_args(int argc, char *argv[], char **directory, int *maxFiles);
void	printHelp();



int main(int argc, char *argv[]){

	// default args
	char 	*directory="../../Inputs/Gutenberg";
	char 	**files;
	struct 	alpha **all;
	int 	noFiles;
	int 	maxFiles = 0;
	int	gutenFiles = 1;

	parse_args(argc, argv, &directory, &maxFiles);
	
	initialiseTFIDF(directory, &all, &files, &noFiles, maxFiles, "../../Inputs/stopwords.txt", gutenFiles);
	if(noFiles > 1){
		printTopXTFIDF(all, files, noFiles, 10);
	}
	else{
		printf("Require two or more documents to calculate tfidf values\n");
	}



	freeTFIDF(&all, &files, noFiles);		
	return 0;
}


void parse_args(int argc, char *argv[], char **directory, int *maxFiles){
	//parse command line arguments
	int opt;
	while((opt=getopt(argc,argv,"d:hm:"))!=-1){
		switch(opt){
			case 'd':
				*directory = optarg;
				break;
			case 'h':
				printHelp();
				exit(1);
			case 'm':
				*maxFiles = atoi(optarg);
				break;
			default:
				printHelp();
				exit(EXIT_FAILURE);
		}
	}
	return;
}


void printHelp(){
	printf("-d [DIR]\tSets directory (default = ../../Gutenberg)\n");
	printf("-h\t\tPrints help\n");
	printf("-m [INT]\tSets max number of files (default: no max)\n");

	return;
}
