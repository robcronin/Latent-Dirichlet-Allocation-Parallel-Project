#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>

// fills a dictionary from a given file
int	checkCopyIndiv(char *filename){
	int i, j;

	FILE *tester;
	if((tester=fopen(filename, "r")) == NULL){
		fprintf(stderr, "*** ERROR ***: Couldn't open %s\n\n", filename);
		return;
	}

	char lineBuffer[1024];
	fgets(lineBuffer, 1024, tester);
	int copyLine = 0;
	while(1){
		if( (fgets(lineBuffer, 1024, tester)) == NULL){
			break;
		}
		if(strncmp(lineBuffer, "rules is very easy", 18) == 0){
			if(strncmp(&(lineBuffer[21]), "You may use this eBook for nearly any purpose", 45) == 0 || strncmp(&(lineBuffer[22]), "You may use this eBook for nearly any purpose") == 0){
				copyLine ++;
			}
			fgets(lineBuffer, 1024, tester);
			if(strncmp(lineBuffer, "such as creation of derivative works, reports, performances and", 63) == 0){
				copyLine ++;
			}
			fgets(lineBuffer, 1024, tester);
			if(strncmp(lineBuffer, "research", 8) == 0){
				copyLine ++;
			}
			break;
		
			
		}
		else if(strncmp(lineBuffer, "for copies of this eBook, complying with the rules is very easy.", 64) == 0){
			if(strncmp( &(lineBuffer[65]), "You", 3) == 0 || strncmp( &(lineBuffer[66]), "You", 3) == 0){
				copyLine ++;
				fgets(lineBuffer, 1024, tester);
				if(strncmp(lineBuffer, "may use this eBook for nearly any purpose such as creation of derivative", 72) == 0){
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "works, reports, performances and research", 41) == 0){
					copyLine ++;	
					}
				}
			}
			break;
		}
		else if(strncmp(lineBuffer, "eBook, complying with the rules is very easy.", 45) == 0){
			if(strncmp(&(lineBuffer[46]), "You may use this eBook", 22) == 0 || strncmp(&(lineBuffer[47]), "You may use this eBook", 22) == 0){
				copyLine ++;
				fgets(lineBuffer, 1024, tester);
				if(strncmp(lineBuffer, "for nearly any purpose such as creation of derivative works, reports,", 69) == 0){
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "performances and research", 25) == 0){
						copyLine ++;
					}
				}
			}
			break;
		}
		else if(strncmp(lineBuffer, "copies of this ebook, complying with the rules is very easy.", 60) == 0){
			if(strncmp(&(lineBuffer[61]), "You may use", 11) == 0 || strncmp(&(lineBuffer[62]), "You may use", 11) == 0){
				copyLine ++;
				fgets(lineBuffer, 1024, tester);
				if(strncmp(lineBuffer, "this ebook for nearly any purpose such as creation of derivative works,", 71) == 0){
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "reports, performances and research", 34) == 0){
						copyLine ++;
					}
				}
			}
			break;
		}
		else if(strncmp(lineBuffer, "this eBook, complying with the rules is very easy.", 50) == 0){
			if(strncmp(&(lineBuffer[51]), "You may use this", 16) == 0 || strncmp(&(lineBuffer[52]), "You may use", 16) == 0){
				if(strncmp(&(lineBuffer[68]), "eBook", 5) == 0 || strncmp(&(lineBuffer[69]), "eBook", 5) == 0){
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "for nearly any purpose such as creation of derivative works,", 60) == 0){
						copyLine ++;
						fgets(lineBuffer, 1024, tester);
						if(strncmp(lineBuffer, "performances and research", 25) == 0){
							copyLine ++;
						}
					}
				
				}
				else{
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "eBook for nearly any purpose such as creation of derivative works,", 66) == 0){
						copyLine ++;
						fgets(lineBuffer, 1024, tester);
						if(strncmp(lineBuffer, "reports, performances and research", 34) == 0){
							copyLine ++;
						}
					}
				}
			}
			break;
		}
		else if(strncmp(lineBuffer, "the rules is very easy.", 23) == 0){
			if(strncmp(&(lineBuffer[24]), "You may use this eBook for nearly any", 37) == 0 || strncmp(&(lineBuffer[25]), "You may use this eBook for nearly any", 37) == 0){
				copyLine ++;
				fgets(lineBuffer, 1024, tester);
				if(strncmp(lineBuffer, "purpose such as creation of derivative works, reports, performances", 67) == 0){
					copyLine ++;
					fgets(lineBuffer, 1024, tester);
					if(strncmp(lineBuffer, "and research", 12) == 0){
						copyLine ++;
					}
				}
			}
			break;
		}
		else if(strncmp(lineBuffer, "very easy.", 10) == 0){
			if(strncmp(&(lineBuffer[11]), "You may use this eBook for nearly any purpose such as", 53) == 0 || strncmp(&(lineBuffer[12]), "You may use this eBook for nearly any purpose such as", 53) == 0){
				copyLine ++;
				fgets(lineBuffer, 1024, tester);
				if(strncmp(lineBuffer, "creation of derivative works, reports, performances and research", 64) == 0){
					copyLine += 2;
				}
			}
			break;
		}
		
	}


	if(copyLine == 3){return 1;}
	else{printf("%s has %d\n", filename, copyLine);}
	return 0;
}

// does TFIDF calculation for files in a given directory
void 	checkCopy(char *directory){
	DIR *FD;
	struct dirent *in_file;

	// opens directory and counts number of files
	if((FD=opendir(directory))==NULL){
		fprintf(stderr, "\t*** ERROR: Couldn't open %s directory\n", directory);
		exit(1);
	}
	
	int noFiles = 0;
	int copyOk = 0;
	// fills a dict for each file and does tf analysis
	while( (in_file = readdir(FD)) != NULL){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			char path[1024];
			sprintf(path, "%s/%s", directory, in_file->d_name);
			
			noFiles ++;

			copyOk += checkCopyIndiv(path);
		}
	}
	closedir(FD);


	// print results
	printf("Of %d files, %d files contain the free to use for reasearch line\n", noFiles, copyOk);

	return;
}


int 	main(int argc, char *argv[]){
	char 	*directory;

	if(argc == 2){
		directory = argv[1];
	}	
	else{
		directory = "../CorpusExtra";
	}
	printf("Checking %s for copyright\n", directory);
	checkCopy(directory);

	return 0;
}
