#include "tfidf.h"


// allocates memory and initialses a dictionary
void create_dict(struct alpha **dict){
	int i;
	*dict=(struct alpha *)malloc(sizeof(struct alpha));
	struct letter_head *temp;
	struct letter_head **tempp;
	tempp=malloc(27*sizeof(struct letter_head *));
	(*dict)->letters=tempp;
	for(i = 0; i < 27; i++){
		temp=malloc(sizeof(struct letter_head));
		(*dict)->letters[i]=temp;
		(*dict)->letters[i]->head=NULL;
	}
	(*dict)->total_count=0;
	(*dict)->uniqueCount=0;

	return;
}

// adds a given word to a dictionary, or adds to the counter if already in
void add_word(struct alpha *dict, char *new_word){
	int verbose = 0;

	// checks if word is empty
	if(strcmp(new_word, "\0")==0 || strlen(new_word) >= 28){
		free(new_word);
		return;
	}
	dict->total_count++;

	// finds index for word
	char let=new_word[0]-97;
	if(let < 0 || let > 25){
		let = 26;
	}
	if(verbose){
		printf("adding %s to letter %d\n", new_word, let);
	}
	struct letter_head *root = dict->letters[(int) let];
	struct wordholder *search = root->head;
	struct wordholder *temp;

	// if that letter index is empty, adds to the root
	if(search == NULL){
		if(verbose){
			printf("adding %s to empty tree %d\n", new_word, let);
		}
		temp=malloc(sizeof(struct wordholder));
		temp->word = new_word;
		temp->count = 1;
		temp->next = NULL;
		root->head = temp;
		dict->uniqueCount++;

	}
	else{
		// adds before the root
		if(strcmp(new_word, search->word) < 0){
			if(verbose){
				printf("Adding %s before root %d\n", new_word, let);
			}
			temp = malloc(sizeof(struct wordholder));
			temp->word = new_word;
			temp->count = 1;
			temp->next = root->head;
			root->head = temp;
			dict->uniqueCount++;
		}
		// root is the word
		else if(strcmp(new_word, search->word) == 0){
			search->count++;
			free(new_word);
		}
		// else is after root
		else{
			// searches for the word
			while(search->next != NULL){
				// if next word is the new word
				if(strcmp(search->next->word, new_word) == 0){
					search->next->count++;
					free(new_word);
					return;
				}
				// else if it is less than, move forward
				else if(strcmp(search->next->word, new_word) < 0){
					search = search->next;
				}
				// else it is greater than, adds into the gap
				else{
					if(verbose){
 						printf("Adding %s inside tree %d\n", new_word, let);
					}
					temp = malloc(sizeof(struct wordholder));
					temp->word = new_word;
					temp->count = 1;
					temp->next = search->next;
					search->next = temp;
					dict->uniqueCount++;
					return;
				}
					
			}
			// if we reach here the word hasn't been added yet, so add to end
			if(verbose){
				printf("Adding %s to end of tree %d\n", new_word, let );
			}
			temp = malloc(sizeof(struct wordholder));
			temp->word = new_word;
			temp->count = 1;
			temp->next = NULL;
			search->next = temp;
			dict->uniqueCount++;
			
		}
	}

	return;
			
}

// checks if a word is in a dict
int	is_present(struct alpha *dict, char *word){
	
	// finds index
	char let = word[0]-97;
	if(let<0 || let>25){
		let=26;
	}
	
	// loops through and searches
	struct wordholder *search = dict->letters[(int) let]->head;
	while(search != NULL && strcmp(search->word, word) <= 0){
		if(strcmp(word, search->word) == 0){
			return 1;
		}
		search = search->next;
	}
	return 0;
}

void 	findTitle(char *filename, char **files, int index){

	FILE *tester;
	if((tester=fopen(filename, "r")) == NULL){
		fprintf(stderr, "*** ERROR ***: Couldn't open %s\n\n", filename);
		exit(1);
	}

	char lineBuffer[1024];

	// finds title of book
	while(fgets(lineBuffer, 1024, tester)!=NULL){
		if(strncmp(lineBuffer, "Title:", 6)==0){
			lineBuffer[strcspn(lineBuffer, "\n")] = '\0';
			strcpy(files[index], &(lineBuffer[7]));
		}
	}

	fclose(tester);

	return;
}

// fills a dictionary from a given file
void	fill_dict(char *filename, struct alpha **dict, struct alpha *stop, int gutenFiles){
	int i, j;

	//printf("Filling dict from %s\n", filename);
	create_dict(dict);
	FILE *tester;
	if((tester=fopen(filename, "r")) == NULL){
		fprintf(stderr, "*** ERROR ***: Couldn't open %s\n\n", filename);
		exit(1);
		return;
	}

	char lineBuffer[1024];

	// finds start and end indices from start of file if inputted manually
	fgets(lineBuffer, 1024, tester);
	fgets(lineBuffer, 1024, tester);
	int start = atoi(lineBuffer);
	fgets(lineBuffer, 1024, tester);
	int end = atoi(lineBuffer);
	

	// else finds them via Gutenberg conventions if requested
	if((start == 0 || end == 0) && gutenFiles == 1){
		rewind(tester);
		start = 1;

		// get first line
		fgets(lineBuffer, 1024, tester);
		int counter = 0;

		// find start of Gutenberg
		while(strncmp(lineBuffer, "*** START OF", 12) != 0 && strncmp(lineBuffer, "***START OF", 11) != 0){
			if((fgets(lineBuffer, 1024, tester)) == NULL){
				printf("\tERROR reached line %d in file %s with no end\n", counter, filename);

				// To avoid seg faults from empty dict
				char *tw;
				tw = malloc(8);
				strcpy(tw, "HOLDER");
				add_word(*dict, tw);
				return;
			}
			start ++;
			counter ++;
		}

		// find end of Gutenberg
		end = start;
		while(strncmp(lineBuffer, "*** END OF", 10) != 0 && strncmp(lineBuffer, "***END OF", 9) != 0){
		//while(strncmp(lineBuffer, "*** END OF", 10) != 0){
			if((fgets(lineBuffer, 1024, tester)) == NULL){
				printf("\tERROR reached line %d in file %s with no stop\n", counter, filename);

				// To avoid seg faults from empty dict
				char *tw;
				tw = malloc(8);
				strcpy(tw, "HOLDER");
				add_word(*dict, tw);
				return;
			}
			end ++;
		}

		start += 5;
		end -= 5;
	}
	// else reads whole file
	else if(start == 0 || end == 0){
		start = 0;
		end = 2147483647;
	}
		
	// reset to start
	rewind(tester);
	for(i = 0; i < start - 1; i++){
		fgets(lineBuffer, 1024, tester);
	}

	char *word;
	
	// loops through whole files
	for(i = 0; i < end - start; i++){
		char temp[1024];
		if(fgets(temp, 1024, tester) == NULL){
			break;
		}

		// if non-empty line
		if(strlen(temp) > 1){
			word = strtok(temp, " ");

			// loops through words in that line
			while(word != NULL && strncmp(word, " ", 1) != 0){
		
				// formats word 
				while(strcspn(word, "\r\n,.!?\"-_:;(){}") == 0 && strlen(word) > 0){
					word++;
				}
				word[strcspn(word, "\r\n,.!?\"-_:;(){}")] = '\0';
				while(word[0] < 0){
					word++;
				}
				while(word[strlen(word)-1] < 0){
					word[strlen(word)-1] = '\0';
				}

				// puts in lower case and fixes ' issue
				for(j = 0; j < (int) strlen(word); j++){
					word[j] = tolower(word[j]);
					// converts ' properly
					if(word[j] < 0){
						word[j] = 39;
						int k;
						for(k = j+1; k < (int) strlen(word); k++){
							word[k]=word[k+2];
						}
					}
				}

				// allocates memory for word and adds to dictionary
				char *tw;
				tw = malloc(strlen(word)+2);
				strcpy(tw, word);
				if(!is_present(stop, tw)){
					add_word(*dict, tw);
				}
				else{
					free(tw);
				}
				word = strtok(NULL, " ");
			}
		}

		/* next line
		if( (fgets(temp, 1024, tester)) == NULL){
			return;
		}*/
	}
	fclose(tester);

	//printf("\tDONE\n");

	return;
}



void	fill_dictUsingMainDict(char *filename, struct alpha **dict, struct alpha *stop, int gutenFiles, struct alpha *phi){
	int i, j;

	//printf("Filling dict from %s\n", filename);
	create_dict(dict);
	FILE *tester;
	if((tester=fopen(filename, "r")) == NULL){
		fprintf(stderr, "*** ERROR ***: Couldn't open %s\n\n", filename);
		exit(1);
		return;
	}

	char lineBuffer[1024];

	// finds start and end indices from start of file if inputted manually
	fgets(lineBuffer, 1024, tester);
	fgets(lineBuffer, 1024, tester);
	int start = atoi(lineBuffer);
	fgets(lineBuffer, 1024, tester);
	int end = atoi(lineBuffer);
	

	// else finds them via Gutenberg conventions if requested
	if((start == 0 || end == 0) && gutenFiles == 1){
		rewind(tester);
		start = 1;

		// get first line
		fgets(lineBuffer, 1024, tester);
		int counter = 0;

		// find start of Gutenberg
		while(strncmp(lineBuffer, "*** START OF", 12) != 0 && strncmp(lineBuffer, "***START OF", 11) != 0){
			if((fgets(lineBuffer, 1024, tester)) == NULL){
				printf("\tERROR reached line %d in file %s with no end\n", counter, filename);

				// To avoid seg faults from empty dict
				char *tw;
				tw = malloc(8);
				strcpy(tw, "HOLDER");
				add_word(*dict, tw);
				return;
			}
			start ++;
			counter ++;
		}

		// find end of Gutenberg
		end = start;
		while(strncmp(lineBuffer, "*** END OF", 10) != 0 && strncmp(lineBuffer, "***END OF", 9) != 0){
		//while(strncmp(lineBuffer, "*** END OF", 10) != 0){
			if((fgets(lineBuffer, 1024, tester)) == NULL){
				printf("\tERROR reached line %d in file %s with no stop\n", counter, filename);

				// To avoid seg faults from empty dict
				char *tw;
				tw = malloc(8);
				strcpy(tw, "HOLDER");
				add_word(*dict, tw);
				return;
			}
			end ++;
		}

		start += 5;
		end -= 5;
	}
	// else reads whole file
	else if(start == 0 || end == 0){
		start = 0;
		end = 2147483647;
	}
		
	// reset to start
	rewind(tester);
	for(i = 0; i < start - 1; i++){
		fgets(lineBuffer, 1024, tester);
	}

	char *word;
	
	// loops through whole files
	for(i = 0; i < end - start; i++){
		char temp[1024];
		if(fgets(temp, 1024, tester) == NULL){
			break;
		}

		// if non-empty line
		if(strlen(temp) > 1){
			word = strtok(temp, " ");

			// loops through words in that line
			while(word != NULL && strncmp(word, " ", 1) != 0){
		
				// formats word 
				while(strcspn(word, "\r\n,.!?\"-_:;(){}") == 0 && strlen(word) > 0){
					word++;
				}
				word[strcspn(word, "\r\n,.!?\"-_:;(){}")] = '\0';
				while(word[0] < 0){
					word++;
				}
				while(word[strlen(word)-1] < 0){
					word[strlen(word)-1] = '\0';
				}

				// puts in lower case and fixes ' issue
				for(j = 0; j < (int) strlen(word); j++){
					word[j] = tolower(word[j]);
					// converts ' properly
					if(word[j] < 0){
						word[j] = 39;
						int k;
						for(k = j+1; k < (int) strlen(word); k++){
							word[k]=word[k+2];
						}
					}
				}

				// allocates memory for word and adds to dictionary
				char *tw;
				tw = malloc(strlen(word)+2);
				strcpy(tw, word);
				if(is_present(stop, tw) == 0 && is_present(phi, tw) == 1){
					add_word(*dict, tw);
				}
				else{
					free(tw);
				}
				word = strtok(NULL, " ");
			}
		}

		/* next line
		if( (fgets(temp, 1024, tester)) == NULL){
			return;
		}*/
	}
	fclose(tester);

	//printf("\tDONE\n");

	return;
}




// calculates total frequencies for a filled dictionary
void	calc_tf(struct alpha *dict){
	
	// if empty, no work to be done
	if(dict == NULL){
		return;
	}
	int i;
	struct wordholder *search;

	// loops through each index
	for(i = 0; i < 27; i++){
		search = dict->letters[i]->head;
			
		// loops through each word in index, calcs tf and moves on
		while(search != NULL){
			search->tf = search->count/(double)dict->total_count;
			search = search->next;
		}
	}
	
	return;
}

// calculates tfidf for a number of documents
void	calc_tfidf(struct alpha **all, int noFiles){
	int i, j, k;
	struct wordholder *search;
	int doc_counter;

	// loops through each file
	for(i = 0; i < noFiles; i++){
		// through each letter index
		for(j = 0; j < 27; j++){
			search=all[i]->letters[j]->head;

			// through each word
			while(search!=NULL){
				doc_counter=1;
	
				// counts number of other docs the word appears in 
				for(k = 0; k < noFiles; k++){
					if(i!=k){
						doc_counter+=is_present(all[k], search->word);
					}
				}

				// does idf and tfidf calculations
				search->idf=log(noFiles/doc_counter);
				search->tfidf=search->tf*search->idf;
				search=search->next;		
			}
		}
	}

	return;
}

// does TFIDF calculation for files in a given directory
void 	initialiseTFIDF(char *directory, struct alpha ***all, char ***files, int *noFiles, int maxFiles, char *stopfile, int gutenFiles){
	int i;
	DIR *FD;
	struct dirent *in_file;
	*noFiles = 0;

	// opens directory and counts number of files
	if((FD=opendir(directory))==NULL){
		fprintf(stderr, "\t*** ERROR: Couldn't open %s directory\n", directory);
		exit(1);
	}
	while( (in_file = readdir(FD)) != NULL && ((*noFiles) < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			(*noFiles)++;
		}
	}
	rewinddir(FD);

	
	// allocates memory
	*files = malloc( (*noFiles) * sizeof(char *));
	*all = malloc( (*noFiles) *sizeof(struct alpha *));
	i=0;

	// creates stop file
	struct alpha *stop;
	if(stopfile != NULL){
		struct alpha *holder;
		create_dict(&holder);
		fill_dict(stopfile, &stop, holder, gutenFiles);
		free_dict(&holder);
	}
	else{
		create_dict(&stop);
	}

	// fills a dict for each file and does tf analysis
	while( (in_file = readdir(FD)) != NULL && (i < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			//(*files)[i]=in_file->d_name;
			(*files)[i] = malloc(1024);
			strcpy( (*files)[i], in_file -> d_name);
			char path[1024];
			sprintf(path, "%s/%s", directory, in_file->d_name);

			fill_dict(path, &((*all)[i]), stop, gutenFiles);
			findTitle(path, *files, i);
			calc_tf((*all)[i]);
			i++;
		}
	}
	free_dict(&stop);
	closedir(FD);


	// does idf and tfidf calculations once all documents added
	calc_tfidf(*all, *noFiles);

	return;
}

// prints the top X words in terms of tfidf
void 	printTopXTFIDF(struct alpha **all, char **files, int noFiles, int no){
	int i, j, k;
	struct wordholder *search;

	
	double *max = malloc(no * sizeof(double));
	char **word = malloc(no * sizeof(char *));
	for(i = 0; i < no; i++){
		word[i] = malloc(1024 * sizeof(char));
		strcpy(word[i], "HOLDER");
	}


	// loops through all files
	for(i = 0; i < noFiles; i++){
		// resets max counters
		for(j = 0; j < no; j++){
			max[j]=0;
		}

		// loops through each letter index
		for(j = 0; j < 27; j++){
			search=all[i]->letters[j]->head;

			// loops through each word
			while(search!=NULL){
				// if the word is in the top X, find its position and update top X
				if(search->tfidf>max[no-1]){
					max[no-1]=search->tfidf;
					strcpy(word[no-1], search->word);
					for(k = no-2; k >= 0; k--){
						if(search->tfidf>max[k]){
							max[k+1]=max[k];
							strcpy(word[k+1], word[k]);
							max[k]=search->tfidf;
							strcpy(word[k], search->word);
						}
						else{break;}
					}
				}
				search=search->next;
			}
		}

		// prints results
		printf("%s\n", files[i]);
		for(j = 0; j < no; j++){
			printf("\t%s with %.10lf\n", word[j], max[j]);
		}
	}

	free(max);
	for(i = 0; i < no; i++){
		free(word[i]);
	}
	free(word);

	return;

}


// prints scores for a given dictionary
void	print_dict(struct alpha *dict, int tfidf){
	
	// checks if dict initialised
	if(dict == NULL){
		fprintf(stderr, "\t*** WARNING: Directory not initialised\n");
		return;
	}
	int i;
	struct wordholder *search;

	// if no tfidf scores requested
	if(!tfidf){
		// loops through the words and prints them with counter
		for(i = 0; i < 27; i++){
			search = dict->letters[i]->head;
			while(search != NULL){
				// basic formatting
				if(strlen(search->word) > 7){
					printf("%s\t%d\n", search->word, search->count);
				}
				else{
					printf("%s\t\t%d\n", search->word, search->count);
				}
				search = search->next;
			}
		}
	}
	// else prints with tfidf scores aswell
	else{
		for(i = 0; i < 27; i++){
			search = dict->letters[i]->head;
			while(search != NULL){
				// basic formatting
				if(strlen(search->word) > 7){
					printf("%s\t%d\t%lf\t%lf\t%lf\n", search->word, search->count, search->tf, search->idf, search->tfidf);
				}
				else{
					printf("%s\t\t%d\t%lf\t%lf\t%lf\n", search->word, search->count, search->tf, search->idf, search->tfidf);
				}
				search = search->next;
			}
		}
	}

	return;
}

// frees each word in a dict and then the dict
void	free_dict(struct alpha **dict){
	int i;
	struct wordholder *temp1;
	struct wordholder *temp2;

	// loops through each letter index
	for(i = 0; i < 27; i++){
		temp1 = (*dict)->letters[i]->head;
		// loops through words
		while(temp1 != NULL){
			temp2 = temp1->next;	// stores pointer to next word
			free(temp1->word);	// frees actual word
			free(temp1);		// frees holder
			temp1 = temp2;		// updates pointer
		}
		free((*dict)->letters[i]);	// frees letter holder
	}
	// frees pointers to letters, then the dict
	free((*dict)->letters);
	free(*dict);

	return;
}

// frees all memory from TFIDF analysis		
void 	freeTFIDF(struct alpha ***all, char *** files, int noFiles){
	int i;
	// frees each dict and the files array
	for( i = 0; i < noFiles; i++){
		free_dict( &((*all)[i]) );
		free( (*files)[i] );
	}
	free(*all);
	free(*files);

	return;
}


// frees all memory from TFIDF analysis		
void 	freeLocalTFIDF(struct alpha ***all, int noFiles){
	int i;
	// frees each dict and the files array
	for( i = 0; i < noFiles; i++){
		free_dict( &((*all)[i]) );
	}
	free(*all);

	return;
}


// PERHAPS CREATE NEW WORDHOLDER FOR BETA
void 	createBeta(struct alpha **all, int noFiles, struct alpha *beta){
	int i, j;
	struct wordholder *search;

	// add all words
	for(i = 0; i < noFiles; i++){
		for(j = 0; j < 27; j++){
			search = all[i]->letters[j]->head;
			while(search != NULL){
				char *tw;
				tw = malloc(strlen(search->word)+2);
				strcpy(tw, search->word);
				add_word(beta, tw);			
				search = search->next;
			}
		}
	}

	return;
}
			
// free Beta
void	freeBeta(struct alpha **dict){
	int i;
	struct wordholder *temp1;
	struct wordholder *temp2;

	// loops through each letter index
	for(i = 0; i < 27; i++){
		temp1 = (*dict)->letters[i]->head;
		// loops through words
		while(temp1 != NULL){
			temp2 = temp1->next;	// stores pointer to next word
			free(temp1->word);	// frees actual word
			free(temp1);		// frees holder
			temp1 = temp2;		// updates pointer
		}
		free((*dict)->letters[i]);	// frees letter holder
	}
	// frees pointers to letters, then the dict
	free((*dict)->letters);
	free(*dict);

	return;
}



// PERHAPS CREATE NEW WORDHOLDER FOR BETA
void 	createBetaTopics(struct alpha **all, int noFiles, struct alpha *beta, int noTopics){
	int i, j;
	struct wordholder *search;

	// add all words
	for(i = 0; i < noFiles; i++){
		for(j = 0; j < 27; j++){
			search = all[i]->letters[j]->head;
			while(search != NULL){
				char *tw;
				tw = malloc(strlen(search->word)+2);
				strcpy(tw, search->word);
				add_word(beta, tw);			
				search = search->next;
			}
		}
	}

	// assign memory for matrix
	for(i = 0; i < 27; i++){
		search = beta->letters[i]->head;
		while(search != NULL){
			//double *temp = malloc(noTopics*sizeof(double));
			(search->matRow) = malloc(noTopics*sizeof(double));
			(search->wordCountPerTopic) = malloc(noTopics*sizeof(int));
			for(j = 0; j < noTopics; j++){
				search->wordCountPerTopic[j] = 0;
			}
			//search->matRow = temp;
			search = search->next;
		}
	}

	return;
}



void	freeBetaTopics(struct alpha **dict){
	int i;
	struct wordholder *temp1;
	struct wordholder *temp2;

	// loops through each letter index
	for(i = 0; i < 27; i++){
		temp1 = (*dict)->letters[i]->head;
		// loops through words
		while(temp1 != NULL){
			temp2 = temp1->next;	// stores pointer to next word
			free(temp1->word);	// frees actual word
			free(temp1->matRow);	// frees Matrix part			// only difference
			free(temp1->wordCountPerTopic);
			free(temp1);		// frees holder
			temp1 = temp2;		// updates pointer
		}
		free((*dict)->letters[i]);	// frees letter holder
	}
	// frees pointers to letters, then the dict
	free((*dict)->letters);
	free(*dict);

	return;
}


void	findNoFiles(char *directory, int *noFiles, int maxFiles){
	DIR *FD;
	struct dirent *in_file;
	*noFiles = 0;

	// opens directory and counts number of files
	if((FD=opendir(directory))==NULL){
		fprintf(stderr, "\t*** ERROR: Couldn't open %s directory\n", directory);
		exit(1);
	}
	while( (in_file = readdir(FD)) != NULL && ((*noFiles) < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			(*noFiles)++;
		}
	}
	closedir(FD);

	return;
}



void	findPathNames(char *directory, char ***pathNames, int *noFiles, int maxFiles){
	int i;
	DIR *FD;
	struct dirent *in_file;
	*noFiles = 0;

	// opens directory and counts number of files
	if((FD=opendir(directory))==NULL){
		fprintf(stderr, "\t*** ERROR: Couldn't open %s directory\n", directory);
		exit(1);
	}
	while( (in_file = readdir(FD)) != NULL && ((*noFiles) < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			(*noFiles)++;
		}
	}
	rewinddir(FD);

	
	// allocates memory
	char *temp = malloc( (*noFiles) * 1024 * sizeof(char));
	*pathNames = malloc( (*noFiles) * sizeof(char *));
	for(i = 0; i < *noFiles; i++){
		(*pathNames)[i] = &(temp[1024*i]);
	}
	i=0;

	// finds Path Names
	while( (in_file = readdir(FD)) != NULL && (i < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			//(*pathNames)[i] = malloc(1024);
			sprintf((*pathNames)[i], "%s/%s", directory, in_file -> d_name); 

			i++;
		}
	}
	closedir(FD);

	return;
}



void	findPathNamesSorted(char *directory, char ***pathNames, int *noFiles, int maxFiles, int size){
	int i;
	DIR *FD;
	struct dirent *in_file;
	*noFiles = 0;

	// opens directory and counts number of files
	if((FD=opendir(directory))==NULL){
		fprintf(stderr, "\t*** ERROR: Couldn't open %s directory\n", directory);
		exit(1);
	}
	while( (in_file = readdir(FD)) != NULL && ((*noFiles) < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			(*noFiles)++;
		}
	}
	rewinddir(FD);

	
	// allocates memory
	char *temp = malloc( (*noFiles) * 1024 * sizeof(char));
	*pathNames = malloc( (*noFiles) * sizeof(char *));

	char *temp2 = malloc( (*noFiles) * 1024 * sizeof(char));
	char **pathHolder = malloc( (*noFiles) * sizeof(char *));
	for(i = 0; i < *noFiles; i++){
		(*pathNames)[i] = &(temp[1024*i]);
		pathHolder[i] = &(temp2[1024*i]);
	}
	i=0;

	// finds Path Names
	while( (in_file = readdir(FD)) != NULL && (i < maxFiles || maxFiles == 0)){
		if(strncmp(in_file->d_name, ".", 1)!=0){
			//(*pathNames)[i] = malloc(1024);
			sprintf(pathHolder[i], "%s/%s", directory, in_file -> d_name); 
			i++;
		}
	}

	// get file sizes
	struct stat fs;
	int	**filesizeMatrix;
	allocMatInt(&filesizeMatrix, (*noFiles), 2);
	for(i = 0; i < (*noFiles); i++){
		filesizeMatrix[i][0] = i;	
		stat(pathHolder[i], &fs);
		filesizeMatrix[i][1] = fs.st_size;
	}
	quicksortMat(filesizeMatrix, (*noFiles), 2, 1);

	int	filesEach = (*noFiles) / size;
	int	noRemainders = (*noFiles) % size; 
	int	noLoops = filesEach / 2;
	int	extraLoop = filesEach % 2;
	int	*startPoint = malloc( (*noFiles) * sizeof(int));
	for(i = 0; i < noRemainders; i++){
		startPoint[i] = i*(filesEach+1);
	}
	for(i = noRemainders; i < size; i++){
		startPoint[i] = i*filesEach + noRemainders;
	}

	int j;
	int	fileIndex = 0;
	int 	nextFile = (*noFiles) - 1;
	for(i = 0; i < noLoops; i++){
		// loop forward
		for(j = 0; j < size; j++){
			sprintf( (*pathNames)[ startPoint[j] + fileIndex], "%s", pathHolder[filesizeMatrix[nextFile][0]]);
			nextFile --;
		}
		fileIndex ++;
		// then back
		for(j = size - 1; j >= 0; j--){
			sprintf( (*pathNames)[ startPoint[j] + fileIndex], "%s", pathHolder[filesizeMatrix[nextFile][0]]);
			nextFile --;
		}
		fileIndex ++;
	}
	// extra loop forward if necessary
	if(extraLoop){
		// loop forward
		for(j = 0; j < size; j++){
			sprintf( (*pathNames)[ startPoint[j] + fileIndex], "%s", pathHolder[filesizeMatrix[nextFile][0]]);
			nextFile --;
		}
		fileIndex ++;
	}
	// then the leftovers
	for(j = 0; j < noRemainders; j++){
		sprintf( (*pathNames)[ startPoint[j] + fileIndex], "%s", pathHolder[filesizeMatrix[nextFile][0]]);
		nextFile --;
	}



	closedir(FD);
	free(temp2);

	return;
}



void 	freePathNames(char **pathNames){
	free(pathNames[0]);
	free(pathNames);

	return;

}



void 	initialiseLocalTFIDF(char **pathNames, struct alpha ***all, int locNoFiles, char *stopfile, int startIndex, int gutenFiles){
	int i;


	
	// allocates memory
	*all = malloc( locNoFiles *sizeof(struct alpha *));

	// creates stop file
	struct alpha *stop;
	if(stopfile != NULL){
		struct alpha *holder;
		create_dict(&holder);
		fill_dict(stopfile, &stop, holder, gutenFiles);
		free_dict(&holder);
	}
	else{
		create_dict(&stop);
	}

	// fills a dict for each file
	
	for(i = 0; i < locNoFiles; i++){
		fill_dict(pathNames[startIndex + i], &((*all)[i]), stop, gutenFiles);
	}


	free_dict(&stop);


	// does idf and tfidf calculations once all documents added
//	calc_tfidf(*all, *noFiles);

	return;
}

void 	initialiseLocalTFIDFDict(char **pathNames, struct alpha ***all, int locNoFiles, char *stopfile, int startIndex, int gutenFiles, struct alpha *phi){
	int i;


	
	// allocates memory
	*all = malloc( locNoFiles *sizeof(struct alpha *));

	// creates stop file
	struct alpha *stop;
	if(stopfile != NULL){
		struct alpha *holder;
		create_dict(&holder);
		fill_dict(stopfile, &stop, holder, gutenFiles);
		free_dict(&holder);
	}
	else{
		create_dict(&stop);
	}

	// fills a dict for each file
	
	for(i = 0; i < locNoFiles; i++){
		fill_dictUsingMainDict(pathNames[startIndex + i], &((*all)[i]), stop, gutenFiles, phi);
	}


	free_dict(&stop);


	// does idf and tfidf calculations once all documents added
//	calc_tfidf(*all, *noFiles);

	return;
}

void	allocPhiBuffer(char ***phiBuffer, int sendCount, int allocLength){
	char *temp;
	int i;
	
	
	temp = malloc(allocLength * sendCount * sizeof(char));
	*phiBuffer = malloc( sendCount * sizeof(char *));

	for(i = 0; i < sendCount; i++){
		(*phiBuffer)[i] = &(temp[allocLength*i]);
	}

	return;

}

void	freePhiBuffer(char **phiBuffer){
	free(phiBuffer[0]);
	free(phiBuffer);

	return;

}

void	fillPhiBuffer(struct alpha *phi, char **phiBuffer, int allocLength){
	int i;
	int index = 0;


	struct wordholder *search;
	for(i = 0; i < 27; i++){
		search = phi -> letters[i] -> head;
		while(search != NULL){
			strncpy( phiBuffer[index] , search -> word, allocLength-2);
			index ++;
			search = search -> next;
			
		}
	}

	return;
}

void	fillBetaFromBuffer(struct alpha *phi, char **phiBuffer, int sendCount){
	int i;
	
	for(i = 0; i < sendCount; i++){
		char *tw;
		tw = malloc(strlen(phiBuffer[i])+2);
		strcpy(tw, phiBuffer[i]);
		add_word(phi, tw);			
	}
	
	return;
}


