// functions for running tfidf analysis on a directory of documents

#ifndef TFIDF
#define TFIDF


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <sys/stat.h>
#include "../general/matrix.h"
#include "../general/quicksort.h"


// holder for a word with its tfidf stats and pointer to next word
struct wordholder{
	char 	*word;
	struct 	wordholder *next;
	int	count;
	double 	tf;
	double 	idf;
	double 	tfidf;
	double 	*matRow;
	int 	*wordCountPerTopic;		// n_k,w
};

// pointer to the first word for a letter index
struct letter_head{
	struct 	wordholder *head;
};

// stores all pointers for each letter index and a total word count
struct alpha{
	struct 	letter_head **letters;
	int 	total_count;
	int 	uniqueCount;
};


void 	initialiseTFIDF(char *directory, struct alpha ***all, char ***files, int *noFiles, int maxFiles, char *stopfile, int gutenFiles);
void 	printTopXTFIDF(struct alpha **all, char **files, int noFiles, int no);
void	print_dict(struct alpha *dict, int tfidf);
void 	freeTFIDF(struct alpha ***all, char *** files, int noFiles);
void    freeLocalTFIDF(struct alpha ***all, int noFiles);




void 	create_dict(struct alpha **dict);
//void 	add_word(struct alpha *dict, char *new_word);
int	is_present(struct alpha *dict, char *word);
void	fill_dict(char *filename, struct alpha **dict, struct alpha *stop, int gutenFiles);
//void	findTitle(char *filename, char **files, int index); 
//void	calc_tf(struct alpha *dict);
//void	calc_tfidf(struct alpha **all, int noFiles);
void	free_dict(struct alpha **dict);

void 	createBeta(struct alpha **all, int noFiles, struct alpha *beta);
void	freeBeta(struct alpha **dict);
void 	createBetaTopics(struct alpha **all, int noFiles, struct alpha *beta, int noTopics);
void	freeBetaTopics(struct alpha **dict);





void	findNoFiles(char *directory, int *noFiles, int maxFiles);
void	findPathNames(char *directory, char ***pathNames, int *noFiles, int maxFiles);
void	findPathNamesSorted(char *directory, char ***pathNames, int *noFiles, int maxFiles, int size);
void    freePathNames(char **pathNames);
void	initialiseLocalTFIDF(char **pathNames, struct alpha ***all, int locNoFiles, char *stopfile, int startIndex, int gutenFiles);
void	initialiseLocalTFIDFDict(char **pathNames, struct alpha ***all, int locNoFiles, char *stopfile, int startIndex, int gutenFiles, struct alpha *phi);
void    allocPhiBuffer(char ***phiBuffer, int sendCount, int allocLength);
void    freePhiBuffer(char **phiBuffer);
void    fillPhiBuffer(struct alpha *phi, char **phiBuffer, int allocLength);
void    fillBetaFromBuffer(struct alpha *phi, char **phiBuffer, int sendCount);


#endif
