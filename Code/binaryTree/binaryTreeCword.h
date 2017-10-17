/*	Binary Search Tree to store Cword information
 *		standart bst functions (add, rebalance)
 *		no delete required
 *		Nodes refer to word indexes and store word counts per topic
 */



#ifndef BSTCWORD
#define BSTCWORD



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../general/matrix.h"


struct nodeCword{
	int 	index;
	struct 	nodeCword *left;
	struct	nodeCword *right;
	int 	bal;

	// ADD HERE
	int 	*topicArray;
};



struct treeRootCword{
	struct 	nodeCword *head;
	int	size;

	// ADD HERE
	int 	noTopics;
	int	*topicCountTotal;
};


// change args
void 	createTreeCword(struct treeRootCword **root, int noTopics);

void	freeTreeCword(struct treeRootCword **root);
void	addValueCword(struct treeRootCword *root, int index);
int	findHeightCword(struct treeRootCword *root);
int	findSizeCword(struct treeRootCword *root);
void	rebalanceTreeCword(struct treeRootCword *root);
void	printTreeCword(struct treeRootCword *root);
struct nodeCword * searchTreeCword(struct treeRootCword *root, int index);
int	countTopicsCword(struct treeRootCword *root);


#ifdef MPI
#include <mpi.h>
void	sendForPrintCword(struct treeRootCword *root, MPI_Comm mainComm);
#endif


#endif
