/*      Binary Search Tree to store Inverted Index information
 *      	standart bst functions (add, rebalance)
 *      	no delete required
 *		Nodes refere to word indexes and store their location in a corpus of documents
 */



#ifndef BSTINV
#define BSTINV


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../general/matrix.h"


struct nodeInv{
	int 	index;
	struct 	nodeInv *left;
	struct	nodeInv *right;
	int 	bal;

	// ADD HERE
	int	*invIndex;
	int	occurrences;
};



struct treeRootInv{
	struct 	nodeInv *head;
	int	size;

	// ADD HERE
	int 	noFiles;
};


// change args
void 	createTreeInv(struct treeRootInv **root, int noFiles);

void	freeTreeInv(struct treeRootInv **root);
void	addValueInv(struct treeRootInv *root, int index, int doc, int location, int occ);
int	findHeightInv(struct treeRootInv *root);
int	findSizeInv(struct treeRootInv *root);
void	rebalanceTreeInv(struct treeRootInv *root);
void	printTreeInv(struct treeRootInv *root);
void	fillOccMat(int **occMat, struct treeRootInv *root, int *wordLocation);
void	fillOccMatMatrix(int **occMat, struct treeRootInv *root, int **wordLocationMat);
struct 	nodeInv * searchTreeInv(struct treeRootInv *root, int index);


#endif
