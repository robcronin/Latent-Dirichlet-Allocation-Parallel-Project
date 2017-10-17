#include "binaryTreeCword.h"

// creates a tree
// CHANGE ARGS
void 	createTreeCword(struct treeRootCword **root, int noTopics){
	*root = malloc(sizeof(struct treeRootCword));
	(*root) -> head = NULL;
	(*root) -> size = 0;

	// CHANGE HERE
	(*root) -> noTopics = noTopics;
	callocArrInt( &((*root) -> topicCountTotal), noTopics);

	/*
	int i;
	for(i = 0; i < noTopics; i++){
		printf("%d  ", (*root)->topicCountTotal[i]);
	}
	printf("\n");
	sleep(2);
	*/


	return;
}


// allocates node and initialises values
void	allocNodeCword(struct nodeCword **newNode, int index, struct treeRootCword *root){
	*newNode = malloc(sizeof(struct nodeCword));
	(*newNode) -> index = index;
	(*newNode) -> left = NULL;
	(*newNode) -> right = NULL;
	(*newNode) -> bal = 0;

	// CHANGE HERE
	callocArrInt( &((*newNode) -> topicArray), root -> noTopics );

	/*
	int i;
	for(i = 0; i < root->noTopics; i++){
		printf("%d  ", (*newNode)->topicArray[i]);
	}
	printf("\n");
	*/

	return;
}


// frees an individual node	
void	freeNodeCword(struct nodeCword **node){
	// CHANGE HERE
	free( (*node) -> topicArray);

	free(*node);
	
	return;
}


// recursively frees nodes
void	recursiveFreeCword(struct nodeCword **node){
	if(*node == NULL){
		return;
	}
	recursiveFreeCword( &((*node) -> left));
	recursiveFreeCword( &((*node) -> right));
	freeNodeCword(node);

	return;
}



// frees the Tree
void	freeTreeCword(struct treeRootCword **root){
	if( *root == NULL){
		return;
	}
	else{
		recursiveFreeCword(&( (*root) -> head));
	}
	free(*root);
	
	return;
}



// recursively looks for place to add a node
void	recursiveAddCword(struct nodeCword **current, int index, struct treeRootCword *root){
	// if NULL, then found location for new node
	if(*current == NULL){
		struct nodeCword *newNode;
		allocNodeCword(&newNode, index, root);
		*current = newNode;
		root -> size ++;
	}

	// move left if new value less
	else if( (*current) -> index > index){
		recursiveAddCword( &((*current) -> left), index, root);
	}
	
	// move right if greater
	else if( (*current) -> index < index){
		recursiveAddCword( &((*current) -> right), index, root);
	}
	
	// if equal then no need to add
	return;
}



// adds a value to the tree
void	addValueCword(struct treeRootCword *root, int index){
	
	// if tree uninitialised, then return error
	if(root == NULL){
		printf("*** ERROR ***\nUse of uninitialised tree\n\n");
		exit(1);
	}
	
	// else add from the head
	else{
		recursiveAddCword( &(root -> head), index, root);
	}

	return;
}



// recursively finds the height
int 	recursiveHeightCword(struct nodeCword *current){
	if(current == NULL){
		return 0;
	}
	int left = recursiveHeightCword(current -> left);
	int right = recursiveHeightCword(current -> right);
	
	if(left > right){
		return left + 1;
	}
	return right + 1;
}



// finds the height of a tree
int	findHeightCword(struct treeRootCword *root){

	if(root == NULL){
		return 0;
	}

	return recursiveHeightCword(root->head);
}



// shifts clockwise
void shiftClockCword(struct nodeCword **topNodePointer){
	// creates holder for current top node
	struct nodeCword *current = *topNodePointer;
	if(current == NULL){
		return;
	}

	// sets the top node to the node on the left
	*topNodePointer = current -> left;

	// sets the old top's left to the the new top's right
	current -> left = (*topNodePointer) -> right;

	// sets the new top's left to be the old tp
	(*topNodePointer) -> right = current;

	// marks as potentially unbalanced
	current -> bal = 0;
	(*topNodePointer) -> bal = 0;
	
	return;
}



// shifts counter clockwise from a given node
void	shiftCounterCword(struct nodeCword **topNodePointer){
	// creates holder for current top node
	struct nodeCword *current = *topNodePointer;
	if(current == NULL){
		return;
	}

	// sets the top node to the node on the right
	*topNodePointer = current -> right;

	// sets the old top's right to the new top's left
	current -> right = (*topNodePointer) -> left;

	// sets the new top's left to be the old top
	(*topNodePointer) -> left = current;

	// marks as potentially unbalanced
	current -> bal = 0;
	(*topNodePointer) -> bal = 0;

	return;
}



// marks every node affected by a shift as unbalanced
void	recursiveMarkUnbalancedCword(struct nodeCword *current, int index){
	if(current == NULL){
		return;
	}
	current -> bal = 0;
	if(index < current -> index){
		recursiveMarkUnbalancedCword(current -> left, index);
	}
	else if(index > current -> index){
		recursiveMarkUnbalancedCword(current -> right, index);
	}
	return;
}



// recursively finds the size
void	recursiveSizeCword(struct nodeCword *current, int *size){
	if(current == NULL){
		return;
	}
	(*size) ++;
	recursiveSizeCword(current -> left, size);
	recursiveSizeCword(current -> right, size);

	return;
}



// finds the size of a tree
int	findSizeCword(struct treeRootCword *root){
	int size = 0;
	recursiveSizeCword(root -> head, &size);
	return size;
}



// recursively rebalances the tree
void	recursiveRebalanceCword(struct nodeCword **topNodePointer, struct treeRootCword *root){
	struct nodeCword *current = *topNodePointer;
	
	// if reach NULL or the node is marked as balanced then skip
	if(current == NULL || current -> bal == 1){
		return;
	} 

	int left = recursiveHeightCword(current -> left);
	int right = recursiveHeightCword(current -> right);
	
	// if left is longer then shift clockwise
	while(left - right > 1){
		shiftClockCword(topNodePointer);
		recursiveMarkUnbalancedCword(root -> head, (*topNodePointer) -> index);
		left = recursiveHeightCword((*topNodePointer) -> left);
		right = recursiveHeightCword((*topNodePointer) -> right);
	}

	// if right is longer then shift counterclockwise
	while(right - left > 1){
		shiftCounterCword(topNodePointer);
		recursiveMarkUnbalancedCword(root -> head, (*topNodePointer) -> index);
		left = recursiveHeightCword((*topNodePointer) -> left);
		right = recursiveHeightCword((*topNodePointer) -> right);
	}

	// mark as balanced
	(*topNodePointer) -> bal = 1;

	// recursively balance rest of tree
	recursiveRebalanceCword(&((*topNodePointer) -> left), root);
	recursiveRebalanceCword(&((*topNodePointer) -> right), root);

	return;
}



// rebalances a tree
void	rebalanceTreeCword(struct treeRootCword *root){
	if(root == NULL){
		return;
	}

	// finds current height and targets optimal height (log_2(size))
	int height = findHeightCword(root);
	int target = log(root -> size) / log(2) + 1;
	int counter = 0, newHeight;

	// loops until it is optimal or a max loop of 100 per height change
	while(target < height && counter < 100){
		counter ++;
		recursiveRebalanceCword(&(root -> head), root);
		newHeight = findHeightCword(root);
		if(newHeight != height){
			counter = 0;
			height = newHeight;
			//printf("new height is %d\n", newHeight);
		}
	}

	return;
}



// recursively prints
void	recursivePrintCword(struct nodeCword *current){
	if(current == NULL){
		return;
	}
	printf("%d\n", current -> index);
	recursivePrintCword(current -> left);
	recursivePrintCword(current -> right);

	return;
}



// prints a tree of less than 10 size
void	printTreeCword(struct treeRootCword *root){
	if(root == NULL || root -> size > 10){
		return;
	}
	else{
		recursivePrintCword(root->head);
		printf("\n");
	}
	
	return;
}



struct nodeCword * recursiveSearchCword(struct nodeCword *current, int index){
	if(current == NULL){
		return NULL;
	}
	else if(index < current -> index){
		return	recursiveSearchCword(current -> left, index);
	}
	else if(index > current -> index){
		return	recursiveSearchCword(current -> right, index);
	}
	else{
		return current;
	}
	
	return NULL;
}


struct nodeCword * searchTreeCword(struct treeRootCword *root, int index){
	if(root == NULL){
		return NULL;
	}

	return recursiveSearchCword(root -> head, index);
}



void	recursiveCountCword(struct nodeCword *current, int *total, struct treeRootCword *root){
	if(current == NULL){
		return;
	}

	int i;
	for(i = 0; i < root -> noTopics; i++){
		(*total) += current -> topicArray[i];
	}
	recursiveCountCword(current -> left, total, root);
	recursiveCountCword(current -> right, total, root);

	return;
}



int	countTopicsCword(struct treeRootCword *root){
	int total = 0;
	recursiveCountCword(root -> head, &total, root);
	
	return total;
}



#ifdef MPI
void	recursiveSendForPrintCword(struct nodeCword *current, int noTopics, MPI_Comm mainComm){
	if(current == NULL){
		return;
	}

	MPI_Send(&(current -> index), 1, MPI_INT, 0, 7, mainComm);
	MPI_Send(current -> topicArray, noTopics, MPI_INT, 0, 7, mainComm);

	recursiveSendForPrintCword(current -> left, noTopics, mainComm);
	recursiveSendForPrintCword(current -> right, noTopics, mainComm);

	return;
}


void	sendForPrintCword(struct treeRootCword *root, MPI_Comm mainComm){
	recursiveSendForPrintCword(root -> head, root -> noTopics, mainComm);	

	// send finish signal
	int end = -1;
	MPI_Send(&end, 1, MPI_INT, 0, 7, mainComm);
	
	return;
}

#endif
//OLD FUNCTIONS
/*
void	simpleRecursiveRebalance(struct nodeCword **topNodePointer){
	struct nodeCword *current = *topNodePointer;
	
	if(current == NULL){
		return;
	} 
	int left = recursiveHeightCword(current -> left);
	int right = recursiveHeightCword(current -> right);
	
	while(left - right > 1){
		shiftClockCword(topNodePointer);
		left = recursiveHeightCword((*topNodePointer) -> left);
		right = recursiveHeightCword((*topNodePointer) -> right);
	}
	while(right - left > 1){
		shiftCounterCword(topNodePointer);
		left = recursiveHeightCword((*topNodePointer) -> left);
		right = recursiveHeightCword((*topNodePointer) -> right);
	}

	simpleRecursiveRebalance(&((*topNodePointer) -> left));
	simpleRecursiveRebalance(&((*topNodePointer) -> right));

	return;
}



void	simpleRebalanceTree(struct treeRootCword *root){
	if(root == NULL){
		return;
	}
	simpleRecursiveRebalance(&(root -> head));

	return;
}


void	recursiveMarkUnbalancedCword(struct nodeCword *current){
	if(current == NULL){
		return;
	}
	current -> bal = 0;
	recursiveMarkUnbalancedCword(current -> left);
	recursiveMarkUnbalancedCword(current -> right);

	return;
}

*/
