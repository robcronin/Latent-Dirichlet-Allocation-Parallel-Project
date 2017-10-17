#include "binaryTreeInv.h"

// creates a tree
// CHANGE ARGS
void 	createTreeInv(struct treeRootInv **root, int noFiles){
	*root = malloc(sizeof(struct treeRootInv));
	(*root) -> head = NULL;
	(*root) -> size = 0;

	// CHANGE HERE
	(*root) -> noFiles = noFiles;

	return;
}


// allocates node and initialises values
void	allocNodeInv(struct nodeInv **newNode, int index, struct treeRootInv *root){
	*newNode = malloc(sizeof(struct nodeInv));
	(*newNode) -> index = index;
	(*newNode) -> left = NULL;
	(*newNode) -> right = NULL;
	(*newNode) -> bal = 0;

	// CHANGE HERE
	(*newNode) -> invIndex = malloc(root -> noFiles * sizeof(int));
	int i;
	for(i = 0; i < root -> noFiles; i++){
		(*newNode) -> invIndex[i] = -1;
	}

	return;
}


// frees an individual node	
void	freeNodeInv(struct nodeInv **node){
	// CHANGE HERE
	free( (*node) -> invIndex);

	free(*node);
	
	return;
}


// recursively frees nodes
void	recursiveFreeInv(struct nodeInv **node){
	if(*node == NULL){
		return;
	}
	recursiveFreeInv( &((*node) -> left));
	recursiveFreeInv( &((*node) -> right));
	freeNodeInv(node);

	return;
}



// frees the Tree
void	freeTreeInv(struct treeRootInv **root){
	if( *root == NULL){
		return;
	}
	else{
		recursiveFreeInv(&( (*root) -> head));
	}
	free(*root);
	
	return;
}


// recursively looks for place to add a node
void	recursiveAddInv(struct nodeInv **current, int index, struct treeRootInv *root, int doc, int location, int occ){
	// if NULL, then found location for new node
	if(*current == NULL){
		struct nodeInv *newNode;
		allocNodeInv(&newNode, index, root);
		newNode -> occurrences = occ;
		newNode -> invIndex[doc] = location;
		*current = newNode;
		root -> size ++;
	}

	// move left if new value less
	else if( (*current) -> index > index){
		recursiveAddInv( &((*current) -> left), index, root, doc, location, occ);
	}
	
	// move right if greater
	else if( (*current) -> index < index){
		recursiveAddInv( &((*current) -> right), index, root, doc, location, occ);
	}
	
	// if equal then no need to add
	else{
		(*current) -> occurrences += occ;
		(*current) -> invIndex[doc] = location;
	}
	return;
}



// adds a value to the tree
void	addValueInv(struct treeRootInv *root, int index, int doc, int location, int occ){
	
	// if tree uninitialised, then return error
	if(root == NULL){
		printf("*** ERROR ***\nUse of uninitialised tree\n\n");
		exit(1);
	}
	
	// else add from the head
	else{
		recursiveAddInv( &(root -> head), index, root, doc, location, occ);
	}

	return;
}



// recursively finds the height
int 	recursiveHeightInv(struct nodeInv *current){
	if(current == NULL){
		return 0;
	}
	int left = recursiveHeightInv(current -> left);
	int right = recursiveHeightInv(current -> right);
	
	if(left > right){
		return left + 1;
	}
	return right + 1;
}



// finds the height of a tree
int	findHeightInv(struct treeRootInv *root){

	if(root == NULL){
		return 0;
	}

	return recursiveHeightInv(root->head);
}



// shifts clockwise
void shiftClockInv(struct nodeInv **topNodePointer){
	// creates holder for current top node
	struct nodeInv *current = *topNodePointer;
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
void	shiftCounterInv(struct nodeInv **topNodePointer){
	// creates holder for current top node
	struct nodeInv *current = *topNodePointer;
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
void	recursiveMarkUnbalancedInv(struct nodeInv *current, int index){
	if(current == NULL){
		return;
	}
	current -> bal = 0;
	if(index < current -> index){
		recursiveMarkUnbalancedInv(current -> left, index);
	}
	else if(index > current -> index){
		recursiveMarkUnbalancedInv(current -> right, index);
	}
	return;
}



// recursively finds the size
void	recursiveSizeInv(struct nodeInv *current, int *size){
	if(current == NULL){
		return;
	}
	(*size) ++;
	recursiveSizeInv(current -> left, size);
	recursiveSizeInv(current -> right, size);

	return;
}



// finds the size of a tree
int	findSizeInv(struct treeRootInv *root){
	int size = 0;
	recursiveSizeInv(root -> head, &size);
	return size;
}



// recursively rebalances the tree
void	recursiveRebalanceInv(struct nodeInv **topNodePointer, struct treeRootInv *root){
	struct nodeInv *current = *topNodePointer;
	
	// if reach NULL or the node is marked as balanced then skip
	if(current == NULL || current -> bal == 1){
		return;
	} 

	int left = recursiveHeightInv(current -> left);
	int right = recursiveHeightInv(current -> right);
	
	// if left is longer then shift clockwise
	while(left - right > 1){
		shiftClockInv(topNodePointer);
		recursiveMarkUnbalancedInv(root -> head, (*topNodePointer) -> index);
		left = recursiveHeightInv((*topNodePointer) -> left);
		right = recursiveHeightInv((*topNodePointer) -> right);
	}

	// if right is longer then shift counterclockwise
	while(right - left > 1){
		shiftCounterInv(topNodePointer);
		recursiveMarkUnbalancedInv(root -> head, (*topNodePointer) -> index);
		left = recursiveHeightInv((*topNodePointer) -> left);
		right = recursiveHeightInv((*topNodePointer) -> right);
	}

	// mark as balanced
	(*topNodePointer) -> bal = 1;

	// recursively balance rest of tree
	recursiveRebalanceInv(&((*topNodePointer) -> left), root);
	recursiveRebalanceInv(&((*topNodePointer) -> right), root);

	return;
}



// rebalances a tree
void	rebalanceTreeInv(struct treeRootInv *root){
	if(root == NULL){
		return;
	}

	// finds current height and targets optimal height (log_2(size))
	int height = findHeightInv(root);
	int target = log(root -> size) / log(2) + 1;
	int counter = 0, newHeight;

	// loops until it is optimal or a max loop of 100 per height change
	while(target < height && counter < 100){
		counter ++;
		recursiveRebalanceInv(&(root -> head), root);
		newHeight = findHeightInv(root);
		if(newHeight != height){
			counter = 0;
			height = newHeight;
			//printf("new height is %d\n", newHeight);
		}
	}

	return;
}



// recursively prints
void	recursivePrintInv(struct nodeInv *current){
	if(current == NULL){
		return;
	}
	printf("%d\n", current -> index);
	recursivePrintInv(current -> left);
	recursivePrintInv(current -> right);

	return;
}



// prints a tree of less than 10 size
void	printTreeInv(struct treeRootInv *root){
	if(root == NULL || root -> size > 10){
		return;
	}
	else{
		recursivePrintInv(root->head);
		printf("\n");
	}
	
	return;
}


void	recursiveFillOcc(int **occMat, struct nodeInv *current, int *counter, int *wordLocation){
	if(current == NULL){
		return;
	}
	occMat[*counter][0] = current -> index;
	occMat[*counter][1] = current -> occurrences;
	occMat[*counter][2] = wordLocation[current -> index];
	(*counter) ++;	

	recursiveFillOcc(occMat, current -> left, counter, wordLocation);
	recursiveFillOcc(occMat, current -> right, counter, wordLocation);

	return;
}



void	fillOccMat(int **occMat, struct treeRootInv *root, int *wordLocation){
	int counter = 0;
	recursiveFillOcc(occMat, root -> head, &counter, wordLocation);

	#ifdef DEBUG
	printf("counter is %d, size is %d\n", counter, root -> size);
	#endif
		
	return;
}


void	recursiveFillOccMatrix(int **occMat, struct nodeInv *current, int *counter, int **wordLocationMat){
	if(current == NULL){
		return;
	}
	occMat[*counter][0] = current -> index;
	occMat[*counter][1] = current -> occurrences;
	occMat[*counter][2] = wordLocationMat[current -> index][0];
	(*counter) ++;	

	recursiveFillOccMatrix(occMat, current -> left, counter, wordLocationMat);
	recursiveFillOccMatrix(occMat, current -> right, counter, wordLocationMat);

	return;
}



void	fillOccMatMatrix(int **occMat, struct treeRootInv *root, int **wordLocationMat){
	int counter = 0;
	recursiveFillOccMatrix(occMat, root -> head, &counter, wordLocationMat);

	#ifdef DEBUG
	printf("counter is %d, size is %d\n", counter, root -> size);
	#endif
		
	return;
}


struct nodeInv * recursiveSearchInv(struct nodeInv *current, int index){
	if(current == NULL){
		return NULL;
	}
	else if(index < current -> index){
		return	recursiveSearchInv(current -> left, index);
	}
	else if(index > current -> index){
		return	recursiveSearchInv(current -> right, index);
	}
	else{
		return current;
	}
	
	return NULL;
}


struct nodeInv * searchTreeInv(struct treeRootInv *root, int index){
	if(root == NULL){
		return NULL;
	}

	return recursiveSearchInv(root -> head, index);
}



#ifdef MAIN
// main for testing
int main(){
	
	struct treeRootInv *root;
	createTreeInv(&root);
	srand(123);
	int i;
	for(i = 0; i < 100000; i++){
		int rando = rand()%10000000;
		addValueInv(root, rando);
	}
	printf("size of tree is %d or %d\n", root -> size, findSizeInv(root));
	printf("height of tree is %d\n", findHeightInv(root));

	printTreeInv(root);

	rebalanceTreeInv(root);

	printf("size of tree is %d or %d\n", root -> size, findSizeInv(root));
	printf("height of tree is %d\n", findHeightInv(root));
	
	printTreeInv(root);

	freeTreeInv(&root);
	
	return 0;
}
#endif

//OLD FUNCTIONS
/*
void	simpleRecursiveRebalance(struct nodeInv **topNodePointer){
	struct nodeInv *current = *topNodePointer;
	
	if(current == NULL){
		return;
	} 
	int left = recursiveHeightInv(current -> left);
	int right = recursiveHeightInv(current -> right);
	
	while(left - right > 1){
		shiftClockInv(topNodePointer);
		left = recursiveHeightInv((*topNodePointer) -> left);
		right = recursiveHeightInv((*topNodePointer) -> right);
	}
	while(right - left > 1){
		shiftCounterInv(topNodePointer);
		left = recursiveHeightInv((*topNodePointer) -> left);
		right = recursiveHeightInv((*topNodePointer) -> right);
	}

	simpleRecursiveRebalance(&((*topNodePointer) -> left));
	simpleRecursiveRebalance(&((*topNodePointer) -> right));

	return;
}



void	simpleRebalanceTree(struct treeRootInv *root){
	if(root == NULL){
		return;
	}
	simpleRecursiveRebalance(&(root -> head));

	return;
}


void	recursiveMarkUnbalancedInv(struct nodeInv *current){
	if(current == NULL){
		return;
	}
	current -> bal = 0;
	recursiveMarkUnbalancedInv(current -> left);
	recursiveMarkUnbalancedInv(current -> right);

	return;
}

*/
