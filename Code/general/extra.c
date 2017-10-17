#include "extra.h"

// returns min integer
int  minInt(int a, int b){
	if(a < b){
		return a;
	}
		
	return b;
}



// returns min double
double  minDouble(double a, double b){
	if(a < b){
		return a;
	}
		
	return b;
}



// returns max integer
int  maxInt(int a, int b){
	if(a > b){
		return a;
	}
		
	return b;
}



// returns max double
double  maxDouble(double a, double b){
	if(a > b){
		return a;
	}
		
	return b;
}



// prints time taken
void	printTime(struct timeval t1, struct timeval t2, char *message){
	double elapsedTime = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec)/1000000.0;
	if(message == NULL){
		printf("Time Taken:\t\t%.6es\n", elapsedTime);
	}
	else{
		int size = strlen(message);
		if(size <= 7){
			printf("Time Taken for %s:\t\t\t\t%.6es\n", message, elapsedTime);
		}
		else if(size <= 15){
			printf("Time Taken for %s:\t\t\t%.6es\n", message, elapsedTime);
		}
		else if(size <= 23){
			printf("Time Taken for %s:\t\t%.6es\n", message, elapsedTime);
		}
		else{
			printf("Time Taken for %s:\t%.6es\n", message, elapsedTime);
		}
	}

	return;

}



// prints numbered debug flags
void	flag(int num){
	printf("FLAG %d\n", num);

	return;
}



#ifdef MPI
// prints flags for MPI
void	mflag(int num, MPI_Comm comm){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank == 0){
		printf("%02d: MPI FLAG %d\n", rank, num);
	}

	return;
}



// Starts timing and puts up a barrier before if requested
void	startTimeMPI(struct timeval *t1, MPI_Comm comm){
	#ifdef ACCTIME
	MPI_Barrier(comm);
	#endif
	
	gettimeofday(t1, NULL);

	return;

}



void	endTimeMPI(struct timeval t1, char *message, int rank, MPI_Comm comm){
	#ifdef ACCTIME
	MPI_Barrier(comm);
	#endif
	
	struct timeval t2;
	gettimeofday(&t2, NULL);

	#ifdef ACCTIME
	MPI_Barrier(comm);
	#endif
	if(rank == 0){
		printTime(t1, t2, message);	
	}

	return;
}
#endif
