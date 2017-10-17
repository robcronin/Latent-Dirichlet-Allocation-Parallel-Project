/*	Various functions that may be reusable in various functions
 *		min/maxes
 *		debug flags
 *		time information
 */


#ifndef EXTRA
#define EXTRA


#include <sys/time.h>
#include <stdio.h>
#include <string.h>


int 	minInt(int a, int b);
double 	minDouble(double a, double b);
int 	maxInt(int a, int b);
double 	maxDouble(double a, double b);
void	printTime(struct timeval t1, struct timeval t2, char *message);
void	flag(int num);


#ifdef MPI
	#include <mpi.h>
	void	mflag(int num, MPI_Comm comm);
	void	startTimeMPI(struct timeval *t1, MPI_Comm comm);
	void	endTimeMPI(struct timeval t1, char *message, int rank, MPI_Comm comm);
#endif


#endif
