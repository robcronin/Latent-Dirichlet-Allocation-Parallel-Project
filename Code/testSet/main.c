#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
//#include "probDist.h"

int main(){
	int i, j;

	srand48(time(NULL));
	int 	noDocs = 100;
	int	noTopics = 10;
	int	wordsPerDoc = 10000;
	int	delTopics = 5;
	char 	*directory = "testCorpus/";

	// create random distribution for the topics
	double 	*probs = malloc(noTopics * sizeof(double));
	double 	*docProbs = malloc(noTopics * sizeof(double));
	double probTotal = 0;
	for(i = 0; i < noTopics; i++){
		probs[i] = drand48();
		probTotal += probs[i];
	}	
	for(i = 0; i < noTopics; i++){
		probs[i] /= probTotal;
		printf("Topic %d:\t%lf\n", i, probs[i]);
	}


	int 	*counters = malloc(noTopics * sizeof(int));
	for(i = 0; i < noTopics; i++){
		counters[i] = 0;
	}

	FILE	*docptr;
	int d;
	for(d = 0; d < noDocs; d++){
		// create file
		char *filename = malloc(1024);
		sprintf(filename, "%sdoc%d.txt", directory, d);
		docptr = fopen(filename, "w");

		fprintf(docptr, "\n\n\n");

		// find topic probs
		for(i = 0; i < noTopics; i++){
			docProbs[i] = probs[i];
		}

		// delete some topics for each doc
		for(i = 0; i < delTopics; i++){
			int del = noTopics*drand48();
			docProbs[del] = 0;
		}

		// fix probs
		probTotal = 0;
		for(i = 0; i < noTopics; i++){
			probTotal += docProbs[i];
		}	
		for(i = 0; i < noTopics; i++){
			docProbs[i] /= probTotal;
			//printf("Doc %d, Topic %d, Prob %lf\n", d, i, docProbs[i]);
		}

		int num;
		int topic;
		double sampler, counter;
		for(j = 0; j < wordsPerDoc; j++){
			// find topic
			topic = 0;
			counter = docProbs[0];
			sampler = drand48();
			//printf("sampler is %lf, counter %lf\n", sampler, counter);
			while(counter < sampler){// && topic < noTopics - 1){
				topic ++;
				counter += docProbs[topic];
			}
			//printf("Topic is %d\n", topic);

			counters[topic] ++;

			num = 10*drand48() + topic*10;
			fprintf(docptr, "%d\n", num);
		}
		fclose(docptr);
		free(filename);
	}	

	printf("\n");
	double act;
	for(i = 0; i < noTopics; i++){
		act = (double) counters[i] / (wordsPerDoc * noDocs);
		printf("Topic %d Act:\t%lf\n", i, act);
	}
	
	return 0;
}


