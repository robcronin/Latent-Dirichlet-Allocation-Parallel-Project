SP = #scorep --user --static
CC = $(SP)  mpicc
CFLAGS = -DMPI -DACCTIME #-DDEBUG -g -W -Wall
LDFLAGS = -lm -lpthread

objects = 	main.o \
		plda+.o \
		gibbsplda+.o \
		initTopics.o \
		pthreadFunctions.o \
		../Gibbs-LDA-MPI/gibbsLDA-MPI.o \
		../Gibbs-LDA2/gibbsLDA2.o \
		../tfidf/tfidf.o \
		../general/matrix.o \
		../general/extra.o \
		../binaryTree/binaryTreeCword.o \
		../binaryTree/binaryTreeInv.o \
		../general/quicksort.o \
#		../linkedList/linkedList.o \

executable = gibbsPLDA+.out

$(executable): $(objects)
	rm -f main.o
	$(CC) $(CFLAGS) -c -o main.o main.c
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

threads: 
	rm -f main.o
	make -f Makefile.Threads

clean :
	$(RM) $(objects) $(executable) gibbsPLDA+Threads.out

test: $(executable)
	mpirun -n 4 ./$(executable) -s 1500547989 -p -i 30 -m 10

test2: $(executable)
	mpirun -n 4 ./$(executable) -s 1500547989 -p -i 100

rtest: $(executable)
	mpirun -n 4 ./$(executable) -s 1500547989 -pr -i 100

mtest: $(executable)
	mpirun -n 32 ./$(executable) -s 1500547989 -pr -i 1000 -d ../../bigGutenberg

