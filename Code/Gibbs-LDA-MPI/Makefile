CC = mpicc
CFLAGS = -W -Wall -g -DMPI -DACCTIME #-DDEBUG
LDFLAGS = -lm

objects = 	main.o \
		gibbsLDA-MPI.o \
		../Gibbs-LDA2/gibbsLDA2.o \
		../tfidf/tfidf.o \
		../general/matrix.o \
		../general/extra.o \
		../general/quicksort.o \

executable = gibbsLDA-MPI.out

$(executable): $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

clean :
	$(RM) $(objects) $(executable)

test: $(executable)
	mpirun -n 4 ./$(executable) -s 1500547989 -p

rtest: $(executable)
	mpirun -n 4 ./$(executable) -s 1500547989 -pr -i 100

