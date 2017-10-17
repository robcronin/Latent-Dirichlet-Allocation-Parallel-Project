CC = gcc
CFLAGS = -W -Wall -g
LDFLAGS = -lm

objects = 	main.o \
		gibbsLDA2.o \
		../tfidf/tfidf.o \
		../general/matrix.o \
		../general/extra.o \
		../general/quicksort.o \

executable = gibbsLDA2.out

$(executable): $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

clean :
	$(RM) $(objects) $(executable)

test: $(executable)
	./$(executable) -d ../../Inputs/Gutenberg -s 1500547989
