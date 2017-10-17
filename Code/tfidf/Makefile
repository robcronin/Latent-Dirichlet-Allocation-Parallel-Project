CC = gcc
CFLAGS = -g -W -Wall
LDFLAGS = -lm

objects = main.o \
		tfidf.o \
		../general/matrix.o \
		../general/quicksort.o \

executable = tfidf.out

$(executable): $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

clean :
	$(RM) $(objects) $(executable)

test: $(executable)
	./$(executable)
