CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -lm

objects = main.o 
executable = createCorpus.out

$(executable): $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

clean :
	$(RM) *.o $(executable)

test: $(executable)
	./$(executable)
