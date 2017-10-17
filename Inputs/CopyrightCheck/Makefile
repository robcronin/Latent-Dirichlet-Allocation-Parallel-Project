CC = gcc
CFLAGS = 
LDFLAGS =

objects = copyright.o
executable = copyrightChecker

$(executable): $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

.PHONY: clean test

clean :
	$(RM) *.o $(executable)

test: $(executable)
	./$(executable)
