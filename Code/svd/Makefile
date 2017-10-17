CC = gcc
CFLAGS = -g -W -Wall
LDFLAGS = -lm

objects = main.o ../general/matrix.o jacobiEigen.o svd.o ../general/extra.o
executable = svd.out

$(executable): $(objects)
	$(CC) $(LDFLAGS) $(objects) -o $@

.PHONY: clean test

clean :
	$(RM) $(objects) $(executable) 

test: $(executable)
	./$(executable) -v
