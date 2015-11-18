# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99  -pedantic -Wall -Wextra -g
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJ = main.o 


main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


clean:
	rm -f *.o
	rm -f *~
	rm main
