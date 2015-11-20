# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99  -pedantic -Wall -Wextra -g
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJ = timempi.o 
OBJ2 = timempi2.o

all: timempi timempi2

timempi: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

timempi2: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm timempi
	rm timempi2
