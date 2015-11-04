# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 
LFLAGS = $(CFLAGS)
LIBS   = -lm
CPPFLAGS = -DCOLS -fopenmp


OPENMP = partdiff-openmp.o askparams.o 


# Rule to create *.o from *.c
.c.o:
	$(CC) $(CPPFLAGS) -c $(CFLAGS) $*.c


partdiff-openmp-zeile: $(OPENMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OPENMP) $(LIBS)

