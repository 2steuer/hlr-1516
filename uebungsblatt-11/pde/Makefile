# Common definitions
#CC = gcc
MPI = vtcc -vt:cc mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS-par = partdiff-par.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(MPI) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-par

partdiff-par: $(OBJS-par) Makefile
	$(MPI) $(LFLAGS) -o $@ $(OBJS-par) $(LIBS)

clean:
	$(RM) *.o *~ partdiff-par

