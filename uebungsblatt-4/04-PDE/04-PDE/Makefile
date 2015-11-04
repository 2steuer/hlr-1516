# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

SEQ = partdiff-seq.o askparams.o 
OPENMP = partdiff-openmp.o askparams.o 


# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-seq partdiff-openmp-zeile partdiff-openmp-spalte partdiff-openmp-element

partdiff-seq: $(OPENMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OPENMP) $(LIBS)

clean:
	$(RM) *.o *~ partdiff-seq partdiff-openmp-zeile partdiff-openmp-spalte partdiff-openmp-element

partdiff-seq.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

partdiff-openmp-zeile:LFLAGS += -DCOLS -fopenmp
partdiff-openmp-zeile: 	$(OPENMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OPENMP) $(LIBS)
				

partdiff-openmp-spalte: LFLAGS += -DROWS -fopenmp
partdiff-openmp-spalte: 	$(OPENMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OPENMP) $(LIBS)

partdiff-openmp-element: LFLAGS += -DELEMENTS -fopenmp
partdiff-openmp-element: 	$(OPENMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OPENMP) $(LIBS)

