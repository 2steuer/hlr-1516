# Common definitions
CC = gcc


	# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2
LDFLAGS = $(CFLAGS)
LIBS   = -lm
OPENMP =  partdiff-openmp-spalte partdiff-openmp-zeile partdiff-openmp-element partdiff-seq


all: $(OPENMP) 


partdiff-openmp-spalte: CFLAGS += -DCOLS  -fopenmp
partdiff-openmp-spalten-obj.o: CFLAGS += -DCOLS  -fopenmp

partdiff-openmp-zeile: CFLAGS += -DROWS  -fopenmp
partdiff-openmp-zeile-obj.o: CFLAGS += -DROWS  -fopenmp

partdiff-openmp-element: CFLAGS += -DELEMENTS -fopenmp
partdiff-openmp-element-obj.o: CFLAGS += -DELEMENTS -fopenmp

partdiff-seq.o: partdiff-openmp.c
	$(CC) $(CFLAGS) -c $< -o $@

%-obj.o: partdiff-openmp.c
	$(CC) $(CFLAGS) -c $< -o $@
	
askparams.o: askparams.c
	$(CC) $(CFLAGS) -c $< -o $@


$(OPENMP): % : %-obj.o askparams.o
	$(CC) $(LDFLAGS) -o $@ $@-obj.o askparams.o $(LIBS)

clean:
	$(RM) $(OPENMP) *.o *~
