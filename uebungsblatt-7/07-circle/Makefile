CC = mpicc
# CC = vtcc -vt:cc mpicc #with vampir trace

all: circle circle-blockierend 

circle: circle.c Makefile
	$(CC) -std=c99 -ggdb -Wall -Wextra circle.c -o circle

circle-blockierend: circle-blockierend.c Makefile
	$(CC) -std=c99 -ggdb -Wall -Wextra circle-blockierend.c -o circle-blockierend

clean:
	rm -f circle circle-blockierend
