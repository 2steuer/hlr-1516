all: pointer map

pointer: pointer.o
	gcc -o pointer pointer.o
	
map: map.o
	gcc -o map map.o
		
%.o: %.c
	gcc -c $<

clean:
	rm -f map
	rm -f pointer
	rm -f *.o
	rm -f *~