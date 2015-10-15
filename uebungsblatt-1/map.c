#include <stdio.h>

// Definieren Sie ein enum cardd
//N -> 0001
//E -> 0010
//S -> 0100
//W -> 1000
typedef enum {
	N = 1, 
	E = 2, 
	S = 4, 
	W = 8
} cardd;

//array 3x3
int map[3][3];


// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
	if ((x >= 0 && x <= 2)&&(y >= 0 && y <= 2))
	{
		map[x][y] = dir;
	}
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
	int x, y;
	
	for (x = 0; x < 3; x++) 
	{
		for (y = 0; y < 3; y++)
		{
			switch(map[x][y]) 
			{
				case 1: printf("N"); break;
				case 2: printf(" E"); break;
				case 4: printf("S"); break;
				case 8: printf("W "); break;
				
				case 9: printf("NW"); break;  	// 0001|1000 = 1001 -> 9
				case 12: printf("SW"); break;	// 1000|0100 = 1100 -> 12
				
				case 3: printf("NE"); break;  // 0001|0010 = 0011 -> 3 (brauchen wir eigentlich nicht, wird überschrieben)
				case 5: printf("NE"); break;  	// wegen 	set_dir(0, 2, N|S); an (0, 2) muss "NE" stehen. 0001|0100 = 0101 -> 5

				case 6: printf("SE"); break;	// 0100|0010 = 0110 -> 6 (brauchen wir eigentlich nicht, wird überschrieben)
				case 10: printf("SE"); break; 	// wegen 	set_dir(2, 2, E|W); an (2, 2) muss "SE" stehen. 0010|1000 = 1010 -> 10
				
				default: printf("0"); break;
			}
			
			if (y < 2){	
				printf("  ");
			}
		}
		printf("\n");
	}
}

int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);		
	set_dir(0, 2, N|S);		//hier wird (0, 2) nocheinmal überschrieben
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W); 	//s.o.

	show_map();

	return 0;
}
