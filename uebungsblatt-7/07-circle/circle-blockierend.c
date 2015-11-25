#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <mpi.h>


#define MASTER 0
#define LAST (nprocs - 1)

#define BEFORE 0
#define AFTER 1

int rank;
int nprocs;

MPI_Status status;

int*
init (int N)
{
	int* buf = malloc(sizeof(int) * N);

	//first place for length of array
	buf[0] = N;

	//+rank, because sometimes prozesses on the same node had the same seed.
	srand(time(NULL)+rank);

	//rest for rand numbers
	for (int i = 1; i < N; i++)
	{
		buf[i] = rand() % 25; //do not modify %25
	}

	return buf;
}


int*
circle (int* buf, int count)
{ 
	int next, previous;

	int new_buf[count];

	//only one process
	if(nprocs == 1)
	{
		return buf;
	}

	//Master process
	else if(rank == MASTER)
	{
		next = rank + 1;
		previous = LAST;
	}

	//last process
	else if(rank == LAST)
	{
		next = MASTER;
		previous = rank - 1;
	}

  	//the other processes
	else
	{
	    next = rank + 1;
	    previous = rank - 1;
 	}
	
	if(rank%2 == 0)
	{
		MPI_Send(buf, count, MPI_INT, next, 0, MPI_COMM_WORLD);
		MPI_Recv(&new_buf, count, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
	}
	else 
	{
		MPI_Recv(&new_buf, count, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
		MPI_Send(buf, count, MPI_INT, next, 0, MPI_COMM_WORLD);		
	}

	memcpy(buf, new_buf, (sizeof(int)*count));

  	return buf;

}


//time, 0 before, 1 after
void 
ausgabe(int* buf, int buf_length, int time)
{
	//Ausgabe "Before" von Master
	if(rank == MASTER)
	{ 
		if(time == BEFORE)
		{
			printf("\nBefore\n");
		}
		else if (time == AFTER)
		{
			printf("\nAfter\n");
		}

		//printf for MASTER, an buf[0] steht die Länge des Arrays
		for (int j = 1; j < buf[0]; j++)
		{
			printf ("rank %d: %d\n", MASTER, buf[j]);
			fflush(stdout);
		}

		int print_buf[buf_length];

		//printf für die restlichen Prozesse, an print_buf[0] steht die Länge des Arrays
		for (int i = 1; i < nprocs; i++)
		{ 
		    //MPI_RECV(buf,count,datatype,source,tag,comm,status)
		    MPI_Recv(&print_buf, buf_length, MPI_INT, i, 0,  MPI_COMM_WORLD, &status);

		    for (int j = 1; j < print_buf[0]; j++)
		    {
		    	printf ("rank %d: %d\n", i, print_buf[j]);
		    	fflush(stdout);
		    }
		}  
	}
	else
	{
		MPI_Send(buf, buf_length, MPI_INT, MASTER, 0,  MPI_COMM_WORLD);
	}
}


int
main (int argc, char** argv)
{
	MPI_Init(NULL, NULL);

	char arg[256];
  	int N;
  	int* buf;

  
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 2)
	{
		printf("Arguments error \n");
		return EXIT_FAILURE;
	}

	sscanf(argv[1], "%s", arg);

	//array length
	N = atoi(arg);

	if (N < nprocs)
	{
		printf("Number of processes should be at least N\n");
		return EXIT_FAILURE;
	}

	//aufteilen, zuerst wird N - (N % nprocs)  gleich verteilt, dann wird der rest auf die Prozesse verteilt
	int rest = N % nprocs;
	N = (N - rest)/nprocs;

	//Maximale Länge aller benutzten Arrays
	int maxbuf_length = N + 2;

	//rest wird verteilt
	if(rank < rest)
	{
		 N = N + 1;
	}

	//N Anzahl der zufälligen Werte
	//+1 Feld für die Länge des Arrays
	buf = init(N+1);

	int first_value;
	//send first value of MASTER to LAST
	if(rank == MASTER)
	{
		first_value = buf[1];
		MPI_Send(&first_value, 1, MPI_INT, LAST, 0, MPI_COMM_WORLD);
	}
	if(rank == LAST)
	{
		MPI_Recv(&first_value, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
	}

	//Ausgabe "Before"
	ausgabe(buf, maxbuf_length, BEFORE);


	int termination = 1;

	//Start circling
	while(termination)
	{	
		buf = circle(buf, maxbuf_length);

		//check for termination condition
		if(rank == LAST)
		{
		    if (buf[1] == first_value)
		    {
		        termination = 0;
		    }


		}
		//wait for all processes to circle
		MPI_Barrier(MPI_COMM_WORLD);
		//Broadcast termination status
		MPI_Bcast(&termination, 1, MPI_INT, LAST, MPI_COMM_WORLD);
	}

	//Ausgabe "After"
	ausgabe(buf, maxbuf_length, AFTER);

	free(buf);
	MPI_Finalize();

	return EXIT_SUCCESS;
}
