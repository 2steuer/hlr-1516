/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#define MASTER 0
#define LAST (size-1)

#define NUM_COLS (N+1)

#define FIRST_ROW 1
#define LAST_ROW numberOfRows

#define NEXT_FROW (numberOfRows + 1)
#define PREVIOUS_LROW 0

#define NEXT_RANK (rank+1)
#define PREVIOUS_RANK (rank-1)

#define PRECISION_NOT_REACHED 0
#define PRECISION_REACHED 1
#define STOP_CALCULATION 2

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>

#include "partdiff-par.h"


struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  numberOfRows;	  /* number of rows per process 				    */
	uint64_t  num_matrices;   /* number of matrices                             */
	int 	  from;			  /* global number of first row */
	int 	  to; 			  /* global number of last row */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
		/* exit program */
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
	//+ 2 rows for calculation
	uint64_t const numberOfRows= arguments->numberOfRows + 2;

	arguments->M = allocateMemory(arguments->num_matrices * numberOfRows * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(numberOfRows * sizeof(double*));

		for (j = 0; j < numberOfRows; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * numberOfRows * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, int rank, int size)
{
	uint64_t g, i, j; /*  local variables for loops   */
	int iAdjusted;                               

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	//+ 2 rows for calculation
	uint64_t const numberOfRows= arguments->numberOfRows + 2;
	int const from = arguments->from;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < numberOfRows; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do), iAdjust hast to be */
	/* calculated in relation to the global row numbers.										 */
	if (options->inf_func == FUNC_F0)
	{
		if (rank == MASTER)
		{
			for (g = 0; g < arguments->num_matrices; g++)
			{
				for (i = 0; i < N; i++)
				{	
					Matrix[g][0][i] = 1.0 - (h * i);
				}

				Matrix[g][0][N] = 0.0;
			}
		}

		if (rank == LAST)
		{
			for (g = 0; g < arguments->num_matrices; g++)
			{
				for (i = 0; i < N; i++)
				{
					Matrix[g][numberOfRows-1][i] = h * i;
				}

				Matrix[g][numberOfRows-1][0] = 0.0;
			}
		}

		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 1; i < numberOfRows; i++)
			{	
				iAdjusted = (i + from-1);
				Matrix[g][i][0] = 1.0 - (h * iAdjusted);
				Matrix[g][i][N] = h * iAdjusted;
			}	
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation     										*/
/* ************************************************************************ */
static
void
calculateGS (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, int rank, int size)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;
	int const numberOfRows = arguments->numberOfRows;
	int const from = arguments->from;

	int term_iteration = options->term_iteration;
	int iteration_count = 1;
	int precisionStateCurrent = PRECISION_NOT_REACHED;
	int precisionStatePrevious = (rank == MASTER) ? PRECISION_REACHED : PRECISION_NOT_REACHED;
	int stopSignal = 0;
	int cnt = 0;


	MPI_Status status;
	//MPI_Request request;

	int stop_tag = 2;
	int stop_flag = 0;
	MPI_Status stop_status;

	double pih = 0.0;
	double fpisin = 0.0;
	double globalresiduum = 0;



	/* initialize m1 and m2 */
	m1 = 0;
	m2 = 0;


	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* if there is only one process no message exchange is needed,		*/ 
		/* if rank > N - 1 no message exchange needed						*/
		if(size > 1 && rank < (N-1))
		{
			
			/* if not MASTER 
				Recieve previous last row (MASTER does not need it)*/
			if(rank != MASTER)
			{		
				MPI_Recv(Matrix_In[PREVIOUS_LROW], NUM_COLS, MPI_DOUBLE, PREVIOUS_RANK, 0, MPI_COMM_WORLD, &status);

				/* get the precision state of the previous process */
				if(options->termination == TERM_PREC)
				{
					MPI_Recv(&precisionStatePrevious, 1, MPI_INT, PREVIOUS_RANK, 0, MPI_COMM_WORLD, &status);
					if(precisionStatePrevious == STOP_CALCULATION)
					{
					printf("PREV State STOP_CALCULATION received\n, i am rank: %i", rank);

					}
				}
			}	

			/* if 	process is not LAST (LAST does not need it)
						or
					if iteraction_count != 1 (In the first Iteration, the processes dont need the next first row)	

						Recieve next first row 				*/	
			if(rank != LAST && iteration_count != 1)
			{	
				MPI_Recv(Matrix_In[NEXT_FROW], NUM_COLS, MPI_DOUBLE, NEXT_RANK, 0, MPI_COMM_WORLD, &status);
			}
		}	
		
		/* over all rows */
		for (i = 1; i <= numberOfRows; i++)
		{
			double fpisin_i = 0.0;
			

			if (options->inf_func == FUNC_FPISIN)
			{
				//sin richtig berechnen!
				fpisin_i = fpisin * sin(pih * (double)(i+from-1));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}

		}		

		/* if there is only one process no message exchange is needed,		*/ 
		/* if rank > N - 1 no message exchange needed						*/
		if(size > 1 && rank < (N-1))
		{		
			if (options->termination == TERM_PREC)
			{
				/* the current precision state is switched to PRECISION_REACHED when
						the maxresiduum is small as requested
				        AND the previous process already has reached the requested maxresiduum */
				if (maxresiduum < options->term_precision && (precisionStatePrevious == PRECISION_REACHED && precisionStateCurrent == PRECISION_NOT_REACHED))
				{
					precisionStateCurrent = PRECISION_REACHED;
				}

				/* if the previous process has the precision state STOP_CALCULATION, it gets passed to the current one */
				if(precisionStatePrevious == STOP_CALCULATION)
				{
					precisionStateCurrent = STOP_CALCULATION;
					printf("rank:%i  has now STOP_CALCULATION\n", rank);

				}

				if(rank == MASTER && stopSignal == 1)
				{
					precisionStateCurrent = STOP_CALCULATION;
					printf("Master has now STOP_CALCULATION\n");
				}
			}

			/* Send the first row, not required for MASTER process and after the last iteration */
			if((rank != MASTER) && (iteration_count != (int)options->term_iteration) && (precisionStateCurrent != STOP_CALCULATION))
			{				
				MPI_Send(Matrix_Out[FIRST_ROW], NUM_COLS, MPI_DOUBLE, PREVIOUS_RANK, 0, MPI_COMM_WORLD);	
			}

			/* Waiting for the LAST process to receive a stop signal, but only one time
			   This is nonblocking, so the final iteration count may fluctuate */
			if(options->termination == TERM_PREC && rank == MASTER && cnt == 0)
			{
				MPI_Iprobe(LAST, stop_tag, MPI_COMM_WORLD, &stop_flag, &stop_status);
				if(stop_flag)
				{
					printf("Recieving stopSignal\n");
					MPI_Recv(&stopSignal, 1, MPI_INT, LAST, stop_tag, MPI_COMM_WORLD, &stop_status);
					printf("stopSignal is now: %i\n", stopSignal);

					cnt++;
				}

				/* only one time */
			}
			

			if(rank != LAST)
			{	
				/* Send last row to next process, LAST has no following process */
				MPI_Send(Matrix_Out[LAST_ROW], NUM_COLS, MPI_DOUBLE, NEXT_RANK, 0, MPI_COMM_WORLD);

				/* update current precision state at next process */
				if(options->termination == TERM_PREC)
				{
					MPI_Send(&precisionStateCurrent, 1, MPI_INT, NEXT_RANK, 0, MPI_COMM_WORLD);
					if(precisionStateCurrent == STOP_CALCULATION)
					{
						printf("rank: %i  has send current state: STOP_CALCULATION\n", rank);
					}
				}	

			}
			
			/* If the last process reach the precision state PRECISION_REACHED, all processes have reached this state.
			   It sends now a stop signal to the MASTER process and afterwards this state is getting passed to following processes. */
			if(options->termination == TERM_PREC && rank == LAST && precisionStateCurrent == PRECISION_REACHED && stopSignal == 0)
			{
				stopSignal = 1;

				MPI_Send(&stopSignal, 1, MPI_INT, MASTER, stop_tag, MPI_COMM_WORLD);
				printf("Sending stopSignal\n");
			}

		}

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;
		
		results->stat_iteration++;

		/* if TERM_ITER, Broadcast maxresiduum in last iteration. */
		if(term_iteration == 1)
		{
			MPI_Allreduce(&maxresiduum, &globalresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			
			results->stat_precision = globalresiduum;
		}

		/* if precision state is STOP_CALCULATION, set term_iteration = 0 for stopping, allreduce maxresiduum*/
		if(options->termination == TERM_PREC && (precisionStateCurrent == STOP_CALCULATION || rank >= (N-1)))
		{
			term_iteration = 0;

			MPI_Allreduce(&maxresiduum, &globalresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			results->stat_precision = globalresiduum;
		} 

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC && size == 1)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
				results->stat_precision = maxresiduum;
			}
		}
		
		if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}

		iteration_count++;
	}

	results->m = m2;
}
/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculateJacobi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, int rank, int size)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;
	int const numberOfRows = arguments->numberOfRows;
	int const from = arguments->from;

	int term_iteration = options->term_iteration;

	MPI_Status status;

	double pih = 0.0;
	double fpisin = 0.0;
	double globalresiduum = 0;

	/* initialize m1 and m2 */
	m1 = 0;
	m2 = 1;


	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		/* if there is only one process no message exchange is needed,		*/ 
		/* if rank > N - 1 no message exchange needed						*/
		if(size > 1 && rank < (N-1))
		{
			//MASTER does not need the previous last row
			if(rank != MASTER)
			{	
				MPI_Send(Matrix_In[FIRST_ROW], NUM_COLS, MPI_DOUBLE, PREVIOUS_RANK, 0, MPI_COMM_WORLD);
				MPI_Recv(Matrix_In[PREVIOUS_LROW], NUM_COLS, MPI_DOUBLE, PREVIOUS_RANK, 0, MPI_COMM_WORLD, &status);		
			}

			//LAST does not need the next first row
			if(rank != LAST)
			{	
				MPI_Recv(Matrix_In[NEXT_FROW], NUM_COLS, MPI_DOUBLE, NEXT_RANK, 0, MPI_COMM_WORLD, &status);	
				MPI_Send(Matrix_In[LAST_ROW], NUM_COLS, MPI_DOUBLE, NEXT_RANK, 0, MPI_COMM_WORLD);
			}
		}

		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i <= numberOfRows; i++)
		{
			double fpisin_i = 0.0;
			

			if (options->inf_func == FUNC_FPISIN)
			{
				//sin richtig berechnen!
				fpisin_i = fpisin * sin(pih * (double)(i+from-1));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}		

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;
		
		//reduce maxresiduum, broadcast afterwards
		MPI_Allreduce(&maxresiduum, &globalresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		results->stat_iteration++;
		results->stat_precision = globalresiduum;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (globalresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
	
	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options, int size)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauss-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("Processes used:     %i\n", size);
	printf("\n");
}


/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrixMPI (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
         printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}



/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{	
	int rank;
	int size;
	int rest;
	int from, to;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	
	/* Parameter nur einmal abfragen */
	if(rank == 0)
	{
		AskParams(&options, argc, argv);   
	}
	MPI_Bcast(&options, (sizeof(options)), MPI_BYTE, MASTER, MPI_COMM_WORLD);
        	
	initVariables(&arguments, &results, &options);          

	//  Damit allocation + initialization richtig läuft, wird für GS size = 1 gesetzt 
	// if(options.method == METH_GAUSS_SEIDEL)
	// {
	// 	size = 1;
	// }

	/* Aufteilen bis auf rest */
	int N_part = arguments.N;
	int lines = N_part - 1;
	rest = lines % size;
	N_part = (lines - rest) / size;

	/* globale zeilennummer berechnen, hier wird der rest beachtet						*/
	/* offset ist (rank + 1) für rank < rest, steigt also linear mit steigendem rang 	*/
	if(rank < rest)
	{
		from = N_part * rank 		+ rank + 1;
		to = N_part * (rank + 1) 	+ (rank + 1);
	}
	/* offset hier ist rest also die der maximale offset von oben */
	else
	{
		from = N_part * rank 		+ rest + 1;
		to = N_part * (rank + 1) 	+ rest ;
	}
	arguments.to = to;
	arguments.from = from;
	
	
	/* at least we only need N - 1 processes for calculation */
	if((unsigned int)size > (arguments.N -1))
	{
		size = (arguments.N - 1);

		if(rank == MASTER )
		{
			printf("\nWarning, you are using more processes than rows.\n This can slow down the calculation process! \n\n");
		}

	} 

	//calculate Number of Rows
	arguments.numberOfRows = ((to - from + 1) > 0 ) ? (to - from + 1) : 0;

	allocateMatrices(&arguments);        
	initMatrices(&arguments, &options, rank, size);            

	gettimeofday(&start_time, NULL);                   /*  start timer         */


	if (options.method == METH_JACOBI )
	{
		calculateJacobi(&arguments, &results, &options, rank, size);   
	}
	else
	{	
		calculateGS(&arguments, &results, &options, rank, size);   
	}

	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	/* only once */
	if(rank == MASTER)
	{
		displayStatistics(&arguments, &results, &options, size);
	} 


	DisplayMatrixMPI(&arguments, &results, &options, rank, size, from, to);


	freeMatrices(&arguments);                                       /*  free memory     */

	MPI_Finalize();
	return 0;
}
