#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#define MASTER 0

int
main()
{	
	int  numtasks, taskid, source, i, tag;
	//char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	tag = 1;


	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);


	char buffer_time[20];
	char buffer_hostname[45];
	char buffer_microsec[7];
	int microsec;
	struct tm* tm;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	
	microsec = tv.tv_usec;
	tm = localtime(&tv.tv_sec);
	snprintf(buffer_microsec, 7, "%07d", microsec);

	strftime(buffer_time, 20, "%Y-%m-%d %H:%M:%S", tm);

	gethostname(buffer_hostname, sizeof buffer_hostname);
	strcat(buffer_hostname, ": ");
	strcat(buffer_hostname, buffer_time);
	strcat(buffer_hostname,".");
	strcat(buffer_hostname, buffer_microsec);
	
	if (taskid == MASTER)
	{
		//own Hostname : timestamp
    	printf("%s\n", buffer_hostname);

    	char formatted_string[45];
		//wait for results from other tasks
		for (i=1; i<numtasks; i++)
		{
			source = i;
			//MPI_RECV(buf,count,datatype,source,tag,comm,status)
			MPI_Recv(formatted_string, 45, MPI_CHAR, source, tag,  MPI_COMM_WORLD, &status);
			printf("%s\n", formatted_string);
		}
	}

	if (taskid > MASTER)
	{
		//MPI_SEND(buf,count,datatype,dest,tag,comm)
		MPI_Send(buffer_hostname, 45, MPI_CHAR, MASTER, tag, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	printf("Rang %i beendet jetzt!\n", taskid);

	MPI_Finalize();
}
