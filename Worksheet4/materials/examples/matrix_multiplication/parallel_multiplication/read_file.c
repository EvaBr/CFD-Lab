#include <stdio.h>
#include <stdlib.h>
#include "mmult_def.h"

int *read_file(char *filename, int *psize, int hor_vert)
{
	FILE *fp;
	int i, j, size, *matrix;
	char size_buffer[16], *buffer;
	

	//open file
	if((fp = fopen(filename, "r")) < 0)
	{
		printf("\nfile error\n");
		exit(-1);
	}
	
	//read the size of the matrix
	fgets(size_buffer, 100, fp);
	size = atoi(size_buffer);

	//dynamically allocate the new matrix and the buffer
	matrix = malloc(size * size * sizeof(int));
	buffer = malloc((size * 5 + 2) * sizeof(char));
	
	//read the rest of the file and write it formated into the two-dimensional matrix-array
	if(hor_vert == MAT_NORM)	//write the data into ram in the normal order
	{
		for(i=0; i < size; i++)		//rows
		{
			fgets(buffer, size * 8, fp);
			*(matrix + (i * size)) = atoi((char *) strtok(buffer, " "));

			for(j=1; j < size; j++)		//columns
			{
				*(matrix + (j + i * size)) = atoi((char *) strtok(NULL, " "));
			}
		}
	}
	else	//write the data into the ram 90 degr turned 
	{
		for(i=0; i < size; i++)		//rows
		{
			fgets(buffer, size * 8, fp);
			*(matrix + i) = atoi((char *) strtok(buffer, " "));
			
			for(j=1; j < size; j++)		//columns
			{
				*(matrix + (j * size + i)) = atoi((char *) strtok(NULL, " "));
			}
		}
	}
	
	*psize = size;
	free(buffer);
	fclose(fp);
	return matrix;
}
