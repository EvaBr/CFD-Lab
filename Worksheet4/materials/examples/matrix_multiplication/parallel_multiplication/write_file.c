#include <stdlib.h>
#include <stdio.h>
#include "mmult_def.h"


int write_file(char *filename, int *matrix, int size)
{
	FILE *fp;
	int i, j;
	
	fp = fopen (filename, "w");
	fprintf(fp, "%d\n", size);
	for(i=0; i < size; i++)
	{
		for(j=0; j < size; j++)
			fprintf(fp, "%d ", *(matrix + j + size * i));
		fprintf(fp, "\n");
	}      
	fclose(fp);
	
	return 0;
}
