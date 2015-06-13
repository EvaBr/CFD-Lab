#include <stdlib.h>
#include <stdio.h>

/* Short function that generates two matrices of size <size>-by-<size> and
 * writes them to files <filename1> and <filename2>. */

int main (int argc, char **argv)
{
	FILE *fp;
	int i, j, size;
	

	if(argc != 4)
	{
		printf("\nusage:\n%s <filename1> <filename2> <size>\n\n", argv[0]);
		exit(-1);
	}

	size = atoi(argv[3]);
	

	//open first file and write formated random numbers
	fp = fopen (argv[1], "w");
	fprintf(fp, "%d\n", size);
	for(i=0; i < size; i++)
	{
		for(j=0; j < size; j++)
			fprintf(fp, "%d ", rand() / 1000000);
		fprintf(fp, "\n");
	}      
	fclose(fp);
	
	
	//open second file and write formated random numbers
	fp = fopen (argv[2], "w");
        fprintf(fp, "%d\n", size);
	for(i=0; i < size; i++)
	{
		for(j=0; j < size; j++)
		fprintf(fp, "%d ", rand() / 1000000);
		fprintf(fp, "\n");
	}
	fclose(fp);
			
	
	return 0;
}
