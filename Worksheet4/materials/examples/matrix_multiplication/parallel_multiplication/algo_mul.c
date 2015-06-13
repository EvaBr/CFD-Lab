#include <stdlib.h>
#include <stdio.h>

int *algo_mul(int *row, int *col, int num_row, int num_col, int matrix_size)
{
	int i, j, k;
	int *fract_matrix;

	fract_matrix = calloc(num_col * num_row, sizeof(int));


	
	for (k=0; k < num_row; k++)
	{
		for (j=0; j < num_col; j++)
		{
			for (i=0; i < matrix_size; i++)
			{
				*(fract_matrix + j + k * num_col) += row[i + k * matrix_size] * col[i + j * matrix_size];
			}
		}
	}
/*	
	printf("\nalgo_mul:\n");
	for(i=0; i < num_col * num_row; i++)
		printf("%d ", fract_matrix[i]);
	printf("\n\n");
*/	
	return fract_matrix;
}
