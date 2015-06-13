#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "mmult_def.h"


int gen_offset(int, int*, int);


void slave()
{
	int i, j, myrank, position, mpi_size, matrix_size, root_info[INFO_SIZE], *numofrows;
	MPI_Status stat;
	int *matrix_col, *matrix_row, *matrix_row_old, *res_fract, *comp_matrix;
	int num_row, num_row_old, num_col, successor, predecessor, offset;


	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	//setup
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	numofrows = calloc(mpi_size, sizeof(int));
	MPI_Bcast(numofrows, mpi_size, MPI_INT, 0, MPI_COMM_WORLD);
	

	//receive information (mt 0)
	MPI_Recv(root_info, INFO_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	num_row = num_col = root_info[2];
	successor = root_info[1];
	predecessor = root_info[0];

	
	//allocate memory for the whole resulting data block
	comp_matrix = calloc(matrix_size * num_col, sizeof(int));
	
	//allocate memory for the single received rows and columns
	matrix_col = malloc(num_col * matrix_size * sizeof(int));
	matrix_row = malloc(num_row * matrix_size * sizeof(int));
	MPI_Recv(matrix_col, num_col * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

	

	MPI_Recv(&num_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	MPI_Recv(matrix_row, num_row * matrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

	
	//generate the first offset for comp_matrix
	position = myrank;
	offset = gen_offset(position, numofrows, num_col);
	//calculate the first fraction of the column from the initial row-data
	res_fract = algo_mul(matrix_row, matrix_col, num_row, num_col, matrix_size);
	//and insert it into the column on its correct position
	for(i=0; i < num_col * num_row; i++)
	{
		*(comp_matrix + offset + i) = *(res_fract + i);
	}
	free(res_fract);
				
	
	for(i=1; i < mpi_size; i++)
	{
		//save old data in dummys
		num_row_old = num_row;
		matrix_row_old = matrix_row;
		//receicve the number of rows that will arrive in the next step
		//this is necessary for an odd number of processes
        	MPI_Recv(&num_row, 1, MPI_INT, successor, 0, MPI_COMM_WORLD, &stat);
		//allocate new memory for matrix_row
		matrix_row = malloc(num_row * matrix_size * sizeof(int));
		//receive the new row-data from the successor
	        MPI_Recv(matrix_row, num_row * matrix_size, MPI_INT, successor, 0, MPI_COMM_WORLD, &stat);
		//send old num_row and matrix_row to the predecessor
		if(myrank > 1)
		{
			MPI_Send(&num_row_old, 1, MPI_INT, predecessor, 0, MPI_COMM_WORLD);
			MPI_Send(matrix_row_old, num_row_old * matrix_size, MPI_INT, predecessor, 0, MPI_COMM_WORLD);
		}
		free(matrix_row_old);
		//calculate the data-fraction from the received rows and the static columns
		res_fract = algo_mul(matrix_row, matrix_col, num_row, num_col, matrix_size);
		if(position < mpi_size-1)
			position++;
		else
			position = 0;
		//insert the calculated data into the column-block
		offset = gen_offset(position, numofrows, num_col);
		for(j=0;  j < num_col * num_row; j++)
			*(comp_matrix + offset + j) = *(res_fract + j);
		free(res_fract);
	}
	
	//finally send the completed column-block back to the master
	MPI_Send(comp_matrix, matrix_size * num_col, MPI_INT, 0, 0, MPI_COMM_WORLD);
	
	free(matrix_row);
	free(comp_matrix);
}


//small sub-function to generate the correct offset
int gen_offset(int position, int *numofrows, int numofcol)
{
	int i, offset=0;

	for (i=0; i < position; i++)
		offset += numofrows[i] * numofcol;
		
	return offset;
}
