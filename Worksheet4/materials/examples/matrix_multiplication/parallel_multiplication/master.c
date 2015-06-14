#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "mmult_def.h"


int * master(int *m1, int *m2, int matrix_size)
{
	int i, j, k, mpi_size, myrank, offset=0;
	int send_info[INFO_SIZE], *numofcol, *numofrow, d1, d2, offset_row, offset_col;
	int *res_fract, *rm;
	MPI_Status stat;

        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
		
	rm = calloc(matrix_size * matrix_size, sizeof(int));
	//how many rows/columns gets every node?
	//**************************************
	//algorithm to calculate matrices with a line-count that cannot be (evenly) divided by the number of processes
	numofcol = malloc(mpi_size * sizeof(int));	//allocate array of the size of the processes - every field represents one process
	numofrow = malloc(mpi_size * sizeof(int));
	d1 = matrix_size / mpi_size;	//how much gets every node?
	d2 = matrix_size % mpi_size;	//rest - how many nodes get one more row/column?
	for(i = mpi_size-1; i >= 0; i--)
	{
		numofrow[i] = numofcol[i] = d1;
		if(d2-- > 0)
			numofcol[i]++;
	}
	
	//***setup***
	//broadcast the size
	MPI_Bcast(&matrix_size, 1, MPI_INT, myrank, MPI_COMM_WORLD);
	//and the array of the dynamic row-numbers
	MPI_Bcast(numofrow, mpi_size, MPI_INT, myrank, MPI_COMM_WORLD);

	//send predecessor, successor and number of columns
	for(i=1; i < mpi_size; i++)
	{
		if(i == 0)
			send_info[0] = mpi_size;	//presuc[0] is predecessor
		else
			send_info[0] = i-1;
		if(i == mpi_size-1)
			send_info[1] = 0;		//presuc[1] is successor
		else
			send_info[1] = i+1;
		send_info[2] = numofcol[i];
		MPI_Send(send_info, INFO_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	
	//send the initial blocks of columns and rows to each node
	offset = numofcol[0] * matrix_size;
	for(i=1; i < mpi_size; i++)
	{
		MPI_Send(&m2[offset], numofcol[i] * matrix_size, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&numofrow[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&m1[offset], numofrow[i] * matrix_size, MPI_INT, i, 0, MPI_COMM_WORLD);
		offset += numofcol[i] * matrix_size;
	}
	res_fract = algo_mul(m1, m2, numofrow[0], numofcol[0], matrix_size);
	
	offset_col = 0;
	i=0;	
        for(j=0; j < numofrow[0]; j++)
        {
		for(k=0; k < numofcol[0]; k++)
		{
			*(rm + offset_col + k + j * matrix_size) = res_fract[k + j * numofcol[0]];
		}
	}
	free(res_fract);

	//perform the calculation dynamically (dependend on the number of processes) and distributed
	offset_row = 0;
	for(i=1; i < mpi_size; i++)
	{
		//first send the row-data to the predecessor
		MPI_Send(&numofrow[i-1], 1, MPI_INT, mpi_size-1, 0, MPI_COMM_WORLD); 
		MPI_Send(&m1[offset_row], numofrow[i-1] * matrix_size, MPI_INT, mpi_size-1, 0, MPI_COMM_WORLD);
		//set the correct offset to write the data into the resulting matrix
		offset_row += numofrow[i-1] * matrix_size;
		//now perform your own calculation
		res_fract = algo_mul(m1 + offset_row, m2, numofrow[i], numofcol[0], matrix_size);
		//and write the resulting data into the resulting matrix
		for(k=0; k < numofrow[i]; k++)
			for(j=0; j < numofcol[0]; j++)
				*(rm + offset_row  + k * matrix_size + j) = *(res_fract + j + k * numofcol[0]);
		//finally free the allocated data of 
		free(res_fract);
	}
	
	//now collect the resulting columns from each slave and store it to the right position of the resulting matrix
	offset_row = numofrow[0];
	for(i=1; i < mpi_size; i++)
	{
		res_fract = malloc(numofcol[i] * matrix_size * sizeof(int));
		MPI_Recv(res_fract, numofcol[i] * matrix_size, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);
                for(k=0; k < matrix_size; k++)
                        for(j=0; j < numofcol[i]; j++)
                                *(rm + offset_row  + k * matrix_size + j) = *(res_fract + j + k * numofcol[0]);
		offset_row += numofrow[i];
	}
	return rm;
}
