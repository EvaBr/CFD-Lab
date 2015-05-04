#include <stdio.h>
#include <stdlib.h>
#include "helper.h"

double **matrix( int nrl, int nrh, int ncl, int nch )
{
   int i;
   int nrow = nrh - nrl + 1;    /* compute number of lines */
   int ncol = nch - ncl + 1;    /* compute number of columns */
   
   double **pArray  = (double **) malloc((size_t)( nrow * sizeof(double*)) );
   double  *pMatrix = (double *)  malloc((size_t)( nrow * ncol * sizeof( double )));

   //if( pArray  == 0)  ERROR("Storage cannot be allocated");
   //if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
{
   double **pArray  = m + nrl;
   double  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}



double mmax( double **U, int imax, int jmax)
{   int i,j;
    double maxij = U[0][5];
    for( i=0; i<imax+1; i++){
        for( j=5; j<=jmax; j++){
            if (fabs(U[i][j])>maxij){
                maxij = fabs(U[i][j]);
            }
        }
    }
    return maxij;
}



void main(){
	int imax = 2, jmax = 10;
	double **M;
	int i, j;

	M = matrix(0, imax, 5, jmax);

/**	for (i=0; i<imax+1; i++){
                for (j=1; j<=jmax; j++){
                        printf("%f ", M[i][j]);
                }
		printf("\n");
        } **/

//	printf("\n\n");
	for (i=0; i<imax+1; i++){
		for (j=5; j<=jmax; j++){
			M[i][j] = (double) rand()/ 1357.591 - 1000000;
		}
	}

	for (i=0; i<imax+1; i++){
                for (j=5; j<=jmax; j++){
                        printf("%f ", M[i][j]);
                }
		printf("\n");
        }


	double maxi = mmax(M, imax, jmax);
	printf("\n maximum is %f \n", maxi);
	
	free_matrix(M, 0, imax, 5, jmax);
}
