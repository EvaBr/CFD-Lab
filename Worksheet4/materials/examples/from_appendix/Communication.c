#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

const double PI = 3.14;                        /* the value to be broadcast */

int main(int argc, char* argv[]) {
   
   int i, myrank, nproc;
   int NeighborID_li, NeighborID_re;
   int recv_NeighborID_li, recv_NeighborID_re;
   double pi; 
   MPI_Status status;


   /* initialisation */


   MPI_Init( &argc, &argv );                    /* execute n processes      */
   MPI_Comm_size( MPI_COMM_WORLD, &nproc );     /* asking for the number of processes  */
   MPI_Comm_rank( MPI_COMM_WORLD, &myrank );    /* asking for the local process id   */




   /* determine the neighbors */

   if ( 0 < myrank )
      NeighborID_li = myrank - 1;
   else 
      NeighborID_li = MPI_PROC_NULL;             /* no neighbor */

   if ( (nproc - 1) > myrank )
      NeighborID_re = myrank + 1; 
   else
      NeighborID_re = MPI_PROC_NULL;             /* no neighbor */

   
   /* broadcast of pi from process 0 to all  */

   if ( 0 == myrank )
      pi = PI;
   else
      pi = 0.0;

   printf("Proc %2d : Before the broadcast pi = %e\n", myrank, pi);

   MPI_Bcast( &pi, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );

   printf("Proc %2d : After the broadcast pi = %e\n", myrank, pi);   


   /* send the ID to left/right neighbors (conventional)
      think over how the following send/receive procedure can be implemented with 
      the help of a single MPI_Sendrecv command */

   /* ID of a non-existent neighbor has the value MPI_PROC_NULL  */

   MPI_Send( &myrank, 1, MPI_INT, NeighborID_li, 1, MPI_COMM_WORLD );
   MPI_Recv( &recv_NeighborID_re, 1, MPI_INT, NeighborID_re, 1, 
             MPI_COMM_WORLD, &status );

   if ( MPI_PROC_NULL != NeighborID_re )
      printf("Proc %2d : ID right neighbor = %2d\n", myrank, recv_NeighborID_re);
   else
      printf("Proc %2d : ID right neighbor = MPI_PROC_NULL\n", myrank);

   MPI_Send( &myrank, 1, MPI_INT, NeighborID_re, 2, MPI_COMM_WORLD );
   MPI_Recv( &recv_NeighborID_li, 1, MPI_INT, NeighborID_li, 2,
             MPI_COMM_WORLD, &status );

   if ( MPI_PROC_NULL != NeighborID_li )
      printf("Proc %2d : ID left neighbor = %2d\n", myrank, recv_NeighborID_li);
   else
      printf("Proc %2d : ID left neighbor = MPI_PROC_NULL\n", myrank);
  
 
   /* Epilog */
   fflush(stdout);                       /* write out all output   */
   fflush(stderr);
   MPI_Barrier( MPI_COMM_WORLD );        /* synchronize all processes */
   MPI_Finalize();                       /* end the MPI session */
 
   return 0;
}
