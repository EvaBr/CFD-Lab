//#include <mpi.h>
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void initialiseBuffers( double **sendBuffer, double **readBuffer, int *subdomain);
/* initialises buffers for one process */


void extractionX ( double **sendBuffer, int rank, int *subdomain);
void extractionY ( double **sendBuffer, int rank, int *subdomain);
void extractionZ ( double **sendBuffer, int rank, int *subdomain);
/* functions for extracting the pdfs that are to be streamed into neighbouring process' region */
