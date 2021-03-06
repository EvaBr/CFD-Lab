//#include <mpi.h>
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include "../../helper.h"
#include "../../LBDefinitions.h"

void Program_Message ( char *txt );
/* produces a stderr text output  */



void Programm_Sync ( char *txt );
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop ( char *txt );
/* all processes will produce a text output, be synchronized and finished */

//void initialiseBuffers ( double **sendBuffer, double **readBuffer, int *subdomain, int* flagField );
/* initialises buffers for one process */


void extractionXleft	( double **sendBuffer, double *collideField, int *subdomain );
void extractionXright	( double **sendBuffer, double *collideField, int *subdomain );
void extractionYfront	( double **sendBuffer, double *collideField, int *subdomain );
void extractionYback	( double **sendBuffer, double *collideField, int *subdomain );
void extractionZtop	( double **sendBuffer, double *collideField, int *subdomain );
void extractionZbottom	( double **sendBuffer, double *collideField, int *subdomain );
/* functions for extracting the pdfs that are to be streamed into neighbouring process' region */

void swapXleft		( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
void swapXright		( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
void swapYfront		( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
void swapYback		( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
void swapZtop 		( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
void swapZbottom	( double **sendBuffer, double **readBuffer, int *subdomain, int *proc, int rank);
/* funtions for the pdfs to send to the neighbor's readbuffer, and receive from neighbor's sendbuffer*/


void injectionXleft	( double **readBuffer, double *collideField, int *subdomain );
void injectionXright	( double **readBuffer, double *collideField, int *subdomain );
void injectionYfront	( double **readBuffer, double *collideField, int *subdomain );
void injectionYback	( double **readBuffer, double *collideField, int *subdomain );
void injectionZtop	( double **readBuffer, double *collideField, int *subdomain );
void injectionZbottom	( double **readBuffer, double *collideField, int *subdomain );
/* functions for extracting the pdfs that are to be streamed into neighbouring process' region */
