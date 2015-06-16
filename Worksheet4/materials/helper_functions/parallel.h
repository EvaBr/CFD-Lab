//#include <mpi.h>
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void Program_Message ( char *txt );
/* produces a stderr text output  */



void Programm_Sync ( char *txt );
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop ( char *txt );
/* all processes will produce a text output, be synchronized and finished */

void initialiseBuffers ( double **sendBuffer, double **readBuffer, int *subdomain );
/* initialises buffers for one process */


void extractionXleft	( double **sendBuffer, double *collideField, int *subdomain );
void extractionXright	( double **sendBuffer, double *collideField, int *subdomain );
void extractionYfront	( double **sendBuffer, double *collideField, int *subdomain );
void extractionYback	( double **sendBuffer, double *collideField, int *subdomain );
void extractionZtop	( double **sendBuffer, double *collideField, int *subdomain );
void extractionZbottom	( double **sendBuffer, double *collideField, int *subdomain );
/* functions for extracting the pdfs that are to be streamed into neighbouring process' region */

void injectionXleft	( double **readBuffer, double *collideField, int *subdomain );
void injectionXright	( double **readBuffer, double *collideField, int *subdomain );
void injectionYfront	( double **readBuffer, double *collideField, int *subdomain );
void injectionYback	( double **readBuffer, double *collideField, int *subdomain );
void injectionZtop	( double **readBuffer, double *collideField, int *subdomain );
void injectionZbottom	( double **readBuffer, double *collideField, int *subdomain );
/* functions for extracting the pdfs that are to be streamed into neighbouring process' region */
