#ifndef __HELPER_H__
#define __HELPER_H__

/* includefiles */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

#ifdef PI
#undef PI
#endif

#define FREE_ARG char*

#define FLUID 0
#define AIR 1
#define NO_SLIP 2
#define FREE_SLIP 3
#define OUTFLOW 4
#define INFLOW 5
#define MOVING_WALL 6

//helper flags for easier checking in boundary functions
#define B_N 8
#define B_O 32
#define B_S 4
#define B_W 16
#define B_U 1
#define B_D 2

#define B_NO 40
#define B_NW 24
#define B_NU 9
#define B_ND 10
#define B_SO 36
#define B_SW 20
#define B_SU 5
#define B_SD 6
#define B_OU 33
#define B_WU 17
#define B_OD 34
#define B_WD 18

#define B_NOU 41
#define B_NOD 42
#define B_NWU 25
#define B_NWD 26
#define B_SOU 37
#define B_SOD 38
#define B_SWU 21
#define B_SWD 22




//#define C_P 32 //additional bit; case of given pressure


/**
 * Maximum length of input lines
 */
#define MAX_LINE_LENGTH 1024

/**
 * Stores the last timer value
 */
extern clock_t last_timer_reset;


int min( int a, int b);
int max( int a, int b);
double fmin( double a, double b);
double fmax( double a, double b);
double mmax(double **M, int imax, int jmax);
double tmax(double ***M, int imax, int jmax, int kmax);
/*Just a helper function, to have less writing*/
int isfluid(int flag);
int pow2 (int en, int dva);
int getbit (int wall);
int interior (int flag);
int getcelltype (int flags); //i, int j, int k, int ***Flag);


/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define ERROR(s)    errhandler( __LINE__, __FILE__, s)

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define ERROUT stdout

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
void  errhandler( int nLine, const char *szFile, const char *szString );


/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_INT( szFileName, VarName)    read_int   ( szFileName, #VarName, &(VarName) )

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_DOUBLE( szFileName, VarName) read_double( szFileName, #VarName, &(VarName) )

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_STRING( szFileName, VarName) read_string( szFileName, #VarName,  (VarName) )

void read_string( const char* szFilename, const char* szName, char*  sValue);
void read_int   ( const char* szFilename, const char* szName, int*    nValue);
void read_double( const char* szFilename, const char* szName, double*  Value);




void write_matrix2( const char* szDebug,       /* filename */
			int	timeStepNumber,
		   double ***m,		       /* matrix */
		   int nrl, int nrh, int ncl, int nch, int nll, int nlh);
/**
 * Writing matrices to a file.
 * -----------------------------------------------------------------------
 * write_matrix(...) wites a matrice to a file
 * the file has the following format
 *
 *    -----------------------------------------
 *    |  xlength          |  float  |  ASCII  |
 *    ----------------------------------------|
 *    |  ylength          |  float  |  ASCII  |
 *    ----------------------------------------|
 *    |  nrl              |  int    |  ASCII  |
 *    ----------------------------------------|
 *    |  nrh              |  int    |  ASCII  |    1. call of the
 *    ----------------------------------------|
 *    |  ncl              |  int    |  ASCII  |
 *    ----------------------------------------|
 *    |  nch              |  int    |  ASCII  |    1. call of the
 *    ----------------------------------------|    function with
 *    |  m[nrl][ncl]      |  float  |  binaer |    bFirst == 1
 *    ----------------------------------------|
 *    |  m[nrl][ncl+1]    |  float  |  binaer |
 *    ----------------------------------------|
 *    |                  .                    |
 *                       .
 *    |                  .                    |
 *    -----------------------------------------
 *    |  m[nrh][nch]      |  float  |  binary |
 *    -----------------------------------------------------------------
 *    |  m[nrl][ncl]      |  float  |  binary |
 *    ----------------------------------------|
 *    |  m[nrl][ncl+1]    |  float  |  binary |     2. call with
 *    ----------------------------------------|     bFirst == 0
 *    |                  .                    |
 *                       .
 *    |                  .                    |
 *    -----------------------------------------
 *    |  m[nrh][nch]      |  float  |  binary |
 *    ------------------------------------------------------------------
 *
 * @param szFileName          name of the file
 * @param m                   matrix
 * @param nrl                 first column
 * @param nrh  		          last column
 * @param ncl                 first row
 * @param nch                 last row
 * @param xlength             size of the geometry in x-direction
 * @param ylength             size of the geometry in y-direction
 * @param xlength             size of the geometry in x-direction
 * @param fFirst              0 == append, else overwrite
 */
void write_matrix(
  const char* szFileName,
  double **m,
  int nrl,
  int nrh,
  int ncl,
  int nch,
  double xlength,
  double ylength,
  int fFirst
);

/**
 * @param szFileName    filehandle
 * @param m             matrix
 * @param nrl           first column
 * @param nrh           last column
 * @param ncl           first row
 * @param nch           last row
 */
void read_matrix( const char* szFileName,	               /* filehandle */
		  double **m,		       /* matrix */
		  int nrl,		       /* first column */
		  int nrh,		       /* last column */
		  int ncl,		       /* first row */
		  int nch );                   /* last row */


/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
double **matrix( int nrl, int nrh, int ncl, int nch );
double ***matrix2( int nrl, int nrh, int ncl, int nch, int nll, int nlh );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch );
void free_matrix2( double ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a);
void init_matrix2( double ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh, double a);

/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
int  **imatrix( int nrl, int nrh, int ncl, int nch );
int  ***imatrix2( int nrl, int nrh, int ncl, int nch, int nll, int nlh );

/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch );
void free_imatrix2( int ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void init_imatrix2( int ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh, int a);
void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a);


/**
 * reads in a ASCII pgm-file and returns the colour information in a two-dimensional integer array.
 * At this, a boundary layer around the image is additionally stored and initialised with 0.
 */
int **read_pgm(const char *filename);


/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMPOUT stdout

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_POSITION() fprintf( DUMPOUT, "%s:%d Dumpposition \n", __FILE__, __LINE__ )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_MESSAGE(s) fprintf( DUMPOUT, "%s:%d %s\n",            __FILE__, __LINE__, s  )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_INT(n)     fprintf( DUMPOUT, "%s:%d %s = %d\n", __FILE__, __LINE__, #n, n )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_DOUBLE(d)  fprintf( DUMPOUT, "%s:%d %s = %f\n", __FILE__, __LINE__, #d, d )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_STRING(s)  fprintf( DUMPOUT, "%s:%d %s = %s\n", __FILE__, __LINE__, #s, s )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define RESET_TIMER()   last_timer_reset = clock()

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_TIMER()    fprintf( DUMPOUT, "%s:%d Timer: %f\n", __FILE__, __LINE__, (float)(clock()-last_timer_reset)/(float)CLOCKS_PER_SEC )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_MATRIX_TO_FILE( m, nrl, nrh, ncl, nch, xlength, ylength) \
        {  \
           static nCount = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #m); \
           write_matrix( szFileName, m, nrl, nrh, ncl, nch, xlength, ylength, nCount == 0); \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_INT_TO_FILE(n) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #n); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%d\n", nCount, n ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_DOUBLE_TO_FILE(d) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #d); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%f\n", nCount, d ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_STRING_TO_FILE(s) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #s); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%s\n", nCount, s ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

#endif
