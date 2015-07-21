#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include "helper.h"

/* ----------------------------------------------------------------------- */
/*                             auxiliary functions                         */
/* ----------------------------------------------------------------------- */


inline int createflag(int** Pic,int i,int j,int k,int imax,int jmax,int kmax){
	int temp, temp2;
	const int tarr[] = {1, 2, 0, 3, 4, 7, 8};
	temp = tarr[Pic[i][j + k*(jmax+2)]]*pow2(2, 12) + min(Pic[i+1][j + k*(jmax+2)]+1, 3)*pow2(2, 10) + min(Pic[i-1][j + k*(jmax+2)]+1, 3)*pow2(2, 8) +   min(Pic[i][j + k*(jmax+2)-1]+1, 3)*pow2(2, 6) +  min(Pic[i][j + k*(jmax+2) + 1]+1, 3)*pow2(2, 4) + min(Pic[i][j + (k-1)*(jmax+2)]+1, 3)*4 + min(Pic[i][j + (k+1)*(jmax+2)]+1, 3);   //use min() bcs obstacle neighbours will have numbrs 3-6 in the pic, but they should be flagged with (11)_2 = 3 in the flag field.
	//check for forbidden cells:
	//if ( ((temp > pow(2, 12)*3) || (temp < pow(2, 12))) /*so it is bound.*/ & (E,W water / N,S water / D,U water) ) { error }
	if (Pic[i][j + k*(jmax+2)]>1 && (min(Pic[i][j+k*(jmax+2)+1]+Pic[i][j+k*(jmax+2)-1], Pic[i-1][j+ k*(jmax+2)]+Pic[i+1][j+ k*(jmax+2)])==0  || !(Pic[i][j+(k-1)*(jmax+2)]+Pic[i][j+(k+1)*(jmax+2)])))  {
		ERROR("Invalid geometry! Forbidden boundary cell found.\n");
	}
	temp2 = getboundarytype(temp)%16;
	if (( ((temp>>12)&15) == 8 ) && (temp2==5 || temp2==6 || temp2==9 || temp2==10)){ //tryin to set moving wall for a B_??? cell. not allowed!
		ERROR("Invalid geometry! Forbidden boundary cell found.\n");
	}
	return temp;
}



inline int getboundarytype (int flag){
	//flag = flag;
	//int flags = Flag[i][j][k];
	//int isboundary = (flags > pow(2, 12)*3) || (flags < pow(2, 10));//check if this is really boundary cell <- for now, this is assumed true.
	int flags = (~(2730&flag))&getwallbit(0); // & (10|10|10|10|10|10) - check where is water, and get just the important bits (00 where water, 10 where b or air)
	flags = ((flags&2) >> 1) + ((flags >> 2)&2) + ((flags >> 3)&4) + ((flags >> 4)&8) + ((flags >> 5)&16) + ((flags >> 6)&32);
	//TODO: when doing free surfaces, this might need to be extended for the cases of water/air cells, not just boundary cells. for now, we dont even need check for it being a boundary cell. (well do this in a loop in boundary.c)
	//e.g.:
	//flags = ((getwallbit(0)/3) & flags); // & (01|01|01|01|01|01) - check where is air, and get just the important bits (00 where air, 01 where b or water)

	//flags = flags&getwallbit(isboundary);  //add the cond. of being a boundarycell
	return flags;
}

inline int getsurfacetype (int flag,int *ui,int *vi,int *wi,int *nx,int*ny, int*nz,int *xdb,int*ydb, int*zdb,int*num){
	                     //Flag        ,&ui,    &vi,    &wi,    &nx,   &ny,    &nz,    &mx,   &my,    &mz,   &num
	int type = 0;
	*nx = 0;
	*ny = 0;
	*nz = 0;
	*ui = 0;
	*vi = 0;
	*wi = 0;

	*num = 0;

	if(issurface(flag)){
		if(getbit(flag,10)==0 && getbit(flag,11)==1){
			*ui = 1;
			*nx = 1;

			type|= B_O;
			*num=*num+1;
		}
	    if(getbit(flag,8)==0 && getbit(flag,9)==1){
			*ui =  0;
			if(*nx == 0)
			   *nx = -1;
			else{
				*xdb = 1;
				*nx = 0;
			}

			type|= B_W;
			*num=*num+1;
		}

		if(getbit(flag,6)==0 && getbit(flag,7)==1){
			*vi = 1;
			*ny = 1;

			type|= B_N;
			*num=*num+1;
		}
		if(getbit(flag,4)==0 && getbit(flag,5)==1){
			*vi =  0;

			if(*ny==0)
				*ny = -1;
			else{
				*ydb = 1;
				*ny = 0;
			}

			type|= B_S;
			*num=*num+1;
		}
		if(getbit(flag,0)==0 && getbit(flag,1)==1){ /*todo ok? */
			*wi = 1;
			*nz = 1;
			type|= B_U;
			*num=*num+1;
		}
		if(getbit(flag,2)==0 && getbit(flag,3)==1){
			*wi =  0;
			if(*nz==0)
				*nz = -1;
			else{
				*zdb = 1;
				*nz = 0;
			};

			type|= B_D;
			*num=*num+1;
		}
	}

	return type;
}


inline int isfluid(int flag) {
	if(isboundary(flag)==0){
		int ret = pow(2, 12);
		return ( (flag & (3*ret)) == ret );
	}
	return 0;
}




inline int isboundary(int flag) {
	return (flag > pow2(2, 12)*3) || (flag < pow2(2, 12));
}


inline int isempty(int flag) {
	return !isfluid(flag) && !isboundary(flag);
}


inline void setcelltype(int*flag,int type){
	if(type == FLUID){
		*flag = *flag&(~C_FLUID_2);
		*flag = *flag|C_FLUID_1;
	}
	else if(type==AIR){
		*flag = *flag&(~C_AIR_2);
		*flag = *flag|C_AIR_1;
	}
}

inline int getcelltype(int flag){
	int temp = flag >> 12;
	return (temp >> 2)*2 + temp%2 + ((temp&1)!=((temp>>1)&1))*5 + 2;
}

inline int emptyneighbor(int flag){
	flag = (~flag)&B_ALL;
	if ((flag&1365) != 0){
		return 1;
	}
	return 0;
}

inline int nonfluidneighbor(int flag){
	flag = flag&B_ALL;
	if((flag&(~1365)) !=0){
		return 1;
	}
	return 0;
}


inline int issurface(int flag){
	if(isfluid(flag)){
		return  emptyneighbor(flag);
	}
	return 0;
}


inline void setbit(int *number,int pos){
	*number |= 1 << pos;
}

inline void clearbit(int *number,int pos){
	*number &= ~(1 << pos);
}

inline void changebit(int* number,int pos, int bit){
	*number ^= (-bit ^ *number) & (1 << pos);
}


inline int getbit(int number,int pos){
	return (number >> pos) & 1;
}





/*helper function for getting the right bit represent. of edges.*/
inline int getwallbit (int wall) {
	const int tarr[] = {0, 2, 0, 3, 4, 7, 8};  //first two should never be used; we set the first one to 0 for convenience when only pow(..) need to be computed
	return (tarr[wall]*pow2(2, 12)+3*(pow2(2,10)+pow2(2,8)+pow2(2,6)+16+4+1)); //this represents a cell that has obstacles all round, and is itself obstacle of sort 'wall'
}

inline int interior (int flag) {
	//int big = Flag[i][j][k];
	int big = (flag > pow2(2, 12)*3) || (flag < pow2(2, 12)); //check the cell itself is b
	int gb = getwallbit(0);
	return (big && ((flag&gb)==gb)); //check also all neighb. are b
}







inline int min( int a, int b)
{
	if( a < b ) return a;
	return b;
}

inline int max( int a, int b)
{
	if( a > b ) return a;
	return b;
}

inline double fmin( double a, double b)
{
	if( a < b ) return a;
	return b;
}

inline double fmax( double a, double b)
{
	if( a > b ) return a;
	return b;
}

inline int pow2 (int en, int dva) {
	int ret = pow(en, dva);
	return ret;
}

inline double mmax( double **U, int imax, int jmax)
{
	double maxij = 0;
	for( int i=0; i<=imax+1; i++){
		for( int j=0; j<=jmax+1; j++){
			if (fabs(U[i][j])>maxij){
				maxij = fabs(U[i][j]);
			}
		}
	}
	return maxij;
}



double tmax( double ***U, int imax, int jmax, int kmax) //added function for getting max of a tensor
{
	double max = 0;
	for( int i=0; i<=imax+1; i++){
		for( int j=0; j<=jmax+1; j++){
			for( int k=0; k<=kmax+1; k++){
				if (fabs(U[i][j][k])>max){
					max = fabs(U[i][j][k]);
				}
			}
		}
	}
	return max;
}






/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int( const int n1, const int n2 )
{
	if( n1 < n2 ) return n1;
	return n2;
}


/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
	int err = errno;

	fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
	fprintf( ERROUT, "\n" );

	/* if an error within the c-library occured, an error code can be   */
	/* found in the global variable err                                 */
	if( err != 0 )
	{
		fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
		fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
	}
	exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
		{ char szTmp[80]; \
		if( nLine ) \
		sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
		else \
		sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
		ERROR( szTmp ); \
		}


/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char* find_string( const char* szFileName, const char *szVarName )
{
	int nLine = 0;
	int i;
	FILE *fh = NULL;

	static char szBuffer[MAX_LINE_LENGTH];	/* containes the line read  */
	/* from the datafile        */

	char* szLine = szBuffer;
	char* szValue = NULL;
	char* szName = NULL;

	/* open file */
	fh = fopen( szFileName, "rt" );
	if( fh == 0 )
		READ_ERROR("Could not open file", szVarName, szFileName, 0);

	/* searching */
	while( ! feof(fh) )
	{
		fgets( szLine, MAX_LINE_LENGTH, fh );
		++nLine;

		/* remove comments */
		for( i = 0; i < strlen(szLine); i++)
			if( szLine[i] == '#' )
			{
				szLine[i] = '\0'; /* Stringende setzen */
				break;
			}

		/* remove empty lines */
		while( isspace( (int)*szLine ) && *szLine) ++szLine;
		if( strlen( szLine ) == 0) continue;

		/* now, the name can be extracted */
		szName = szLine;
		szValue = szLine;
		while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;

		/* is the value for the respective name missing? */
		if( *szValue == '\n' || strlen( szValue) == 0)
			READ_ERROR("wrong format", szName, szFileName, nLine);

		*szValue = 0;		/* complete szName! at the right place */
		++szValue;

		/* read next line if the correct name wasn't found */
		if( strcmp( szVarName, szName)) continue;

		/* remove all leading blnkets and tabs from the value string  */
		while( isspace( (int)*szValue) ) ++szValue;
		if( *szValue == '\n' || strlen( szValue) == 0)
			READ_ERROR("wrong format", szName, szFileName, nLine);

		fclose(fh);
		return szValue;
	}

	READ_ERROR("variable not found", szVarName, szFileName, nLine);

	return NULL;		/* dummy to satisfy the compiler  */
}

void read_string( const char* szFileName, const char* szVarName, char*   pVariable)
{
	char* szValue = NULL;	/* string containg the read variable value */

	if( szVarName  == 0 )  ERROR("null pointer given as variable name" );
	if( szFileName == 0 )  ERROR("null pointer given as filename" );
	if( pVariable  == 0 )  ERROR("null pointer given as variable" );

	if( szVarName[0] == '*' )
		szValue = find_string( szFileName, szVarName +1 );
	else
		szValue = find_string( szFileName, szVarName );

	if( sscanf( szValue, "%s", pVariable) == 0)
		READ_ERROR("wrong format", szVarName, szFileName,0);

	printf( "File: %s\t\t%s%s= %s\n", szFileName,
			szVarName,
			&("               "[min_int( strlen(szVarName), 15)]),
			pVariable );
}

void read_int( const char* szFileName, const char* szVarName, int* pVariable)
{
	char* szValue = NULL;	/* string containing the read variable value */

	if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
	if( szFileName == 0 )  ERROR("null pointer given as filename" );
	if( pVariable  == 0 )  ERROR("null pointer given as variable" );

	if( szVarName[0] == '*' )
		szValue = find_string( szFileName, szVarName +1 );
	else
		szValue = find_string( szFileName, szVarName );

	if( sscanf( szValue, "%d", pVariable) == 0)
		READ_ERROR("wrong format", szVarName, szFileName, 0);

	printf( "File: %s\t\t%s%s= %d\n", szFileName,
			szVarName,
			&("               "[min_int( strlen(szVarName), 15)]),
			*pVariable );
}

void read_double( const char* szFileName, const char* szVarName, double* pVariable)
{
	char* szValue = NULL;	/* String mit dem eingelesenen Variablenwert */

	if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
	if( szFileName == 0 )  ERROR("null pointer given as filename" );
	if( pVariable  == 0 )  ERROR("null pointer given as variable" );

	if( szVarName[0] == '*' )
		szValue = find_string( szFileName, szVarName +1 );
	else
		szValue = find_string( szFileName, szVarName );

	if( sscanf( szValue, "%lf", pVariable) == 0)
		READ_ERROR("wrong format", szVarName, szFileName, 0);

	printf( "File: %s\t\t%s%s= %f\n", szFileName,
			szVarName,
			&("               "[min_int( strlen(szVarName), 15)]),
			*pVariable );
}


/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */




void write_matrix2( const char* szDebug,       /* filename */
		int	timeStepNumber,
		double ***m,		       /* matrix */
		int nrl, int nrh, int ncl, int nch, int nll, int nlh)
{
	int i, j,k;
	FILE * fh = 0;


	char szFileName[80];
	if(timeStepNumber<1000){ ///media/norbert/940CB6150CB5F27A/Documents/
		sprintf( szFileName, "simulation/debug/d%i_%s", timeStepNumber,szDebug );
	}
	else
	{
		timeStepNumber = timeStepNumber-1000; ///media/norbert/940CB6150CB5F27A/Documents
		sprintf( szFileName, "simulation/debug/P/d%i_%s", timeStepNumber,szDebug );
	}


	fh = fopen( szFileName, "w");	/* overwrite file/write new file */
	if( fh == NULL )			/* opening failed ? */
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
		ERROR( szBuff );
	}

	/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */

	for( k = nll; k <= nlh; k++){
		for( j = ncl; j <= nch; j++){
			for( i = nrl; i <= nrh; i++){
				fprintf( fh, "%f ", m[i][j][k] );
			}
			fprintf( fh, "\n" );
		}
		fprintf( fh, "-----------------------------------------------------------\n" );
	}

	if( fclose(fh) )
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

}



void write_matrix( const char* szFileName,       /* filename */
		double **m,		       /* matrix */
		int nrl,		       /* first column */
		int nrh,		       /* last column */
		int ncl,		       /* first row */
		int nch,		       /* last row */
		double xlength,	       /* size of the geometry in */
		/* x-direction */
		double ylength,	       /* size of the geometry in */
		/* y-direction  */
		int fFirst ) 	       /* 0 == append, else overwrite*/
{
	int i, j;
	FILE * fh = 0;
	int nSize = (nrh-nrl+1) * (nch-ncl+1);
	float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
	int k = 0;

	if( fFirst )				/* first call of the function ? */
	{
		fh = fopen( szFileName, "w");	/* overwrite file/write new file */
		if( fh == NULL )			/* opening failed ? */
		{
			char szBuff[80];
			sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
			ERROR( szBuff );
		}

		/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
	}
	else
	{
		fh = fopen( szFileName ,"a");	/* append to the file */
		if( fh == NULL )			/* opening failed ? */
		{
			char szBuff[80];
			sprintf( szBuff, "Outputfile %s cannot be opened", szFileName );
			ERROR( szBuff );
		}
	}

	for( j = ncl; j <= nch; j++)
		for( i = nrl; i <= nrh; i++)
			tmp[k++] = (float)m[i][j];

	fwrite( tmp, sizeof(float), nSize, fh);

	if( fclose(fh) )
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

	free( tmp );
}


void read_matrix( const char* szFileName,       /* filename */
		double **m,		       /* matrix */
		int nrl,		       /* first column */
		int nrh,		       /* last column */
		int ncl,		       /* first row */
		int nch		       /* last row */
)
{
	int i, j;
	FILE * fh = 0;
	int nSize = (nrh-nrl+1) * (nch-ncl+1);
	float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
	int k = 0;

	fh = fopen( szFileName, "r");	/* overwrite file/write new file */
	if( fh == NULL )			/* opening failed ? */
	{
		char szBuff[80];
		sprintf( szBuff, "Can not read file %s !!!", szFileName );
		ERROR( szBuff );
	}


	fread( tmp, sizeof(float), nSize, fh);

	for( j = ncl; j <= nch; j++)
		for( i = nrl; i <= nrh; i++)
			m[i][j]=tmp[k++];

	if( fclose(fh) )
	{
		char szBuff[80];
		/*orig bug:
       sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
		sprintf( szBuff, "Inputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

	free( tmp );
}



/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */


double ***matrix2( int nrl, int nrh, int ncl, int nch, int nll, int nlh ){
	double ***pArray;
	int i,j;

	pArray=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double **));
	if (pArray== 0) ERROR("Storage cannot be allocated");
	pArray -= nrl;

	for(i=nrl;i<=nrh;i++) {
		pArray[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double *));
		if (pArray[i]== 0) ERROR("Storage cannot be allocated");
		pArray[i] -= ncl;
		for(j=ncl;j<=nch;j++) {
			pArray[i][j]=(double *) malloc((unsigned) (nlh-nll+1)*sizeof(double));
			if (pArray[i][j] == 0) ERROR("Storage cannot be allocated");
			pArray[i][j] -= nll;
		}

	}
	return pArray;
}


/*  allocates storage for a matrix                                         */
double **matrix( int nrl, int nrh, int ncl, int nch )
{
	int i;
	int nrow = nrh - nrl + 1;	/* compute number of lines */
	int ncol = nch - ncl + 1;	/* compute number of columns */

	double **pArray  = (double **) malloc((size_t)( nrow * sizeof(double*)) );
	double  *pMatrix = (double *)  malloc((size_t)( nrow * ncol * sizeof( double )));

	if( pArray  == 0)  ERROR("Storage cannot be allocated");
	if( pMatrix == 0)  ERROR("Storage cannot be allocated");

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

/* deallocates the storage of a matrix  */
void free_matrix2( double ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh ){
	short i,j;
	for(i=nrh;i>= nrl;i--)
	{
		for(j=nch;j>= ncl;j--){

			free(m[i][j]+nll);
		}
		free( (m[i]+ncl));
	}
	free( (m+nrl));
}

void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
{
	double **pArray  = m + nrl;
	double  *pMatrix = m[nrl]+ncl;

	free( pMatrix );
	free( pArray );
}

void init_matrix2( double ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh, double a){
	int i,j,k;
	for( i = nrl; i <= nrh; i++)
		for( j = ncl; j <= nch; j++)
			for( k = nll; k <= nlh; k++)
				m[i][j][k] = a;
}

void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a)
{
	int i,j;
	for( i = nrl; i <= nrh; i++)
		for( j = ncl; j <= nch; j++)
			m[i][j] = a;
}


int  ***imatrix2( int nrl, int nrh, int ncl, int nch, int nll, int nlh ){
	int ***pArray;
	int i,j;
	pArray=(int ***) malloc((unsigned) (nrh-nrl+1)*sizeof(int **));
	if (pArray== 0) ERROR("Storage cannot be allocated");
	pArray -= nrl;
	for(i=nrl;i<=nrh;i++) {
		pArray[i]=(int **) malloc((unsigned) (nch-ncl+1)*sizeof(int *));
		if (pArray[i]== 0) ERROR("Storage cannot be allocated");
		pArray[i] -= ncl;
		for(j=ncl;j<=nch;j++) {
			pArray[i][j]=(int *) malloc((unsigned) (nlh-nll+1)*sizeof(int));
			if (pArray[i][j] == 0) ERROR("Storage cannot be allocated");
			pArray[i][j] -= nll;}

	}
	return pArray;
}

/* allocates storage for a matrix */
int **imatrix( int nrl, int nrh, int ncl, int nch )
{
	int i;

	int nrow = nrh - nrl + 1;	/* compute number of rows */
	int ncol = nch - ncl + 1;	/* compute number of columns */

	int **pArray  = (int **) malloc((size_t)( nrow * sizeof( int* )) );
	int  *pMatrix = (int *)  malloc((size_t)( nrow * ncol * sizeof( int )));


	if( pArray  == 0)  ERROR("Storage cannot be allocated");
	if( pMatrix == 0)  ERROR("Storage cannot be allocated");

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

void write_flag_imatrix( const char* szDebug,       /* filename */
		int	timeStepNumber,
		int ***m,		       /* matrix */
		int nrl, int nrh, int ncl, int nch, int nll, int nlh)
{
	int i, j,k;
	FILE * fh = 0;


	char szFileName[80];
///media/norbert/940CB6150CB5F27A/Documents/
	sprintf( szFileName, "simulation/debug/d%i_%s", timeStepNumber,szDebug );


	fh = fopen( szFileName, "w");	/* overwrite file/write new file */
	if( fh == NULL )			/* opening failed ? */
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
		ERROR( szBuff );
	}

	/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */

	for( k = nll; k <= nlh; k++){
		for( j = ncl; j <= nch; j++){
			for( i = nrl; i <= nrh; i++){
				//fprintf( fh, "%d ",getcelltype( m[i][j][k]) );
				if(!isfluid(m[i][j][k])){
					fprintf( fh, "%d ",getboundarytype( m[i][j][k]) );
				}
				else{
					fprintf( fh, "F " );
				}
			}
			fprintf( fh, "\n" );
		}
		fprintf( fh, "-----------------------------------------------------------\n" );
	}

	if( fclose(fh) )
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

}










void write_imatrix2( const char* szDebug,       /* filename */
		int	timeStepNumber,
		int ***m,		       /* matrix */
		int nrl, int nrh, int ncl, int nch, int nll, int nlh)
{
	int i, j,k;
	FILE * fh = 0;


	char szFileName[80];
///media/norbert/940CB6150CB5F27A/Documents/
	sprintf( szFileName, "simulation/debug/d%i_%s", timeStepNumber,szDebug );


	fh = fopen( szFileName, "w");	/* overwrite file/write new file */
	if( fh == NULL )			/* opening failed ? */
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
		ERROR( szBuff );
	}

	/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */

	for( k = nll; k <= nlh; k++){
		for( j = ncl; j <= nch; j++){
			for( i = nrl; i <= nrh; i++){
				fprintf( fh, "%d ", m[i][j][k] );
			}
			fprintf( fh, "\n" );
		}
		fprintf( fh, "-----------------------------------------------------------\n" );
	}

	if( fclose(fh) )
	{
		char szBuff[80];
		sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
		ERROR( szBuff );
	};

}





void free_imatrix2( int ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh ){
	short i,j;
	for(i=nrh;i>= nrl;i--)
	{
		for(j=nch;j>= ncl;j--){

			free(m[i][j]+nll);
		}
		free( (m[i]+ncl));
	}
	free( (m+nrl));
}

/* deallocates the storage of a matrix  */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch )
{
	int **pArray  = m + nrl;
	int  *pMatrix = m[nrl]+ncl;

	free( pMatrix );
	free( pArray );
}

void init_imatrix2( int ***m, int nrl, int nrh, int ncl, int nch, int nll, int nlh, int a){
	int i,j,k;
	for( i = nrl; i <= nrh; i++)
		for( j = ncl; j <= nch; j++)
			for( k = nll; k <= nlh; k++)
				m[i][j][k] = a;
}

void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a)
{
	int i,j;
	for( i = nrl; i <= nrh; i++)
		for( j = ncl; j <= nch; j++)
			m[i][j] = a;
}


int **read_pgm(const char *filename)
{
	FILE *input = NULL;
	char line[1024];
	int levels;
	int xsize, ysize;
	int i1, j1;
	int **pic = NULL;


	if ((input=fopen(filename,"rb"))==0)
	{
		char szBuff[80];
		sprintf( szBuff, "Can not read file %s!", filename );
		ERROR( szBuff );
	}

	/* check for the right "magic number" */
	if ( fread(line,1,3,input)!=3 )
	{
		fclose(input);
		ERROR("Error: Wrong Magic field!");
	}

	/* skip the comments */
	do
		fgets(line,sizeof line,input);
	while(*line=='#');

	/* read the width and height */
	sscanf(line,"%d %d\n",&xsize,&ysize);

	printf("Image set size: %d x %d\n", xsize,ysize);

	/* read # of gray levels */
	fgets(line,sizeof line,input);
	sscanf(line,"%d\n",&levels);

	/* allocate memory for image */
	pic = imatrix(0,xsize-1,0,ysize-1);
	printf("Image initialised...\n");

	/* read pixel row by row */
	for(j1=0; j1 <= ysize-1; j1++)
	{
		for (i1=0; i1 <= xsize-1; i1++)
		{
			int byte;
			fscanf(input, "%d", &byte);

			if (byte==EOF)
			{
				fclose(input);
				ERROR("Read failed");
			}
			else
			{
				pic[i1][j1] = byte;
				//printf("%d,%d: %d\n", i1, ysize+1-j1, min(byte, 1));
			}
		}
	}
	/*
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][0] = 0;
    }
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][ysize+1] = 0;
    }
    for (j1 = 0; j1 < ysize+2; j1++)
    {
        pic[0][j1] = 0;
        pic[xsize+1][j1] = 0;
    }
	 */

	/* close file */
	fclose(input);

	return pic;
}
