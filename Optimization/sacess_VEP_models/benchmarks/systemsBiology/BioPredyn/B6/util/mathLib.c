
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <global.h>
#include <mathLib.h>

double
Trunc( double x, int n ) {
    int power = pow( 10, n );
    return floor( x * power ) / power;
}

double
AbsDiff( double x, double y ) {
    return abs( x - y );
}

double
SqrDiff( double x, double y ) {
    return ( x - y ) * ( x - y );
}

/*********************************************************************
 ** to malloc memories                                              **
 ** AllocChar1D(size): size*char                                 **
 ** to realloc memories                                             **
 ** ReallocChar1D(s, size)                                           **
 ** to free memories                                                **
 ** SharefreeM1c(s)                                                 **
 *********************************************************************/
char *
AllocChar1D( int size ) {
    char *s = NULL;
    if( ( s = ( char * ) calloc( size, sizeof( char ) ) ) == NULL ) {
        printf( "char * malloc error!\n" );
        exit( 1 );
    }
    return s;
}

char *
ReallocChar1D( char *s, int size ) {
    if( !s ) {
        s = AllocChar1D( size );
        return s;
    }
    if( ( s = ( char * ) realloc( s, size * sizeof( char ) ) ) == NULL ) {
        printf( "char * realloc error!" );
        exit( 1 );
    }
    return s;
}

void
FreeChar1D( char *s ) {
    if( s )
        free( ( void * ) s );
    return;
}

char **
AllocChar2D( int size1, int size2 ) {
    int i;
    char **s = NULL;

    if( ( s = ( char ** ) calloc( size1, sizeof( char * ) ) ) == NULL ) {
        printf( "char ** malloc error!\n" );
        exit( 1 );
    }
    if( size2 > 0 ) {
        for( i = 0; i < size1; i++ )
            s[i] = AllocChar1D( size2 );
    }
    return s;
}

char **
ReallocChar2D( char **s, int size1, int size2 ) {
    int i;

    if( ( s = ( char ** ) realloc( s, size1 * sizeof( char * ) ) ) == NULL ) {
        printf( "char ** realloc error!" );
        exit( 1 );
    }
    if( size2 > 0 ) {
        for( i = 0; i < size1; i++ )
            s[i] = ReallocChar1D( s[i], size2 );
    }
    return s;
}

void
FreeChar2D( char **s, int size ) {
    int i;

    if( !s )
        return;
    for( i = 0; i < size; i++ )
        FreeChar1D( s[i] );
    free( ( void * ) s );

    return;
}

int *
AllocInt1D( int size ) {
    int *s = NULL;
    if( ( s = ( int * ) calloc( size, sizeof( int ) ) ) == NULL ) {
        printf( "int * malloc error!\n" );
        exit( 1 );
    }
    return s;
}

void
FreeInt1D( int *s ) {
    if( s )
        free( ( void * ) s );
    return;
}

double *
AllocDbl1D( int size ) {
    double *s = NULL;
    if( ( s = ( double * ) calloc( size, sizeof( double ) ) ) == NULL ) {
        printf( "double * malloc error!\n" );
        exit( 1 );
    }
    return s;
}

void
FreeDbl1D( double *s ) {
    if( s )
        free( ( void * ) s );
    return;
}

double **
AllocDbl2D( int size1, int size2 ) {
    int i;
    double **s = NULL;

    if( ( s = ( double ** ) calloc( size1, sizeof( double * ) ) ) == NULL ) {
        printf( "double ** malloc error!\n" );
        exit( 1 );
    }
    if( size2 > 0 ) {
        for( i = 0; i < size1; i++ )
            s[i] = AllocDbl1D( size2 );
    }

    return s;
}

void
FreeDbl2D( double **s, int size ) {
    int i;

    if( !s )
        return;
    for( i = 0; i < size; i++ )
        FreeDbl1D( s[i] );
    free( ( void * ) s );

    return;
}

double ***
AllocDbl3D( int size1, int size2, int size3 ) {
    int i;
    double ***s = NULL;

    if( ( s = ( double *** ) calloc( size1, sizeof( double ** ) ) ) == NULL ) {
        printf( "double *** malloc error!\n" );
        exit( 1 );
    }
    if( size2 > 0 ) {
        for( i = 0; i < size1; i++ )
            s[i] = AllocDbl2D( size2, size3 );
    }

    return s;
}

void
FreeDbl3D( double ***s, int size1, int size2 ) {
    int i;

    if( !s )
        return;
    for( i = 0; i < size1; i++ )
        FreeDbl2D( s[i], size2 );
    free( ( void * ) s );

    return;
}

int
IsZero( double x ) {
    if( fabs( x ) <= ZeroValue )
        return BoolTrue;
    else
        return BoolFalse;
}
