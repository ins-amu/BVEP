/**
 *                                                               
 *   @file error.c                                               
 *                                                               
 *****************************************************************
 *                                                               
 *   written by JR, modified by Yoginho                          
 *   thanks to Marcel Wolf for explaining varargs() to me (JJ)   
 *                                                               
 *****************************************************************
 *                                                               
 *   contains error and warning functions plus a function that   
 *   prints out a help or usage message                          
 *                                                               
 *****************************************************************
 *                                                               
 * Copyright (C) 1989-2003 John Reinitz                          
 * the full GPL copyright notice can be found in lsa.c           
 *                                                               
 */

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

const int MAX_ARGS = 25;        /* max number of error() and warning() arguments */

/** ERROR HANDLING FUNCTIONS 
 *   The two routines 'error' and 'warning' print error messages with value.
 *   'error'then exits, while 'warning' returns to the calling function.   
 *                                                                         
 *   Both functions are called analogous to the way you would call printf. 
 *   The only conversion specs that are legal are %c, %d, %f, %g and %s.   
 *   No modifiers (h,l or float precision or size specifiers are allowed). 
 *                                                                         
 *   No % sign means just print msg and exit.                              
 *                                                                         
 */
void
error( const char *format, ... ) {
    int i;                      /* loop counter */
    char *msg;                  /* same as format but with appended /n */
    char *msg_ptr;              /* used to parse format string */
    char **msg_array;           /* holds different parts of the parsed msg */
    char *msg_types;            /* holds types of the messages to be printed */
    size_t length;              /* length of the format string */
    int arg_count = 0;          /* number of additional arguments */
    va_list arg_ptr;            /* pointer to additional arguments */

    /* initialize variable argument pointer */

    va_start( arg_ptr, format );        /* initialize the argument pointer */

    /* allocate memory for message arrays */

    msg_array = ( char ** ) calloc( MAX_ARGS, sizeof( char * ) );
    for( i = 0; i < MAX_ARGS; i++ )
        msg_array[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    msg_types = ( char * ) calloc( MAX_ARGS, sizeof( char ) );

    /* copy the format string to msg and attach a newline; note that we allo-  *
     * cate length+5 bytes for the string; I dunno exactly why that is nece-   *
     * ssary but 5 is the minimum number for which Third Degree will not com-  *
     * plain                                                                   */

    length = strlen( format );
    msg = ( char * ) calloc( ( length + 5 ), sizeof( char ) );
    msg = strcpy( msg, format );
    msg = strcat( msg, "\n" );

    /* finds '%' in msg by reverse search; then copies subportions of msg to   *
     * the msg_array and the corresponding types to msg_types; these arrays    *
     * contain the parts of the error message *in reverse*, from end to be-    *
     * ginning; a new /0 is then inserted at the position of the % and the     *
     * search is repeated until no more % can be found in the format string    *
     * (the upper limit of %'s is given by MAX_ARGS in error.h)                */

    while( ( msg_ptr = strrchr( msg, '%' ) ) != NULL ) {
        if( arg_count >= MAX_ARGS ) {
            fprintf( stderr, "error: too many arguments (max. %d)!\n", MAX_ARGS );
            exit( 1 );
        }
        msg_array[arg_count] = strcpy( msg_array[arg_count], msg_ptr );
        msg_types[arg_count] = *( msg_ptr + 1 );
        ++arg_count;
        *msg_ptr = '\0';
    }

    /* then print error message according to format string; note that only     *
     * self-promoting types are allowed for va_arg(); that's why we use int    *
     * for char, since char would get promoted to an int in the process and    *
     * some compilers (like gcc) have a problem with that                      */

    fprintf( stderr, msg );
    if( arg_count == 0 ) {
        exit( 1 );
    } else
        for( i = 1; i <= arg_count; i++ )
            switch ( msg_types[arg_count - i] ) {
            case 's':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, char * ) );
                break;
            case 'c':
            case 'd':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, int ) );
                break;
            case 'f':
            case 'g':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, double ) );
                break;
            default:
                fprintf( stderr, "error: error in message parsing, bailing out!\n" );
                exit( 1 );
            }

    /* clean up and go home */

    va_end( arg_ptr );

    for( i = 0; i < MAX_ARGS; i++ )
        free( msg_array[i] );
    free( msg_array );
    free( msg_types );
    free( msg );

    exit( 1 );

}

/** ERROR HANDLING FUNCTIONS 
 *   The two routines 'error' and 'warning' print error messages with value.
 *   'error'then exits, while 'warning' returns to the calling function.   
 *                                                                         
 *   Both functions are called analogous to the way you would call printf. 
 *   The only conversion specs that are legal are %c, %d, %f, %g and %s.   
 *   No modifiers (h,l or float precision or size specifiers are allowed). 
 *                                                                                                                                                   
 *   No % sign means just print msg and exit.                              
 *                                                                         
 */
void
warning( const char *format, ... ) {
    int i;                      /* loop counter */
    char *msg;                  /* same as format but with appended /n */
    char *msg_ptr;              /* used to parse format string */
    char **msg_array;           /* holds different parts of the parsed msg */
    char *msg_types;            /* holds types of the messages to be printed */
    size_t length;              /* length of the format string */
    int arg_count = 0;          /* number of additional arguments */
    va_list arg_ptr;            /* pointer to additional arguments */

    /* initialize variable argument pointer */

    va_start( arg_ptr, format );        /* initialize the argument pointer */

    /* allocate memory for message arrays */

    msg_array = ( char ** ) calloc( MAX_ARGS, sizeof( char * ) );
    for( i = 0; i < MAX_ARGS; i++ )
        msg_array[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    msg_types = ( char * ) calloc( MAX_ARGS, sizeof( char ) );

    /* copy the format string to msg and attach a newline; note that we allo-  *
     * cate length+5 bytes for the string; I dunno exactly why that is nece-   *
     * ssary but 5 is the minimum number for which Third Degree will not com-  *
     * plain                                                                   */

    length = strlen( format );
    msg = ( char * ) calloc( ( length + 5 ), sizeof( char ) );
    msg = strcpy( msg, format );
    msg = strcat( msg, "\n" );

    /* finds '%' in msg by reverse search; then copies subportions of msg to   *
     * the msg_array and the corresponding types to msg_types; these arrays    *
     * contain the parts of the error message *in reverse*, from end to be-    *
     * ginning; a new /0 is then inserted at the position of the % and the     *
     * search is repeated until no more % can be found in the format string    *
     * (the upper limit of %'s is given by MAX_ARGS in error.h)                */

    while( ( msg_ptr = strrchr( msg, '%' ) ) != NULL ) {
        if( arg_count >= MAX_ARGS ) {
            fprintf( stderr, "warning: too many arguments (max. %d)!\n", MAX_ARGS );
            exit( 1 );
        }
        msg_array[arg_count] = strcpy( msg_array[arg_count], msg_ptr );
        msg_types[arg_count] = *( msg_ptr + 1 );
        ++arg_count;
        *msg_ptr = '\0';
    }

    /* then print error message according to format string; note that only     *
     * self-promoting types are allowed for va_arg(); that's why we use int    *
     * for char, since char would get promoted to an int in the process and    *
     * some compilers (like gcc) have a problem with that                      */

    fprintf( stderr, msg );
    if( arg_count == 0 ) {
        return;
    } else
        for( i = 1; i <= arg_count; i++ )
            switch ( msg_types[arg_count - i] ) {
            case 's':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, char * ) );
                break;
            case 'c':
            case 'd':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, int ) );
                break;
            case 'f':
            case 'g':
                fprintf( stderr, msg_array[arg_count - i], va_arg( arg_ptr, double ) );
                break;
            default:
                fprintf( stderr, "warning: error in message parsing ... bailing out!\n" );
                exit( 1 );
            }

    /* clean up and go home */

    va_end( arg_ptr );

    for( i = 0; i < MAX_ARGS; i++ )
        free( msg_array[i] );
    free( msg_array );
    free( msg_types );
    free( msg );

    return;

}

/** file_error prints an file handling error using perror(); it only 
 *   needs the name of the calling function as an argument                 
 */
void
file_error( const char *call_name ) {
    perror( call_name );
    exit( 1 );
}

/** PrintMsg: prints a string to stderr, then quits with exit status 0;
 *             useful for printing help and usage messages             
 */
void
PrintMsg( const char *msg, int exit_status ) {
    fprintf( stderr, msg );
    exit( exit_status );
}
