/**
 * @file savestate.c                                           
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Three function implementations that read, write and remove 
 * the state file for an annealing run.
 * 
 * The frequency with  
 * which state are saved can be chosen by the command line op-   
 * tion -b (for backup stepsize). The state file is very useful  
 * for the case when long annealing runs have to be interrupted  
 * or crash for some reason or another. The run can then be re-  
 * sumed by indicating the state file as an additional argument  
 * to fly_sa.                                                    
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>             /* getopt stuff */

//#include <error.h>
//#include <moves.h>
#include <random.h>
#include "fly_io.h"


#ifdef MPI
#include <MPI.h>                /* for myid */
#endif


/*** A STATIC VARIABLE *****************************************************/
static char *filename;          /* name of state file */


/*** FUNCTION DEFINITIONS **************************************************/

/**  StateRead: reads Lam statistics, move state and erand state from a 
 *              state file and restores the annealer's state to the same   
 *              state it was in before it got interrupted                  
 *     CAUTION: InitMoves must be called before calling StateRead!         
 */
void
StateRead( char *statefile, Opts * options, double *stats, char* rand, double *delta ) {
    int a;                   /* local loop counter */
    FILE *infile;               /* pointer to state file */
    char* line;

    /* make the state file name static to savestate.c */

    filename = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    filename = strcpy( filename, statefile );
    
    line = ( char * ) calloc( MAX_RECORD, sizeof( char ) );


    /* open the state file and read it */

    infile = fopen( filename, "r" );
    if( !infile )
        file_error( "StateRead" );

    options->argv = fgets( options->argv, MAX_RECORD, infile );

    fscanf( infile, "%s\n", options->inname );
    fscanf( infile, "%s\n", options->outname );
    fscanf( infile, "%s\n", options->derivfunc );
    fscanf( infile, "%s\n", options->solver );
    fscanf( infile, "%d\n", &( a ) );
    fscanf( infile, "%d\n", &( options->log_flag ) );
    fscanf( infile, "%d\n", &( options->time_flag ) );
    fscanf( infile, "%ld\n", &( options->state_write ) );
    fscanf( infile, "%ld\n", &( options->print_freq ) );
    fscanf( infile, "%ld\n", &( options->captions ) );
    fscanf( infile, "%d\n", &( options->olddivstyle ) );
    fscanf( infile, "%d\n", &( options->precision ) );
    fscanf( infile, "%lg\n", &( options->stepsize ) );

    if( options->time_flag ) {
        fscanf( infile, "%lf\n", &( delta[0] ) );
        fscanf( infile, "%lf\n", &( delta[1] ) );
    }

    //Keep in nmand if using this function: I disabled dsfmt because for now I don't use it in SSm and it complains when compiling 
    //32 bit. In the future if we will activate state saving, we may again reactivate dsfmt
    strncpy (rand, "", 1);
    while ( strncmp( ( line = fgets( line, MAX_RECORD, infile ) ), "dsfmt_", 6 ) == 0 ) {
        rand = strcat( rand, line );
    }
    fclose( infile );
    free(line);
}

/**  StateWrite: collects Lam statistics, move state and the state of the 
 *               erand48 random number generator and writes all that into 
 *               the state file, which can then be used to restore the run
 *               in case it gets interrupted                              
 */
void
StateWrite( char *statefile ) {
    if( debug ) {
        printf( "Writing state to statefile %s\n", statefile );
    }
    FILE *outfile;              /* state file pointer */
    Opts *options;              /* command line opts to be saved */
    //char *prand;                /* dSFMT() state to be saved as a string */
    double *delta = 0;              /* wallclock and user time to be saved */

    /* if StateWrite() called for the first time: make filename static */

    if( filename == NULL ) {
        filename = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        filename = strcpy( filename, statefile );
    }

    /* collect the state and the options */

    options = GetOptions(  );
    //prand = GetDSFMTState();

    /* write the answer; now *fully* portable, no binary!!! */

    outfile = fopen( filename, "w" );
    if( !outfile )
        file_error( "StateWrite" );

    fprintf( outfile, "%s", options->argv );

    fprintf( outfile, "%s\n", options->inname );
    fprintf( outfile, "%s\n", options->outname );
    fprintf( outfile, "%s\n", options->derivfunc );
    fprintf( outfile, "%s\n", options->solver );
    fprintf( outfile, "%d\n", options->log_flag );
    fprintf( outfile, "%d\n", options->time_flag );
    fprintf( outfile, "%ld\n", options->state_write );
    fprintf( outfile, "%ld\n", options->print_freq );
    fprintf( outfile, "%ld\n", options->captions );
    fprintf( outfile, "%d\n", options->olddivstyle );
    fprintf( outfile, "%d\n", options->precision );
    fprintf( outfile, "%.16g\n", options->stepsize );

    if( options->time_flag ) {
        fprintf( outfile, "%.3f\n", delta[0] );
        fprintf( outfile, "%.3f\n", delta[1] );
    }

    fclose( outfile );

    free( options->derivfunc );
    free( options->solver );
    free( options );
    if( options->time_flag )
        free( delta );
}

/**  StateRm: removes the state file after the run has been completed; 
 *            unless we're tuning in parallel, only the root node needs to  
 *            delete a state file                                          
 */
void
StateRm( void ) {
    if( remove( filename ) )
        warning( "StateRm: could not delete %s", filename );

    free( filename );
}
