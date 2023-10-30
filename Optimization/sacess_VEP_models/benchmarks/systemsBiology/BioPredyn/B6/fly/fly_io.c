/**
 * @file fly_io.c
 * @author JR, modified by Yoginho,                          
 *   modified by Damjan Cicin-Sain in 2011
 * @contact damjan.cicin@crg.es
 * @date Created on May 27, 2010, 12:07 PM
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 * 
 * @brief Implementation of input/output functions for the SimAnn code.                                                 
 */

#include <string.h>
#include <math.h>
#include <mathLib.h>

#include "fly_io.h"

//
// INPUT
//

/** Eliminate newline characters */
void
chomp( char *s ) {
    while( *s && *s != '\n' && *s != '\r' ) {
        s++;
    }
    // terminate a char string with a zero
    *s = 0;
}

/** @brief ReadDivTimes: read divison times from file */
/** 
     * MITOSIS SCHEDULE: hard-wired cell division tables ***********************                                                                         
     * The division schedule is taken from the $times section in the input file; 
     * If there is none, hard-coded times from maternal.c is taken.
     * 
     * We assume that ndivs is the number of divisions that we want to simulate:
     * 
     * i.e. 
     * 
     * - simulation starts from ccycle 14 means ndivs = 0
     * - simulation starts from ccycle 13 means ndivs = 1
     * - simulation starts from ccycle 12 means ndivs = 2
     * ...
     * 
     * the section entryes are defined as following:
     * 
     * <total_divisions> - this is the maximum number of divisions for that time 
     * schedule. total_divisions >= ndivs
     * 
     * <gastrulation_times> - each entry here corresponds to a gastrulation time
     * assuming that time 0 corresponds to the first division that we are 
     * simulating. 
     * 
     * i.e. 
     * gast_time = 50.0 if starting from ccycle 14 (ndivs = 0)
     * gast_time = 71.1 if starting from ccycle 13 (ndivs = 1)
     * 
     * <division_times> - in each line we have the times of all (total_divisions) 
     * divisions assuming that time 0 corresponds to the start of simulation, and
     * time entries are ordered in descending order of the cell cycles 
     * (cc14, cc13, cc12...).
     * 
     * i.e. 
     * for ndivs = 0 we will read the first line. We start the simulation 
     * from the cycle 14, so we assume that the first simulated cycle (cycle 14) 
     * started at time 0 and no other divisions will happen. So all the times 
     * corresponding to prevous (not simulated) divisions are negative.
     * for ndivs = 1 we will read the second line. We start the simulation 
     * from the cycle 13, so we assume that the first simulated cycle (cycle 13) 
     * started at time 0 and only one other division will happen (cc 13 -> cc 14)
     * at time 21.1. So (again) all the times corresponding to prevous 
     * (not simulated) divisions are negative.
     * 
     *  
     * division_durations - those are the division durations (in seconds) for
     * each division starting from the last one (cc 14).
     *                                                                        
     * Here's an example of the $times section
     *
     * $times
       total_divisions
       6
       gastrulation_times
       50.0 71.1 83.5 93 100.8
       division_times 
       0 -21.1 -33.5 -43.0 -51.8 -57.8
       21.1 0.0 -12.4 -21.9 -30.7 -36.7
       33.5 12.4 0.0 -9.5 -18.3 -24.3
       43.0 21.9  9.5 0.0 -8.8 -14.8
       50.8 29.7 17.3 7.8 -1.0 -7.0
       division_durations
       5.1  3.3  3.0 3.3 3.0 3.0
       $$   
     */   
Times
ReadDivTimes( FILE * fp, TheProblem defs ) {
     
    int n = defs.ndivs;         // ndivs
    int i, j;                   // loop counters
    double buf;                 // temporary variable - double
    char *line;                 // temporary variable - string

    Times times;                // structure containing all the times from the $times section and more ...

    // init
    line = ( char * ) malloc( MAX_RECORD * sizeof( char ) );

    // and find that section
    fp = FindSection( fp, "times" );
    if( !fp ) {
        if( debug ) {
            printf( "ReadDivTimes: cannot locate $times section in input file - using default times\n" );
        }
        times = InitTimes( defs );
    } else {
        fgets( line, MAX_RECORD, fp );                  // advance to next line
        fscanf( fp, "%d", &( times.total_divs ) );      // reading total number of divisions
        fscanf( fp, "%*s\n" );                          // advance to next line

        times.div_times = ( double * ) calloc( n, sizeof( double ) );
        times.div_duration = ( double * ) calloc( n, sizeof( double ) );
        times.full_div_times = ( double * ) calloc( times.total_divs, sizeof( double ) );
        times.full_div_durations = ( double * ) calloc( times.total_divs, sizeof( double ) );

        if( n == 0 ) {
            times.div_times = NULL;
            times.div_duration = NULL;
        }

        for( i = 0; i <= n; i++ ) {
            fscanf( fp, "%lg", &buf );
        }
        times.gast_time = buf;                          // reading gastrulation time

        fgets( line, MAX_RECORD, fp );                  // skip the rest of the line
        fgets( line, MAX_RECORD, fp );                  // advance to next data line

        for( j = 0; j < times.total_divs; j++ ) {
            if( j == n ) {
                for( i = 0; i < times.total_divs; i++ ) {
                    fscanf( fp, "%lg", &( times.full_div_times[i] ) );  // read division times
                    if( i < n ) {
                        times.div_times[i] = times.full_div_times[i];   // fill div_times array
                    }
                }
            }
            fgets( line, MAX_RECORD, fp );              // advance to the right data line
        }

        fgets( line, MAX_RECORD, fp );                  // advance to next line, skipping text
        for( i = 0; i < times.total_divs; i++ ) {
            fscanf( fp, "%lg", &( times.full_div_durations[i] ) );      // read division durations
            if( i < n ) {
                times.div_duration[i] = times.full_div_durations[i];    // fill div_times array
            }
        }
    }

    free( line );
    return times;
}

/** Here we read initial state from the input file (for a gene circuit) */
EqParms
ReadParameters( FILE * fp, TheProblem defs, char *section_title ) {    
    
    EqParms l_parm;                 // local copy of EqParm struct
    double *tempparm, *tempparm1;   // temporary array to read parms

    int i;                      // local loop counter
    int c;                      // used to parse text lines
    int lead_punct;             // flag for leading punctuation
    int linecount = 0;          // keep track of # of lines read
    int Tcount = 0;             // keep track of T lines read
    int Ecount = 0;             // keep track of E lines read

    char *base;                 // pointer to beginning of line string
    char *record;               // string for reading whole line of params

    char **fmt, **fmt1;         // array of format strings for reading params
    char *skip, *skip1;         // string of values to be skipped

    const char read_fmt[] = "%lg";      // read a double
    const char skip_fmt[] = "%*lg ";    // ignore a double

    base = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    skip = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    skip1 = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    fmt = ( char ** ) calloc( defs.ngenes, sizeof( char * ) );
    fmt1 = ( char ** ) calloc( defs.egenes, sizeof( char * ) );

    tempparm = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    tempparm1 = ( double * ) calloc( defs.egenes, sizeof( double ) );

    // create format strings according to the number of genes
    for( i = 0; i < defs.ngenes; i++ ) {
        fmt[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt[i] = strcpy( fmt[i], skip );
        fmt[i] = strcat( fmt[i], read_fmt );
        skip = strcat( skip, skip_fmt );
    }

    // create format strings according to the number of external inputs
    for( i = 0; i < defs.egenes; i++ ) {
        fmt1[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt1[i] = strcpy( fmt1[i], skip1 );
        fmt1[i] = strcat( fmt1[i], read_fmt );
        skip1 = strcat( skip1, skip_fmt );
    }

    // initialize the EqParm struct
    l_parm.R = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.T = ( double * ) calloc( defs.ngenes * defs.ngenes, sizeof( double ) );
    l_parm.E = ( double * ) calloc( defs.ngenes * defs.egenes, sizeof( double ) );
    l_parm.m = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.h = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        l_parm.d = ( double * ) malloc( sizeof( double ) );
    } else {
        l_parm.d = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    }
    l_parm.lambda = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.tau = ( double * ) calloc( defs.ngenes, sizeof( double ) );

    fp = FindSection( fp, section_title );      // find input section
    if( !fp )
        error( "ReadParameters: cannot locate %s section", section_title );

    while( strncmp( ( base = fgets( base, MAX_RECORD, fp ) ), "$$", 2 ) ) {

        record = base;
        lead_punct = 0;

        c = ( int ) *record;

        while( c != '\0' ) {

            if( isdigit( c ) ) {        // line contains data
                record = base;
                // usually read ngenes parameters, but for diff. schedule A or C only read *
                // one d parameter

                if( ( linecount == 5 ) && ( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) ) {
                    if( 1 != sscanf( record, fmt[0], &tempparm[0] ) )
                        error( "ReadParameters: error reading parms" );
                } else if( linecount == 2 ) {
                    for( i = 0; i < defs.egenes; i++ ) {
                        if( 1 != sscanf( record, fmt1[i], &tempparm1[i] ) )
                            error( "ReadParameters: error reading parms" );
                    }
                } else {
                    for( i = 0; i < defs.ngenes; i++ ) {
                        if( 1 != sscanf( record, fmt[i], &tempparm[i] ) )
                            error( "ReadParameters: error reading parms" );
                    }
                }
                switch ( linecount ) {  // copy read parameters into the right array
                case 0: // R
                    for( i = 0; i < defs.ngenes; i++ ) {
                        l_parm.R[i] = tempparm[i];
                    }
                    linecount++;
                    break;
                case 1: // T: keep track of read lines with Tcount
                    for( i = 0; i < defs.ngenes; i++ )
                        l_parm.T[i + Tcount * defs.ngenes]
                            = tempparm[i];
                    Tcount++;
                    if( Tcount == defs.ngenes )
                        linecount++;
                    break;
                case 2: // E: keep track of read lines with Ecount
                    for( i = 0; i < defs.egenes; i++ )
                        l_parm.E[i + Ecount * defs.egenes]
                            = tempparm1[i];
                    Ecount++;
                    if( Ecount == defs.ngenes )
                        linecount++;
                    break;
                case 3: // m
                    for( i = 0; i < defs.ngenes; i++ )
                        l_parm.m[i] = tempparm[i];
                    linecount++;
                    break;
                case 4:
                    for( i = 0; i < defs.ngenes; i++ )  // h
                        l_parm.h[i] = tempparm[i];
                    linecount++;
                    break;
                case 5: // d: consider diffusion schedule
                    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
                        l_parm.d[0] = tempparm[0];
                    } else {
                        for( i = 0; i < defs.ngenes; i++ )
                            l_parm.d[i] = tempparm[i];
                    }
                    linecount++;
                    break;
                case 6: // lambda
                    for( i = 0; i < defs.ngenes; i++ ) {
                        l_parm.lambda[i] = tempparm[i];
                        l_parm.lambda[i] = log( 2. ) / l_parm.lambda[i];
                    } // conversion done here?
                    linecount++;
                    break;
                case 7: // tau: translational/transcription delay
                    for( i = 0; i < defs.ngenes; i++ )
                        l_parm.tau[i] = tempparm[i];
                    linecount++;
                    break;
                default:
                    error( "ReadParameters: too many lines in parameter section" );
                }
                break; // don't do rest of loop anymore!
            } else if( isalpha( c ) ) { // letter means comment
                break;
            } else if( c == '-' ) {     // next two else-ifs for punct
                if( ( ( int ) *( record + 1 ) ) == '.' )
                    record++;
                lead_punct = 1;
                c = ( int ) *( ++record );
            } else if( c == '.' ) {
                lead_punct = 1;
                c = ( int ) *( ++record );
            } else if( ispunct( c ) ) { // other punct means comment
                break;
            } else if( isspace( c ) ) { // ignore leading white space
                if( lead_punct )        // white space after punct means comment
                    break;
                else {
                    c = ( int ) *( ++record );  // get next character in record
                }
            } else {
                error( "ReadParameters: illegal character in %s" );
            }
        }
    }

    free( tempparm );
    free( tempparm1 );
    free( base );
    free( skip );
    free( skip1 );

    for( i = 0; i < defs.ngenes; i++ )
        free( fmt[i] );
    free( fmt );

    for( i = 0; i < defs.egenes; i++ )
        free( fmt1[i] );
    free( fmt1 );

    return l_parm;
}

/** Here we read initial state (and further states) from the x variable
    NOTE: this is a Matlab-specific function 
    IMPORTANT: here x is the list of variable parameters ONLY (parameters for which tweak/mask == 1)
 
 */
EqParms
ReadParametersX( double *x, int *mask, EqParms *iparm, TheProblem defs ) {
        
    EqParms l_parm;             // local copy of EqParm struct
    int i, j, k;                   // local loop counter
    //int len_x;                // length of the input (parameter) array

    // initialize the EqParm struct
    l_parm.R = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.T = ( double * ) calloc( defs.ngenes * defs.ngenes, sizeof( double ) );
    l_parm.E = ( double * ) calloc( defs.ngenes * defs.egenes, sizeof( double ) );
    l_parm.m = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.h = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        l_parm.d = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    } else {
        l_parm.d = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    }
    l_parm.lambda = ( double * ) calloc( defs.ngenes, sizeof( double ) );
    l_parm.tau = ( double * ) calloc( defs.ngenes, sizeof( double ) );

    // len_x = sizeof( x ) / sizeof( double );
    j = 0;
    k = 0;
    for( i = 0; i < defs.ngenes; i++ ) { // R 
        l_parm.R[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->R[i];
        k+=mask[j + i];
    }
    j += defs.ngenes;

    for( i = 0; i < defs.ngenes * defs.ngenes; i++ ) { // T 
        l_parm.T[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->T[i];
        k+=mask[j + i];
    }
    j += defs.ngenes * defs.ngenes;

    for( i = 0; i < defs.egenes * defs.egenes; i++ ) { // E 
        l_parm.E[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->E[i];
        k+=mask[j + i];
    }
    j += defs.egenes * defs.egenes;
        
    for( i = 0; i < defs.ngenes; i++ ) { // m 
        l_parm.m[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->m[i];
        k+=mask[j + i];
    }
    j += defs.ngenes;
 
    for( i = 0; i < defs.ngenes; i++ ) { // h 
        l_parm.h[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->h[i];
        k+=mask[j + i];
    }
    j += defs.ngenes;
    
    // usually read ngenes parameters, but for diff. schedule A or C only read
    // one d parameter
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) { // d          
        l_parm.d[0] = mask[j]*x[k] + (1 - mask[j])*iparm->d[0];
        k+=mask[j];
    } else {
        
        for( i = 0; i < defs.ngenes; i++ ) {
            l_parm.d[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->d[i];
            k+=mask[j + i];
        }
    }
    j += defs.ngenes;

    for( i = 0; i < defs.ngenes; i++ ) { // lambda 
        l_parm.lambda[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->lambda[i];        
        l_parm.lambda[i] = log( 2. ) / l_parm.lambda[i];
        k+=mask[j + i];        
    }
    j += defs.ngenes;

    for( i = 0; i < defs.ngenes; i++ ) { // tau 
        l_parm.tau[i] = mask[j + i]*x[k] + (1 - mask[j + i])*iparm->tau[i];
        k+=mask[j + i];
    }

    return l_parm;
}


/** @brief A function that writes parameters into the data file */
/** WriteParameters: writes the out_parm struct into a new section in the 
 *                    file specified by filename; the new 'eqparms' sec-   
 *                    tion is inserted right after the 'input' section;    
 *                    to achieve this, we need to write to a temporary     
 *                    file which is then renamed to the output file name   
 *              NOTE: lambdas are converted back into protein half lives!! 
 */
void
WriteParameters( char *filename, EqParms * p, char *title, int ndigits, TheProblem defs ) {
    char *temp;                 /* temporary file name */
    char *record;               /* record to be read and written */
    char *record_ptr;           /* pointer used to remember record for 'free' */
    char *saverec;              /* used to save following section title */
    char *shell_cmd;            /* used by 'system' below */

    FILE *outfile;              /* name of output file */
    FILE *tmpfile;              /* name of temporary file */

    temp = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    record = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    saverec = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    record_ptr = record;        /* this is to remember record for 'free' */

    /* open output and temporary file */

    outfile = fopen( filename, "r" );   /* open outfile for reading */
    if( !outfile )              /* sorry for the confusion! */
        error( "WriteParameters: error opening output file" );
    temp = strcpy( temp, "parmXXXXXX" );        /* required by mkstemp() */
    if( mkstemp( temp ) == -1 ) /* get unique name for temp file */
        error( "WriteParameters: error creating temporary file" );

    tmpfile = fopen( temp, "w" );       /* ... and open it for writing */
    if( !tmpfile )
        error( "WriteParameters: error opening temporary file" );

    if( FindSection( outfile, title ) ) {       /* erase section if already present */
        fclose( outfile );      /* this is a little kludgey but */
        KillSection( filename, title ); /* since KillSection needs to be in */
        outfile = fopen( filename, "r" );       /* total control of the file */
    }
    rewind( outfile );

    /* the following two loops look for the appropriate file position to write */
    /* the eqparms section (alternatives are input and eqparms)                */

    if( !strcmp( title, "input" ) ) {
        while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$genotypes", 10 ) )
            fputs( record, tmpfile );
    } else if( !strcmp( title, "eqparms" ) ) {
        while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$input", 6 ) )
            fputs( record, tmpfile );
    }

    fputs( record, tmpfile );

    while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$$", 2 ) )
        fputs( record, tmpfile );

    fputs( record, tmpfile );

    do {
        record = fgets( record, MAX_RECORD, outfile );
        if( !record )
            break;
    } while( strncmp( record, "$", 1 ) );

    fputs( "\n", tmpfile );

    if( record )
        saverec = strcpy( saverec, record );

    /* now we write the eqparms section into the tmpfile */
    PrintParameters( tmpfile, p, title, ndigits, defs );

    fprintf( tmpfile, "\n" );

    /* ... and then write all the rest */
    if( record )
        fputs( saverec, tmpfile );

    while( ( record = fgets( record, MAX_RECORD, outfile ) ) )
        fputs( record, tmpfile );

    fclose( outfile );
    fclose( tmpfile );
    /* rename tmpfile into new file */

    sprintf( shell_cmd, "cp -f %s %s", temp, filename );

    if( -1 == system( shell_cmd ) )
        error( "WriteParameters: error renaming temp file %s" );

    if( remove( temp ) )
        warning( "WriteParameters: temp file %s could not be deleted", temp );

    /* clean up */
    free( temp );
    free( record_ptr );
    free( saverec );
    free( shell_cmd );
}

/** PrintParameters: prints an eqparms section with 'title' to the stream   
 *  indicated by fp 
 */
void
PrintParameters( FILE * fp, EqParms * p, char *title, int ndigits, TheProblem defs ) {
    int i, j;                   /* local loop counters */
    double lambda_tmp;          /* temporary var for lambda */


    fprintf( fp, "$%s\n", title );
    fprintf( fp, "promoter_strengths:\n" );     /* Rs are written here */

    for( i = 0; i < defs.ngenes; i++ ) {
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->R[i] );
    }

    fprintf( fp, "\n" );
    fprintf( fp, "genetic_interconnect_matrix:\n" );    /* Ts written here */

    for( i = 0; i < defs.ngenes; i++ ) {
        for( j = 0; j < defs.ngenes; j++ )
            fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->T[( i * defs.ngenes ) + j] );
        fprintf( fp, "\n" );
    }

    fprintf( fp, "external_input_strengths:\n" );       /* Es written here */

    for( i = 0; i < defs.ngenes; i++ ) {
        for( j = 0; j < defs.egenes; j++ )
            fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->E[( i * defs.egenes ) + j] );
        fprintf( fp, "\n" );
    }

    fprintf( fp, "maternal_connection_strengths:\n" );  /* ms written here */

    for( i = 0; i < defs.ngenes; i++ )
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->m[i] );

    fprintf( fp, "\n" );
    fprintf( fp, "promoter_thresholds:\n" );    /* hs are written here */

    for( i = 0; i < defs.ngenes; i++ )
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->h[i] );

    fprintf( fp, "\n" );
    fprintf( fp, "diffusion_parameter(s):\n" ); /* ds are written here */

    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) )
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->d[0] );
    else
        for( i = 0; i < defs.ngenes; i++ )
            fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->d[i] );

    fprintf( fp, "\n" );
    fprintf( fp, "protein_half_lives:\n" );     /* lambdas are written here */

    for( i = 0; i < defs.ngenes; i++ ) {
        lambda_tmp = log( 2. ) / p->lambda[i];  /* conversion done here */
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, lambda_tmp );
    }

    fprintf( fp, "\n" );
    fprintf( fp, "translational_transcriptional_delays:\n" );   /* taus are written here */

    for( i = 0; i < defs.ngenes; i++ )
        fprintf( fp, "%*.*f ", ndigits + 4, ndigits, p->tau[i] );

    fprintf( fp, "\n$$\n" );
    fflush( fp );

}

/*
void 
WriteScore(char *filename,  double theScore, char *title,  int ndigits)
{
  char   *temp;                                     // temporary file name 
  char   *record;                         // record to be read and written 
  char   *record_ptr;         // pointer used to remember record for 'free' 
  char   *saverec;                 // used to save following section title 
  char   *shell_cmd;                             // used by 'system' below 

  FILE   *outfile;                                  // name of output file 
  FILE   *tmpfile;                               // name of temporary file 


  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  saverec   = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            // this is to remember record for 'free' 

// open output and temporary file 

  outfile = fopen(filename, "r");              // open outfile for reading 
  if ( !outfile )                              // sorry for the confusion! 
    error("WriteScore error opening %s", filename);

  temp = strcpy(temp,"aparmXXXXXX");              // required by mkstemp() 
  if ( mkstemp(temp) == -1 )              // get unique name for temp file 
    error("WriteScore: error creating temporary file name");

  tmpfile = fopen(temp, "w");               // ... and open it for writing 
  if ( !tmpfile )
    error("WriteScore: error opening temporary file %s", temp);

  if ( FindSection(outfile,  title) ) {
    fclose(outfile);                   // erase section if already present 
    KillSection(filename,  title);
    outfile = fopen(filename, "r");
  }
  rewind(outfile);

// the follwoing three loops look for the appropriate file position to     
// write the annealing_output section                                      

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile),
          "$local_search_parameters", 23) )
      fputs(record, tmpfile);
  fputs(record, tmpfile);

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), "$$", 2) )
    fputs(record, tmpfile);
  fputs(record, tmpfile);

  do {
    record = fgets(record, MAX_RECORD, outfile);
    if ( !record ) break;
  } while ( strncmp(record, "$", 1) );

  fputs("\n", tmpfile);

  if ( record )
    saverec = strcpy(saverec, record);

// now we write the eqparms section into the tmpfile 

  fprintf(tmpfile, "$%s\n", title);
  fprintf(tmpfile, "%*.*f\n", ndigits+4, ndigits, theScore);
  fprintf(tmpfile, "$$\n");
  fprintf(tmpfile, "\n");

// ... and then write all the rest 

  if ( record )
    fputs(saverec, tmpfile);

  while ( (record=fgets(record, MAX_RECORD, outfile)) )
    fputs(record, tmpfile);

  fclose(outfile);
  fclose(tmpfile);

// rename tmpfile into new file 

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("WriteAParameters: error renaming temp file %s", temp);

  if ( remove(temp) )
    warning("WriteAParameters: temp file %s could not be deleted", temp);

// clean up 

  free(record_ptr);
  free(saverec);
  free(temp);
  free(shell_cmd);
}
*/

/*FUNCTIONS THAT READ DATA FROM FILE INTO STRUCTS                 *********/

/** ReadTheProblem: reads the problem section of a data file into the 
                     TheProblem struct.                                    
 */
TheProblem
ReadTheProblem( FILE * fp ) {
    TheProblem l_defs;          /* local copy of TheProblem struct */
    fp = FindSection( fp, "problem" );  /* find problem section */
    if( !fp )
        error( "ReadTheProblem: cannot locate problem section" );
    fscanf( fp, "%*s\n" );      /* advance pointer past first text line */
    
    if( 1 != ( fscanf( fp, "%d\n", &l_defs.ngenes ) ) ) /* read # of genes */
        error( "ReadTheProblem: error reading problem section (ngenes)" );
    fscanf( fp, "%*s\n" );      /* advance the pointer an extra line */

    if( 1 != ( fscanf( fp, "%d\n", &l_defs.egenes ) ) ) /* read # of external inputs */
        error( "ReadTheProblem: error reading problem section (egenes)" );

    fscanf( fp, "%*s\n" );      /* advance the pointer past the second text line */

    l_defs.gene_ids = ( char * ) calloc( l_defs.ngenes + 1, sizeof( char ) );
    if( 1 != ( fscanf( fp, "%s\n", l_defs.gene_ids ) ) )        /* read geneID string */
        error( "ReadTheProblem: error reading problem section (gene_ids)" );

    fscanf( fp, "%*s\n" );      /* advance the pointer past the second text line */

    l_defs.egene_ids = ( char * ) calloc( l_defs.egenes + 1, sizeof( char ) );
    if( 1 != ( fscanf( fp, "%s\n", l_defs.egene_ids ) ) )       /* read geneID string */
        error( "ReadTheProblem: error reading problem section (egene_ids)" );

    fscanf( fp, "%*s\n" );      /* next line (ignore comment) */

    if( 1 != ( fscanf( fp, "%d\n", &l_defs.ndivs ) ) )  /* read # of cell divs */
        error( "ReadTheProblem: error reading problem section (ndivs)" );

    fscanf( fp, "%*s\n" );      /* advance the pointer past the third text line */
    if( 1 != ( fscanf( fp, "%d\n", &l_defs.nnucs ) ) )  /* read the max # of nucs (at cycle 14) */
        error( "ReadTheProblem: error reading problem section (nnucs)" );

    fscanf( fp, "%*s\n" );      /* advance the pointer once more */

    if( 1 != ( fscanf( fp, "%c\n", &l_defs.diff_schedule ) ) )  /* read diff sched */
        error( "ReadTheProblem: error reading problem section (diff. schedule)" );
    return l_defs;
}

/** ReadGenotypes: This function reads all the genotypes in a datafile & 
 *                  returns an SList with genotype number and pointers to 
 *                  the corresponding section titles for bias, bcd & facts 
 */
Slist *
ReadGenotypes( FILE * fp, int ngenes ) {
    /* Buffer strings for reading:       */
    char biasbuf[MAX_RECORD];   /* section title of bias section     */
    char factsbuf[MAX_RECORD];  /* section title of data section     */
    char matbuf[MAX_RECORD];    /* section title of bcd section      */
    char histbuf[MAX_RECORD];   /* section title of history section      */
    char extbuf[MAX_RECORD];    /* section title of external input section      */
    char weightsbuf[MAX_RECORD];        /* section title of weights */
    char gtbuf[MAX_RECORD];     /* genotype string                   */

    char *record;               /* pointer to current data record    */

    Slist *current;             /* holds current element of Slist    */
    Slist *first;               /* pointer to first element of Slist */
    Slist *last = NULL;         /* pointer to last element of Slist  */

    int n;
    /*** open the data file and locate genotype section ************************/

    fp = FindSection( fp, "genotypes" );
    if( !fp )
        error( "ReadGenotypes: cannot locate genotypes" );
    if( !( record = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
        error( "ReadGenotypes: error allocating record" );
    first = init_Slist(  );
    current = first;
    /*** read titles of bias, data and bcd sections and genotype number for    *
     *   each genotype *********************************************************/

    while( strncmp( ( record = fgets( record, MAX_RECORD, fp ) ), "$$", 2 ) ) {
        n = sscanf( record, "%s %s %s %s %s %s %s", biasbuf, factsbuf, matbuf, histbuf, extbuf, weightsbuf, gtbuf );
        if( ( n < 6 ) || ( n > 7 ) )
            error( "ReadGenotypes: error reading %s", record );
        else if( n == 6 ) {     // if there are 6 words instead of 7 in the genotype line, we assume that the 
            strcpy( gtbuf, weightsbuf );        // weights section title is missing and that we are reading one of the older 
            strcpy( weightsbuf, "" );   // input files
        }

        /* we eventually want to get rid of the hard wired genotype identifiers    *
         * and replace them by some kind of mechanism by which we can include any  *
         * genes we want in a scoring or annealing run. Think about this!          */
        if( strlen( gtbuf ) != ngenes )
            error( "ReadGenotypes: bad genotype string %s (does not match ngenes)", gtbuf );

        if( !current )
            current = init_Slist(  );
        if( !( current->bias_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating bias_section" );
        current->bias_section = strcpy( current->bias_section, biasbuf );

        if( !( current->fact_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating fact_section" );
        current->fact_section = strcpy( current->fact_section, factsbuf );
        if( !( current->bcd_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating bcd_section" );
        current->bcd_section = strcpy( current->bcd_section, matbuf );

        if( !( current->hist_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating hist_section" );
        current->hist_section = strcpy( current->hist_section, histbuf );

        if( !( current->ext_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating ext_section" );
        current->ext_section = strcpy( current->ext_section, extbuf );

        if( !( current->weights_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating weights_section" );
        current->weights_section = strcpy( current->weights_section, weightsbuf );

        if( !( current->genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGenotypes: error allocating genotype string" );
        current->genotype = strcpy( current->genotype, gtbuf );
        addto_Slist( last, current );

        last = current;
        current = NULL;
    }
    free( record );
    return first;
}

/** ReadBicoid: reads the bcd section of a data file into a linked list; 
 *               also determines maxconc from the bicoid gradient          
 */
Blist *
ReadBicoid( FILE * fp, char *section ) {
    int c;                      /* holds char for parser          */
    int lead_punct;             /* flag for leading punctuation   */

    double maxv = -1.;          /* maximum v (protein conc.)      */

    char *base;                 /* pointer to beginning of string */
    char *record;               /* pointer to string, used as     */
    /* counter                        */
    Blist *current;             /* holds current element of Blist */
    Blist *inlist;              /* holds whole read Blist         */

    if( ( fp = FindSection( fp, section ) ) ) { /* position the fp */
        base = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        current = NULL;
        inlist = NULL;

        /* while loop: reads and processes lines from file until sections ends **** */

        while( strncmp( ( base = fgets( base, MAX_RECORD, fp ) ), "$$", 2 ) ) {

            record = base;      /* base always points to start */
            lead_punct = 0;     /* of string */

            /* for loop: parses and processes each line from the data file ************ */
            //record is a pointer to a string (line of text) while c is the same address, pointing to an int (first word)
            c = ( int ) *record;
            while( c != '\0' ) {

                if( isdigit( c ) ) {    /* number means data */

                    record = base;      /* reset pointer to start of str */
                    current = init_Blist(  );   /* allocate memory for lnkd list */

                    if( 2 != sscanf( record, "%d %lg", &( current->lineage ), &( current->conc ) ) ) {
                        error( "ReadBicoid: error reading %s", base );
                    }

                    if( current->conc > maxv )
                        maxv = current->conc;

                    inlist = addto_Blist( inlist, current );
                    break;

                } else if( isalpha( c ) ) {     /* letter means comment */
                    break;
                } else if( c == '-' ) { /* next two elsifs for punct */
                    if( ( ( int ) *( record + 1 ) ) == '.' )
                        record++;
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( c == '.' ) {
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( ispunct( c ) ) {     /* other punct means comment */
                    break;
                } else if( isspace( c ) ) {     /* ignore leading white space */
                    if( lead_punct )    /* white space after punct means */
                        break;  /* comment */
                    else {
                        c = ( int ) *( ++record );      /* get next character in record */
                    }
                } else {
                    error( "ReadBicoid: illegal character in %s", base );
                }
            }
        }

        if( maxv > 12. )        /* oldstyle or newstyle data? */
            maxconc = 255.;
        else
            maxconc = 12.;

        free( base );
        return inlist;

    } else {

        return NULL;
    }
}

/** ReadData: reads in a data or bias section and puts it in a linked 
 *             list of arrays, one line per array; ndp is used to count
 *             the number of data points in a data file (ndp), which is
 *             used to calculate the root mean square (RMS) if required    
 *                                                                         
 *             ReadData allows comments that start with a letter or punc-  
 *             tuation mark, and treats lines starting with a number or    
 *             .num or -.num as data. It expects to read an int and        
 *             ngenes + 1 data points: lineage, time and ngenes protein    
 *             concentrations. It expects data to be in increasing spatial 
 *             and temporal order.                                         
 */
Dlist *
ReadData( FILE * fp, char *section, int *ndp, TheProblem * defs ) {
    int c;                      /* holds char for parser            */
    int lead_punct;             /* flag for leading punctuation     */
    int i;                      /* loop counter                     */

    char *base;                 /* pointer to beginning of string   */
    char *record;               /* pointer to string, used as       */
    /* counter                          */
    Dlist *current;             /* holds current element of Dlist   */
    Dlist *inlist;              /* holds whole read Dlist           */

    /* the following chars are used for parsing lines of data (record)         */
    /* using sscanf                                                            */

    char *fmt = NULL;           /* used as format for sscanf        */
    char *skip = NULL;          /* skip this stuff!                 */

    const char init_fmt[] = "%*d ";     /* initial format string - * means ignore */
    const char skip_fmt[] = "%*lg ";    /* skip one more float              */
    const char read_fmt[] = "%lg ";     /* read one more float              */

    if( ( fp = FindSection( fp, section ) ) ) { /* position the fp */

        base = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        skip = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

        current = NULL;
        inlist = NULL;

        /* while loop: reads and processes lines from file until sections ends **** */

        while( strncmp( ( base = fgets( base, MAX_RECORD, fp ) ), "$$", 2 ) ) {

            record = base;      /* base always points to start   */
            lead_punct = 0;     /* of string                     */

            /* while loop: parses and processes each line from the data file ************ */

            c = ( int ) *record;

            while( c != '\0' ) {

                if( isdigit( c ) ) {    /* number means data */

                    record = base;      /* reset pointer to start of str */
                    current = init_Dlist( defs->ngenes + 1 );

                    if( 1 != sscanf( record, "%d ", &( current->lineage ) ) )
                        error( "ReadData: error reading %s", base );

                    /* the following loop reads a line of data of variable length into a d-    */
                    /* array, skipping lineage number and all previously read data values      */

                    skip = strcpy( skip, init_fmt );    /* format to be skipped by sscanf */

                    for( i = 0; i < ( defs->ngenes + 1 ); i++ ) {

                        fmt = strcpy( fmt, skip );      /* all this stuff here is to get */
                        fmt = strcat( fmt, read_fmt );  /* the fmt string for sscanf to */
                        skip = strcat( skip, skip_fmt );        /* the correct format */

                        if( 1 != sscanf( record, ( const char * ) fmt, &( current->d[i] ) ) )
                            error( "ReadData: error reading %s", base );

                        /* update number of data points */
                        if( ( i != 0 ) && ( current->d[i] != IGNORE ) ) {
                            ( *ndp )++;
                        }
                    }
                    /* now add this to the lnkd list */
                    inlist = addto_Dlist( inlist, current );
                    break;
                } else if( isalpha( c ) ) {     /* letter means comment */
                    break;
                } else if( c == '-' ) { /* next two elsifs for punct */
                    if( ( ( int ) *( record + 1 ) ) == '.' )
                        record++;
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( c == '.' ) {
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( ispunct( c ) ) {     /* other punct means comment */
                    break;
                } else if( isspace( c ) ) {     /* ignore leading white space */
                    if( lead_punct )    /* white space after punct means */
                        break;  /* comment */
                    else {
                        c = ( int ) *( ++record );      /* get next character in record */
                    }
                } else {
                    error( "ReadData: Illegal character in %s", base );
                }
            }
        }

        free( base );
        free( fmt );
        free( skip );
        //printf("NDP ReadData = %d\n", *ndp);
        return inlist;

    } else {
        return NULL;

    }
}

/** ReadTimes: reads a time table from a file and returns a DArrPtr 
 * FILE FORMAT: one time per line separated with newlines                  
 *        NOTE: max. times tab size is 10000                               
 */
DArrPtr
ReadTimes( char *timefile, Zygote zyg ) {
    FILE *fp;

    int count = 0;
    double gast;

    DArrPtr tabtimes;
    double *tabptr;

    gast = zyg.times.gast_time;

    tabtimes.size = 10000;
    tabtimes.array = ( double * ) calloc( 10000, sizeof( double ) );
    tabptr = tabtimes.array;

    fp = fopen( timefile, "r" );
    if( !fp )
        file_error( "ReadTimes" );

    if( 1 != ( fscanf( fp, "%lg\n", tabptr ) ) )
        error( "ReadTimes: time file %s empty!", timefile );
    tabptr++;
    count++;

    while( ( fscanf( fp, "%lg\n", tabptr ) ) != EOF ) {
        if( ( *tabptr < 0 ) || ( *tabptr <= *( tabptr - 1 ) ) || ( *tabptr > gast ) )
            error( "ReadTimes: invalid time(s) in %s!", timefile );
        tabptr++;
        count++;
    }

    fclose( fp );

    tabtimes.size = count;
    tabtimes.array = ( double * ) realloc( tabtimes.array, count * sizeof( double ) );
    return tabtimes;
}

/** ReadGuts: reads the $gutsdefs section in a data file into an array 
 *             of strings which then need to get parsed                    
 */
char **
ReadGuts( FILE * fp ) {
    char **gutsbuf = NULL;      /* buffer strings for reading */
    char *record = NULL;
    int i;                      /* loop counters */

    if( !( gutsbuf = ( char ** ) calloc( MAX_RECORD, sizeof( char * ) ) ) )
        error( "ReadGuts: error allocating memory for gutsdefs" );

    if( !( record = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
        error( "ReadGuts: error allocating record" );

    /*** locate gutsdefs section ***********************************************/

    fp = FindSection( fp, "gutsdefs" );
    if( !fp )
        error( "ReadGuts: cannot locate gutsdefs" );

    /* Remember the location of the first line and count the number of lines * */

    i = 0;
    while( strncmp( ( record = fgets( record, MAX_RECORD, fp ) ), "$$", 2 ) ) {
        if( !( *( gutsbuf + i ) = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "ReadGuts: error allocating memory for gutsdefs strings" );
        *( gutsbuf + i ) = strcpy( *( gutsbuf + i ), record );
        i++;
    }

    if( !( *gutsbuf ) )
        error( "ReadGuts: gutsdefs section is empty!" );

    free( record );

    return gutsbuf;
}

/** ReadInterpData: reads the history section from a data file into an array 
  *  a dedicated structure
  */
Dlist *
ReadInterpData( FILE * fp, char *section, int num_genes, int *ndp ) {
    int c;                      /* holds char for parser            */
    int lead_punct;             /* flag for leading punctuation     */
    int i;                      /* loop counter                     */

    char *base;                 /* pointer to beginning of string   */
    char *record;               /* pointer to string, used as       */
    /* counter                          */
    Dlist *current;             /* holds current element of Dlist   */
    Dlist *inlist;              /* holds whole read Dlist           */

    /* the following chars are used for parsing lines of data (record)         */
    /* using sscanf                                                            */

    char *fmt = NULL;           /* used as format for sscanf        */
    char *skip = NULL;          /* skip this stuff!                 */

    const char init_fmt[] = "%*d ";     /* initial format string            */
    const char skip_fmt[] = "%*lg ";    /* skip one more float              */
    const char read_fmt[] = "%lg ";     /* skip one more float              */

    if( ( fp = FindSection( fp, section ) ) ) { /* position the fp */
        base = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        skip = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

        current = NULL;
        inlist = NULL;

        /* while loop: reads and processes lines from file until sections ends **** */

        while( strncmp( ( base = fgets( base, MAX_RECORD, fp ) ), "$$", 2 ) ) {

            record = base;      /* base always points to start   */
            lead_punct = 0;     /* of string                     */

            /* while loop: parses and processes each line from the data file ************ */

            c = ( int ) *record;

            while( c != '\0' ) {

                if( isdigit( c ) ) {    /* number means data */

                    record = base;      /* reset pointer to start of str */
                    current = init_Dlist( num_genes + 1 );

                    if( 1 != sscanf( record, "%d ", &( current->lineage ) ) )
                        error( "ReadInterpData: error reading %s", base );

                    /* the following loop reads a line of data of variable length into a d-    */
                    /* array, skipping lineage number and all previously read data values      */

                    skip = strcpy( skip, init_fmt );    /* format to be skipped by sscanf */

                    for( i = 0; i < ( num_genes + 1 ); i++ ) {

                        fmt = strcpy( fmt, skip );      /* all this stuff here is to get */
                        fmt = strcat( fmt, read_fmt );  /* the fmt string for sscanf to */
                        skip = strcat( skip, skip_fmt );        /* the correct format */

                        if( 1 != sscanf( record, ( const char * ) fmt, &( current->d[i] ) ) )
                            error( "ReadInterpData: error reading %s", base );

                        /* update number of data points */
                        if( ( i != 0 ) && ( current->d[i] != IGNORE ) ) {
                            ( *ndp )++;
                        }

                    }
                    /* now add this to the lnkd list */
                    inlist = addto_Dlist( inlist, current );
                    break;
                } else if( isalpha( c ) ) {     /* letter means comment */
                    break;
                } else if( c == '-' ) { /* next two elsifs for punct */
                    if( ( ( int ) *( record + 1 ) ) == '.' )
                        record++;
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( c == '.' ) {
                    lead_punct = 1;
                    c = ( int ) *( ++record );
                } else if( ispunct( c ) ) {     /* other punct means comment */
                    break;
                } else if( isspace( c ) ) {     /* ignore leading white space */
                    if( lead_punct )    /* white space after punct means */
                        break;  /* comment */
                    else {
                        c = ( int ) *( ++record );      /* get next character in record */
                    }
                } else {
                    error( "ReadInterpData: Illegal character in %s", base );
                }
            }
        }

        free( base );
        free( fmt );
        free( skip );
        return inlist;

    } else {

        return NULL;

    }
}


/*** FUNCTIONS THAT READ STUFF INTO STRUCTURES *****************************/

/*** ReadLimits: reads the limits section of a data file and returns the  
 *               approriate SearchSpace struct to the calling function     
 */
// Experimental feature: free limits
// @author Anton Crombach
// @date 2012, August

/** ReadRangeElement:
 * Read in an upper or lower bound until a comma, parenthese or whitespace is
 * encountered. Return number of characters read.
 *
 * @note No interpretation is done, only parsing.
 */
int
ReadRangeElement( char *line, int i, char *element ) {
    
    // find first whitespace char, comma or closing parenthese
    char stop[] = " \t\n,)";
    int j = strcspn( line + i, stop );
    // copy the characters inbetween
    strncpy( element, &( line[ i ] ), j );
    return j+1;
}

/**
 * Chops a line of ranges up in an array of \c Range elements
 */
void
ReadLineOfRanges( char *line, struct Range **ranges, int len_ranges ) {
    
    const char *myNAN = "N/A";
    /* a little parser */
    int len = strlen( line );
    int i = 0, j = 0;
    char **tokens = ( char ** ) calloc( 2 * len_ranges,  sizeof( char * ) );
    while( i < len ) {
        switch( line[ i ] ) {
        case ' ':
        case '\t':
        case '\n':
        case ')':
        case ',':
        case '(':
            // whitespace, just skip
            ++i;
            break;
        case '-':
        case '.':
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case 'N':
            tokens[ j ] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            i += ReadRangeElement( line, i, tokens[ j ] );
            ++j;
            break;
        }
    }
    
    /* all tokenized, now they need to be interpreted */
    for( i = 0; i < len_ranges; ++i ) {
        // lower bound
        int ii = 2 * i;
        if( strcmp( tokens[ ii ], myNAN ) == 0 ) {
            ranges[ i ]->lower = -DBL_MAX;
        } else {
            sscanf( tokens[ ii ], "%lg", &( ranges[ i ]->lower ) );
        }
        // upper bound
        if( strcmp( tokens[ ii+1 ], myNAN ) == 0 ) {
            ranges[ i ]->upper = DBL_MAX;
        } else {
            sscanf( tokens[ ii+1 ], "%lg", &( ranges[ i ]->upper ) );
        }

        free( tokens[ ii ] );
        free( tokens[ ii+1 ] );
    }
    free( tokens );
}

/**
 * Set all limits in a line to (-DBL_MAX, DBL_MAX) in case of NA limits
 */
void
CreateLineOfRangesFromNA(struct Range **ranges, int ncols ) {
    
    int i = 0;
    /* all tokenized, now they need to be interpreted */
    for( i = 0; i < ncols; ++i ) {
        // lower bound
        ranges[ i ]->lower = -DBL_MAX;
        // upper bound
        ranges[ i ]->upper = DBL_MAX;
    }
}

/**
 * Set all limits in a matrix to (-DBL_MAX, DBL_MAX) in case of NA limits 
 */
void
CreateMatrixOfRangesFromNA(struct Range **ranges, int ncols, int nlines ) {
    
    int i = 0, j = 0;
    /* all tokenized, now they need to be interpreted */
    for( j = 0; j < nlines; j++ ) {
        for( i = 0; i < ncols; i++ ) {
            // lower bound
            ranges[ j * ncols + i ]->lower = -DBL_MAX;
            // upper bound
            ranges[ j * ncols + i ]->upper = DBL_MAX;
        }
    }
}

/**
 * Allocate memory for the elements in an array of Range.
 *
 * @todo It's worth finding out if I can make it inline C99-style
 */
Range **
mallocRanges( Range **r, int n ) {
    int i;
    for(  i = 0; i < n; ++i ) {
        r[ i ] = ( Range * ) malloc( sizeof( Range ) );
        if( !r[ i ] )
            error( "'mallocRanges' failed" );
        r[ i ]->lower = -1.0;
        r[ i ]->upper = -1.0;
    }
    return r;
}

/**
 * Read a single item of ranges. Used for promoter's strength, maternal
 * contribution, threshold of expression, protein decay and translational delay
 */
Range **
ReadSingleRanges( FILE *fp, char *record, Range **r, int n ) {
    
    /* advance past first title line */
    fscanf( fp, "%*s\n" );
    /* reserve memory */
    r = ( Range ** ) malloc( n * sizeof( Range * ) );
    r = mallocRanges( r, n );
    /* and read them in */
    record = fgets( record, MAX_RECORD, fp );
    if( !strncmp( record, "N/A", 3 )) {
        // If line starts with N/A we are dealing with an old formatted file; 
        // in that case set all limits to (-DBL_MAX, DBL_MAX)
        CreateLineOfRangesFromNA( r, n );  
    } else {
        ReadLineOfRanges( record, r, n );
    }
    return r;
}

/**
 * Read a compound item of ranges. Used for gap binding, gap hill-coefficients,
 * external binding, and external hill-coefficients
 */
Range **
ReadCompoundRanges( FILE * fp, char *record, Range **r, int n, int e ) {
    int i; 
    /* advance past title line */
    fscanf( fp, "%*s\n" );
    /* reserve memory */
    r = ( Range ** ) malloc( n * e * sizeof( Range * ) );
    r = mallocRanges( r, n * e ); 
    /* and fill it */
    for( i = 0; i < n; ++i ) {
        record = fgets( record, MAX_RECORD, fp );
        if ( !strncmp( record, "N/A", 3 ) ) {
            CreateMatrixOfRangesFromNA( r, n, e );
            break;
        }
        ReadLineOfRanges( record, r + i * e, e );
    }
    return r;
}

/**
 * ReadLimits: reads the limits section of a data file and returns the  
 *             approriate SearchSpace struct to the calling function     
 *
 Experimental feature: free limits
 @author Anton Crombach
 @date 2012, August
*/ 
SearchSpace *
ReadLimits( FILE * fp, TheProblem defs ) {
 
    const char *myNAN = "N/A";
    const int ncols = defs.ngenes;
    const int nrows = defs.ngenes;
    const int egenes = defs.egenes;
    int i; 
    
    /* string for reading whole line of limits */
    char *record;
    record = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
    /* struct that needs to be filled */
    SearchSpace *l_limits;
    l_limits = ( SearchSpace * ) malloc( sizeof( SearchSpace ) );
    
    /* find limits section */
    fp = FindSection( fp, "limits" );
    if( !fp )
        error( "ReadLimitsXXX: cannot locate limits section" );

    /* advance past first title line */
    fscanf( fp, "%*s\n" );
    
    /* read Lamda for penalty */
    record = fgets( record, MAX_RECORD, fp );
    /* check if it is not-a-number */
    if( !strncmp( record, myNAN, 3 ) ) {
        /* we are not using penalties, so all ranges are explicitly given, 
         * hence no NANs allowed in the ranges!
         */
        l_limits->pen_vec = NULL;
    } else {
        /* we may have some penalties somewhere, need to read entire section
         * to establish where and when...
         */
        l_limits->pen_vec = ( double * ) calloc( 2 + ncols + egenes, sizeof( double ) );
        if( 1 != sscanf( record, "%lg", l_limits->pen_vec ) )
            error( "ReadLimitsXXX: error reading Lambda for penalty" );
    }
    
    
    /* read Promoter strengths */
    l_limits->Rlim = ReadSingleRanges( fp, record, l_limits->Rlim, ncols );
    /* read T matrix ranges */
    l_limits->Tlim = ReadCompoundRanges( fp, record, l_limits->Tlim, nrows, ncols );
    /*for( int p = 0; p < ncols; ++p ) {
        for( int q = 0; q < nrows; ++q ) {
            printf( "T lim %.3e, %.3e\n", l_limits->Tlim[ q*ncols + p ]->lower,
                l_limits->Tlim[ q*ncols + p ]->upper );
        }
    }*/
    
    /* read E matrix ranges */
    l_limits->Elim = ReadCompoundRanges( fp, record, l_limits->Elim, nrows, egenes );
    /* read Maternal contribution ranges */
    l_limits->mlim = ReadSingleRanges( fp, record, l_limits->mlim, ncols );
    /* read Threshold of constitutive expression ranges */
    l_limits->hlim = ReadSingleRanges( fp, record, l_limits->hlim, ncols );
    
    /* read Diffusion parameter ranges, some ancient code here */
    /* advance past title line */
    fscanf( fp, "%*s\n" );    
    /* reserve memory */
    l_limits->dlim = ( Range ** ) calloc( ncols, sizeof( Range * ) );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        l_limits->dlim[0] = ( Range * ) malloc( sizeof( Range ) );
    } else {
        l_limits->dlim = mallocRanges( l_limits->dlim, ncols );
    }
    /* and read them in */
    record = fgets( record, MAX_RECORD, fp );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        ReadLineOfRanges( record, l_limits->dlim, 1 );
    } else {
        ReadLineOfRanges( record, l_limits->dlim, ncols );
    }

    /* read Gene product decay parameter ranges */
    l_limits->lambdalim = ReadSingleRanges( fp, record, l_limits->lambdalim,
        ncols );
    for( i = 0; i < ncols; ++i ) {
        /* need to swap lower and upper limit due to conversion to half lives */
        double aux = l_limits->lambdalim[ i ]->lower;
        l_limits->lambdalim[ i ]->lower = log( 2. ) / 
            l_limits->lambdalim[ i ]->upper; 
        l_limits->lambdalim[ i ]->upper = log( 2. ) / aux;
    }        
    /* read Transcriptional/translational delay parameter ranges */
    l_limits->taulim = ReadSingleRanges( fp, record, l_limits->taulim, ncols );
    
    free( record );
    return l_limits;
}
//------------------------------------------------------------------------------
// End of experimental feature
//------------------------------------------------------------------------------

/**
 * Old (obsolete) function to read limits
 */
SearchSpace *
ReadLimits_OLD( FILE * fp, TheProblem defs ) {
    SearchSpace *l_limits;

    char *record;               /* string for reading whole line of limits */

    int i, j;                   /* local loop counter */

    char **fmt1;                /* format strings for reading lower limits */
    char **fmt2;                /* format strings for reading upper limits */
    char **efmt1;               /* format strings for reading lower
                                   limits of external inputs */
    char **efmt2;               /* format strings for reading upper
                                   limits of external inputs */
    char *skip;                 /* string of values to be skipped */
    const char read_fmt[] = "%lg";      /* read a double */
    const char skip_fmt1[] = "%*lg, %*lg) (";   /* ignore lower limit */
    const char skip_fmt2[] = "%*lg) (%*lg, ";   /* ignore upper limit */
    l_limits = ( SearchSpace * ) malloc( sizeof( SearchSpace ) );
    record = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
    skip = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
    fmt1 = ( char ** ) calloc( defs.ngenes, sizeof( char * ) );
    fmt2 = ( char ** ) calloc( defs.ngenes, sizeof( char * ) );
    efmt1 = ( char ** ) calloc( defs.egenes, sizeof( char * ) );
    efmt2 = ( char ** ) calloc( defs.egenes, sizeof( char * ) );
    /* create format strings according to the number of genes */

    skip = strcpy( skip, "(" ); /* first for lower limits */

    for( i = 0; i < defs.ngenes; i++ ) {
        fmt1[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt1[i] = strcpy( fmt1[i], skip );
        fmt1[i] = strcat( fmt1[i], read_fmt );
        skip = strcat( skip, skip_fmt1 );
    }

    skip = strcpy( skip, "(%*lg, " );   /* then for upper limits */

    for( i = 0; i < defs.ngenes; i++ ) {
        fmt2[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        fmt2[i] = strcpy( fmt2[i], skip );
        fmt2[i] = strcat( fmt2[i], read_fmt );
        skip = strcat( skip, skip_fmt2 );
    }

    /* create format strings according to the number of external inputs  -
     * this could have been done with just fmt1 and fmt2, picking their
     * size to be whatever was bigger - ngenes or egenes, however this
     * maintains the seperation of reading in normal things and external
     * inputs and avoid confusion */

    skip = strcpy( skip, "(" ); /* first for lower limits */

    for( i = 0; i < defs.egenes; i++ ) {
        efmt1[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        efmt1[i] = strcpy( efmt1[i], skip );
        efmt1[i] = strcat( efmt1[i], read_fmt );
        skip = strcat( skip, skip_fmt1 );
    }

    skip = strcpy( skip, "(%*lg, " );   /* then for upper limits */

    for( i = 0; i < defs.egenes; i++ ) {
        efmt2[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        efmt2[i] = strcpy( efmt2[i], skip );
        efmt2[i] = strcat( efmt2[i], read_fmt );
        skip = strcat( skip, skip_fmt2 );
    }



    /* find limits section and check if penalty or explicit ranges are used    */

    fp = FindSection( fp, "limits" );   /* find limits section */
    if( !fp )
        error( "ReadLimits: cannot locate limits section" );

    fscanf( fp, "%*s\n" );      /* advance past first title line */

    record = fgets( record, MAX_RECORD, fp );   /* read Lamda for penalty */

    if( !strncmp( record, "N/A", 3 ) ) {
        l_limits->pen_vec = NULL;
    } else {
        l_limits->pen_vec = ( double * ) calloc( 2 + defs.ngenes + defs.egenes, sizeof( double ) );
        if( 1 != sscanf( record, "%lg", l_limits->pen_vec ) )
            error( "ReadLimits: error reading Lambda for penalty" );
    }
    // initialize limits struct according to penalty or not
    // first the stuff that's not dependent on penalty

    l_limits->Rlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );
    for( i = 0; i < defs.ngenes; i++ ) {
        l_limits->Rlim[i] = ( Range * ) malloc( sizeof( Range ) );
    }

    l_limits->dlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );

    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        l_limits->dlim[0] = ( Range * ) malloc( sizeof( Range ) );
    } else {
        for( i = 0; i < defs.ngenes; i++ )
            l_limits->dlim[i] = ( Range * ) malloc( sizeof( Range ) );
    }

    l_limits->lambdalim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );

    for( i = 0; i < defs.ngenes; i++ ) {
        l_limits->lambdalim[i] = ( Range * ) malloc( sizeof( Range ) );
    }

    l_limits->taulim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );

    for( i = 0; i < defs.ngenes; i++ ) {
        l_limits->taulim[i] = ( Range * ) malloc( sizeof( Range ) );
    }

    /* ... and now the stuff that depends on whether we use penalty or not     */
    if( l_limits->pen_vec == 0 ) {

        l_limits->Tlim = ( Range ** ) calloc( defs.ngenes * defs.ngenes, sizeof( Range * ) );
        for( i = 0; i < defs.ngenes; i++ ) {
            for( j = 0; j < defs.ngenes; j++ )
                l_limits->Tlim[( i * defs.ngenes ) + j] = ( Range * ) malloc( sizeof( Range ) );
        }
        l_limits->Elim = ( Range ** ) calloc( defs.ngenes * defs.egenes, sizeof( Range * ) );
        for( i = 0; i < defs.ngenes; i++ ) {
            for( j = 0; j < defs.egenes; j++ )
                l_limits->Elim[( i * defs.egenes ) + j] = ( Range * ) malloc( sizeof( Range ) );
        }
        l_limits->mlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );
        for( i = 0; i < defs.ngenes; i++ ) {
            l_limits->mlim[i] = ( Range * ) malloc( sizeof( Range ) );
        }
        l_limits->hlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );
        for( i = 0; i < defs.ngenes; i++ ) {
            l_limits->hlim[i] = ( Range * ) malloc( sizeof( Range ) );
        }
    } else {
        l_limits->Tlim = NULL;
        l_limits->Elim = NULL;
        l_limits->mlim = NULL;
        l_limits->hlim = NULL;
    }

    /* ... and finally read the actual values from the data file               */

    fscanf( fp, "%*s\n" );      /* advance past title line */

    record = fgets( record, MAX_RECORD, fp );   /* loop to read R limits */
    for( i = 0; i < defs.ngenes; i++ ) {
        if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->Rlim[i]->lower ) )
            error( "ReadLimits: error reading promoter strength (R) limits" );
        if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->Rlim[i]->upper ) )
            error( "ReadLimits: error reading promoter strength (R) limits" );
    }

    fscanf( fp, "%*s\n" );

    /* using explicit limits? */
    if( l_limits->pen_vec == 0 ) {
        for( i = 0; i < defs.ngenes; i++ ) {    /* loops to read T matrix limits */
            record = fgets( record, MAX_RECORD, fp );
            for( j = 0; j < defs.ngenes; j++ ) {
                if( 1 != sscanf( record, ( const char * ) fmt1[j], &l_limits->Tlim[( i * defs.ngenes ) + j]->lower ) )
                    error( "ReadLimits:: error reading T matrix limits" );
                if( 1 != sscanf( record, ( const char * ) fmt2[j], &l_limits->Tlim[( i * defs.ngenes ) + j]->upper ) )
                    error( "ReadLimits:: error reading T matrix limits" );
            }
        }

        fscanf( fp, "%*s\n" );

        for( i = 0; i < defs.ngenes; i++ ) {    /* loops to read E matrix limits */
            record = fgets( record, MAX_RECORD, fp );
            for( j = 0; j < defs.egenes; j++ ) {
                if( 1 != sscanf( record, ( const char * ) efmt1[j], &l_limits->Elim[( i * defs.egenes ) + j]->lower ) )
                    error( "ReadLimits:: error reading E matrix limits" );
                if( 1 != sscanf( record, ( const char * ) efmt2[j], &l_limits->Elim[( i * defs.egenes ) + j]->upper ) )
                    error( "ReadLimits:: error reading E matrix limits" );
            }
        }

        fscanf( fp, "%*s\n" );

        record = fgets( record, MAX_RECORD, fp );       /* loop to read m limits */
        for( i = 0; i < defs.ngenes; i++ ) {
            if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->mlim[i]->lower ) )
                error( "ReadLimits: error reading maternal interconnect (m) limits" );
            if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->mlim[i]->upper ) )
                error( "ReadLimits: error reading maternal interconnect (m) limits" );
        }

        fscanf( fp, "%*s\n" );

        record = fgets( record, MAX_RECORD, fp );       /* loop to read h limits */
        for( i = 0; i < defs.ngenes; i++ ) {
            if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->hlim[i]->lower ) )
                error( "ReadLimits: error reading promoter threshold (h) limits" );
            if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->hlim[i]->upper ) )
                error( "ReadLimits: error reading promoter threshold (h) limits" );
        }

        /* using penalty? -> ignore this part of the limit section */

    } else {
        for( i = 0; i < 7; i++ )
            fscanf( fp, "%*s\n" );
    }
    fscanf( fp, "%*s\n" );      /* diffusion paramter limit(s) are read here */

    record = fgets( record, MAX_RECORD, fp );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        if( 1 != sscanf( record, ( const char * ) fmt1[0], &l_limits->dlim[0]->lower ) )
            error( "ReadLimits: error reading diffusion parameter limit (d)" );
        if( 1 != sscanf( record, ( const char * ) fmt2[0], &l_limits->dlim[0]->upper ) )
            error( "ReadLimits: error reading diffusion parameter limit (d)" );
    } else {
        for( i = 0; i < defs.ngenes; i++ ) {
            if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->dlim[i]->lower ) )
                error( "ReadLimits: error reading diffusion parameter limits (d)" );
            if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->dlim[i]->upper ) )
                error( "ReadLimits: error reading diffusion parameter limits (d)" );
        }
    }

    fscanf( fp, "%*s\n" );

    record = fgets( record, MAX_RECORD, fp );   /* loop to read lambda lims */
    for( i = 0; i < defs.ngenes; i++ ) {
        if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->lambdalim[i]->upper ) )
            error( "ReadLimits: error reading lambda limits" );
        l_limits->lambdalim[i]->upper = log( 2. ) / l_limits->lambdalim[i]->upper;
        if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->lambdalim[i]->lower ) )
            error( "ReadLimits: error reading lambda limits" );
        l_limits->lambdalim[i]->lower = log( 2. ) / l_limits->lambdalim[i]->lower;
    }

    fscanf( fp, "%*s\n" );

    record = fgets( record, MAX_RECORD, fp );   /* loop to read tau lims */
    for( i = 0; i < defs.ngenes; i++ ) {
        if( 1 != sscanf( record, ( const char * ) fmt1[i], &l_limits->taulim[i]->lower ) )
            error( "ReadLimits: error reading tau limits" );
        if( 1 != sscanf( record, ( const char * ) fmt2[i], &l_limits->taulim[i]->upper ) )
            error( "ReadLimits: error reading tau limits" );
    }
    free( record );
    free( skip );
    for( i = 0; i < defs.ngenes; i++ ) {
        free( fmt1[i] );
        free( fmt2[i] );
    }

    for( i = 0; i < defs.egenes; i++ ) {
        free( efmt1[i] );
        free( efmt2[i] );
    }

    free( fmt1 );
    free( fmt2 );
    free( efmt1 );
    free( efmt2 );

    return l_limits;
}

/** ReadTweak: reads the tweak array passed from the calling 
 *              function through the mask array. This array has a value of  
 *              1 or 0 for each parameter in the model and is used by       
 *              Translate to create the array of pointers to the            
 *              parameters-to-be-tweaked. If mask array is NULL, it reads   
 *              that information from the input file                        
 */
Tweak
ReadTweak( FILE * fp, int *mask, TheProblem defs ) {
    Tweak l_tweak;              // local Tweak struct
    int *temptweak, *temptweak1;        // temp arrays to read tweaks

    int i, j;                   // local loop counter
    int c;                      // counter
    int linecount = 0;          // keep track of # of lines read
    int Tcount = 0;             // keep track of T lines read
    int Ecount = 0;             // keep track of E lines read

    char *base;                 // pointer to beginning of line string
    char *record;               // string for reading whole line of params

    char **fmt;                 // array of format strings for reading params
    char **fmt1;                // array of format strings for reading E tweaks
    char *skip, *skip1;         // string of values to be skipped

    const char read_fmt[] = "%d";       /* read an int */
    const char skip_fmt[] = "%*d ";     /* ignore an int */

    // initialize the Tweak struct

    l_tweak.Rtweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );
    l_tweak.Ttweak = ( int * ) calloc( defs.ngenes * defs.ngenes, sizeof( int ) );
    l_tweak.Etweak = ( int * ) calloc( defs.ngenes * defs.egenes, sizeof( int ) );
    l_tweak.mtweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );
    l_tweak.htweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );
    if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
        l_tweak.dtweak = ( int * ) malloc( sizeof( int ) );
    } else {
        l_tweak.dtweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );
    }
    l_tweak.lambdatweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );
    l_tweak.tautweak = ( int * ) calloc( defs.ngenes, sizeof( int ) );

    if( mask != NULL ) {        //reading mask from the mask array
        j = 0;
        for( i = 0; i < defs.ngenes; i++ ) {
            l_tweak.Rtweak[i] = mask[j + i];
        }
        j += defs.ngenes;

        for( i = 0; i < defs.ngenes * defs.ngenes; i++ ) {
            l_tweak.Ttweak[i] = mask[j + i];
        }
        j += defs.ngenes * defs.ngenes;

        for( i = 0; i < defs.ngenes * defs.egenes; i++ ) {
            l_tweak.Etweak[i] = mask[j + i];
        }
        j += defs.ngenes * defs.egenes;

        for( i = 0; i < defs.ngenes; i++ ) {
            l_tweak.mtweak[i] = mask[j + i];
        }
        j += defs.ngenes;

        for( i = 0; i < defs.ngenes; i++ ) {
            l_tweak.htweak[i] = mask[j + i];
        }
        j += defs.ngenes;

        if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
            l_tweak.dtweak[0] = mask[j];
        } else {
            for( i = 0; i < defs.ngenes; i++ ) {
                l_tweak.dtweak[i] = mask[j + i];
            }
        }
        j += defs.ngenes;

        for( i = 0; i < defs.ngenes; i++ ) {
            l_tweak.lambdatweak[i] = mask[j + i];
        }
        j += defs.ngenes;

        for( i = 0; i < defs.ngenes; i++ ) {
            l_tweak.tautweak[i] = mask[j + i];
        }
        j += defs.ngenes;
    } else {                    //reading mask from file

        base = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
        skip = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
        skip1 = ( char * ) calloc( MAX_RECORD, sizeof( char * ) );
        fmt = ( char ** ) calloc( defs.ngenes, sizeof( char * ) );
        fmt1 = ( char ** ) calloc( defs.egenes, sizeof( char * ) );
        temptweak = ( int * ) calloc( defs.ngenes, sizeof( int * ) );
        temptweak1 = ( int * ) calloc( defs.egenes, sizeof( int * ) );

        /* create format strings according to the number of genes */

        for( i = 0; i < defs.ngenes; i++ ) {
            fmt[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            fmt[i] = strcpy( fmt[i], skip );
            fmt[i] = strcat( fmt[i], read_fmt );
            skip = strcat( skip, skip_fmt );
        }

        /* create format strings according to the number of external inputs */

        for( i = 0; i < defs.egenes; i++ ) {
            fmt1[i] = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            fmt1[i] = strcpy( fmt1[i], skip1 );
            fmt1[i] = strcat( fmt1[i], read_fmt );
            skip1 = strcat( skip1, skip_fmt );
        }



        fp = FindSection( fp, "tweak" );        // find tweak section
        if( !fp )
            error( "ReadTweak: could not locate tweak\n" );
        while( strncmp( ( base = fgets( base, MAX_RECORD, fp ) ), "$$", 2 ) ) {
            record = base;
            c = ( int ) *record;
            while( c != '\0' ) {

                if( isdigit( c ) ) {    // line contains data
                    record = base;

                    // usually read ngenes parameters, but for diff. schedule A or C only read
                    // one d parameter
                    if( ( linecount == 5 ) && ( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) ) {
                        if( 1 != sscanf( record, fmt[0], &temptweak[0] ) ) {
                            error( "ReadTweak: error reading tweaks" );
                        }
                    } else if( linecount == 2 ) {
                        for( i = 0; i < defs.egenes; i++ ) {
                            if( 1 != sscanf( record, fmt1[i], &temptweak1[i] ) )
                                error( "ReadTweak: error reading tweak variables" );
                        }
                    } else {
                        for( i = 0; i < defs.ngenes; i++ ) {
                            if( 1 != sscanf( record, fmt[i], &temptweak[i] ) )
                                error( "ReadTweak: error reading tweak variables" );
                        }
                    }
                    switch ( linecount ) {      // copy read parameters into the right array
                    case 0:
                        for( i = 0; i < defs.ngenes; i++ )      // R tweaks
                            l_tweak.Rtweak[i] = temptweak[i];
                        linecount++;
                        break;
                    case 1:    // T tweaks: keep track of read lines with Tcount
                        for( i = 0; i < defs.ngenes; i++ )
                            l_tweak.Ttweak[i + Tcount * defs.ngenes] = temptweak[i];
                        Tcount++;
                        if( Tcount == defs.ngenes )
                            linecount++;
                        break;
                    case 2:    // E tweaks: keep track of read lines with Ecount
                        for( i = 0; i < defs.egenes; i++ )
                            l_tweak.Etweak[i + Ecount * defs.egenes] = temptweak1[i];
                        Ecount++;
                        if( Ecount == defs.ngenes )
                            linecount++;
                        break;
                    case 3:    // m tweaks
                        for( i = 0; i < defs.ngenes; i++ )
                            l_tweak.mtweak[i] = temptweak[i];
                        linecount++;
                        break;
                    case 4:
                        for( i = 0; i < defs.ngenes; i++ )      // h tweaks
                            l_tweak.htweak[i] = temptweak[i];
                        linecount++;
                        break;
                    case 5:    // d tweaks: consider diff. schedule
                        if( ( defs.diff_schedule == 'A' ) || ( defs.diff_schedule == 'C' ) ) {
                            l_tweak.dtweak[0] = temptweak[0];
                        } else {
                            for( i = 0; i < defs.ngenes; i++ )
                                l_tweak.dtweak[i] = temptweak[i];
                        }
                        linecount++;
                        break;
                    case 6:    // lambda tweaks
                        for( i = 0; i < defs.ngenes; i++ )
                            l_tweak.lambdatweak[i] = temptweak[i];
                        linecount++;
                        break;
                    case 7:    // lambda tweaks
                        for( i = 0; i < defs.ngenes; i++ )
                            l_tweak.tautweak[i] = temptweak[i];
                        linecount++;
                        break;
                    default:
                        error( "ReadTweak: too many data lines in tweak section" );
                    }
                    break;      // don't do rest of loop anymore!
                } else if( isspace( c ) ) {     // ignore leading white space
                    c = ( int ) *( ++record );
                } else {        // anything but space or digit means comment
                    break;
                }
            }
        }
        free( temptweak );
        free( temptweak1 );
        free( base );
        free( skip );
        free( skip1 );

        for( i = 0; i < defs.ngenes; i++ )
            free( fmt[i] );
        free( fmt );

        for( i = 0; i < defs.egenes; i++ )
            free( fmt1[i] );
        free( fmt1 );
    }
    return l_tweak;
}


/*
 * This routine reads the distribution parameters  from the input file
 * and stores them in DistP.xx from distributions.h                   
 * LG 03-02: need q for gen visiting distribution input file          
 * LG 05-02: set factors only dependent on q for general visiting     
 * distribution by calling qgt2_init or qlt2_init from distributions.c
 */
/*
DistParms
InitDistribution( FILE * fp ) {
    DistParms DistP;            // variables for distributions 

    // read in the data from file
    // read in header info then output it
    // At the end of file exit from the loop

    // find distribution section
    fp = FindSection( fp, "distribution_parameters" );

    // to be backward compatible set default to exp
    if( !fp ) {
        if( debug ) {
            printf( "ReadTune: no distribution parameters, using exponential.\n" );
        }
        DistP.distribution = 1;
        DistP.q = 1.0;
    } else {
        // input user selected distribution parameters
        // advance past title line no blanks allowed!
        fscanf( fp, "%*s\n" );
        // read distribution stuff
        if( 2 != ( fscanf( fp, "%d %lf\n", &( DistP.distribution ), &( DistP.q ) ) ) )
            error( "ReadTune: error reading distribution stuff" );
        if( DistP.distribution > 11 || DistP.distribution < 1 ) {
            error( "fly_**: distribution must be int [1, 11] \n" );
        } else if( DistP.distribution == 4 || DistP.distribution == 3 ) {
            error( "fly_**: PLEASE use 5 for Lorentz or 10 for normal distribution \n" );
        } else if( DistP.distribution == 6 || DistP.distribution == 9 ) {
            error( "fly_**: 6=poisson or 9=pareto distribution returns positive values--do not use for fly \n" );
        } else if( DistP.distribution == 7 ) {
            // general distribution
            if( DistP.q >= 3.0 || DistP.q <= 1.0 ) {
                error( "tsp_sa: q must be between 1 and 3 \n" );
            } else if( DistP.q == 2.0 ) {
                DistP.distribution = 5;
                // fly needs lorentz, tsp use abs lorentz(4) 
                printf( "fly_**: q=2 is lorentz--setting distribution to 5\n" );
            } else if( DistP.q > 2.0 ) {
                qgt2_init( &DistP );
            } else {
                qlt2_init( &DistP );
            }
        }

        // LG 05-02
        // calculate q dependent factors that
        // do not change for entire run.
        // these live in distributions.h
        // LG 05-02
    }
    return DistP;
}
*/
//OUTPUT
//------

/** PrintBlastoderm: writes the output of the model to a stream specified 
 *                    by the fp file pointer. The Table is a solution of   
 *                    the model as returned by Blastoderm, the id speci-   
 *                    fies the title of the output and ndigits specifies   
 *                    the floating point precision to be printed.          
 *                    PrintBlastoderm adjusts its format automatically to  
 *                    the appropriate number of genes.                     
 */
void
PrintBlastoderm( FILE * fp, NArrPtr table, char *id, int ndigits, Zygote * zyg ) {
    int i, j, k;                /* local loop counters */
    int lineage;                /* lineage number for nucleus */
    int columns = zyg->defs.ngenes;
    /* print title (id) */
    fprintf( fp, "$%s\n", id );
    /* print table with correct lineage numbers (obtained from maternal.c) */
    for( i = 0; i < table.size; i++ ) {
        for( j = 0; j < ( table.array[i].state.size / columns ); j++ ) {
            lineage = GetStartLin( table.array[i].time, zyg->defs, zyg->lin_start, &( zyg->times ) ) + j;
            fprintf( fp, "%5d %9.3f", lineage, table.array[i].time );
            for( k = 0; k < columns; k++ ) 
                fprintf( fp, " %*.*f", ndigits + 5, ndigits, table.array[i].state.array[k + ( j * columns )] );                            
            fprintf( fp, "\n" );
        }
        fprintf( fp, "\n\n" );
    }
    fprintf( fp, "$$\n" );
    fflush( fp );
}


/** WriteVersion: prints the version and the complete command line used 
 *                 to run fly_sa into the $version section of the data  
 *                 file                                                    
 */
void
WriteVersion( char *filename, char *version, char *argvsave ) {
    char *temp;                 /* temporary file name */
    char *record;               /* record to be read and written */
    char *record_ptr;           /* pointer used to remember record for 'free' */
    char *convline;             /* temporarily saves conv version line */
    char *shell_cmd;            /* used by 'system' below */

    FILE *outfile;              /* name of output file */
    FILE *tmpfile;              /* name of temporary file */

    int versionflag = 0;        /* version section already present or not? */



    temp = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    record = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    convline = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    record_ptr = record;        /* this is to remember record for 'free' */

    outfile = fopen( filename, "r" );   /* open outfile for reading */
    if( !outfile )              /* sorry for the confusion! */
        error( "WriteVersion: error opening %s", filename );

    if( FindSection( outfile, "version" ) )     /* version section already there? */
        versionflag = 1;
    rewind( outfile );

    temp = strcpy( temp, "versionXXXXXX" );     /* required by mkstemp() */
    if( mkstemp( temp ) == -1 ) /* get unique name for temp file */
        error( "WriteVersion: error creating temporary file" );

    tmpfile = fopen( temp, "w" );       /* ... and open it for writing */
    if( !tmpfile )
        error( "WriteVersion: error opening temporary file %s", temp );

    if( !versionflag ) {        /* CASE 1: no version section yet */
        fprintf( tmpfile, "$version\n" );       /* -> write new one at beginning */
        fprintf( tmpfile, "%s\n", version );    /* of the data file */
        fprintf( tmpfile, "%s", argvsave );
        fprintf( tmpfile, "$$\n\n" );
        while( ( record = fgets( record, MAX_RECORD, outfile ) ) )
            fputs( record, tmpfile );
        /* CASE 2: version section already present */
    } else {                    /* -> append new lines to existing section */
        while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$version", 8 ) )        /* first write everything before */
            fputs( record, tmpfile );   /* version section */
        fputs( record, tmpfile );

        while( strncmp( record = fgets( record, MAX_RECORD, outfile ), "$$", 2 ) )      /* save existing converter line */
            if( !strncmp( record, "converted", 9 ) )    /* skip all the rest */
                convline = strcpy( convline, record );

        fprintf( tmpfile, "%s\n", version );    /* new lines are written here */
        fprintf( tmpfile, "%s", argvsave );

        if( strlen( convline ) > 0 )    /* conversion line is appended at end */
            fprintf( tmpfile, "%s", convline );

        fprintf( tmpfile, "$$\n" );

        while( ( record = fgets( record, MAX_RECORD, outfile ) ) )
            fputs( record, tmpfile );   /* ... and then write all the rest */
    }

    fclose( outfile );
    fclose( tmpfile );

    /* rename tmpfile into new file */

    sprintf( shell_cmd, "cp -f %s %s", temp, filename );

    if( -1 == system( shell_cmd ) )
        error( "WriteVersion: error renaming temp file %s", temp );

    if( remove( temp ) )
        warning( "WriteVersion: temp file %s could not be deleted", temp );

    /* clean up */
    free( temp );
    free( record_ptr );
    free( convline );
    free( shell_cmd );
}

/** PrintEquil: writes an 'equilibrate_variance' section with 'title' 
 *               to the stream specified by fp                             
 */
void
PrintEquil( FILE * fp, double *equil_var, char *title ) {
    fprintf( fp, "$%s\n", title );
    fprintf( fp, "variance:\n" );
    fprintf( fp, "%g\n", equil_var[0] );
    fprintf( fp, "equilibrate_final_energy:\n" );
    fprintf( fp, "%g\n", equil_var[1] );
    fprintf( fp, "$$\n" );
}

/** PrintTimes: writes two (parallel: three) times sections */
void
PrintTimes( FILE * fp, double *times ) {
    fprintf( fp, "wallclock: %.3f\n", times[0] );
    fprintf( fp, "user:      %.3f\n", times[1] );
}

