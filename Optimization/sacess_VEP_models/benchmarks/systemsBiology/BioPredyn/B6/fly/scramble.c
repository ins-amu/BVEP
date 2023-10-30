/**
 * @file scramble.c                                              
 * @author JR, modified by Yoginho                            
 *                                                                  
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                                 
 * @brief Reads the 'input' or 'eqparms' section of a data file and 
 * returns random parameter values that are within the limits specified by 
 * the $limits section of the data file. 
 *
 * Note that only those parameters will get scrambled   
 * which are indicated as to-be-tweaked in the $tweak section of   
 * the data file.                                                  
 *                                                                 
 * \par A short comment on penalty limits (JJ, Aug 7, 2001):            
 *                                                                 
 * In order to produce random values of T, m and h that are        
 * within reasonable limits when no explicit limits are used for   
 * them, we do the following approximation:                        
 *                                                                 
 * First, we determine the limits of u in g(u) within which        
 * there is no penalty assigned. This is true for all g(u) be-     
 * tween Lambda and (1-Lambda). Therefore, we calculate the in-    
 * verse function of g(u) for those to values to get upper and      
 * lower limits for u (based on ideas by Eric Mjolsness).           
 * Since we need to sum up all T * vmax, m * mmax and h to get     
 * u, we'll compensate by dividing the limits by the sqrt of the   
 * number of genes in the problem. This way we think, we'll get    
 * reasonable limits for single parameters (idea by JR).           
 * All the above happens in the Penalty2Limits function in sco-    
 * re.c. When comparing parameters to limits, don't forget to      
 * multiply Ts with vmax and ms with mmax. This happens in main    
 * of scramble below. hs are compared as they are.                 
 *                                                                 
 * See JJs lab notes for further detail on g(u)-inverse and such.  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

//#include <moves.h>              /* for tweak */
//#include <error.h>
#include <maternal.h>           /* for defs */
#include <score.h>              /* for limits */
#include <zygotic.h>            /* for EqParms */
#include <../util/random.h>


const char *OPTS = ":f:hvw:x:"; /* command line option string */



/*** Help, usage and version messages **************************************/

static const char usage[] =
    "Usage: scramble [-f <float_prec>] [-h] [-v]\n" "                [-w <out_file>] [-x <sect_title>]\n" "                <datafile>\n";

static const char help[] =
    "Usage: scramble [options] <datafile>\n\n"
    "Argument:\n"
    "  <datafile>          data file to be scambled\n\n"
    "Options:\n"
    "  -f <float_prec>     float precision of output is <float_prec>\n"
    "  -h                  prints this help message\n"
    "  -v                  print version and compilation date\n"
    "  -w <out_file>       write output to <out_file> instead of <datafile>\n"
    "  -x <sect_title>     scrambles eq params from section <sect_title>\n\n" "Please report bugs to <yoginho@usa.net>. Thank you!\n";

static const char verstring[] =
    "%s version %s\n" "compiled by:      %s\n" "         on:      %s\n" "      using:      %s\n" "      flags:      %s\n" "       date:      %s at %s\n";


// GLOBAL CONSTANTS

/** The following defines the maximum float precision that is supported by  
 * the code.
 */
const int MAX_PRECISION = 16;
/* the following constant as a score tells the annealer to reject a move,  */
/* no matter what. It had better not be a number that could actually be a  */
/* score.                                                                  */
const double FORBIDDEN_MOVE = DBL_MAX;  /* the biggest possible score, ever */

const int OUT_OF_BOUND = -1;

/** scramble.c main function */
int
main( int argc, char **argv ) {
    int c;                      /* used to parse command line options */
    FILE *fp;                   /* parameter file to be scrambled */

    char *section;              /* section title of section to be scrambled */
    char *outname = NULL;       /* output file name for -o */
    int ndigits = 8;            /* precision of output, default = 8 */

    int i, j;               /* local loop counters */
    long pid;                   /* process ID of scramble: used to seed dSFMT */

    int penaltyflag = 0;

    double out;                 /* temporary storage for random numbers */
    double T_pen;               /* temp var for T * vmax in case penalty is used */
    double m_pen;               /* temp var for m * mmax in case penalty is used */

    char *shell_cmd;            /* used by 'system' below */

    /* external declarations for command line option parsing (unistd.h) */

    extern char *optarg;        /* command line option argument */
    extern int optind;          /* pointer to current element of argv */
    extern int optopt;          /* contain option character upon error */

    Input inp;
    int nrows, ncols, egenes;


    /* initialize section to 'input' as default */

    section = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    section = strcpy( section, "input" );       /* input section is default */

    /* following part parses command line for options and their arguments      */

    optarg = NULL;
    while( ( c = getopt( argc, argv, OPTS ) ) != -1 )
        switch ( c ) {
        case 'f':
            ndigits = atoi( optarg );   /* -f determines float precision */
            if( ndigits > MAX_PRECISION )
                error( "scramble: max. float precision is %d!", MAX_PRECISION );
            break;
        case 'h':              /* -h help option */
            PrintMsg( help, 0 );
            break;
        case 'v':              /* -v prints version number and exit */
            //fprintf(stderr, verstring, *argv, VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
            exit( 0 );
        case 'w':              /* -o sets output file */
            outname = calloc( MAX_RECORD + 1, sizeof( char ) );
            outname = strcpy( outname, optarg );
            break;
        case 'x':
            if( ( strcmp( optarg, "input" ) ) && ( strcmp( optarg, "eqparms" ) ) )
                error( "scramble: invalid section title (%s), must be input or eqparms", optarg );
            section = strcpy( section, optarg );
            break;
        case ':':
            error( "scramble: need an argument for option -%c", optopt );
            break;
        case '?':
        default:
            error( "scramble: unrecognized option -%c", optopt );
        }

    /* error check */

    if( argc - ( optind - 1 ) != 2 )
        PrintMsg( usage, 1 );

    /* now, let's get started */
    fp = fopen( argv[optind], "r" );
    if( !fp )
        file_error( "scramble" );
    inp.zyg = InitZygote( fp, pd, pj, &inp, "input" );
    inp.sco.searchspace = InitLimits( fp, &inp );
    //inp.sco = InitScoring(fp, method, &inp);
    inp.twe = InitTweak( fp, NULL, inp.zyg.defs );

    //Here we read the parameters given by the optimization algorithm
    fclose( fp );

    if( inp.sco.searchspace->pen_vec != 0 ) {   /* using penalty? */
        //printf( "Converting penalty to limits...\n" );
        Penalty2Limits( inp.sco.searchspace, inp.zyg.defs );    /* convert to explicit limits */
        penaltyflag = 1;        /* see also comment above */
    }
    pid = getpid(  );           /* get the process ID for seeding erand */
    
    InitRand(pid);

    /* the following part of the code assigns each parameter-to be tweaked a   */
    /* random value within the limits specified by the limits struct using the */
    /* formula:                                                                */
    /*             new_value = lowerlim + rand_real * (upperlim - lowerlim)    */
    /*                                                                         */
    /* where rand_real is a random number between 0 and 1 (dSFMT())            */
    /* (for the penalty case, see comment above)                               */
    ncols  = inp.zyg.defs.ngenes;
    nrows  = inp.zyg.defs.ngenes;
    egenes = inp.zyg.defs.egenes;


    for( i = 0; i < ncols; ++i ) {
        if( inp.twe.Rtweak[i] == 1 ) {  /* scramble Rs here */
            out = RandomReal();
            
            inp.zyg.parm.R[i] = inp.sco.searchspace->Rlim[i]->lower + out * ( inp.sco.searchspace->Rlim[i]->upper - inp.sco.searchspace->Rlim[i]->lower );
        }
        

        if( inp.twe.mtweak[i] == 1 ) {  /* scramble ms here */
            out = RandomReal();
            if( !penaltyflag ) {
                inp.zyg.parm.m[i] = inp.sco.searchspace->mlim[i]->lower + out * ( inp.sco.searchspace->mlim[i]->upper - inp.sco.searchspace->mlim[i]->lower );
            } else {
                m_pen = inp.sco.searchspace->mlim[i]->lower + out * ( inp.sco.searchspace->mlim[i]->upper - inp.sco.searchspace->mlim[i]->lower );
                inp.zyg.parm.m[i] = m_pen / inp.sco.searchspace->pen_vec[1];    /* m_pen / mmax */
            }
        }

        if( inp.twe.htweak[i] == 1 ) {  /* scramble hs here */
            out = RandomReal();
            inp.zyg.parm.h[i] = inp.sco.searchspace->hlim[i]->lower + out * ( inp.sco.searchspace->hlim[i]->upper - inp.sco.searchspace->hlim[i]->lower );
        }

        if( inp.twe.lambdatweak[i] == 1 ) {     /* scramble lambdas here */
            out = RandomReal();
            inp.zyg.parm.lambda[i] = inp.sco.searchspace->lambdalim[i]->lower +
                out * ( inp.sco.searchspace->lambdalim[i]->upper - inp.sco.searchspace->lambdalim[i]->lower );
        }

        if( inp.twe.tautweak[i] == 1 ) {        /* scramble delays here */
            out = RandomReal();
            inp.zyg.parm.tau[i] = inp.sco.searchspace->taulim[i]->lower +
                out * ( inp.sco.searchspace->taulim[i]->upper - inp.sco.searchspace->taulim[i]->lower );
        }
    }
    
    for( i = 0; i < nrows; ++i ) {

        for( j = 0; j < ncols; j++ ) { /* scramble Ts here */
            if( inp.twe.Ttweak[( i * ncols ) + j] == 1 ) {
                out = RandomReal();
                if( !penaltyflag ) {
                    inp.zyg.parm.T[( i * ncols ) + j] =
                        inp.sco.searchspace->Tlim[( i * ncols ) + j]->lower +
                        out * ( inp.sco.searchspace->Tlim[( i * ncols ) + j]->upper - inp.sco.searchspace->Tlim[( i * ncols ) + j]->lower );
                } else {
                    T_pen = inp.sco.searchspace->Tlim[( i * ncols ) + j]->lower +
                        out * ( inp.sco.searchspace->Tlim[( i * ncols ) + j]->upper - inp.sco.searchspace->Tlim[( i * ncols ) + j]->lower );
                    inp.zyg.parm.T[( i * ncols ) + j] = T_pen / inp.sco.searchspace->pen_vec[j + 2];
                }               /* above is T_pen / vmax[i] */
            }
        }

        for( j = 0; j < egenes; j++ ) { /* scramble Es here */
            if( inp.twe.Etweak[( i * egenes ) + j] == 1 ) {
                out = RandomReal();
                if( !penaltyflag ) {
                    inp.zyg.parm.E[( i * egenes ) + j] =
                        inp.sco.searchspace->Elim[( i * egenes ) + j]->lower +
                        out * ( inp.sco.searchspace->Elim[( i * egenes ) + j]->upper - inp.sco.searchspace->Elim[( i * egenes ) + j]->lower );
                } else {
                    T_pen = inp.sco.searchspace->Elim[( i * egenes ) + j]->lower +
                        out * ( inp.sco.searchspace->Elim[( i * egenes ) + j]->upper - inp.sco.searchspace->Elim[( i * egenes ) + j]->lower );

                    inp.zyg.parm.E[( i * egenes ) + j] = T_pen / inp.sco.searchspace->pen_vec[ncols + j + 2];

                }               /* above is T_pen / vmax[i] */
            }
        }
    }
    /* ds need to be srambled separately, since for diffusion schedules A & C  */
    /* there's just one single d                                               */

    if( ( inp.zyg.defs.diff_schedule == 'A' ) || ( inp.zyg.defs.diff_schedule == 'C' ) ) {
        if( inp.twe.dtweak[0] == 1 ) {
            out = RandomReal();
            inp.zyg.parm.d[0] = inp.sco.searchspace->dlim[0]->lower + out * ( inp.sco.searchspace->dlim[0]->upper - inp.sco.searchspace->dlim[0]->lower );
        }
    } else {
        for( i = 0; i < ncols; i++ ) {
            if( inp.twe.dtweak[i] == 1 ) {
                out = RandomReal();
                inp.zyg.parm.d[i] = inp.sco.searchspace->dlim[i]->lower + out * ( inp.sco.searchspace->dlim[i]->upper - inp.sco.searchspace->dlim[i]->lower );
            }
        }
    }

    if( outname != NULL ) {     /* -o used? */

        shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        sprintf( shell_cmd, "cp -f %s %s", argv[optind], outname );
        if( -1 == system( shell_cmd ) )
            error( "scramble: error creating the output file %s", outname );
        free( shell_cmd );
        WriteParameters( outname, &( inp.zyg.parm ), section, ndigits, inp.zyg.defs );
    } else {
        WriteParameters( argv[optind], &( inp.zyg.parm ), section, ndigits, inp.zyg.defs );
    }

    /* clean up */
    FreeZygote(  );

    free( section );

    free( inp.twe.Rtweak );
    free( inp.twe.Ttweak );
    free( inp.twe.Etweak );
    free( inp.twe.mtweak );
    free( inp.twe.htweak );
    free( inp.twe.dtweak );
    free( inp.twe.lambdatweak );
    free( inp.twe.tautweak );

    for( i = 0; i < ncols; i++ ) {
        free( inp.sco.searchspace->Rlim[i] );
        free( inp.sco.searchspace->mlim[i] );
        free( inp.sco.searchspace->hlim[i] );
        free( inp.sco.searchspace->lambdalim[i] );
        free( inp.sco.searchspace->taulim[i] );
    }
    
    for (i=0; i < nrows; i++) {
        for( j = 0; j < ncols; j++ )
            free( inp.sco.searchspace->Tlim[( i * ncols ) + j] );
        for( j = 0; j < egenes; j++ )
            free( inp.sco.searchspace->Elim[( i * egenes ) + j] );
    }
    
    if( ( inp.zyg.defs.diff_schedule == 'A' ) || ( inp.zyg.defs.diff_schedule == 'C' ) )
        free( inp.sco.searchspace->dlim[0] );
    else
        for( i = 0; i < ncols; i++ )
            free( inp.sco.searchspace->dlim[i] );

    if( inp.sco.searchspace->pen_vec )
        free( inp.sco.searchspace->pen_vec );
    free( inp.sco.searchspace->Rlim );
    free( inp.sco.searchspace->Tlim );
    free( inp.sco.searchspace->Elim );
    free( inp.sco.searchspace->mlim );
    free( inp.sco.searchspace->hlim );
    free( inp.sco.searchspace->dlim );
    free( inp.sco.searchspace->lambdalim );
    free( inp.sco.searchspace->taulim );

    free( inp.sco.searchspace );

    free( inp.zyg.defs.egene_ids );
    free( inp.zyg.defs.gene_ids );
    free( inp.zyg.nnucs );
    free( inp.zyg.full_nnucs );
    free( inp.zyg.full_lin_start );
    free( inp.zyg.lin_start );
    free( inp.zyg.parm.E );
    free( inp.zyg.parm.R );
    free( inp.zyg.parm.T );
    free( inp.zyg.parm.d );
    free( inp.zyg.parm.h );
    free( inp.zyg.parm.lambda );
    free( inp.zyg.parm.m );
    free( inp.zyg.parm.tau );

    for( i = 0; i < inp.zyg.nalleles; i++ ) {
        free( inp.zyg.bias.bt[i].genotype );
        free( inp.zyg.bias.bt[i].ptr.times.array );
        free( inp.zyg.bias.biastype[i].genotype );

        for( j = 0; j < inp.zyg.bias.biastype[i].ptr.bias.size; j++ ) {
            free( inp.zyg.bias.biastype[i].ptr.bias.array[j].state.array );
        }
        free( inp.zyg.bias.biastype[i].ptr.times.array );
        free( inp.zyg.bcdtype[i].genotype );
        for( j = 0; j < inp.zyg.bcdtype[i].ptr.bicoid.size; j++ ) {
            free( inp.zyg.bcdtype[i].ptr.bicoid.array[j].gradient.array );
        }
        free( inp.zyg.bcdtype[i].ptr.bicoid.array );
    }
    free( inp.zyg.bias.bt );
    free( inp.zyg.bias.biastype );
    free( inp.zyg.bcdtype );


    if( outname )
        free( outname );

    return 0;
}
