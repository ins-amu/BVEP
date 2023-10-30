/**
 * @file printscore.c
 * @author JR, modified by Yoginho,
 *   -g option by Yousong Wang in Feb 2002,                           
 *   -a option by Marcel Wolf in Apr 2002
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                                 
 * @brief Prints score and penalty, as well as the root mean   
 * square (RMS) of a gene circuit to STDOUT.                        
 *                                                                 
 * See the dataformatX.X file for further details on the current   
 * data file format (X.X stands for the current code version).     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>             /* for getopt */

//#include <error.h>
#include <integrate.h>
#include <maternal.h>
#include <score.h>
#include <solvers.h>
#include <zygotic.h>


/* *Constants *************************************************************/

const char *OPTS = ":a:Df:g:Ghi:m:opqr:s:vx:";  /* command line option string */


/*** Help, usage and version messages **************************************/

static const char usage[] =
    "Usage: printscore [-a <accuracy>] [-D] [-f <float_prec>] [-g <g(u)>] [-G]\n"
    "                  [-h] [-i <stepsize>] [-m <score_method>] [-o] [-p] \n"
    "                  [-s <solver>] [-v] [-x <sect_title>]\n" 
    "                  <datafile>\n";

static const char help[] =
    "Usage: printscore [options] <datafile>\n\n"
    "Arguments:\n"
    "  <datafile>          data file for which we evaluate score and RMS\n\n"
    "Options:\n"
    "  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
    "  -D                  debugging mode, prints all kinds of debugging info\n"
    "  -f <float_prec>     float precision of output is <float_prec>\n"
    "  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
    "  -G                  gut mode: prints squared diffs for all datapoints\n"
    "  -h                  prints this help message\n"
    "  -m <score_method>   w = wls, o=ols score calculation method\n"
    "  -i <stepsize>       sets ODE solver stepsize (in minutes)\n"
    "  -o                  use oldstyle cell division times (3 div only)\n"
    "  -p                  prints penalty in addition to score and RMS\n"
    "  -s <solver>         choose ODE solver\n"
    "  -v                  print version and compilation date\n"
    "  -x <sect_title>     uses equation paramters from section <sect_title>\n\n" "Please report bugs to <yoginho@usa.net>. Thank you!\n";

static const char verstring[] =
    "%s version %s\n" 
    "compiled by:      %s\n" 
    "         on:      %s\n" 
    "      using:      %s\n" 
    "      flags:      %s\n" 
    "       date:      %s at %s\n";


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

/** printscore main() function */
int
main( int argc, char **argv ) {
    int c;                      /* used to parse command line options */
    char *infile;               /* pointer to input file name */
    FILE *fp;                   /* pointer to input data file */

    char *slogfile;             /* name of solver log file */
    FILE *slog;                 /* solver log file pointer */
    int i, j;

    /* the follwoing few variables are read as arguments to command line opts  */

    int ndigits = 12;           /* precision of output, default = 12 */
    int gutndigits = 6;         /* precision for gut output, def = 6 */
    int penaltyflag = 0;        /* flag for printing penalty */
    int rmsflag = 1;            /* flag for printing root mean square */
    int gutflag = 0;            /* flag for root square diff guts */

    double stepsize = 1.;       /* stepsize for solver */
    double accuracy = 0.001;    /* accuracy for solver */
    int method = 0;             /* 0 for wls, 1 for ols */


    char *section_title;        /* parameter section name */

    /* two format strings */

    char *format;               /* output format string */
    char *precision;            /* precision string for output format string */

    /* output values */

    double chisq = 0.;          /* sum of squared differences for output */

    double ols_chisq = 0;       /* chisq for RMS is always calculated with ols method */

    double penalty = 0.;        /* variable for penalty */

    double rms = 0.;            /* root mean square (RMS) */
    double ms = 0.;             /* mean square (= rms * rms) */
    int ndp = 0;                /* number of datapoints */
    struct rusage begin, end;   /*        structs for measuring time */

    Input inp;
    ScoreOutput out;

    out.score = 1e38;           // start with a very large number
    out.penalty = 0;
    out.size_resid_arr = 0;
    out.jacobian = NULL;
    out.residuals = NULL;

    /* the following lines define a pointers to:                               */
    /*            - pd:    dvdt function, currently only DvdtOrig in zygotic.c */
    /*            - pj:    Jacobian function, in zygotic.c                     */
    /*                                                                         */
    /* NOTE: ps (solver) is declared as global in integrate.h                  */

    void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * );
    void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * );
    /* external declarations for command line option parsing (unistd.h) */

    extern char *optarg;        /* command line option argument */
    extern int optind;          /* pointer to current element of argv */
    extern int optopt;          /* contain option character upon error */

    getrusage( RUSAGE_SELF, &begin );   /*          get start time */

    /* following part sets default values for deriv, Jacobian and solver funcs */

    pd = DvdtOrig;
    //pd = Dvdt_sqrt;
    dd = DvdtDelay;             /* delayed derivative fnuction */
    pj = JacobnOrig;
    ps = Rkck;

    section_title = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    section_title = strcpy( section_title, "eqparms" ); /* default is eqparms */

    /* following part parses command line for options and their arguments      */
    optarg = NULL;
    while( ( c = getopt( argc, argv, OPTS ) ) != -1 )
        switch ( c ) {
        case 'a':
            accuracy = atof( optarg );
            if( accuracy <= 0 )
                error( "fly_sa: accuracy (%g) is too small", accuracy );
            break;
        case 'D':              /* -D runs in debugging mode */
            debug = 1;
            break;
        case 'f':
            ndigits = atoi( optarg );   /* -f determines float precision */
            gutndigits = atoi( optarg );
            if( ndigits < 0 )
                error( "printscore: what exactly would a negative precision be???" );
            if( ndigits > MAX_PRECISION )
                error( "printscore: max. float precision is %d!", MAX_PRECISION );
            break;
        case 'g':              /* -g choose g(u) function */
            pd = DvdtOrig;
            if( !( strcmp( optarg, "s" ) ) )
                gofu = Sqrt;
            else if( !( strcmp( optarg, "t" ) ) )
                gofu = Tanh;
            else if( !( strcmp( optarg, "e" ) ) )
                gofu = Exp;
            else if( !( strcmp( optarg, "h" ) ) )
                gofu = Hvs;
            else if( !( strcmp( optarg, "k" ) ) ) {
                gofu = Kolja;
            } else
                error( "printscore: %s is an invalid g(u), should be e, h, s or t", optarg );
            break;
        case 'G':              /* -G guts mode */
            gutflag = 1;
            break;
        case 'h':              /* -h help option */
            PrintMsg( help, 0 );
            break;
        case 'i':              /* -i sets the stepsize */
            stepsize = atof( optarg );
            if( stepsize < 0 )
                error( "printscore: going backwards? (hint: check your -i)" );
            if( stepsize == 0 )
                error( "printscore: going nowhere? (hint: check your -i)" );
            if( stepsize > MAX_STEPSIZE )
                error( "printscore: stepsize %g too large (max. is %g)", stepsize, MAX_STEPSIZE );
            break;
        case 'm':              /* -m sets the score method: w for wls, o for ols */
            if( !( strcmp( optarg, "w" ) ) )
                method = 0;
            else if( !( strcmp( optarg, "o" ) ) )
                method = 1;
            break;
        case 'o':              /* -o sets old division style (ndivs = 3 only! ) */
            olddivstyle = 1;
            break;
        case 'p':
            penaltyflag = 1;
            break;
        case 'q':
            warning( "printscore: RMS now printed by default, no need for -q" );
            break;
        case 'r':
            error( "printscore: -r is not supported anymore, use -g instead" );
            break;
        case 's':              /* -s sets solver to be used */
            if( !( strcmp( optarg, "a" ) ) )
                ps = Adams;
            else if( !( strcmp( optarg, "bd" ) ) )
                ps = BaDe;
            else if( !( strcmp( optarg, "bs" ) ) )
                ps = BuSt;
            else if( !( strcmp( optarg, "e" ) ) )
                ps = Euler;
            else if( !( strcmp( optarg, "h" ) ) )
                ps = Heun;
            else if( !( strcmp( optarg, "mi" ) ) || !( strcmp( optarg, "m" ) ) )
                ps = Milne;
            else if( !( strcmp( optarg, "me" ) ) )
                ps = Meuler;
            else if( !( strcmp( optarg, "r4" ) ) || !( strcmp( optarg, "r" ) ) )
                ps = Rk4;
            else if( !( strcmp( optarg, "r2" ) ) )
                ps = Rk2;
            else if( !( strcmp( optarg, "rck" ) ) )
                ps = Rkck;
            else if( !( strcmp( optarg, "rf" ) ) )
                ps = Rkf;
            else if( !( strcmp( optarg, "sd" ) ) )
                ps = SoDe;
            else if( !( strcmp( optarg, "kr" ) ) )
                ps = Krylov;
            /*
               else if (!(strcmp(optarg, "bnd")))
               ps = Band;
             */
            else
                error( "printscore: bad solver (%s), use: a,bd,bs,e,h,kr,mi,me,r{2,4,ck,f}", optarg );
            break;
        case 'v':              /* -v prints version number */
            //fprintf(stderr, verstring, *argv, VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
            exit( 0 );
        case 'x':
            if( ( strcmp( optarg, "input" ) ) && ( strcmp( optarg, "eqparms" ) ) && ( strcmp( optarg, "parameters" ) ) )
                error( "printscore: invalid section title (%s)", optarg );
            section_title = strcpy( section_title, optarg );
            break;
        case ':':
            error( "printscore: need an argument for option -%c", optopt );
            break;
        case '?':
        default:
            error( "printscore: unrecognized option -%c", optopt );
        }

    /* error check */

    if( ( argc - ( optind - 1 ) ) != 2 )
        PrintMsg( usage, 1 );

    /* dynamic allocation of output format strings */

    precision = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    format = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    /* initialize guts */

    if( gutflag )
        SetGuts( gutflag, gutndigits );

    /* let's get started and open data file here */
    infile = argv[optind];
    fp = fopen( infile, "r" );
    if( !infile )
        file_error( "printscore" );

    if( fp == NULL )
        printf( "FP NULLPOINTER file=%s \n", infile );

    /* if debugging: initialize solver log file */

    if( debug ) {

        slogfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        sprintf( slogfile, "%s.slog", argv[optind] );

        slog = fopen( slogfile, "w" );  /* delete existing slog file */
        fclose( slog );

        slog = fopen( slogfile, "a" );  /* now keep open for appending */
    } else { // suppressing compiler warnings
        slogfile = NULL;
        slog = NULL;
    }

    /* Initialization code here: InitZygote() initializes everything needed    *
     * for running the model, InitScoring() everything for scoring; InitStep-  *
     * size sets solver stepsize and accuracy in score.c to pass it on to      *
     * Blastoderm                                                              */
    //printf("InitZyg...\n");

    inp.zyg = InitZygote( fp, pd, pj, &inp, section_title );

    inp.sco = InitScoring( fp, method, &inp );
    //printf("...ok!\nInitHis...");
    inp.his = InitHistory( fp, &inp );  //It fills the polations vector
    //printf("...ok!\nInitExtinp...");
    inp.ext = InitExternalInputs( fp, &inp );
    //printf("...ok!\nInitStepsize...");
    inp.ste = InitStepsize( stepsize, accuracy, slog, infile );
    // read the list of parameters to be tweaked
    //printf("...ok!\nInitTweak...");
    inp.twe = InitTweak( fp, NULL, inp.zyg.defs );
    //printf("...ok!\nInitTransl...");
    inp.tra = Translate( &inp );
    // array of pointers to parameters and ranges
    //ParamList *ptab = inp.tra.array;
    // initialize distribution stuff
    //printf("...ok!\nInitDist...");

    inp.dis = InitDistribution( fp );

    //printf("ReadParameters() done\n");
    inp.lparm = CopyParm( inp.zyg.parm, &( inp.zyg.defs ) );
    //printf("Initlparm() done\n");
    fclose( fp );


    /* Scoring happens here (and also printing of guts if -G) */

    for( i = 0; i < 1; ++i ) {
        if( i > 0 ) {
            free( out.residuals );
            out.residuals = NULL;
        }
        Score( &inp, &out, 0 );
        chisq = out.score + out.penalty;
        inp.sco.method = 1;
        Score( &inp, &out, 0 ); //in order to get ols score for RMS calculation
        ols_chisq = out.score;  //ols score for calculating rms - without penalty
        inp.sco.method = method;
    }
    /* exit here if we've printed guts already */

    /*if ( gutflag ) {              
       FreeHistory(inp.zyg.nalleles, inp.his);
       FreeExternalInputs(inp.zyg.nalleles, inp.ext);
       free(precision);
       free(format);
       free(section_title);
       return 0;
       }
     */

    /* next few lines format & print the output to the appropriate precision   */

    sprintf( precision, "%d", ndigits );

    format = strcpy( format, " chisq = %." );
    format = strcat( format, precision );
    format = strcat( format, "f" );

    printf( format, chisq );
    if( penaltyflag ) {         /* in case of -p, print penalty too */
        penalty = GetPenalty( &inp, inp.sco.searchspace );
        if( penalty == -1 ) {
            printf( "    penalty = -1.000" );
        } else {
            format = strcpy( format, "     penalty = %." );
            format = strcat( format, precision );
            format = strcat( format, "f" );
            printf( format, penalty );
        }
    }
    if( rmsflag ) {             /* print rms */
        ndp = inp.zyg.ndp;
        ms = ( ols_chisq ) / ( double ) ndp;
        rms = sqrt( ms );
        format = strcpy( format, "     rms = %." );
        format = strcat( format, precision );
        format = strcat( format, "f" );
        printf( format, rms );
    }
    printf( "\n" );
    /* clean up before you go home... */
    getrusage( RUSAGE_SELF, &end );     /* get end time */
    printf( "# Printscore ran for %.13f seconds\n", tvsub( end, begin ) );

    for( i = 0; i < inp.zyg.nalleles; i++ ) {
        free( inp.sco.facts.tt[i].genotype );
        free( inp.sco.facts.tt[i].ptr.times.array );
        free( inp.sco.facts.facttype[i].genotype );
        free( inp.sco.facts.facttype[i].ptr.times.array );
        for( j = 0; j < inp.sco.facts.facttype[i].ptr.facts->size; j++ ) {
            free( inp.sco.facts.facttype[i].ptr.facts->record[j].array );
        }
        free( inp.sco.facts.facttype[i].ptr.facts->record );
        free( inp.sco.facts.facttype[i].ptr.facts );

        if( method == 0 ) {
            free( inp.sco.weights.weighttype[i].genotype );
            free( inp.sco.weights.weighttype[i].ptr.times.array );
            for( j = 0; j < inp.sco.weights.weighttype[i].ptr.facts->size; j++ ) {
                free( inp.sco.weights.weighttype[i].ptr.facts->record[j].array );
            }
            free( inp.sco.weights.weighttype[i].ptr.facts->record );
            free( inp.sco.weights.weighttype[i].ptr.facts );
        }

    }
    free( inp.sco.facts.tt );
    free( inp.sco.facts.facttype );
    free( inp.sco.weights.weighttype );
    if( ( inp.zyg.defs.diff_schedule == 'A' ) || ( inp.zyg.defs.diff_schedule == 'C' ) ) {
        free( inp.sco.searchspace->dlim[0] );
    } else {
        for( i = 0; i < inp.zyg.defs.ngenes; i++ )
            free( inp.sco.searchspace->dlim[i] );
    }

    for( i = 0; i < inp.zyg.defs.ngenes; i++ ) {
        free( inp.sco.searchspace->Rlim[i] );
        free( inp.sco.searchspace->taulim[i] );
        free( inp.sco.searchspace->lambdalim[i] );
    }
    free( inp.sco.searchspace->dlim );
    free( inp.sco.searchspace->Rlim );
    free( inp.sco.searchspace->lambdalim );
    free( inp.sco.searchspace->taulim );
    free( inp.sco.searchspace->pen_vec );
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

    free( inp.twe.Etweak );
    free( inp.twe.Rtweak );
    free( inp.twe.Ttweak );
    free( inp.twe.dtweak );
    free( inp.twe.htweak );
    free( inp.twe.lambdatweak );
    free( inp.twe.mtweak );
    free( inp.twe.tautweak );

    free( inp.tra.array );

    FreeMutant( inp.lparm );
    FreeHistory( inp.zyg.nalleles, inp.his );
    FreeExternalInputs( inp.zyg.nalleles, inp.ext );
    FreeZygote(  );

    free( precision );
    free( format );
    free( section_title );
    free( out.residuals );
    if( debug )
        free( slogfile );

    return 0;
}
