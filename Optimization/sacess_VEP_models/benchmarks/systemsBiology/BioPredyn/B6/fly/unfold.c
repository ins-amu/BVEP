/**
 * @file unfold.c                                                
 * @version 9.3.2 (obsolete numbering)                                                 
 * @author JR, modified by Yoginho, 
 * \c -g option by Yousong Wang in Feb 2002,                           
 * \c -a option by Marcel Wolf in Apr 2002,                            
 * \c -G (guts) option by Manu in June 2002                           
 *                                                                 
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                                 
 * @brief Runs the fly model once (unfolding of development). 
 *
 * It uses the parameters and the maternal contributions stored for 
 * the specified genotype in the data file specified by a filename.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>             /* for getopt */
#include <time.h>               /* for time calculation */
#include <sys/resource.h>       /* for time calculation */

//#include <error.h>
#include <integrate.h>
#include <maternal.h>
#include <solvers.h>
#include <zygotic.h>
#include <score.h>
#include <fly_io.h>
//#include <moves.h>


/*** Constants *************************************************************/

const char *OPTS = ":a:Df:g:Ghi:j:op:r:s:t:vx:z:";      /* cmd line opt string */

/*** Help, usage and version messages **************************************/

static const char usage[] =
    "Usage: unfold [-a <accuracy>] [-D] [-f <float_prec>] [-g <g(u)>] [-G]\n"
    "              [-h] [-i <stepsize>] [-j <timefile>] [-o] [-p <pstep>]\n"
    "              [-s <solver>] [-t <time>] [-v] [-x <sect_title>]\n" 
    "              [-z <gast_time>]\n" 
    "              <datafile> [<genotype>]\n";

static const char help[] =
    "Usage: unfold [options] <datafile> [<genotype>]\n\n"
    "Arguments:\n"
    "  <datafile>          data file which contains equation parameters\n"
    "  <genotype>          genotype string (W, R, S or T for each gene)\n\n"
    "Options:\n"
    "  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
    "  -D                  debugging mode, prints all kinds of debugging info\n"
    "  -f <float_prec>     float precision of output is <float_prec>\n"
    "  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
    "  -G                  prints guts instead of model output\n"
    "  -h                  prints this help message\n"
    "  -i <stepsize>       sets ODE solver stepsize (in minutes)\n"
    "  -j <timefile>       reads output fimes from <timefile>\n"
    "  -o                  use oldstyle cell division times (3 div only)\n"
    "  -p <pstep>          prints output every <pstep> minutes\n"
    "  -s <solver>         choose ODE solver\n"
    "  -t <time>           prints output for <time>\n"
    "  -v                  print version and compilation date\n"
    "  -x <sect_title>     uses equation paramters from section <sect_title>\n\n"
    "  -z <gast_time>      set custom gastrulation time (max. 10'000'000 (min));\n"
    "                      custom gastrulation MUST be later than the normal one\n" "Please report bugs to <yoginho@usa.net>. Thank you!\n";

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

/** unfold.c main function */
int
main( int argc, char **argv ) {
    char *infile;               /* pointer to input file name */
    FILE *fp;                   /* pointer to input data file */

    int c;                      /* used to parse command line options */
    int j, i, ii;               /* loop counter */
    int numguts;                /* number of guts elements calculated */

    char *dumpfile;             /* name of debugging model output */
    FILE *dumpptr;              /* pointer for dumping raw model output */

    char *slogfile;             /* name of solver log file */
    FILE *slog = NULL;          /* solver log file pointer */

    /* the follwoing few variables are read as arguments to command line opts */

    int ndigits = 6;            /* precision of output, default = 6 */
    int guts = 0;               /* flag for printing guts */

    double stepsize = 1.;       /* stepsize for solver */
    double accuracy = 0.001;    /* accuracy for solvers */
    double p_stepsize = 0.;     /* output stepsize */

    double time = -999999999.;  /* time for -t option */
    char *timefile = NULL;      /* file for -j option */

    char *section_title;        /* parameter section name */

    /* stuff used as input/output for blastoderm */

    NArrPtr answer;             /* model output is stored in this */
    NArrPtr outtab;             /* times to be output */
    NArrPtr goutput;            /* output of guts is stored in this */
    DArrPtr tt;                 /* array with times for which there's data */

    int genindex = 0;           /* genotype index */
    char *genotype = NULL;      /* genotype string */
    Slist *geno, *curr;         /* We will read in the genotypes to
                                   set facts in score.c */
    char **gutsdefs = NULL;     /* array of strings to hold guts defs */
    char *gutstitle = NULL;
    int length;                 /* length of genotype string */
    struct rusage begin, end;   /*        structs for measuring time */

    DataTable **interp_dat;     /* for providing
                                   history to blastoderm */
    double *theta_discons;
    int theta_discons_size;
    double *temp_divtable;
    double *temp_durations;
    InterpObject *polation;
    InterpObject *extinp_polation;

    Input inp;

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

    /* following part sets default values for deriv, Jacobian and solver funcs */
    /* and section title                                                       */

    getrusage( RUSAGE_SELF, &begin );   /*          get start time */

    pd = DvdtOrig;
    //pd = Dvdt_sqrt;
    dd = DvdtDelay;             /* delayed derivative fnuction */
    pj = JacobnOrig;
    ps = Rk4;

    section_title = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    section_title = strcpy( section_title, "eqparms" ); /* default is eqparms */

    /* following part parses command line for options and their arguments      */
    /* modified from original getopt manpage                                   */

    optarg = NULL;
    while( ( c = getopt( argc, argv, OPTS ) ) != -1 )
        switch ( c ) {
        case 'a':
            accuracy = atof( optarg );
            if( accuracy <= 0 )
                error( "unfold: accuracy (%g) is too small", accuracy );
            break;
        case 'D':              /* -D runs in debugging mode */
            debug = 1;
            break;
        case 'f':
            ndigits = atoi( optarg );   /* -f determines float precision */
            if( ndigits < 0 )
                error( "unfold: what exactly would a negative precision be???" );
            if( ndigits > MAX_PRECISION )
                error( "unfold: max. float precision is %d!", MAX_PRECISION );
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
                error( "unfold: %s is an invalid g(u), should be e, h, s or t", optarg );
            break;
        case 'G':              /* -G guts options */
            guts = 1;
            break;
        case 'h':              /* -h help option */
            PrintMsg( help, 0 );
            break;
        case 'i':              /* -i sets stepsize for the solver */
            stepsize = atof( optarg );
            if( stepsize < 0 )
                error( "unfold: going backwards? (hint: check your -i)" );
            if( stepsize == 0 )
                error( "unfold: going nowhere? (hint: check your -i)" );
            if( stepsize > MAX_STEPSIZE )
                error( "unfold: stepsize %g too large (max. is %g)", stepsize, MAX_STEPSIZE );
            break;
        case 'j':
            if( p_stepsize > 0. )
                error( "unfold: can't use -p and -j together!" );
            if( time != -999999999. )
                error( "unfold: can't use -t and -j together!" );
            timefile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            timefile = strcpy( timefile, optarg );
            break;
        case 'o':              /* -o sets old division style (ndivs = 3 only! ) */
            olddivstyle = 1;
            break;
        case 'p':              /* -p sets print stepsize for output */
            if( timefile )
                error( "unfold: can't use -j and -p together!" );
            if( time != -999999999. )
                error( "unfold: can't use -t and -p together!" );
            p_stepsize = atof( optarg );
            if( p_stepsize < 0.001 )
                error( "unfold: output stepsize (%g) too small (min 0.001)", p_stepsize );
            break;
        case 'r':
            error( "unfold: -r is not supported anymore, use -g instead" );
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
            /*else if (!(strcmp(optarg, "bnd")))
               ps = Band; */
            else
                error( "unfold: invalid solver (%s), use: a,bd,bs,e,h,kr,mi,me,r{2,4,ck,f}", optarg );
            break;
        case 't':
            if( timefile )
                error( "unfold: can't use -j and -t together!" );
            if( p_stepsize > 0. )
                error( "unfold: can't use -p and -t together!" );
            time = atof( optarg );
            if( ( time < 0 ) && ( time != -999999999 ) )
                error( "unfold: the time (%g) doesn't make sense", time );
            break;
        case 'v':              /* -v prints version message */
            //fprintf(stderr, verstring, *argv, VERS, USR, MACHINE, COMPILER, FLAGS, __DATE__, __TIME__);
            exit( 0 );
        case 'x':
            if( ( strcmp( optarg, "input" ) ) && ( strcmp( optarg, "eqparms" ) ) && ( strcmp( optarg, "parameters" ) ) )
                error( "unfold: invalid section title (%s)", optarg );
            section_title = strcpy( section_title, optarg );
            break;
        case 'z':
            custom_gast = atof( optarg );
            if( custom_gast < 0. )
                error( "unfold: gastrulation time must be positive" );
            if( custom_gast > 10000000. )
                error( "unfold: gastrulation time must be smaller than 10'000'000" );
            break;
        case ':':
            error( "unfold: need an argument for option -%c", optopt );
        case '?':
        default:
            error( "unfold: unrecognized option -%c", optopt );
        }

    /* error check */

    if( stepsize > p_stepsize && p_stepsize != 0 )
        error( "unfold: print-stepsize (%g) smaller than stepsize (%g)!", p_stepsize, stepsize );

    if( ( argc - ( optind - 1 ) ) < 2 || ( argc - ( optind - 1 ) ) > 4 )
        PrintMsg( usage, 1 );

    /* let's get started and open the data file here */

    infile = argv[optind];
    fp = fopen( infile, "r" );
    if( !fp )
        file_error( "unfold" );

    /* Initialization code here */

    if( ( optind + 1 ) < argc )
        genindex = atoi( argv[optind + 1] );

    if( ( optind + 2 ) < argc ) {
        genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        genotype = strcpy( genotype, argv[optind + 2] );        /* genotype string */
    }

    inp.zyg = InitZygote( fp, pd, pj, &inp, section_title );
    inp.sco.facts = InitFacts( fp, &inp );      /* initializes facts */

    //printf("...ok!\nInitHis...");
    inp.his = InitHistory( fp, &inp );  //It fills the polations vector
    //printf("...ok!\nInitExtinp...");
    inp.ext = InitExternalInputs( fp, &inp );
    //printf("...ok!\nInitStepsize...");
    inp.ste = InitStepsize( stepsize, accuracy, slog, infile );
    // read the list of parameters to be tweaked
    inp.lparm = CopyParm( inp.zyg.parm, &( inp.zyg.defs ) );

    /* initialize genotype if necessary, otherwise check for errors */
    if( !( genotype ) ) {
        genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        for( i = 0; i < inp.zyg.defs.ngenes; i++ )
            genotype = strcat( genotype, "W" ); /* construct wt genotype string */
    } else {
        length = ( int ) strlen( genotype );
        if( length != inp.zyg.defs.ngenes )
            error( "unfold: genotype length doesn't match number of genes (%d)", inp.zyg.defs.ngenes );
        for( i = 0; i < inp.zyg.defs.ngenes; i++ ) {
            c = ( int ) *( genotype + i );
            if( ( c != 87 ) && ( ( c < 82 ) || ( c > 84 ) ) )
                error( "unfold: genotype string can only contain R, S, T or W" );
        }
    }

    /* Now we will read and traverse the genotypes list and find the
       facts section title for the genotype and set the facts list in
       integrate.c */
    geno = ReadGenotypes( fp, inp.zyg.defs.ngenes );
    curr = geno;

    if( !( genindex < count_Slist( geno ) ) )
        error( "Unfold: genindex more than last allele!\n" );

    i = 0;
    while( i++ < genindex )
        curr = curr->next;

    if( !( interp_dat = ( DataTable ** ) calloc( 1, sizeof( DataTable * ) ) ) )
        error( "Unfold: could not allocate interp_dat struct" );

    if( !( polation = ( InterpObject * ) calloc( 1, sizeof( InterpObject ) ) ) )
        error( "Unfold: could not allocate polation struct" );

    GetInterp( fp, curr->hist_section, &inp, inp.zyg.defs.ngenes, interp_dat );

    DoInterp( *interp_dat, polation, inp.zyg.defs.ngenes, &inp.zyg );
    FreeFacts( *interp_dat );

    theta_discons = Get_Theta_Discons( &theta_discons_size, &( inp.zyg ) );


    polation->fact_discons = ( double * ) realloc( polation->fact_discons, ( polation->fact_discons_size + theta_discons_size ) * sizeof( double ) );

    for( ii = 0; ii < theta_discons_size; ii++ )
        polation->fact_discons[polation->fact_discons_size + ii] = theta_discons[ii];

    polation->fact_discons_size += theta_discons_size;
    free( theta_discons );
    if( inp.zyg.defs.ndivs > 0 ) {
        if( !( temp_divtable = inp.zyg.times.div_times ) )
            error( "Unfold: error getting temp_div_times" );
        if( !( temp_durations = inp.zyg.times.div_duration ) )
            error( "Unfold: error getting division temp_durations" );

        for( ii = 0; ii < inp.zyg.defs.ndivs; ii++ ) {
            polation->fact_discons = ( double * ) realloc( polation->fact_discons, ( polation->fact_discons_size + 4 ) * sizeof( double ) );

            polation->fact_discons[polation->fact_discons_size] = temp_divtable[ii];
            polation->fact_discons[polation->fact_discons_size + 1] = temp_divtable[ii] + EPSILON;
            polation->fact_discons[polation->fact_discons_size + 2] = temp_divtable[ii] - temp_durations[ii];
            polation->fact_discons[polation->fact_discons_size + 3] = temp_divtable[ii] - temp_durations[ii] + EPSILON;

            polation->fact_discons_size += 4;
        }
    }

    qsort( ( void * ) polation->fact_discons, polation->fact_discons_size, sizeof( double ), ( int ( * )( const void *, const void * ) ) compare );

    /* Below we prepare the interpolant for external inputs */

    if( !( extinp_polation = ( InterpObject * ) calloc( 1, sizeof( InterpObject ) ) ) )
        error( "Unfold: could not allocate extinp_polation struct" );

    GetInterp( fp, curr->ext_section, &inp, inp.zyg.defs.egenes, interp_dat );

    DoInterp( *interp_dat, extinp_polation, inp.zyg.defs.egenes, &( inp.zyg ) );

    FreeFacts( *interp_dat );

    free_Slist( geno );


    /* initialize debugging file names */

    if( debug ) {

        slogfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        dumpfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

        sprintf( slogfile, "%s.slog", infile );
        sprintf( dumpfile, "%s.%s.%d.uout", infile, genotype, genindex );

        slog = fopen( slogfile, "w" );  /* delete existing slog file */
        fclose( slog );

        slog = fopen( slogfile, "a" );  /* now keep open for appending */
    } else {
        
        slogfile = NULL;
        dumpfile = NULL;
        slog = NULL;
    }

    /* the following establishes the tabulated times according to either:      */
    /* -t (output single time step)                                            */
    /* -j (read times from file)                                               */
    /* -p (stepsize set by argument)                                           */
    /* default behavior is producing output for gastrulation time only         */

    if( time != -999999999. ) {
        if( time > inp.zyg.times.gast_time )
            error( "unfold: time (%g) larger than gastrulation time!", time );
        tt.size = 1;
        tt.array = ( double * ) calloc( 1, sizeof( double ) );
        tt.array[0] = time;
    } else if( timefile != NULL ) {
        tt = ReadTimes( timefile, inp.zyg );
        free( timefile );
    } else if( p_stepsize != 0. ) {
        tt = MakeTable( p_stepsize, &( inp.zyg ) );
    } else {
        tt.size = 1;
        tt.array = ( double * ) calloc( 1, sizeof( double ) );
        tt.array[0] = inp.zyg.times.gast_time;
    }
    free( inp.sco.facts.tt[genindex].ptr.times.array );
    inp.sco.facts.tt[genindex].ptr.times = tt;
    /* Run the model... */

    for( i = 0; i < 1; ++i ) {
        if( i > 0 ) {
            FreeSolution( &answer );
        }
        answer = Blastoderm( genindex, genotype, &inp, slog );
    }
    /* if debugging: print out the innards of the model to unfold.out */

    if( debug ) {
        dumpptr = fopen( dumpfile, "w" );
        if( !dumpptr ) {
            perror( "unfold" );
            exit( 1 );
        }
        PrintBlastoderm( dumpptr, answer, "debug_output", MAX_PRECISION, &( inp.zyg ) );
        fclose( dumpptr );
    }

    /* strip output of anything that's not in tt */
    

    outtab = ConvertAnswer( answer, tt );
 //   printf("HAS TO BE AROUND 241.96 AND IS %lg\n\n", outtab.array[3].state.array[1]); //TESTING FOR BERTA

    //printf("outtabsize %d, answersize %d, ttsize %d\n", outtab.size, answer.size, tt.size);

    FreeMutant( inp.lparm );
    /* code for printing guts... first read gutsdefs section */
    if( guts ) {

        fp = fopen( infile, "r" );
        if( !fp ) {
            perror( "unfold" );
            exit( 1 );
        }
        gutsdefs = ReadGuts( fp );

        /* Make sure that interp_info in integrate.c is pointing at the
         * external inputs InterpObject */

        inp.ext = extinp_polation;

        gutstitle = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

        /* then calculate guts and print them */
        Zygote fakezyg;
        char *title;
        fakezyg.lin_start = inp.zyg.lin_start;
        fakezyg.defs = inp.zyg.defs;
        fakezyg.full_nnucs = NULL;
        fakezyg.nnucs = NULL;
        fakezyg.bcdtype = NULL;
        fakezyg.times = inp.zyg.times;
        for( i = 0; *( gutsdefs + i ); i++ ) {
            title = *( gutsdefs + i );
            if( ( numguts = CalcGuts( genindex, genotype, polation, extinp_polation, outtab, &goutput, title, &inp ) ) ) {
                gutstitle = strcpy( gutstitle, "guts_for_" );
                fakezyg.defs.ngenes = numguts;
                PrintBlastoderm( stdout, goutput, strcat( gutstitle, title ), ndigits, &fakezyg );
                FreeSolution( &goutput );
            }
            free( *( gutsdefs + i ) );
        }
        fclose( fp );
        free( gutsdefs );
        free( gutstitle );
        free( fakezyg.lin_start );
        free( fakezyg.full_nnucs );
        free( fakezyg.nnucs );
        free( fakezyg.bcdtype );

        /* code for printing model output */

    } else {
        PrintBlastoderm( stdout, outtab, "output\n", ndigits, &( inp.zyg ) );
    }

    /* ... and then clean up before you go home. */

    getrusage( RUSAGE_SELF, &end );     /*    get end time */
    printf( "# Unfold ran for %.13f seconds\n", tvsub( end, begin ) );

    FreeHistory( inp.zyg.nalleles, inp.his );
    FreeSolution( &answer );
    FreeExternalInputs( inp.zyg.nalleles, inp.ext );
    FreeZygote(  );
    free( extinp_polation );
    free( polation );
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
    }
    free( inp.sco.facts.tt );
    free( inp.sco.facts.facttype );
    //free(inp.sco.weights.weighttype);

    free( inp.zyg.defs.egene_ids );
    free( inp.zyg.defs.gene_ids );
    free( inp.zyg.nnucs );
    free( inp.zyg.full_nnucs );
    free( inp.zyg.full_lin_start );
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
    free( interp_dat );

    for( i = 0; i < outtab.size; i++ ) {
        free( outtab.array[i].state.array );
    }
    free( outtab.array );
    free( section_title );

    free( genotype );

    if( debug ) {
        fclose( slog );
        free( dumpfile );
        free( slogfile );
    }
    return 0;
}
