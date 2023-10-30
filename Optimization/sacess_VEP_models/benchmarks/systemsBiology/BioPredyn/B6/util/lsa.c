/**
 *                                                                 
 *   @file lsa.c                                                   
 *                                                                 
 *******************************************************************
 *   Credits:                                                      
 *                                                                 
 *   Originally written by Jimmy Lam and Dan Greening              
 *   Adaptation for continuous problems and original implemen-     
 *   tation of the parallel algorithm by John Reinitz              
 *   Tuning, equilibration and quenchit mode by King-Wai Chu       
 *   Partly rewritten, extended & documented by Johannes Jaeger    
 *                                                                 
 *******************************************************************
 *                                                                 
 *   see ../doc/prog_man.ps for details of how this works          
 *                                                                 
 *******************************************************************
 *                                                                 
 * IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL    
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND     
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!    
 *                                                                 
 *******************************************************************
 *                                                                 
 * Copyright (C) 1989-2003 John Reinitz                            
 *                                                                 
 * This program is free software; you can redistribute it and/or   
 * modify it under the terms of the GNU General Public License as  
 * published by the Free Software Foundation; either version 2 of  
 * the License, or (at your option) any later version.             
 *                                                                 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   
 * GNU General Public License for more details.                    
 *                                                                 
 * You should have received a copy of the GNU General Public       
 * License along with this program; if not, write to the           
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,     
 * Boston, MA  02111-1307, U.S.A.                                  
 *                                                                 
 */


#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <float.h>              /* for double limits */
#include <limits.h>             /* for integer limits */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>          /* these two are for times() */
#include <sys/times.h>
#include <time.h>               /* this is for time() */
#include <unistd.h>             /* for command line option stuff and access() */

/* NOTE: do not ever dare to include moves.h in here or in any of the hea- */
/*       ders below; lsa.c must remain truly problem independent           */

#include <sa.h>
//#include <error.h>
#include <random.h>
#include <distributions.h>

#ifdef MPI
#include <mpi.h>
#include <MPI.h>
#endif



/* STATIC VARIABLES ********************************************************/

static Files files;
/* some energy-related variables *******************************************/

static double energy;           /* current energy */

static double S;                /* current inverse energy */
static double dS;               /* delta S: change in S during move */
static double S_0;              /* the initial inverse temperature */

static double exp_arg;          /* the exponent of the Metropolis criterion (-dE/T) */

/* Lam stats stuff: estimators, stats and acceptance ratios ****************/

static double mean;             /* mean energy, collected from tau last steps */
static double vari;             /* energy variance, collected from tau last steps */

static double estimate_mean;    /* Lam estimator for mean energy */
static double estimate_sd;      /* Lam estimator for energy standard deviation */

static int success;             /* number of successful moves */
#ifdef MPI
static int l_success;           /* local number of successful moves */
#endif

static double alpha;            /* the third term of the Lam schedule formula */
static double acc_ratio;        /* average acceptance ratio for all parameters */

DistParms distp;                /*variable set by InitDistribution() */

#ifdef MPI
/* parallel code: variables for local Lam stats for tuning */

/* local weights: there are two sets of local estimators for the mean ******
 * energy, since we require a local Lam estimator of the mean energy with  *
 * a small weight for the calculation of cross-correlation of processors   *
 * (for estimating the lower bound of M), whereas we need a local Lam esti-*
 * mator for the mean energy with a large weight for the calculation of    *
 * the variance of local means; these two different weights for local Lam  *
 * stats have been determined experimentally by King-Wai Chu (although     *
 * it's not mentioned in his thesis); it works for all practical purposes, *
 * but in the future, a better way of sampling local statistics will be    *
 * needed (probably as part of a general theory of parallel Lam simulated  *
 * annealing)                                                              *
 *                                                                         *
 * variables for estimators for the lower bound of M end with _l           *
 * variables for estimators for the upper bound of M end with _u           *
 *                                                                         *
 * note that the l_vari is calculated using the upper bound estimator for  *
 * the mean energy (l_estimate_mean_u); while tuning, the upper bound Lam  *
 * statistics are written to the .llog files                               *
 *                                                                         *
 * note that we don't need seperate estimators for the standard since the  *
 * only thing we do with the local sd estimators is writing them to local  *
 * .llog files (only the upper bound estimators get written there)         */

static double l_mean;           /* local mean energy, from proc_tau last steps */
static double l_vari;           /* local energy variance, from proc_tau last steps */

static double l_estimate_mean_l;        /* estimator for mean for lower bound */
static double l_estimate_mean_u;        /* estimator for mean for upper bound */
static double l_estimate_sd;    /* Lam estimator for sd for upper bound */

static double l_alpha;          /* the third term of the Lam schedule formula */
static double l_acc_ratio;      /* average acceptance ratio for all parameters */

#endif

/* Lam stats stuff: variables for calculating the estimators  **************/

/* mean estimator */

static double w_a;              /* w_a is the weight for the mean */

static double usyy;             /* these parameters store intermediate results */
static double usxy;             /* for the updating formulas for A and B */
static double usxx;             /* see Lam & Delosme, 1988b, p10 */
static double usx;
static double usy;
static double usum;

static double A;                /* A and B are the parameters for the rational */
static double B;                /* function for the estimation of the mean */

#ifdef MPI

/* parallel code: variables for calculating local mean estimator for lower *
 *                bound of M (see also comment on local weights above)     */

static double l_w_a_l;          /* l_w_a_l: weight for the mean for the lower bound */

static double l_usyy_l;         /* these parameters store intermediate results */
static double l_usxy_l;         /* for the updating formulas for l_A_l and l_B_l */
static double l_usxx_l;
static double l_usx_l;
static double l_usy_l;
static double l_usum_l;

static double l_A_l;            /* l_A_l and l_B_l: parameters for the rational func */
static double l_B_l;            /* for the estimation of the mean for the lower bound */

/* parallel code: variables for calculating local mean estimator for upper *
 *                bound of M (see also comment on local weights above)     */

static double l_w_a_u;          /* l_w_a_u: weight for the mean for the upper bound */

static double l_usyy_u;         /* these parameters store intermediate results */
static double l_usxy_u;         /* for the updating formulas for l_A_u and l_B_u */
static double l_usxx_u;
static double l_usx_u;
static double l_usy_u;
static double l_usum_u;

static double l_A_u;            /* A and B are the parameters for the rational */
static double l_B_u;            /* function for the estimation of the mean */

#endif

/* sd estimator */

static double w_b;              /* w_b is the weight for the standard deviation */

static double vsyy;             /* these parameters store intermediate results */
static double vsxy;             /* for the updating formulas for D and E */
static double vsxx;             /* see Lam & Delosme, 1988b, p10 */
static double vsx;
static double vsy;
static double vsum;

static double D;                /* D and E are the parameters for the rational */
static double E;                /* function for the estimation of the standard deviation */

#ifdef MPI

/* parallel code: variables for local sd estimator (see also comment on    *
 *                local weights above)                                     */

static double l_w_b;            /* l_w_b: local weight for sd */

static double l_vsyy;           /* these parameters store intermediate results */
static double l_vsxy;           /* for the updating formulas for l_D and l_E */
static double l_vsxx;
static double l_vsx;
static double l_vsy;
static double l_vsum;

static double l_D;              /* l_D_l and l_E_l: parameters for the rational func */
static double l_E;              /* for the estimation of the sd for the lower bound */

#endif

/* Lam stats stuff: variables related to tau *******************************/

static double Tau;              /* double version of tau to calculate mean and vari */
static int proc_tau;            /* proc_tau = tau                     in serial   */
/* proc_tau = tau / (# of processors) in parallel */
static long count_tau;          /* how many times we did tau (or proc_tau) moves */

/* the actual number of moves for collecting initial statistics ************/

static int proc_init;           /* number of initial moves */

#ifdef MPI
/* an array needed for mixing in parallel code *****************************/

static int *dance_partner;      /* stores dance partners for each node */
#endif

/* stuff used by Frozen ****************************************************/

static double old_mean;         /* old mean as stored by Frozen */
static int counter;             /* counter used by Frozen */

/* skip tells the annealer how often it should update S ********************
 * implemented by Lam for increasing performance, i.e. currently we only   *
 * update the temperature every 10 moves (possibly obsolete now)           */

static int skip = -1;

/* flag used by Landscape generation ****************************************/
/*      Set by InitLandscape called from xxx_sa.c****************************/
static int landscape = 0;

/* vars used for equilibration runs ****************************************/
/* note: these need to be static to be passed on to the file that writes   */
/*       them to wherever they need to be written to (see GetEquil())      */

static ChuParam equil_param;    /* equilibration parameter struct */

static double fix_T_avg = 0.0;  /* overall energy average at fixed temp */
static double fix_T_var = 0.0;  /* energy variance at fixed temperature */

#ifdef MPI
//MPI constants

const int MAX_MIX = 10000;      /* max number of mixes during a tuning run */
/* set this to a lower number if you run out of memory */
const int GROUP_SIZE = 10;      /* group size for calculating upper bound */
const int STOP_TUNE_CNT = 20;   /* stop tune count */
const int STOP_TUNE_CRIT = 0.05;        /* tuning stop criterion */
const int LSTAT_LENGTH = 1;     /* length of Lam msg array when annealing */
const int LSTAT_LENGTH_TUNE = 28;       /* length of Lam msg array when tuning */



static double pfix_T_avg = 0.0; /* global energy average at fixed temp */
static double pfix_T_var = 0.0; /* global energy variance at fixed temp */

/* vars used for tuning runs ***********************************************/

/* general variables that determine tune- and sub_tune_interval */

static int covar_sample;        /* # of moves over which we sample tuning stat */
static int sample_size;         /* so many samples per processor */
static int tune_interval;       /* total # of moves to be sampled for tuning */
static int sub_tune_interval;   /* at end of this we write tune stats */

/* various counter variables and an array to save counts */

static int count_sample = 0;    /* how many samples did we collect? */
static int count_tune = 0;      /* how many times did we do sub_tune? */
static int count_mix = 0;       /* counts the times we've been mixing */

static int moves_tune = 0;      /* move counter: reset every tune_interval */

static int *tau_count;          /* # of tau's we've done for whole tune_interval */

/* arrays for sampling stats for the lower bound */

static double *dev;             /* standard deviations for a sub_tune_interval */
static double *tot_dev;         /* gathers all local standard deviations */
static double *coll_dev;        /* collects local standard deviations */
static double *cross_correl;    /* cross-correlations for lower bound */

/* arrays for sampling stats for the upper bound */

static double *means;           /* used to save local means for upper bound */
static double *tot_means;       /* gathers all local means */
static double *coll_means;      /* collects local means */
static double *var_means;       /* variance of local means for lower bound */
static int *midpoints;          /* midpoints of groups for upper bound */

/* a variable for the tuning stop criterion */

static int stop_tune_count;     /* stop tune count */
#endif

/* variables used for timing */

/* these variables are used to evaluate real time using time() */
static double start;            /* wallclock time before run */
static double finish;           /* wallclock time after run */

/* these structs are used to evaluate user time using times() or MPI_Wtime */
static struct tms *cpu_start;   /* user time before run */
static struct tms *cpu_finish;  /* user time after run */

ScoreOutput out;



// GLOBAL CONSTANTS

/** The following defines the maximum float precision that is supported by  
 * the code.
 */
extern const int MAX_PRECISION;
/* the following constant as a score tells the annealer to reject a move,  */
/* no matter what. It had better not be a number that could actually be a  */
/* score.                                                                  */
extern const double FORBIDDEN_MOVE;  /* the biggest possible score, ever */

const int OUT_OF_BOUND = -1;

//GLOBAL SA CONSTANTS

const int MIN_DELTA = -100;     /* minimum exponent for Metropolis criterion */
/* provides a minimum probability for really bad moves */



/* MAIN HERE ***************************************************************/

/** main function for the Simulated Annealing optimization code */
int
main( int argc, char **argv ) {
    double *delta;              /* used to store elapsed times */
    proc_id = -1;

#ifdef MPI
    // MPI initialization steps
    MPI_Init( &argc, &argv );   // initializes the MPI execution environment
    MPI_Comm_size( MPI_COMM_WORLD, &nnodes );   // number of processors?
    MPI_Comm_rank( MPI_COMM_WORLD, &myid );     // ID of local processor?

    stop_tune_count = STOP_TUNE_CNT;    /* stop tune count */
    proc_id = myid;

    printf( "#%d: Starting Simulated Annealing\n", proc_id );
#else
    printf( "Starting Simulated Annealing\n" );
#endif



    /* code for timing: wallclock and user times */
    cpu_start = ( struct tms * ) malloc( sizeof( struct tms ) );        /* user time */
    cpu_finish = ( struct tms * ) malloc( sizeof( struct tms ) );
    times( cpu_start );

#ifdef MPI
    start = MPI_Wtime(  );      /* returns wallclock time on the calling node */
#else
    start = time( NULL );       /* returns wallclock time since EPOCH (1/1/1970) */
#endif
    /* initialize cost function and move state, do initial moves (or restore   */
    /* annealing state if restart)
     *    */
    Initialize( argc, argv );
    /* the following is for non-equlibration runs and equilibration runs that  */
    /* have not yet settled to their equilibrium temperature                   */
    if( ( bench != 1 ) && ( ( equil != 1 ) || ( 1.0 / S > equil_param.end_T ) ) ) {
        Loop(  );
    }

    /* there's an alternative Loop for equlibration runs at stable temperature */
    int fixloopcounter = 0;

    if( equil == 1 ) {
        printf( "%d: Starting fixloop %d\n", proc_id, fixloopcounter );
        FixTLoop(  );
        fixloopcounter++;
    }

    /* code for timing */

    if( time_flag ) {
        delta = GetTimes(  );   /* calculates times to be printed */
#ifdef MPI
        if( myid == 0 )
#endif
            WriteTimes( delta );        /* then write them out */
        free( delta );
    }

    /* clean up MPI and return */

#ifdef MPI
    MPI_Finalize(  );           /* terminates MPI execution environment */
#endif

    return 0;

}



/*** INITIALIZING FUNCTIONS ************************************************/

/**  Initialize: calls ParseCommandLine first; then does either initial    
 *               randomization and collecting Lam stats or restores state  
 *               of the annealer as saved in the state file                
 */
void
Initialize( int argc, char **argv ) {
    debug = 0;
    int opt_index;              /* pointer to current argument of command line */
    int stateflag = 0;          /* state file or not? */
#ifdef MPI
    int flagsum;                /* used for state file check below */
#endif
    double initial_temp;        /* initial temperature for annealer */

    //initialize output structure

    out.score = 1e38;           // start with a very large number
    out.penalty = 0;
    out.size_resid_arr = 0;
    out.jacobian = NULL;
    out.residuals = NULL;

    /* allocate memory for static file names */

    files.inputfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files.statefile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    /* allocate memory for static Lam parameters and dance partners */

    state = ( NucStateType * ) malloc( sizeof( NucStateType ) );

#ifdef MPI
    dance_partner = ( int * ) calloc( nnodes, sizeof( int ) );
#endif
    /* parse the command line and return index to input file name; then we     *
     * initialize the static file names for inputfile and state file; these    *
     * two are needed here, since they are required for reading the state file *
     * whereas all other file names depend on the output file name and need to *
     * be initialized after calling RestoreState() below                       */
    opt_index = ParseCommandLine( argc, argv );
    strcpy( files.inputfile, argv[opt_index] );

    /* state files: used for the case that a run terminates or crashes unex-   *
     * pectedly; we can then restore the state of the run *precisely* as it    *
     * was before the crash by restarting it from the state file               */

#ifdef MPI
    if( nnodes > 1 ) {
        if( nnodes <= 10 )
            sprintf( files.statefile, "%s_%d.state", files.inputfile, myid );
        else if( ( nnodes > 10 ) && ( nnodes <= 100 ) )
            sprintf( files.statefile, "%s_%02d.state", files.inputfile, myid );
        else if( ( nnodes > 100 ) && ( nnodes <= 1000 ) )
            sprintf( files.statefile, "%s_%03d.state", files.inputfile, myid );
        else if( ( nnodes > 1000 ) && ( nnodes <= 10000 ) )
            sprintf( files.statefile, "%s_%04d.state", files.inputfile, myid );
        else if( ( nnodes > 10000 ) && ( nnodes <= 100000 ) )
            sprintf( files.statefile, "%s_%05d.state", files.inputfile, myid );
        else
            error( "Initialize: can't open more than 100'000 state files" );
    } else
#endif
        sprintf( files.statefile, "%s.state", files.inputfile );

    /* check if a state file exists (access() is in unistd.h) */

    if( 0 == access( files.statefile, F_OK ) )
        stateflag = 1;

#ifdef MPI
    /* parallel code: make sure that all state files are present */

    MPI_Allreduce( &stateflag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( ( flagsum > 0 ) && ( flagsum != nnodes ) && ( stateflag == 0 ) )
        error( "Initialize: state file for process %d is missing" );
#endif

    /* first get Lam parameters, initial temp and energy and initialize S_0 */
        initial_temp = MoveSA( state, &distp, &out, &files, 1, 0 );

        //&energy is now the score value in the output structure
        //last two bits decide weather this is a initialization loop and if we want the Jacobian matrix to be calculated
        energy = out.score + out.penalty;
        S_0 = 1. / initial_temp;
    /* initialize those static file names that depend on the output file name */
    InitFilenames( &files );
#ifdef MPI
    /* note that for parallel code both tau and init must be divided by nnodes */
    /* and we need to account for the case when tau isn't divisible by nnodes  */

    if( ( state->tune.tau % nnodes ) != 0 )
        error( "fly_sa: the number of processors (%d) must be a divisor of tau (%d)", nnodes, state->tune.tau );
    proc_tau = state->tune.tau / nnodes;        /* local copy of tau */

    //printf("#%d: proc_tau = %d = %d/%d\n", proc_id, proc_tau, state->tune.tau, nnodes);

    if( ( state->tune.initial_moves % nnodes ) != 0 )
        error( "fly_sa: number of init moves must be divisible by nnodes (%d)", nnodes );
    proc_init = state->tune.initial_moves / nnodes;     /* # of initial moves */
#else
    proc_tau = state->tune.tau; /* static copy to tau */
    proc_init = state->tune.initial_moves;      /* # of initial moves */
#endif
    Tau = ( double ) state->tune.tau;   /* double version to tau */
    /* for calculating estimators */

    /* if we're not restarting: do the initial moves for randomizing and ga-   *
     * thering initial statistics                                              */
    if( !stateflag ) {
        InitialLoop(  );
    }
    /* write first .log entry and write first statefile right after init; note *
     * that equilibration runs are short and therefore don't need state files  *
     * which would be rather complicated because of all the stats collected    *
     * during equilibration; for the exact same reasons we don't write state   *
     * files for tuning either; in case of a restart, we just restore the .log */

    if( !equil && !bench && !nofile_flag ) {
        if( !stateflag ) {
            WriteLog(  );
#ifdef MPI
            if( !tuning )
#endif
                StateWrite( files.statefile );
        } else
            RestoreLog(  );
    }
#ifdef MPI
    /* if we are in tuning mode: initialize/restore tuning structs */

    if( tuning )
        InitTuning(  );

#endif

}

/**  InitFilenames: initializes static file names that depend on the out- 
 *                  put file name                  
 */
void
InitFilenames( Files * files ) {

    /* allocate memory for static file names */

    files->logfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
#ifdef MPI
    files->l_logfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files->lbfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files->ubfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files->mbfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
#endif

    /* if we're not using -w: output = inputfile */

    if( !files->outputfile )
        files->outputfile = files->inputfile;

#ifdef MPI
    if( myid == 0 ) {
#endif
        /* the global .log file: stores iterations, temperature, change in tempe-  *
         * rature, global means and standard deviation, Lam estimators for mean    *
         * and sd as well as global acceptance ratios over an annealing run        */

        sprintf( files->logfile, "%s.log", files->outputfile );

#ifdef MPI
        /* the following files are only needed for tuning */

        /* the lower_bound file: contains cross-correlations for a tune_interval   *
         * averaged over all the past tune_intervals used to calculate the lower   *
         * bound for the mixing_interval M                                         */

        sprintf( files->lbfile, "%s.lb", files->outputfile );

        /* the upper_bound file: for saving variance of local means for a tune_    *
         * interval averaged over all the past tune_intervals; used to calculate   *
         * the upper bound for the mixing_interval M when tuning                   */

        sprintf( files->ubfile, "%s.ub", files->outputfile );

        /* the mix_bound file: contains the history of both upper and lower bound  *
         * for M over all past tune_intervals; used to check for the convergence   *
         * of the estimate for the upper bound of M during a tuning run; also sto- *
         * res the cross-corrlation and variance of local means for both bounds    */

        sprintf( files->mbfile, "%s.mb", files->outputfile );
    }

    /* the local .llog file: used to store iterations, temperature, tempera-   *
     * ture change, local mean and stamdard deviation, local Lam estimators    *
     * for mean and sd (for the upper bound) and local acceptance ratios;      *
     * this is only needed when tuning, otherwise we just write one global log */

    if( nnodes <= 10 )
        sprintf( files->l_logfile, "%s_%d.llog", files->outputfile, myid );
    else if( ( nnodes > 10 ) && ( nnodes <= 100 ) )
        sprintf( files->l_logfile, "%s_%02d.llog", files->outputfile, myid );
    else if( ( nnodes > 100 ) && ( nnodes <= 1000 ) )
        sprintf( files->l_logfile, "%s_%03d.llog", files->outputfile, myid );
    else if( ( nnodes > 1000 ) && ( nnodes <= 10000 ) )
        sprintf( files->l_logfile, "%s_%04d.llog", files->outputfile, myid );
    else if( ( nnodes > 10000 ) && ( nnodes <= 100000 ) )
        sprintf( files->l_logfile, "%s_%05d.llog", files->outputfile, myid );
    else
        error( "Initialize: can't open more than 100'000 llog files" );
#endif
}

/**  InitialLoop: performs the two sets of initial moves: 
 *                   1. randomizing moves (not parallelized)               
 *                   2. loop for initial collection of statistics          
 */
void
InitialLoop( void ) {
    int i;                      /* local loop counter */

    double energy_change;       /* change of energy during move */

#ifdef MPI
    double total[2];            /* will hold mean & vari when pooling stats */
    double tmptotal[2];         /* temporary array for above */

    long nodesuccess[2];        /* array for success and summed init moves */
    long tmpsuccess[2];         /* temporary array for above */
#endif


    /* randomize initial state; throw out results; DO NOT PARALLELIZE! */
    for( i = 0; i < state->tune.initial_moves; i++ ) {
        /* make a move: will either return the energy change or FORBIDDEN_MOVE */
        energy_change = GenerateMove( &files, &distp, &out );
        if( ( i % 20000 == 0 ) && ( i > 0 ) ) {
            printf( "#%d Initial Loop Move = %d\n", proc_id, i );
        }


        /* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
         * that's why we want to prevent overflows here (hence the 'if')           */

        if( energy_change != FORBIDDEN_MOVE )
            exp_arg = -S_0 * energy_change;

        /* MIN_DELTA provides a min. probability with which any move is accepted */

        if( exp_arg <= MIN_DELTA )
            exp_arg = MIN_DELTA;

        /* below, we apply the Metropolis criterion to accept or reject a move */

        if( energy_change == FORBIDDEN_MOVE ) {
            RejectMove(  );
        } else if( ( energy_change <= 0.0 ) || ( exp( exp_arg ) > RandomReal(  ) ) ) {
            energy += energy_change;
            AcceptMove(  );
        } else {
            RejectMove(  );
        }

    }                           /* end randomize initial state */

    /* set all stats to zero, collection starts below */

    mean = 0.0;
    vari = 0.0;
    success = 0;

    /* loop to collect initial statistics; this one is parallelized */
    for( i = 0; i < proc_init; i++ ) {
        /* make a move: will either return the energy change or FORBIDDEN_MOVE */
        energy_change = GenerateMove( &files, &distp, &out );
        /* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
         * that's why we want to prevent overflows here (hence the 'if')           */
        if( energy_change != FORBIDDEN_MOVE )
            exp_arg = -S_0 * energy_change;

        /* MIN_DELTA provides a min. probability with which any move is accepted */

        if( exp_arg <= MIN_DELTA )
            exp_arg = MIN_DELTA;

        /* below, we apply the Metropolis criterion to accept or reject a move */
        if( energy_change == FORBIDDEN_MOVE ) {
            RejectMove(  );

        } else if( ( energy_change <= 0.0 ) || ( exp( exp_arg ) > RandomReal(  ) ) ) {
            energy += energy_change;
            AcceptMove(  );
            success++;
        } else {
            RejectMove(  );
        }

        /* collect stats */
        mean += energy;
        vari += ( energy * energy );
    }
#ifdef MPI

    /* collect local stats if tuning */

    if( tuning && nnodes > 1 ) {
        l_mean = mean;
        l_vari = vari;
        l_success = success;
    }

    /* parallel code: pool initial statistics from all nodes */

    total[0] = mean;
    total[1] = vari;
    nodesuccess[0] = ( long ) success;
    nodesuccess[1] = ( long ) i;        /* sum of i's must equal initial moves */

    /* the tmp arrays are used to send messages */

    tmptotal[0] = total[0];
    tmptotal[1] = total[1];
    tmpsuccess[0] = nodesuccess[0];
    tmpsuccess[1] = nodesuccess[1];

    /* stats from all nodes are summed up here */

    MPI_Allreduce( tmptotal, total, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( tmpsuccess, nodesuccess, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

    success = ( int ) nodesuccess[0];   /* success is now global success! */

    /* sanity check: have we done the correct number of initial moves? */

    if( nodesuccess[1] - state->tune.initial_moves )
        error( "InitialLoop: initial moves was %d?!\n", nodesuccess[1] );

    /* mean and variance are now summed over all nodes */

    mean = total[0];
    vari = total[1];

    /* local stats are calculated here (only if tuning) */

    if( tuning && nnodes > 1 ) {
        l_mean /= ( double ) proc_init;
        l_vari = l_vari / ( ( double ) proc_init ) - l_mean * l_mean;
        l_acc_ratio = ( ( double ) l_success ) / ( ( double ) proc_init );
    }
#endif
    /* global stats are calculated here */
    mean /= ( double ) state->tune.initial_moves;
    vari = vari / ( ( double ) state->tune.initial_moves ) - mean * mean;
    acc_ratio = ( ( double ) success ) / ( ( double ) state->tune.initial_moves );

    /* initialize Lam parameters used for calculating Lam estimators */
    InitializeParameter(  );
}

/** InitializeParameter: initializes variables for Lam annealing: this is  
 *                        executed after doing the initial steps;          
 *                        the local parameters only need to be set for tu- 
 *                        ning code                                        
 */
void
InitializeParameter( void ) {
    double d;                   /* d is used to store intermediate results below */

#ifdef MPI
    double l_d;                 /* l_d is used to store intermediate results */
#endif


    /* 1. set global parameters (serial and parallel code) ******************** */

    /* set estimators to stats collected during initializing phase of run */

    estimate_sd = sqrt( vari );
    estimate_mean = mean;

    /* initialize A,B,D,E according to Lam & Delosme, 1988b, p10 */

    A = estimate_sd * estimate_sd / ( estimate_mean * estimate_mean );
    B = ( 1.0 / estimate_mean ) - ( A * S_0 );
    D = estimate_sd / estimate_mean;
    E = ( 1.0 / estimate_sd ) - ( D * S_0 );

    /* initialize these intermediate variables for updating funcs for A,B,D,E */

    usum = vsum = 1.0;
    usxy = usxx = usx = 0.0;
    usy = 1.0 / estimate_mean;
    usyy = usy * usy;
    vsxy = vsxx = vsx = 0.0;
    vsy = 1.0 / estimate_sd;
    vsyy = vsy * vsy;

    /* set the initial temperature and the initial delta S */

    S = S_0;
    dS = 0.5 / estimate_sd;     /* keep--may not need--based on s_0=0 */

    /* alpha is the third term of the main Lam schedule formula */

    d = ( 1.0 - acc_ratio ) / ( 2.0 - acc_ratio );
    alpha = 4.0 * acc_ratio * d * d;

#ifdef MPI
    /* 2. set local parameters *************************************************
     * these are only needed when tuning...                                    *
     * note that the two sets of local estimators used for tuning only differ  *
     * in their weights, so they can actually be initialized the same way      */

    if( tuning && nnodes > 1 ) {

        /* set estimators to stats collected during initializing phase of run */

        l_estimate_sd = sqrt( l_vari );
        l_estimate_mean_l = l_estimate_mean_u = l_mean;

        /* initialize local A,B,D,E */

        l_A_l = l_A_u = l_estimate_sd * l_estimate_sd / ( l_estimate_mean_u * l_estimate_mean_u );
        l_B_l = l_B_u = ( 1.0 / l_estimate_mean_u ) - ( l_A_u * S_0 );

        l_D = l_estimate_sd / l_estimate_mean_u;
        l_E = ( 1.0 / l_estimate_sd ) - ( l_D * S_0 );

        /* initialize these intermediate variables for updating funcs for A,B,D,E */

        l_usum_l = l_usum_u = l_vsum = 1.0;

        l_usxy_l = l_usxy_u = l_usxx_l = l_usxx_u = l_usx_l = l_usx_u = 0.0;
        l_usy_l = l_usy_u = 1.0 / l_estimate_mean_u;
        l_usyy_l = l_usyy_u = l_usy_u * l_usy_u;

        l_vsxy = l_vsxx = l_vsx = 0.0;
        l_vsy = 1.0 / l_estimate_sd;
        l_vsyy = l_vsy * l_vsy;

        /* alpha is the third term of the main Lam schedule formula */

        l_d = ( 1.0 - l_acc_ratio ) / ( 2.0 - l_acc_ratio );
        l_alpha = 4.0 * l_acc_ratio * l_d * l_d;

    }
#endif

    /* weights determine how the estimators are sampled for times before tau */

    InitializeWeights(  );

}

/**  InitializeWeights: initialize weights a and b for calculating Lam 
 *                      estimators; these weights are computed from the    
 *                      lambda memory length products                      
 */
void
InitializeWeights( void ) {
    FILE *logptr;
#ifdef MPI
    FILE *l_logptr;
#endif

    /* w_a is the weight for the mean */

    w_a = state->tune.lambda_mem_length_u / state->tune.lambda;
    w_a = 1.0 - state->tune.tau / w_a;
    if( w_a < 0.0 )
        w_a = 0.0;

    /* w_b is the weight for the standard deviation */

    w_b = state->tune.lambda_mem_length_v / state->tune.lambda;
    w_b = 1.0 - state->tune.tau / w_b;
    if( w_b < 0.0 )
        w_b = 0.0;

#ifdef MPI
    /* local weights: there are two sets of local statistics, since we require *
     * a local Lam estimator of the mean energy with a small weight for the    *
     * calculation of cross-correlation of processors (for estimating the      *
     * lower bound of M), whereas we need a local Lam estimator for the mean   *
     * energy with a large weight for the calculation of the variance of local *
     * means; these two different weights for local Lam stats have been deter- *
     * mined experimentally by King-Wai Chu (although it's not mentioned in    *
     * his thesis); it works for all practical purposes, but in the future, a  *
     * better way of sampling local statistics will be needed (probably as     *
     * part of a general theory of parallel Lam simulated annealing)           */

    if( tuning && nnodes > 1 ) {

        /* l_w_a_{l|u} are the local weights for the mean */

        l_w_a_l = w_a / ( double ) nnodes;
        l_w_a_u = w_a;

        /* l_w_b is the local weight for the standard deviation */

        l_w_b = w_b;

    }
#endif

    if( !equil && !bench && !nofile_flag ) {
#ifdef MPI
        if( myid == 0 ) {
#endif
            logptr = fopen( files.logfile, "w" );
            if( !logptr )
                file_error( "InitializeWeights" );
            fprintf( logptr, "InitializeWeights:  w_a = %g w_b = %g\n", w_a, w_b );
            fclose( logptr );
#ifdef MPI
        }

        if( write_llog && tuning && nnodes > 1 ) {
            l_logptr = fopen( files.l_logfile, "w" );
            if( !l_logptr )
                file_error( "InitializeWeights" );
            fprintf( l_logptr, "InitializeWeights:  l_w_a = %g l_w_b = %g\n", l_w_a_u, l_w_b );
            fclose( l_logptr );
        }

        if( myid == 0 )
#endif
            if( log_flag )
                printf( "InitializeWeights:  w_a = %g w_b = %g\n", w_a, w_b );
    }
}



#ifdef MPI

/** InitTuning: sets up/restores structs and variables for tuning runs */
void
InitTuning( void ) {
    int i;                      /* loop counter */
    FILE *mbptr;                /* pointer for mix_bound file */


    /* error check */

    if( nnodes <= 1 )
        error( "fly_sa: tuning does not make sense on one processor" );

    if( covar_index > state->tune.mix_interval )
        error( "fly_sa: you can't sample over more than the whole mix interval" );

    /* set size of sample interval */

    covar_sample = covar_index * ( int ) Tau;

    /* set size of sample interval per processor */
    /* by the way: sample_size corresponds to covar_index * proc_tau */

    if( ( covar_sample % nnodes ) != 0 )
        error( "fly_sa: covar_sample (%d) not divisible by nnodes", covar_sample );

    sample_size = covar_sample / nnodes;

    /* tune_interval: how many covar_samples per mix_interval? */

    if( ( state->tune.mix_interval % covar_index ) != 0 )
        error( "fly_sa: mix interval (%d) not divisible by covar_index", state->tune.mix_interval );

    tune_interval = state->tune.mix_interval / covar_index;

    /* size of every tune interval in between writing tuning stats */

    if( write_tune_stat > tune_interval )
        error( "fly_sa: freq of writing tune stats (%d) must be smaller than %d", write_tune_stat, tune_interval );

    if( ( tune_interval % write_tune_stat ) != 0 )
        error( "fly_sa: tune_interval (%d) not divisible by write_tune_stat (%d)", tune_interval, write_tune_stat );

    sub_tune_interval = tune_interval / write_tune_stat;

    /* allocate memory for various tuning-specific arrays */

    dev = ( double * ) calloc( sub_tune_interval * sample_size, sizeof( double ) );
    tot_dev = ( double * ) calloc( sub_tune_interval * covar_sample, sizeof( double ) );
    coll_dev = ( double * ) calloc( covar_sample, sizeof( double ) );

    means = ( double * ) calloc( sub_tune_interval * sample_size, sizeof( double ) );
    tot_means = ( double * ) calloc( sub_tune_interval * covar_sample, sizeof( double ) );
    coll_means = ( double * ) calloc( covar_sample, sizeof( double ) );

    tau_count = ( int * ) calloc( tune_interval, sizeof( int ) );
    cross_correl = ( double * ) calloc( tune_interval, sizeof( double ) );
    var_means = ( double * ) calloc( tune_interval, sizeof( double ) );

    midpoints = ( int * ) calloc( MAX_MIX, sizeof( int ) );

    /* initialize the arrays used to collect tuning information */

    for( i = 0; i < tune_interval; i++ ) {
        tau_count[i] = ( i + 1 ) * covar_index;
        cross_correl[i] = 0.;
        var_means[i] = 0.;
    }

    /* initialize the mb file */

    if( myid == 0 ) {
        mbptr = fopen( files.mbfile, "w" );
        if( !mbptr )
            file_error( "InitTuning (writing mb captions)" );
        fprintf( mbptr, "# mix_count  M_lower    cross_corr" );
        fprintf( mbptr, "  M_upper    var_l_mean\n\n" );
        fclose( mbptr );
    }
}
#endif





/*** MAIN LOOP AND UPDATE FUNCTIONS ****************************************/

/**  Loop: loops (making moves, updating stats etc.) until the system is 
 *         considered frozen according to the stop criterion               
 */
void
Loop( void ) {
    int i;                      /* local loop counter */
    double energy_change;       /* local Delta E */
    double d;                   /* difference between energy and estimated mean */


    /* quenchit mode: set temperature to (approximately) zero immediately */

    if( quenchit )
        S = DBL_MAX;

    /* loop till the end of the universe (or till the stop criterion applies) */
    int loopcounter = 0;
    while( 1 ) {
        /*if (loopcounter % 100000 == 0)
           printf("#%d: Loop %d\n", proc_id, loopcounter);
         */
        loopcounter++;
        /* reset statistics */

        mean = 0.0;
        vari = 0.0;
        success = 0;

#ifdef MPI
        l_mean = 0.0;
        l_vari = 0.0;
        l_success = 0;
#endif
        /* do proc_tau moves here */
        for( i = 0; i < proc_tau; i++ ) {

            /* make a move: will either return the energy change or FORBIDDEN_MOVE */
            energy_change = GenerateMove( &files, &distp, &out );
            /* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
             * that's why we want to prevent overflows here (hence the 'if'); we also  *
             * need to avoid overflows with quenchit (where S is (almost) infinite!)   */

            if( !quenchit && ( energy_change != FORBIDDEN_MOVE ) )
                exp_arg = -S * energy_change;

            /* MIN_DELTA provides a min. probability with which any move is accepted */

            if( ( exp_arg <= MIN_DELTA ) )
                exp_arg = MIN_DELTA;

            /* below, we apply the Metropolis criterion to accept or reject a move; in *
             * quenchit mode, only lower energies are accepted                         */

            if( energy_change == FORBIDDEN_MOVE ) {
                RejectMove(  );
            } else if( ( energy_change <= 0.0 ) || ( ( !quenchit ) && ( exp( exp_arg ) > RandomReal(  ) ) ) ) {
                energy += energy_change;
                AcceptMove(  );
                success++;
#ifdef MPI
                if( tuning && nnodes > 1 )
                    l_success++;
#endif
            } else {
                RejectMove(  );
            }

            /* update statistics */

            mean += energy;
            d = energy - estimate_mean;
            vari += d * d;
#ifdef MPI
            /* if tuning: calculate the local mean and variance */

            if( tuning && nnodes > 1 ) {

                l_mean += energy;
                d = energy - l_estimate_mean_u;
                l_vari += d * d;

                /* Collect samples for tuning: local estimated deviations (for calculating *
                 * cross-correlation between processors for the lower bound of M) and lo-  *
                 * cal estimated means (for calculating the variance of those local means  *
                 * between processors for the upper bound of M); here, we also need to     *
                 * keep track of the moves we've done in the current sub_tune_interval and *
                 * of the number of samples we've collected in the tune_interval           */

                dev[moves_tune] = energy - l_estimate_mean_l;
                means[moves_tune] = l_estimate_mean_u;
                moves_tune++;
                if( ( moves_tune % sample_size ) == 0 )
                    count_sample++;
            }
#endif

            /* update temperature every 'skip' steps; this was put here by Jimmy Lam,  *
             * probably to save computation time on his old Spark; I guess it's obso-  *
             * lete by now, but since it doesn't seem to do any harm and saves us some *
             * time, we left it in here                                                */

            if( !quenchit )
                if( --skip <= 0 )
                    UpdateS(  );
        }                       /* this is the end of the proc_tau loop */

        /* have done tau moves here: update the 'tau' counter */

        count_tau++;

        /* calculate mean, variance and acc_ratio for the last tau steps; i is     *
         * passed as an argument for checking if all local moves add up to Tau     */

        UpdateStats( i );

        /* check if the stop criterion applies: annealing and tuning runs (that    *
         * aren't stopped by the tuning stop criterion) leave the loop here; equi- *
         * libration runs exit below                                               */

        if( Frozen(  ) && !equil ) {
            FinalMove(  );
            return;
        }

        /* equilibration: if we have reached equilibration temperature, we exit    *
         * Loop() and continue the run at a fixed temperature with FixTLoop()      */
        if( ( equil == 1 ) && ( 1.0 / S <= equil_param.end_T ) ) {
            S = 1. / equil_param.end_T;
            return;
        }
        /* update Lam stats: estimators for mean, sd and alpha from acc_ratio (we  *
         * don't need this in quenchit mode since the temperature is fixed to 0)   */

        if( !quenchit )
            UpdateParameter(  );
#ifdef MPI

        /* tuning code: first update local Lam estimators */

        if( tuning && nnodes > 1 ) {
            UpdateLParameter(  );

            /* Do tuning every sub_tune_interval, only after first mix */

            if( ( count_sample % sub_tune_interval ) == 0 && count_sample > 0 ) {
                if( count_mix > 0 ) {
                    DoTuning(  );
                } else {
                    moves_tune = 0;
                }
            }

            /* At the end of each tune interval: ***************************************
             * 1. Root process writes tuning every tune_interval, only after first mix *
             *    NOTE: this needs to get done BEFORE mixing, otherwise count_mix      *
             *          isn't right and you'll get a division by zero in WriteTuning() *
             * 2. If we've done more that 'stop_tune_count' mixes, we check if we can  *
             *    stop tuning run (runs can be forced to continue using the -S option) *
             ***************************************************************************/

            if( ( count_sample % tune_interval ) == 0 && count_sample > 0 ) {

                if( count_mix > 0 )
                    WriteTuning(  );

                if( ( count_mix >= stop_tune_count ) && auto_stop_tune )
                    if( StopTuning(  ) ) {
                        free( dev );
                        free( tot_dev );
                        free( coll_dev );
                        free( means );
                        free( tot_means );
                        free( coll_means );
                        free( tau_count );
                        free( cross_correl );
                        free( var_means );
                        free( midpoints );
                        FinalMove(  );
                        return; /* exit the loop here if finished tuning */
                    }
            }
        }

        /* at each mix_interval: do some mixing */

        if( count_tau % state->tune.mix_interval == 0 ) {
            DoMix(  );
        }
#endif

        /* write the log every print_freq * tau (not proc_tau!) */

#ifdef MPI
        if( ( count_tau % ( print_freq * nnodes ) == 0 ) && !equil && !nofile_flag )
#else
        if( ( count_tau % print_freq == 0 ) && !equil && !nofile_flag )
#endif
            WriteLog(  );
        /* the state file gets written here every state_write * tau */

#ifdef MPI
        if( ( count_tau % ( state_write * nnodes ) == 0 ) && !equil && !tuning && !nofile_flag ) {
#else
        if( ( count_tau % state_write == 0 ) && !equil && !nofile_flag ) {
#endif
            StateWrite( files.inputfile );
        }

    }                           /* this is the end of the while(1) loop */

}                               /* this is the end of Loop */

/** UpdateS: update inverse temperature S at every Sskip step */
void
UpdateS( void ) {
    register double d;          /* used to store intermediate results */

    S += dS;                    /* here, inverse temperature is updated by dS */

    /* we need to update Lam parameters here, since S has changed; A, B, C and */
    /* D get updated in UpdateParameters                                       */

    estimate_mean = 1.0 / ( A * S + B );        /* for temperature updating formulas */
    estimate_sd = 1.0 / ( D * S + E );  /* see Lam & Delosme, 1988b, p10 */

#ifdef MPI

    /* do the same for local Lam parameters (for both lower and upper bounds) */

    if( tuning && nnodes > 1 ) {
        l_estimate_mean_l = 1.0 / ( l_A_l * S + l_B_l );
        l_estimate_mean_u = 1.0 / ( l_A_u * S + l_B_u );
        l_estimate_sd = 1.0 / ( l_D * S + l_E );
    }
#endif

    d = S * estimate_sd;        /* intermediate for the specific heat */

    /* following lines implement the main Lam schedule formula */

    dS = state->tune.lambda * alpha / ( d * d * estimate_sd );
    dS *= state->tune.update_S_skip;    /* ... we have to muliply by skip */

#ifdef MPI
    dS *= Tau / ( ( double ) proc_tau );
#endif

    /* reset skip */

    skip = state->tune.update_S_skip;

}

/** UpdateStats: updates mean, variance and acc_ratio after tau moves
 *                it needs i to do sanity check in parallel code     
 */
void
UpdateStats( int i ) {
#ifdef MPI
    double total[2];            /* will hold mean & vari when pooling stats */
    double tmptotal[2];         /* temporary array for above */

    long nodesuccess[2];        /* array for success and summed moves */
    long tmpsuccess[2];         /* temporary array for above */


    /* parallel code: pool statistics from all nodes */

    total[0] = mean;
    total[1] = vari;

    nodesuccess[0] = ( long ) success;
    nodesuccess[1] = ( long ) i;        /* sum of i's must equal total tau moves */

    tmptotal[0] = total[0];     /* the tmp arrays are used to send messages */
    tmptotal[1] = total[1];

    tmpsuccess[0] = nodesuccess[0];
    tmpsuccess[1] = nodesuccess[1];

    /* stats from all nodes are summed up here */

    MPI_Allreduce( tmptotal, total, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( tmpsuccess, nodesuccess, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

    total[0] /= Tau;            /* need to divide mean and vari by Tau steps */
    total[1] /= Tau;

    success = ( int ) nodesuccess[0];   /* success is now global success! */

    /* sanity check: have we done the correct number of initial moves? */
    if( nodesuccess[1] - state->tune.tau )
        error( "Loop: total moves was %d ?!\n", nodesuccess[1] );

    mean = total[0];            /* mean and variance are now summed over all nodes */
    vari = total[1];

    /* local stats are updated below */

    if( tuning && nnodes > 1 ) {
        l_mean /= ( double ) proc_tau;
        l_vari /= ( double ) proc_tau;
        l_acc_ratio = ( ( double ) l_success ) / ( ( double ) proc_tau );
    }
#else

    mean /= Tau;                /* collect some statistics */
    vari /= Tau;

#endif

    acc_ratio = ( ( double ) success ) / Tau;   /* update acceptance ratio */

}

/**  UpdateParameter: update parameters A, B, D and E and the estimators 
 *                    for mean and standard deviation for the current S    
 */
void
UpdateParameter( void ) {
    register double d;          /* d is used to store intermediate results below */


    /* this part of the code updates the estimator for the mean */

    d = 1.0 / mean;

    /* first: multiply all intermediate vars by weights */

    usyy *= w_a;
    usxy *= w_a;
    usy *= w_a;
    usx *= w_a;
    usxx *= w_a;
    usum *= w_a;

    /* then: update all intermediate vars */

    usyy += d * d;
    usxy += S * d;
    usy += d;
    usx += S;
    usxx += S * S;
    usum += 1.0;

    /* ... and use intermediate vars to update A and B ... */

    A = ( usum * usxy - usx * usy ) / ( usum * usxx - usx * usx );
    B = ( usy - A * usx ) / usum;

    /* ... which are then used to update the estimator for the mean */

    estimate_mean = 1.0 / ( A * S + B );


    /* this part of the code updates the estimator for the standard deviation */

    if( vari > 0.0 ) {

        d = 1.0 / sqrt( vari );

        /* first: multiply all intermediate vars by weights */

        vsyy *= w_b;
        vsxy *= w_b;
        vsy *= w_b;
        vsx *= w_b;
        vsxx *= w_b;
        vsum *= w_b;

        /* then: update all intermediate vars */

        vsyy += d * d;
        vsxy += S * d;
        vsy += d;
        vsx += S;
        vsxx += S * S;
        vsum += 1.0;

        /* ... and use intermediate vars to update D and E ... */

        D = ( vsum * vsxy - vsx * vsy ) / ( vsum * vsxx - vsx * vsx );
        E = ( vsy - D * vsx ) / vsum;
    }

    /* ... which are then used to update the estimator for the std dev */

    estimate_sd = 1.0 / ( D * S + E );

    /* alpha corresponds to the third term in the main Lam schedule formula,   *
     * which is a measure of how efficiently the space state is samples; this  *
     * term is at a maximum for acc_ratio = 0.44, see Lam & Delosme, 1988b, p1 */

    d = ( 1.0 - acc_ratio ) / ( 2.0 - acc_ratio );
    alpha = 4.0 * acc_ratio * d * d;

}



#ifdef MPI

/**  UpdateLParameter: update local parameters l_A_{l|u}, l_B_{l|u}, l_D   *
 *                     and l_E and the local estimators for mean and stan- *
 *                     dard deviation for both upper and lower bounds of M *
 *                     when tuning                                         *
 */
void
UpdateLParameter( void ) {
    register double d;          /* d is used to store intermediate results below */

    /* update the estimator for the local mean */

    d = 1.0 / l_mean;

    /* first: multiply all intermediate vars by weights */
    /* lower bound variables: */

    l_usyy_l *= l_w_a_l;
    l_usxy_l *= l_w_a_l;
    l_usy_l *= l_w_a_l;
    l_usx_l *= l_w_a_l;
    l_usxx_l *= l_w_a_l;
    l_usum_l *= l_w_a_l;

    /* upper bound variables: */

    l_usyy_u *= l_w_a_u;
    l_usxy_u *= l_w_a_u;
    l_usy_u *= l_w_a_u;
    l_usx_u *= l_w_a_u;
    l_usxx_u *= l_w_a_u;
    l_usum_u *= l_w_a_u;

    /* then: update all intermediate vars */
    /* lower bound variables: */

    l_usyy_l += d * d;
    l_usxy_l += S * d;
    l_usy_l += d;
    l_usx_l += S;
    l_usxx_l += S * S;
    l_usum_l += 1.0;

    /* upper bound variables: */

    l_usyy_u += d * d;
    l_usxy_u += S * d;
    l_usy_u += d;
    l_usx_u += S;
    l_usxx_u += S * S;
    l_usum_u += 1.0;

    /* ... and use intermediate vars to update l_A_{l|u} and l_B_{l|u} ... */

    /* lower bound variables: */

    l_A_l = ( l_usum_l * l_usxy_l - l_usx_l * l_usy_l ) / ( l_usum_l * l_usxx_l - l_usx_l * l_usx_l );
    l_B_l = ( l_usy_l - l_A_l * l_usx_l ) / l_usum_l;

    /* upper bound variables: */

    l_A_u = ( l_usum_u * l_usxy_u - l_usx_u * l_usy_u ) / ( l_usum_u * l_usxx_u - l_usx_u * l_usx_u );
    l_B_u = ( l_usy_u - l_A_u * l_usx_u ) / l_usum_u;

    /* ... which are then used to update the local estimators for the mean */

    l_estimate_mean_l = 1.0 / ( l_A_l * S + l_B_l );
    l_estimate_mean_u = 1.0 / ( l_A_u * S + l_B_u );

    /* update the local estimator for the standard deviation */

    if( l_vari > 0.0 ) {

        d = 1.0 / sqrt( l_vari );

        /* first: multiply all intermediate vars by weights */

        l_vsyy *= l_w_b;
        l_vsxy *= l_w_b;
        l_vsy *= l_w_b;
        l_vsx *= l_w_b;
        l_vsxx *= l_w_b;
        l_vsum *= l_w_b;

        /* then: update all intermediate vars */

        l_vsyy += d * d;
        l_vsxy += S * d;
        l_vsy += d;
        l_vsx += S;
        l_vsxx += S * S;
        l_vsum += 1.0;

        /* ... and use intermediate vars to update l_D and l_E ... */

        l_D = ( l_vsum * l_vsxy - l_vsx * l_vsy ) / ( l_vsum * l_vsxx - l_vsx * l_vsx );
        l_E = ( l_vsy - l_D * l_vsx ) / l_vsum;
    }

    /* ... which are then used to update the local estimator for the std dev */

    l_estimate_sd = 1.0 / ( l_D * S + l_E );

    /* alpha corresponds to the third term in the main Lam schedule formula,   *
     * which is a measure of how efficiently the space state is samples; this  *
     * term is at a maximum for acc_ratio = 0.44, see Lam & Delosme, 1988b, p1 */

    d = ( 1.0 - l_acc_ratio ) / ( 2.0 - l_acc_ratio );
    l_alpha = 4.0 * l_acc_ratio * d * d;

}
#endif

/**  Frozen: returns TRUE if frozen, FALSE otherwise; 'frozen' is defined  
 *           by 'freeze_count', 'stop_flag' and 'criterion' (see also sa.h 
 *           for a more extensive comment on this)                         
 */
int
Frozen( void ) {
    double delta;
    delta = 0.0;

    if( stop_flag == proportional_freeze )
        delta = ( mean - old_mean ) / mean;
    else if( stop_flag == absolute_freeze )
        delta = mean - old_mean;
    else if( stop_flag == absolute_energy )
        delta = mean;

    if( delta <= 0. )
        delta = -delta;

    if( delta <= state->tune.criterion ) {
        counter++;
        if( counter > 3 ) {
            printf( "DELTA %lg < CRITERION %lg %d\n", delta, state->tune.criterion, counter );
        }
    } else {
        counter = 0;
        old_mean = mean;
    }
    return ( counter >= state->tune.freeze_count );
}



#ifdef MPI
/*** MIXING ****************************************************************/

/**  DoMix: does the mixing; sends move state and local Lam stats to the 
 *          dance partner(s)                                               
 */
void
DoMix( void ) {
    int i;                      /* loop counter */

    /* variables needed for evaluating the dance partners; note that the dance *
     * partner array is static to lsa.c, since it's also needed by tuning code */

    double prob;                /* probability of choosing a node's energy upon mixing */
    double norm;                /* sum used to normalize probabilities */
    double theirprob = 0.0;     /* probability of the dance partner */
    double psum;                /* sum of probabilities for a certain dance partner */
    double *node_prob;          /* all probabilities for all nodes */

    /* buffers for sending/receiving Lam stats at mix */

    double *sendbuf;
    double *recvbuf;

    /* buffers for send/recv move state from/to move(s).c (incl. sizes) */

    int lsize;
    int dsize;

    long *recv_longbuf;
    double *recv_doublebuf;

    long *send_longbuf;
    double *send_doublebuf;

    /* MPI status and handle (request) arrays */

    MPI_Status *status;         /* status array of the MPI_Waitall function */
    MPI_Request *request;       /* handle array for receiving messages */

    recvbuf = NULL;
    recv_longbuf = NULL;
    recv_doublebuf = NULL;

    /* allocate probability array and arrays for MPI_Waitall below */

    node_prob = ( double * ) calloc( nnodes, sizeof( double ) );

    request = ( MPI_Request * ) calloc( 3, sizeof( MPI_Request ) );
    status = ( MPI_Status * ) calloc( 3, sizeof( MPI_Status ) );

    /* initialize probability & dance partner arrays */

    for( i = 0; i < nnodes; i++ ) {
        dance_partner[i] = 0;   /* static to lsa.c since needed for tuning */
        node_prob[i] = 0.;
    }

    /* update mix counter (used by tuning code only) */

    if( tuning )
        count_mix++;            /* counts the number of mixings we have done */

    /* calculate probabilities for accepting local state upon mixing; we sub-  *
     * tract estimate_mean (which is the same in all processes from the total  *
     * energy to avoid overflows; see also Chu (2001, p.41) for details        */

    prob = exp( ( estimate_mean - energy ) * S );

    if( prob < DBL_MIN )        /* set some limits to prob to avoid under-/overflows */
        prob = DBL_MIN;
    else if( prob >= HUGE_VAL )
        prob = DBL_MAX / nnodes;

    /* sum up probabilities from all nodes for normalization, then normalize */

    MPI_Allreduce( &prob, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    prob /= norm;

    /* gather local probabilities from all nodes */

    MPI_Allgather( &prob, 1, MPI_DOUBLE, node_prob, 1, MPI_DOUBLE, MPI_COMM_WORLD );

    /* theirprob determines the dance partner we choose */

    if( nnodes > 1 )
        theirprob = RandomReal(  );
    else if( nnodes == 1 )
        theirprob = 0.0;        /* guarantee same Markov proc. as serial case */
    else
        error( "DoMix: you can't compute on %d nodes!", nnodes );

    /* sum up probabilities to determine dance partner */

    psum = 0.;
    for( i = 0; i < nnodes; i++ ) {
        psum += node_prob[i];
        if( psum > theirprob )
            break;
    }

    /* pool dance partners */
    MPI_Allgather( &i, 1, MPI_INT, dance_partner, 1, MPI_INT, MPI_COMM_WORLD );
    /* get move state from move(s).c and collect local Lam stats for sending */

    MakeStateMsg( &send_longbuf, &lsize, &send_doublebuf, &dsize );
    MakeLamMsg( &sendbuf );

    /* if I'm not dancing with myself: receive new state and Lam stats         */

    if( dance_partner[myid] != myid ) {

        /* allocate receive buffers for message; receive buffers need to be of the *
         * same length as the sending buffers                                      */

        recv_longbuf = ( long * ) calloc( lsize, sizeof( long ) );
        recv_doublebuf = ( double * ) calloc( dsize, sizeof( double ) );
        if( tuning && nnodes > 1 )
            recvbuf = ( double * ) calloc( LSTAT_LENGTH_TUNE, sizeof( double ) );
        else
            recvbuf = ( double * ) calloc( LSTAT_LENGTH, sizeof( double ) );

        /* MPI_Irecv is a non-blocking receive which allows a process that is      *
         * waiting to receive to send out its message while waiting i.e. we post   *
         * the receive now, check for reception after send (using MPI_Waitall)     */

        if( dance_partner[myid] > nnodes ) {
            printf( "%d WARNING: my dance_partner is %d > nnodes(%d)\n", myid, dance_partner[myid], nnodes );
        }

        MPI_Irecv( recv_doublebuf, dsize, MPI_DOUBLE, dance_partner[myid], dance_partner[myid], MPI_COMM_WORLD, &request[0] );
        MPI_Irecv( recv_longbuf, lsize, MPI_LONG, dance_partner[myid], dance_partner[myid], MPI_COMM_WORLD, &request[1] );

        /* the length of the Lam statistics message depends on if we're tuning or  *
         * not (local stats are only sent when tuning)                             */

        if( tuning && nnodes > 1 )
            MPI_Irecv( recvbuf, LSTAT_LENGTH_TUNE, MPI_DOUBLE, dance_partner[myid], dance_partner[myid], MPI_COMM_WORLD, &request[2] );
        else
            MPI_Irecv( recvbuf, LSTAT_LENGTH, MPI_DOUBLE, dance_partner[myid], dance_partner[myid], MPI_COMM_WORLD, &request[2] );

    }

    /* send messages to dance partners, if requested */

    for( i = 0; i < nnodes; i++ )
        if( ( dance_partner[i] == myid ) && ( i != myid ) ) {
            MPI_Send( send_doublebuf, dsize, MPI_DOUBLE, i, myid, MPI_COMM_WORLD );
            MPI_Send( send_longbuf, lsize, MPI_LONG, i, myid, MPI_COMM_WORLD );
            if( tuning && nnodes > 1 )
                MPI_Send( sendbuf, LSTAT_LENGTH_TUNE, MPI_DOUBLE, i, myid, MPI_COMM_WORLD );
            else
                MPI_Send( sendbuf, LSTAT_LENGTH, MPI_DOUBLE, i, myid, MPI_COMM_WORLD );
        }

    /* if I'm not dancing with myself, we need a new state; MPI_Waitall will   *
     * collect the three messages we need to receive; then we install the move *
     * state in move(s).c and the Lam stats in lsa.c                           */

    if( dance_partner[myid] != myid ) {

        MPI_Waitall( 3, request, status );

        AcceptStateMsg( recv_longbuf, recv_doublebuf );
        AcceptLamMsg( recvbuf );

    }

    /* clean up message buffers and MPI arrays ... */

    free( send_longbuf );
    free( send_doublebuf );
    free( sendbuf );

    free( request );
    free( status );

    /* ... and the probability array */

    free( node_prob );

}

/**  MakeLamMsg: packages local Lam stats into send buffer */
void
MakeLamMsg( double **sendbuf ) {
    if( tuning && nnodes > 1 )
        *sendbuf = ( double * ) calloc( LSTAT_LENGTH_TUNE, sizeof( double ) );
    else
        *sendbuf = ( double * ) calloc( LSTAT_LENGTH, sizeof( double ) );

    ( *sendbuf )[0] = energy;

    if( tuning && nnodes > 1 ) {

        ( *sendbuf )[1] = l_estimate_mean_l;
        ( *sendbuf )[2] = l_estimate_mean_u;
        ( *sendbuf )[3] = l_estimate_sd;

        ( *sendbuf )[4] = l_usyy_l;
        ( *sendbuf )[5] = l_usxy_l;
        ( *sendbuf )[6] = l_usy_l;
        ( *sendbuf )[7] = l_usx_l;
        ( *sendbuf )[8] = l_usxx_l;
        ( *sendbuf )[9] = l_usum_l;

        ( *sendbuf )[10] = l_usyy_u;
        ( *sendbuf )[11] = l_usxy_u;
        ( *sendbuf )[12] = l_usy_u;
        ( *sendbuf )[13] = l_usx_u;
        ( *sendbuf )[14] = l_usxx_u;
        ( *sendbuf )[15] = l_usum_u;

        ( *sendbuf )[16] = l_A_l;
        ( *sendbuf )[17] = l_B_l;

        ( *sendbuf )[18] = l_A_u;
        ( *sendbuf )[19] = l_B_u;

        ( *sendbuf )[20] = l_vsyy;
        ( *sendbuf )[21] = l_vsxy;
        ( *sendbuf )[22] = l_vsy;
        ( *sendbuf )[23] = l_vsx;
        ( *sendbuf )[24] = l_vsxx;
        ( *sendbuf )[25] = l_vsum;

        ( *sendbuf )[26] = l_D;
        ( *sendbuf )[27] = l_E;

    }
}

/**  AcceptLamMsg: receives new energy and Lam stats upon mixing */
void
AcceptLamMsg( double *recvbuf ) {
    energy = recvbuf[0];

    if( tuning && nnodes > 1 ) {

        l_estimate_mean_l = recvbuf[1];
        l_estimate_mean_u = recvbuf[2];
        l_estimate_sd = recvbuf[3];

        l_usyy_l = recvbuf[4];
        l_usxy_l = recvbuf[5];
        l_usy_l = recvbuf[6];
        l_usx_l = recvbuf[7];
        l_usxx_l = recvbuf[8];
        l_usum_l = recvbuf[9];

        l_usyy_u = recvbuf[10];
        l_usxy_u = recvbuf[11];
        l_usy_u = recvbuf[12];
        l_usx_u = recvbuf[13];
        l_usxx_u = recvbuf[14];
        l_usum_u = recvbuf[15];

        l_A_l = recvbuf[16];
        l_B_l = recvbuf[17];

        l_A_u = recvbuf[18];
        l_B_u = recvbuf[19];

        l_vsyy = recvbuf[20];
        l_vsxy = recvbuf[21];
        l_vsy = recvbuf[22];
        l_vsx = recvbuf[23];
        l_vsxx = recvbuf[24];
        l_vsum = recvbuf[25];

        l_D = recvbuf[26];
        l_E = recvbuf[27];

    }

    free( recvbuf );

}





/*** TUNING CODE ***********************************************************/

/**  DoTuning: calculates the cross-correlation (for lower bound) and the  
 *             variance of local means (for upper bound) for a sub_tune_   
 *             interval; the results are added to arrays, which span the   
 *             whole tuning_interval; when these values get written to the 
 *             'bound' files, they will be divided by the number of mixes  
 *             we have done already to average them out (see also King-Wai 
 *             Chu's thesis, p. 63)                                        
 */
void
DoTuning( void ) {
    int i, j, k, l;             /* loop counters */

    /* variables for calculating cross-correlation for lower bound */

    double var;                 /* variance over sample interval */
    double correl = 0.;         /* temp variable for calculating cross-correlation */
    int ncorrel = 0;            /* counter for cross-correlations we've calculated */

    /* variables for calculating variance of local means for upper bound */

    double avg_l_mean = 0.;     /* the average local mean across all processors */
    double var_l_mean = 0.;     /* temp variable for the var of the local means */



    /* For the lower bound of M: ***********************************************
     * here we calculate a intermediate product for the calculation of the     *
     * cross-correlation between processors; what we end up having in each     *
     * element of dev is the local dev / sqrt(var); we can then multi-         *
     * ply this for two different processors to obtain the cross correlation   *
     * below (according to King-Wai's thesis, pp. 55-56)                       *
     ***************************************************************************/

    for( j = 0; j < sub_tune_interval; j++ ) {
        var = 0.;
        for( k = 0; k < sample_size; k++ )
            var += dev[j * sample_size + k] * dev[j * sample_size + k];
        var /= sample_size;
        for( k = 0; k < sample_size; k++ )
            if( var > 0.0 )
                dev[j * sample_size + k] /= sqrt( var );
            else
                dev[j * sample_size + k] /= DBL_MIN;    /* not quite right, but close */
    }

    /* distribute devs and means for this sub_tune_interval among all          *
     * processors                                                              *
    */

    MPI_Allgather( dev, sample_size * sub_tune_interval, MPI_DOUBLE, tot_dev, sample_size * sub_tune_interval, MPI_DOUBLE, MPI_COMM_WORLD );

    MPI_Allgather( means, sample_size * sub_tune_interval, MPI_DOUBLE, tot_means, sample_size * sub_tune_interval, MPI_DOUBLE, MPI_COMM_WORLD );

    /* the following stuff is done for each sample over the sub_tune_inteval */

    for( i = 0; i < sub_tune_interval; i++ ) {

        /* this loop takes all sample_size devs from all processors for the        *
         * current sample within sub_tune_interval and rearranges them into the    *
         * coll arrays, which we use below to calculate cross correlation and      *
         * variances of the local means                                            *
         */

        for( j = 0; j < nnodes; j++ )
            for( k = 0; k < sample_size; k++ ) {
                coll_dev[j * sample_size + k] = tot_dev[j * sub_tune_interval * sample_size + i * sample_size + k];
                coll_means[j * sample_size + k] = tot_means[j * sub_tune_interval * sample_size + i * sample_size + k];
            }


        /* For the lower bound of M: ***********************************************
         * calculate cross-correlation for those processors that had the same dan- *
         * ce partner at the last mix; most of the work has already been done for  *
         * us in a loop above, where we calculate the intermediate product that    *
         * corresponds to dev / sqrt(var); here, we only need to multiply such     *
         * intermediate products for processors that shared the same state after   *
         * the previous mix; then save the results in the cross_correl array       *
         */

        correl = 0.;
        ncorrel = 0;
        for( j = 1; j < nnodes; j++ )   /* loop through all combinations of */
            for( k = 0; k < j; k++ )    /* processors */
                if( dance_partner[j] == dance_partner[k] )
                    for( l = 0; l < sample_size; l++ ) {
                        correl += coll_dev[j * sample_size + l] * coll_dev[k * sample_size + l];
                        ncorrel++;
                    }

        if( ncorrel > 0 )       /* only divide if anything was calculated */
            correl /= ( double ) ncorrel;

        cross_correl[count_tune * sub_tune_interval + i] += correl;

        /* For the upper bound of M: ***********************************************
         * here we calculate the variance of the local mean for the current sample *
         * interval; avg_l_mean is the average local mean from all processors over *
         * the covar_sample; it is then used to calculate the variance of the      *
         * local means (var_l_mean) over the sample interval, which is then stored *
         * in the var_means array                                                  *
         */

        avg_l_mean = 0.;
        var_l_mean = 0.;

        for( j = 0; j < covar_sample; j++ )
            avg_l_mean += coll_means[j];
        avg_l_mean /= ( double ) covar_sample;

        for( j = 0; j < covar_sample; j++ ) {
            if( avg_l_mean == 0.0 )
                coll_means[j] = ( coll_means[j] - DBL_MIN ) / DBL_MIN;
            else
                coll_means[j] = ( coll_means[j] - avg_l_mean ) / avg_l_mean;
            var_l_mean += coll_means[j] * coll_means[j];
        }
        var_l_mean /= ( double ) covar_sample;

        var_means[count_tune * sub_tune_interval + i] += var_l_mean;

    }

    /* we're done with this sub_tune_interval: on to the next one */

    count_tune++;

    /* if we're at a tune_interval: reset the count_tune counter */

    if( count_tune == write_tune_stat )
        count_tune = 0;

    /* restart the moves counter */

    moves_tune = 0;

}

/**  WriteTuning: every tune_interval (which basically corresponds to the 
 *                mix_interval), we write the cross-correlations and vari- 
 *                ances of the local mean to various files; the way this   
 *                is done is a little bizarre: the two .lb and .ub files   
 *                are averaged over all previous tune_intervals and then   
 *                overwritten every time, the .mb file on the other hand   
 *                contains the history of the upper and lower bounds over  
 *                past tune_intervals; here are the details in the words   
 *                of the Great Dr. Chu:                                    
 *                When tuning finishes, there is a file called <infile>.lb 
 *                (Fig. 5.1 in thesis) and a file called <infile>.ub (Fig. 
 *                5.2). The history of the lower and upper bound values    
 *                during tuning are in <infile>.mb (Fig. 5.3). The average 
 *                of the last pair of lower and upper bound values is the  
 *                M from tuning.                                           
 */
void
WriteTuning( void ) {
    int i, j;                   /* loop counters */

    /* file pointers */

    FILE *lbptr;                /* pointer for lower_bound file */
    FILE *ubptr;                /* pointer for upper_bound file */
    FILE *mbptr;                /* pointer for mix_bound file */

    lbptr = NULL;
    ubptr = NULL;
    mbptr = NULL;

    /* variables for lower bound */

    double long_avg;            /* average cross_correl over whole mix_interval */

    /* variables for upper bound */

    double min;                 /* these are for storing the minimum and */
    double max;                 /* maximum variances of local means over a mix_interval */
    int group_size = GROUP_SIZE;/* group for evaluating upper bound */
    int midpoint = 0;           /* midpoint of a group */
    double sum;                 /* summed variances of a group */
    double avg = 0.0;           /* average of a group */
    
    /* both upper and lower bound information goes into the mix_bound file     *
     * that's why we keep it open during both writing lower and upper bounds   */

    if( myid == 0 ) {
        mbptr = fopen( files.mbfile, "a" );
        if( !mbptr )
            file_error( "WriteTuning" );

        /* first write cross-correlation for lower bounds */

        lbptr = fopen( files.lbfile, "w" );
        if( !lbptr )
            file_error( "WriteTuning" );

        fprintf( lbptr, "# tau no.    cross_corr\n\n" );

        for( i = 0; i < tune_interval; i++ )
            fprintf( lbptr, "   %6d   %11.8f\n", tau_count[i], cross_correl[i] / ( double ) count_mix );

        fclose( lbptr );


        /* lower bound: write sample-count for first cross-correlation which is    *
         * below the average cross-correlation over the whole tune_interval to     *
         * the mix-bound file (see King-Wai's thesis, p. 58)                       */

        long_avg = 0.;
        for( i = 0; i < tune_interval; i++ )
            long_avg += cross_correl[i];
        long_avg /= ( double ) tune_interval;

        for( i = 0; i < tune_interval; i++ )
            if( cross_correl[i] <= long_avg ) {
                fprintf( mbptr, "       %4d   %6d   %11.8f", count_mix, tau_count[i], cross_correl[i] / ( double ) count_mix );
                break;
            }

        /* write variance of local means for upper bounds */

        ubptr = fopen( files.ubfile, "w" );
        if( !ubptr )
            file_error( "WriteTuning" );

        fprintf( ubptr, "# tau no.    var_l_mean\n\n" );
    }

    /* here we evaluate min and max value of var_means over the whole interval */

    min = var_means[0] / ( double ) count_mix;
    max = -DBL_MIN;
    for( i = 0; i < tune_interval; i++ ) {
        if( ( var_means[i] / ( double ) count_mix ) > max )
            max = var_means[i] / ( double ) count_mix;
        if( myid == 0 )
            fprintf( ubptr, "   %6d   %11.8f\n", tau_count[i], var_means[i] / ( double ) count_mix );
    }

    if( myid == 0 )
        fclose( ubptr );

    /* upper bounds: evaluate group of variances that differs from the first   *
     * group of variances by more than 7% (see King-Wai's thesis, p. 61)       */

    group_size /= covar_index;  /* set group size for upper bound */
    if( group_size < 1 )
        group_size = 1;

    for( i = 0; i < tune_interval - group_size; i++ ) {
        midpoint = i + group_size / 2;
        sum = 0.;
        for( j = 0; j < group_size; j++ )
            sum += var_means[i + j];
        avg = sum / ( ( double ) group_size * ( double ) count_mix );

        if( fabs( ( avg - min ) / ( max - min ) ) >= 0.07 )
            break;
    }

    midpoints[count_mix - 1] = midpoint;

    if( myid == 0 ) {
        fprintf( mbptr, "   %6d   %11.8f\n", midpoint, avg );
        fclose( mbptr );
    }

    /* reset the sample counter */

    count_sample = 0;

}

/**  StopTuning: is to tuning runs what Frozen() is to a normal annealing 
 *               run: it basically checks if the tuning stop criterion    
 *               applies and returns true if that's the case              
 */
int
StopTuning( void ) {
    int i;                      /* loop counter */

    FILE *logptr;               /* pointer to .log file */

    double tol_tune = STOP_TUNE_CRIT;   /* tuning stop criterion */
    double avg;                 /* temp variable for the average */

    /* calculate average upper bound for the last 'stop_tune_count' mixes */

    avg = 0.;
    for( i = count_mix - stop_tune_count; i < count_mix; i++ )
        avg += ( double ) midpoints[i];
    avg /= ( double ) stop_tune_count;

    /* check if any upper bound is further than tol_tune away from that avg */

    if( avg == 0 )
        error( "StopTuning: average midpoint was zero!!?!\n" );
    else
        for( i = count_mix - stop_tune_count; i < count_mix; i++ )
            if( ( fabs( ( double ) midpoints[i] - avg ) / avg ) >= tol_tune )
                return 0;       /* if yes, return 'false' */

    /* if we're done: write a message into the log file and quit */

    if( myid == 0 ) {

        logptr = fopen( files.logfile, "a" );   /* write comment to global .log file */
        if( !logptr )
            file_error( "StopTuning" );

        fprintf( logptr, "Tuning stops before the end of an annealing run.\n" );
        fprintf( logptr, "Therefore, the score and iterations will not be\n" );
        fprintf( logptr, "the true final score and iterations.\n" );

        fclose( logptr );
    }

    return 1;

}
#endif





/*** FUNCTIONS REQUIRED FOR EQUILIBRATION RUNS ONLY ************************/

/**  FixTLoop: loops (making moves, updating stats etc.) keeping the temp- 
 *             erature fixed for an equilibration run                      
 */
void
FixTLoop( void ) {
    int i;                      /* loop counters */
#ifndef MPI
    int j;
#endif

    char *varfile;              /* global variance file name */
    FILE *varptr;               /* global variance file pointer */
    varptr = NULL;

#ifdef MPI
    char *lvarfile;             /* local variance file name */
    FILE *lvarptr;              /* local variance file pointer */
#else
    char *acfile;               /* autocorrelation file name */
    FILE *acptr;                /* autocorrelation file pointer */
#endif

    double energy_change;       /* local Delta E */

    /* array(s) to store energies */

    double *energy_storage;     /* array used to store intermediate energies */
#ifdef MPI
    double *tot_energy;         /* array for total energy pooled from all nodes */
#endif

    /* below are vars that are used to calculate statistics at equilibrium */

    /* first some variables to calculate average energy and variances */

    double tmp1 = 0.0;          /* three temporary variables */
    double tmp2 = 0.0;          /* for calculating variances */
    double var_sum = 0.0;

    double instant_avg;         /* used for calculating the inst. average energy */

    /* below some variables used to calculate autocorrelation */

#ifndef MPI                     /* autocorrelation is only calculated in serial */

    double *gamma;              /* array for covariances */
    double *rho;                /* array of autocorrelations */

    /* array of intervals for calculating autocorrelations */
    int h[56] = { 0,
        1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 20, 30, 40, 50, 60, 70, 80, 90,
        100, 200, 300, 400, 500, 600, 700, 800, 900,
        1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
        10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
        100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
        1000000
    };

    int num_of_h = 56;          /* number of elements in h array */
    int ntemp;                  /* temporary var for max. step minus h */

#endif

#ifdef MPI
    /* allocate parallel-specific file names */

    varfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    lvarfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    /* allocate the energy arrays: we need as many elements as there are mix-  *
     * ing events during the sampling loop                                     */

    energy_storage = ( double * ) calloc( equil_param.fix_T_step / ( state->tune.mix_interval * proc_tau ) + 1, sizeof( double ) );
    tot_energy = ( double * ) calloc( equil_param.fix_T_step / ( state->tune.mix_interval * proc_tau ) + 1, sizeof( double ) );

#else

    /* allocate serial-specific file names */

    varfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    acfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    /* allocate the energy array */

    energy_storage = ( double * ) calloc( equil_param.fix_T_step + 1, sizeof( double ) );

    /* allocate autocorrelation arrays */

    gamma = ( double * ) calloc( num_of_h, sizeof( double ) );
    rho = ( double * ) calloc( num_of_h, sizeof( double ) );

#endif

#ifdef MPI

    if( myid == 0 )
        sprintf( varfile, "%s_%d.var", files.outputfile, state->tune.mix_interval );

    sprintf( lvarfile, "%s_%d_%d.lvar", files.outputfile, myid, state->tune.mix_interval );

#else

    sprintf( varfile, "%s.var", files.outputfile );
    sprintf( acfile, "%s.ac", files.outputfile );

#endif

#ifdef MPI

    count_mix = 0;              /* we need to reset this counter, since it's needed below */

#endif

    /* Equilibrate and throw away statistics */

    for( i = 0; i < equil_param.fix_T_skip; i++ ) {

        /* randomize initial state; throw out results; DO NOT PARALLELIZE! */

        energy_change = GenerateMove( &files, &distp, &out );

        /* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
         * that's why we want to prevent overflows here (hence the 'if')           */

        if( energy_change != FORBIDDEN_MOVE )
            exp_arg = -S * energy_change;

        /* MIN_DELTA provides a min. probability with which any move is accepted */

        if( exp_arg <= MIN_DELTA )
            exp_arg = MIN_DELTA;

        /* below, we apply the Metropolis criterion to accept or reject a move */

        if( energy_change == FORBIDDEN_MOVE ) {
            RejectMove(  );
        } else if( ( energy_change <= 0.0 ) || ( exp( exp_arg ) > RandomReal(  ) ) ) {
            energy += energy_change;
            AcceptMove(  );
        } else {
            RejectMove(  );
        }

    }

    /* Equilibrate and gather statistics */

    for( i = 0; i <= equil_param.fix_T_step; i++ ) {

        /* make a move: will either return the energy change or FORBIDDEN_MOVE */

        energy_change = GenerateMove( &files, &distp, &out );

        /* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
         * that's why we want to prevent overflows here (hence the 'if')           */

        if( energy_change != FORBIDDEN_MOVE )
            exp_arg = -S * energy_change;

        /* MIN_DELTA provides a min. probability with which any move is accepted */

        if( exp_arg <= MIN_DELTA )
            exp_arg = MIN_DELTA;

        /* below, we apply the Metropolis criterion to accept or reject a move */

        if( energy_change == FORBIDDEN_MOVE ) {
            RejectMove(  );
        } else if( ( energy_change <= 0.0 ) || ( exp( exp_arg ) > RandomReal(  ) ) ) {
            energy += energy_change;
            AcceptMove(  );
            if( landscape ) {
                WriteLandscape( files.landscapefile, i, energy_change );
            }
        } else {
            RejectMove(  );
        }

#ifdef MPI

        /* parallel code: we only collect energies every mix_interval; we mix at   *
         * exactly the same interval as if we would be doing a normal annealing    *
         * run, therefore we need to multiply by proc_tau below                    */

        if( i % ( state->tune.mix_interval * proc_tau ) == 0 ) {

            energy_storage[count_mix] = energy;
            count_mix++;

            /* here we mix; since temperature does not change, we don't need to mix    *
             * Lam stats but only current energies and move state                      */

            DoMix(  );

        }
#else

        /* serial code: collect energies at every step */

        energy_storage[i] = energy;

#endif

    }                           /* end of loop to gather statistics */



#ifdef MPI
    /* parallel code: pool all energies from all nodes into tot_energy array   *
     * and average them                                                        */

    MPI_Allreduce( energy_storage, tot_energy, ( count_mix - 1 ), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    for( i = 0; i < ( count_mix - 1 ); i++ )
        tot_energy[i] /= nnodes;

    /* ... then calculate statistics and write them to files */

    if( myid == 0 ) {           /* this is the global .var file */
        varptr = fopen( varfile, "w" );
        if( !varptr )
            file_error( "FixTLoop" );
    }

    lvarptr = fopen( lvarfile, "w" );   /* this is the local .lvar file */
    if( !lvarptr )
        file_error( "FixTLoop" );

#else

    /* open files to write into */

    varptr = fopen( varfile, "w" );     /* this is the .var file */
    if( !varptr )
        file_error( "FixTLoop" );

    acptr = fopen( acfile, "w" );       /* this is the .ac file */
    if( !acptr )
        file_error( "FixTLoop" );

#endif

    /* first calculate local average and variances */

    /* calculate overall energy average */

#ifdef MPI
    for( i = 0; i < ( count_mix - 1 ); i++ )
#else
    for( i = 0; i <= equil_param.fix_T_step; i++ )
#endif
        fix_T_avg += energy_storage[i];
    fix_T_avg /= equil_param.fix_T_step;

    /* initialize variables */

    tmp1 = ( energy_storage[0] - fix_T_avg ) * ( energy_storage[0] - fix_T_avg );
    tmp2 = energy_storage[0] - fix_T_avg;
    instant_avg = energy_storage[0];

    /* print captions */
#ifdef MPI
    fprintf( lvarptr, "# nmixes     variance     " );
    fprintf( lvarptr, "inst. avg.   overall avg.\n\n" );
#else
    fprintf( varptr, "# nsteps     variance     " );
    fprintf( varptr, "inst. avg.   overall avg.\n\n" );
#endif

    /* calculate and print local variances here */

#ifdef MPI
    for( i = 1; i < ( count_mix - 1 ); i++ ) {
#else
    for( i = 1; i <= equil_param.fix_T_step; i++ ) {
#endif

        instant_avg += energy_storage[i];

        tmp1 += ( energy_storage[i] - fix_T_avg ) * ( energy_storage[i] - fix_T_avg );
        tmp2 += ( energy_storage[i] - fix_T_avg );
        var_sum += ( tmp1 - tmp2 * tmp2 / ( i + 1 ) ) / i;

        fix_T_var = var_sum / ( i + 1 );

        /* print variances and averages every 1000 steps or when finished */

#ifdef MPI
        if( !( i % 10 ) ) {
            fprintf( lvarptr, "%8d %12.5E   %12.5E   %12.5E\n", i, fix_T_var, instant_avg / ( i + 1 ), fix_T_avg );
            fflush( lvarptr );
        }
#else
        if( !( i % 1000 ) || ( i == equil_param.fix_T_step ) ) {
            fprintf( varptr, "%8d %12.5E   %12.5E   %12.5E\n", i, fix_T_var, instant_avg / ( i + 1 ), fix_T_avg );
            fflush( varptr );
        }
#endif

    }

#ifdef MPI
    /* then calculate global average and variances here */

    /* calculate global overall energy average */

    var_sum = 0.0;
    for( i = 0; i < ( count_mix - 1 ); i++ )
        pfix_T_avg += tot_energy[i];
    pfix_T_avg /= ( count_mix - 1 );

    /* initialize variables */

    tmp1 = ( tot_energy[0] - pfix_T_avg ) * ( tot_energy[0] - pfix_T_avg );
    tmp2 = tot_energy[0] - pfix_T_avg;
    instant_avg = tot_energy[0];

    /* print captions */

    if( myid == 0 ) {
        fprintf( varptr, "# nmixes     variance     " );
        fprintf( varptr, "inst. avg.   overall avg.\n\n" );
    }

    /* calculate and print global variances here */

    for( i = 1; i < ( count_mix - 1 ); i++ ) {

        instant_avg += tot_energy[i];

        tmp1 += ( tot_energy[i] - pfix_T_avg ) * ( tot_energy[i] - pfix_T_avg );
        tmp2 += ( tot_energy[i] - pfix_T_avg );
        var_sum += ( tmp1 - tmp2 * tmp2 / ( i + 1 ) ) / i;

        pfix_T_var = var_sum / ( i + 1 );

        if( myid == 0 ) {
            if( !( i % 10 ) && ( myid == 0 ) ) {
                fprintf( varptr, "%8d %12.5E   %12.5E   %12.5E\n", i, pfix_T_var, instant_avg / ( i + 1 ), pfix_T_avg );
                fflush( varptr );
            }
        }
    }

#else

    /* calculate auto-correlation below: this is important for evaluating the  */
    /* appropriate number of initial steps for a problem (see Kingwai's thesis */
    /* section 2.5.3 (pp. 23 - 24)                                             */
    /* gamma is the covariance for a certain interval (taken from the h array) */
    /* rho is the auto-correlation for the same interval                       */

    fprintf( acptr, "#      h           rho          gamma\n\n" );

    for( i = 0; i < num_of_h; i++ ) {   /* calculate autocorrelation for all h's */

        ntemp = equil_param.fix_T_step - h[i];

        /* note: if h[i] is bigger than fix_T_step, gamma and rho will simply be 0 */

        gamma[i] = 0.0;
        for( j = 0; j <= ntemp; j++ )
            gamma[i] += ( energy_storage[j + h[i]] - fix_T_avg ) * ( energy_storage[j] - fix_T_avg );
        gamma[i] /= equil_param.fix_T_step;

        rho[i] = gamma[i] / gamma[0];
        fprintf( acptr, "%8d   % 10.8f   %12.5E\n", h[i], rho[i], gamma[i] );

    }

#endif

    /* clean up and go home ... */

#ifdef MPI
    if( myid == 0 )
#endif
        fclose( varptr );
    free( varfile );
#ifdef MPI
    fclose( lvarptr );
    free( lvarfile );
#else
    fclose( acptr );
    free( energy_storage );
    free( acfile );
    free( gamma );
    free( rho );
#endif

    FinalMove(  );
    return;

}

/** GetEquil: returns the results of an equilibration run */
void
GetEquil( double *equil_var ) {
#ifdef MPI
    equil_var[0] = pfix_T_var;
    equil_var[1] = pfix_T_avg;
#else
    equil_var[0] = fix_T_var;
    equil_var[1] = fix_T_avg;
#endif
}

/** SetEquilibrate: simply makes the equil_param struct static to lsa.c */
void
SetEquilibrate( ChuParam ep ) {
    equil_param = ep;
}




/*** FUNCTIONS USED TO COMMUNICATE WITH OTHER FILES  ***********************/


 /** InitLandscape: sets flag for printing landscape output and acceptance 
 *                 landscape and initializes the landscape file names      
 *                 called from xxx_sa.c to make filenames static and       
 *                 set landscape flag                                      
 */
void
InitLandscape( int value, char *file ) {

    const char *suffix = ".landscape";  /* landscape file in equilibrate */

    /* sets the landscape and acceptance landscape file names static to lsa.c */

    files.landscapefile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    files.landscapefile = strcpy( files.landscapefile, file );
    files.landscapefile = strcat( files.landscapefile, suffix );

    landscape = value;


}

/**  GetLamstats: returns Lam statistics in an array of doubles; used to   
 *                store Lam statistics in a state file                     
 */
double *
GetLamstats( void ) {
    double *stats;

    stats = ( double * ) calloc( 31, sizeof( double ) );

    stats[0] = ( double ) counter;

    stats[1] = old_mean;
    stats[2] = energy;

    stats[3] = mean;
    stats[4] = vari;

    stats[5] = estimate_mean;
    stats[6] = estimate_sd;

    stats[7] = S;
    stats[8] = dS;
    stats[9] = S_0;

    stats[10] = alpha;
    stats[11] = acc_ratio;

    stats[12] = w_b;
    stats[13] = vsyy;
    stats[14] = vsxy;
    stats[15] = vsxx;
    stats[16] = vsx;
    stats[17] = vsy;
    stats[18] = vsum;
    stats[19] = D;
    stats[20] = E;

    stats[21] = w_a;
    stats[22] = usyy;
    stats[23] = usxy;
    stats[24] = usxx;
    stats[25] = usx;
    stats[26] = usy;
    stats[27] = usum;
    stats[28] = A;
    stats[29] = B;

    stats[30] = ( double ) count_tau;

    return ( stats );
}

/**  GetTimes: returns a two-element array with the current wallclock and  
 *             user time to be saved in the state file; for parallel code  
 *             we average the times for all processes                      
 */
double *
GetTimes( void ) {
    double *delta;
    double clk_tck = ( double ) sysconf( _SC_CLK_TCK ); // system-specific clock tick duration
#ifdef MPI
    double temp;
#endif

    delta = ( double * ) calloc( 2, sizeof( double ) );
    // measure user time
    times( cpu_finish );
    // then wallclock time
#ifdef MPI
    finish = MPI_Wtime(  );

    temp = finish - start;
    MPI_Allreduce( &temp, &delta[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    delta[0] /= nnodes;

    temp = ( double ) ( cpu_finish->tms_utime - cpu_start->tms_utime ) / clk_tck;
    MPI_Allreduce( &temp, &delta[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    delta[1] /= nnodes;

#else
    finish = time( NULL );

    delta[0] = finish - start;
    delta[1] = ( cpu_finish->tms_utime - cpu_start->tms_utime ) / clk_tck;
#endif

    return delta;
}

/**  SetOutname: sets the output filename in lsa.c; this is necessary to   
 *               have the diverse .log and .ac and .mb etc files have the  
 *               name of the output file if -w is chosen                   
 */
void
SetOutname( char *outname ) {
    files.outputfile = outname;
}





/*** FUNCTIONS TO RESTORE THINGS IN LSA.C AFTER A RESTART ******************/

/**  RestoreLamstats: restores static Lam statistics in lsa.c from an      
 *                    array of doubles; used to restore runs from a state  
 *                    file.                                                
 */
void
RestoreLamstats( double *stats ) {
    counter = ( int ) rint( stats[0] );

    old_mean = stats[1];
    energy = stats[2];
    mean = stats[3];
    vari = stats[4];

    estimate_mean = stats[5];
    estimate_sd = stats[6];

    S = stats[7];
    dS = stats[8];
    S_0 = stats[9];
    alpha = stats[10];
    acc_ratio = stats[11];

    w_b = stats[12];
    vsyy = stats[13];
    vsxy = stats[14];
    vsxx = stats[15];
    vsx = stats[16];
    vsy = stats[17];
    vsum = stats[18];
    D = stats[19];
    E = stats[20];

    w_a = stats[21];
    usyy = stats[22];
    usxy = stats[23];
    usxx = stats[24];
    usx = stats[25];
    usy = stats[26];
    usum = stats[27];
    A = stats[28];
    B = stats[29];

    count_tau = ( long ) rint( stats[30] );

    free( stats );
}

/**  RestoreTimes: restores the wallclock and user times if -t is used */
void
RestoreTimes( double *delta ) {
    double clk_tck = ( double ) CLOCKS_PER_SEC; /* system-specific clock tick duration */

#ifdef MPI
    start = MPI_Wtime(  );
#else
    start = time( NULL );
#endif
    start -= delta[0];

    times( cpu_start );
    cpu_start->tms_utime -= ( delta[1] * clk_tck );
}

/**  RestoreLog: restores the .log (and the .llog files) upon restart */
void
RestoreLog( void ) {
    char *shell_cmd;            /* used by 'system' below */
    char *outfile;              /* temporary output file name */
    outfile = NULL;

    FILE *logptr;               /* .log file pointer */
    FILE *outptr;

#ifdef MPI
    char *l_outfile;            /* temporary output file name */
    l_outfile = NULL;
    
    FILE *l_logptr;             /* local .llog file pointer */
    FILE *l_outptr;
#endif

    char *logline;              /* array of read buffers */
    long saved_count_tau;       /* count_tau as read from the .log file */
    long max_saved_count;       /* last count_tau that was saved in .log file */
    long i;                     /* loop counter */


    /* this is the last line we've written into the .log file */

    max_saved_count = state->tune.initial_moves + proc_init + count_tau * proc_tau;

    logline = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    /* get temporary file name for output file */

#ifdef MPI
    if( myid == 0 ) {
#endif
        outfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        outfile = strcpy( outfile, "logXXXXXX" );       /* required by mkstemp() */
        if( mkstemp( outfile ) == -1 )  /* get unique name for temp file */
            error( "RestoreLog: error creating temporary (log) file" );
#ifdef MPI
    }

    if( tuning && nnodes > 1 ) {
        l_outfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
        l_outfile = strcpy( l_outfile, "llogXXXXXX" );  /* required by mkstemp() */
        if( mkstemp( l_outfile ) == -1 )        /* get unique name for temp file */
            error( "RestoreLog: error creating temporary (llog) file" );
    }

    if( myid == 0 ) {
#endif

        /* restore the global .log file */

        saved_count_tau = -1;   /* reset the counter */

        /* open the log file and a temporary output file */

        logptr = fopen( files.logfile, "r" );
        if( !logptr )
            file_error( "RestoreLog (at open log file for reading)" );

        outptr = fopen( outfile, "w" );
        if( !outptr )
            file_error( "RestoreLog (at open temp log file for writing)" );

        /* read and write the first few title and caption lines */

        for( i = 0; i < 4; i++ ) {
            if( NULL == fgets( logline, MAX_RECORD, logptr ) )
                error( "RestoreLog: error reading log file captions" );
            fprintf( outptr, "%s", logline );
        }

        /* read and write the actual log lines till we are at current time */

        while( ( saved_count_tau < max_saved_count ) && ( NULL != fgets( logline, MAX_RECORD, logptr ) ) ) {
            if( 1 != sscanf( logline, "%ld", &saved_count_tau ) )
                error( "RestoreLog: error reading saved_count_tau (after %d)", saved_count_tau );
            fprintf( outptr, "%s", logline );
        }

        fclose( logptr );
        fclose( outptr );

        /* rename tmpfile into new file */

        sprintf( shell_cmd, "cp -f %s %s", outfile, files.logfile );

        if( -1 == system( shell_cmd ) )
            error( "RestoreLog: error renaming temp file %s", outfile );

        if( remove( outfile ) )
            warning( "RestoreLog: temp file %s could not be deleted", outfile );

#ifdef MPI
    }

    /* now do the same stuff for the local .llog files */

    if( tuning && nnodes > 1 ) {

        saved_count_tau = -1;   /* reset the counter */

        /*open the log file and a temporary output file */

        l_logptr = fopen( files.l_logfile, "r" );
        if( !l_logptr )
            file_error( "RestoreLog (at open llog file for reading)" );

        l_outptr = fopen( l_outfile, "w" );
        if( !l_outptr )
            file_error( "RestoreLog (at open temp llog file for writing)" );

        /* read and write the first few title and caption lines */

        for( i = 0; i < 4; i++ ) {
            if( NULL == fgets( logline, MAX_RECORD, l_logptr ) )
                error( "RestoreLog: error reading llog file captions" );
            fprintf( l_outptr, "%s", logline );
        }

        /* read and write the actual llog lines till we are at current time */

        while( ( saved_count_tau < max_saved_count ) && ( NULL != fgets( logline, MAX_RECORD, l_logptr ) ) ) {
            if( 1 != sscanf( logline, "%ld", &saved_count_tau ) )
                error( "RestoreLog: error reading saved_count_tau (after %d)", saved_count_tau );
            fprintf( l_outptr, "%s", logline );
        }

        fclose( l_logptr );
        fclose( l_outptr );

        /* rename tmpfile into new file */

        sprintf( shell_cmd, "cp -f %s %s", l_outfile, files.l_logfile );

        if( -1 == system( shell_cmd ) )
            error( "RestoreLog: error renaming temp file %s", l_outfile );

        if( remove( l_outfile ) )
            warning( "RestoreLog: temp file %s could not be deleted", l_outfile );
    }
#endif

}





/*** FUNCTIONS WHICH WRITE LOG FILES ***************************************/

/**  WriteLandscape: write iterations, temperature, dS/S, energy,        
 *                   delta_energy,  mean, std deviation, estimate_mean,  
 *                   estimate_sd, acceptance ratio                       
 *              to look at the landscape of the problem Will be either   
 *              the entire landscape or only the accepted landscape      
 */
void
WriteLandscape( char *landfile, int iteration, double delta_energy ) {
    const char *format = "  %9d %14.6f  %10.6e %16.6f %16.6f %5.2f \n";

    FILE *landptr;

    landptr = fopen( landfile, "a" );   /* first  open appropriate landscape file */
    if( !landptr )
        file_error( "WriteLandscape" );
    fprintf( landptr, format, iteration, 1.0 / S, 1.0 / dS, energy, delta_energy, acc_ratio );
    fclose( landptr );

}

/** WriteLog: writes things like mean and variation, Lam estimators, dS,   
 *             alpha and acceptance ratio to the log files and to stdout   
 *             (if -l is chosen).                                          
 */
void
WriteLog( void ) {
    FILE *logptr;               /* file pointer for global log file */

#ifdef MPI
    FILE *l_logptr;             /* file pointer for local log file (.llog) */

    if( myid == 0 ) {
#endif
        logptr = fopen( files.logfile, "a" );   /* first write to the global .log file */
        if( !logptr )
            file_error( "WriteLog" );
        PrintLog( logptr, 0 );
        fclose( logptr );

#ifdef MPI
    }

    if( write_llog && tuning && nnodes > 1 ) {
        l_logptr = fopen( files.l_logfile, "a" );       /* then do the same for .llog file */
        if( !l_logptr )
            file_error( "WriteLog" );
        PrintLog( l_logptr, 1 );
        fclose( l_logptr );
    }

    if( myid == 0 ) {
#endif

        if( log_flag ) {        /* display log to the screen? */
            PrintLog( stdout, 0 );
            fflush( stdout );
        }
#ifdef MPI
    }
#endif

}

/** PrintLog: actually prints the log to wherever it needs to be printed */
void
PrintLog( FILE * outptr, int local_flag ) {
    const char *format = "  %9d %14.6f  %10.6e %16.6f %16.6f %16.6f %16.6f %5.2f %8.5f\n";

    if( count_tau % ( print_freq * captions ) == 0 ) {
        fprintf( outptr, "\n iterations              T          dS/S            meanE" );
        fprintf( outptr, "              sdE         (e)meanE           (e)sdE" );
        fprintf( outptr, "   acc    alpha\n\n" );
    }
    /* print data */
#ifdef MPI
    if( local_flag ) {
        fprintf( outptr, format,
                 ( state->tune.initial_moves + proc_init + count_tau * proc_tau ),
                 1.0 / S, dS / S, l_mean, sqrt( l_vari ), l_estimate_mean_u, l_estimate_sd, l_acc_ratio, l_alpha );
    } else {
#endif
        fprintf( outptr, format,
                 ( state->tune.initial_moves + proc_init + count_tau * proc_tau ),
                 1.0 / S, dS / S, mean, sqrt( vari ), estimate_mean, estimate_sd, acc_ratio, alpha );
#ifdef MPI
    }
#endif
}

