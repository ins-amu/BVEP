/**
 *                                                               
 *   @file sa.h                                                  
 *                                                               
 *****************************************************************
 *                                                               
 *   written by Jimmy Lam and Dan Greening                       
 *   modified by John Reinitz and Johannes Jaeger                
 *                                                               
 *****************************************************************
 *                                                               
 *   see ../doc/prog_man.ps for details of how this works        
 *                                                               
 *****************************************************************
 *                                                               
 * IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL  
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND   
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!  
 *                                                               
 *****************************************************************
 *                                                               
 * The following came with lsa.c, obtained from Dan Greening,    
 * who got it from Jimmy Lam. It is probably copyright Jimmy Lam 
 * 1988, with modifications by Greening. We (JR & JJ) have       
 * mostly modified it by adding some comments and restructuring  
 * the code to make it more easily understandable. JR has added  
 * the criterion for continous search spaces and separated Lam   
 * statistics from move acceptance statistics. King-Wai Chu has  
 * written the equilibration code.                               
 *                                                               
 *****************************************************************
 *                                                               
 * NOTE: this header only contains prototypes for functions used 
 *       for serial annealing code; all prototypes that are spe- 
 *       cific to parallel code need to go into MPI.h            
 *                                                               
 */

#ifndef SA_INCLUDED
#define SA_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include <global.h>
#endif


//#ifndef DISTRIBUTIONS_INCLUDED
//#include <distributions.h>      /* DistP.variables and prototypes */
//#endif







/*** CONSTANTS *************************************************************/

extern const int MIN_DELTA;     /* minimum exponent for Metropolis criterion */
/* provides a minimum probability for really bad moves */




/*** STRUCTS AND ENUMS *****************************************************/

/* These are the Lam parameters: *******************************************
 * (see King-Wai Chu's thesis for a detailed explanation)                  *
 *                                                                         *
 * lambda:              overall annealing accuracy                         *
 * lambda_mem_length_u: weighting factor for estimating mean energy        *
 * lambda_mem_length_v: weighting factor for estimating standard deviation *
 * initial_moves:       number of initial moves to randomize initial state *
 *                      of the system and gather initial statistics        *
 * tau:                 number of moves between updating statistics and    *
 *                      Lam parameters                                     *
 * freeze_count:        number of times system needs to be 'Frozen()' be-  *
 *                      fore we stop annealing                             *
 * update_S_skip:       number of iterations between changing the inverse  *
 *                      temperature; exists for increasing code efficiency *
 *                      and is probably obsolete now                       *
 * control:             OBSOLETE, used to be the control parameter for the *
 *                      move generation; acceptance statistics now live in *
 *                      move(s).c, since they are problem-specific         *
 * criterion:           defines the limits of what counts as 'Frozen()';   *
 *                      depends on the stop_flag (see below), i.e. it is   *
 *                      either a fixed energy (for discrete problems), a   *
 *                      limit to the change in absolute energy or the lim- *
 *                      it of proportional change in energy (good for con- *
 *                      tinous search spaces)                              *
 * mix_interval:        the mixing interval for parallel code in 'tau'     *
 *                      units, i.e. 100 means mix every 100 tau            *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * NOTE: the stopping criterion is actually defined by the combination of  *
 *       'freeze_count' and 'criterion' below plus the stop_flag which     *
 *       defines the type of the criterion (absolute or proportional);     *
 *       usually when we talk about criterion (e.g. Chu, 2001) we mean     *
 *       only the value of criterion (kappa) though, since 'freeze_count'  *
 *       and the stop_flag should always be constant for specific problems *
 *                                                                         *
 ***************************************************************************/
/*
typedef struct Files {
    char *inputfile;            // name of the input file 
    char *outputfile;           // name of the output file 
    char *statefile;            // name of the state file 
    char *logfile;              // name of the global .log file 
    char *landscapefile;        // filename of landscape file 
#ifdef MPI
    char *l_logfile;            // name of the local .llog file 

    // Files for tuning 

    char *lbfile;               // name of .lb file (for lower_bound) 
    char *ubfile;               // name of .ub file (for upper_bound) 
    char *mbfile;               // name of .mb file (for mix_bound) 
#endif
} Files;*/

typedef struct {
    double lambda;
    double lambda_mem_length_u;
    double lambda_mem_length_v;
    int initial_moves;
    int tau;
    int freeze_count;
    int update_S_skip;
    double control;
    double criterion;
//#ifdef MPI
    int mix_interval;
//#endif   

    /* These were marked "Application program must set these." in the code     */
    /* from Greening/Lam; only progname is used at the moment; we kept them    */
    /* mainly for historical reasons (and to write funny things into tunename) */

    char progname[128];         /* name of my prog */
    FILE *tunefile;             /* NULL */

} SAType;

/** This dummy is provided to allow SA to incorporate any move space; it's  
 * kinda obsolete now, but I'm too lazy to rewrite the parts of the code   
 * that use it                                                             
 */
typedef struct {
    SAType tune;
} NucStateType, *NucStatePtr;

/** this is the equilibration parameter struct for equilibration runs; it's 
 * dedicated to the legendary Chu, a rather nocturnal creature who creates 
 * bizarre code structures using arcane variable names in great redundancy;
 * although rarely seen these days, the Chu had a considerable influence   
 * on the evolution of this code; alas, most of his workings have been un- 
 * done, most of his beautifully intertwined loops have vanished and most  
 * of the redundancy has been sacrificed to the gods of parallel computing 
 * on the great altar of the Emacs-Cmd-K; up to this day, some people still
 * lament that although the code might be leaner and cleaner than in the   
 * old days, some of its soul and charm was lost in those two months of the
 * Great and Awful Cleaning that got rid of most of the traces the Chu had 
 * left around these files; Chu: may the goto be (gone) with you and may   
 * you find a well-paid job that does not involve writing code that other  
 * people ever, ever have to use!                                          
 */
typedef struct {
    double end_T;               /* equilibration temperature */
    int fix_T_skip;             /* number of steps without collecting stats */
    int fix_T_step;             /* number of steps for collecting stats */
} ChuParam;

/** Flag for type of stopping criterion (added by JR) 
 *                                                    
 * 'criterion' sets the limit within which the following must lie for
 * 'freeze_count' number of times:                                   
 *                                                                   
 * proportional freeze: (mean - old_mean)/mean                       
 * absolute freeze:      mean - old_mean                             
 * absolute energy:      mean                                        
 */
typedef enum StopStyle {
    proportional_freeze,
    absolute_freeze,
    absolute_energy
} StopStyle;





/*** GLOBALS ***************************************************************/

NucStatePtr state;              // global annealing parameter struct 

StopStyle stop_flag;            // type of stop criterion (see above) 

int time_flag;                  // flag for timing the code 
int log_flag;                   // flag for displaying log to stdout 
int nofile_flag;                /* flag for not writing .state or .log files */

long state_write;               /* frequency for writing state files (in tau) */
long print_freq;                /* frequency for printing to log (in tau) */
long captions;                  /* option for printing freqency of log captions */

/* Special annealing modes: ************************************************
 *                                                                         * 
 * the following are flags that should be set by cmd line opts for some    *
 * special annealing modes that we use for testing or gathering stats:     *
 *                                                                         *
 * - bench:    does no loop, i.e. only runs through initial steps; this    *
 *             can be useful to do solver benchmarks in -B mode            *
 * - equil:    runs at a fixed temperature to collect equilibrium stats    *
 * - quenchit: means lower the temperature to zero *instantly* after fini- *
 *             shing initialize; this is not implemented for parallel code *
 *             since it doesn't take a long time to finish                 *
 *                                                                         *
 ***************************************************************************/

int bench;
int equil;
int quenchit;




/*** FUNCTION PROTOTYPES ***************************************************/

/* problem independent functions that live in lsa.c ************************/

/* Initializing functions */


/**  Initialize: calls ParseCommandLine first; then does either initial    
 *               randomization and collecting Lam stats or restores state  
 *               of the annealer as saved in the state file                
 */
void Initialize( int argc, char **argv );

/**  InitFilenames: initializes static file names that depend on the out- 
 *                  put file name (i.e. this *must* be called after we     
 *                  have called RestoreState())                            
 */
void InitFilenames( Files * files );

/**  InitialLoop: performs the two sets of initial moves:                  
 *                   1. randomizing moves (not parallelized)               
 *                   2. loop for initial collection of statistics          
 */
void InitialLoop( void );

/** InitializeParameter: initializes variables for Lam annealing: this is  
 *                        executed after doing the initial steps;          
 *                        the local parameters only need to be set for tu- 
 *                        ning code                                        
 */
void InitializeParameter( void );

/**  InitializeWeights: initialize weights a and b for calculating Lam 
 *                      estimators; these weights are computed from the    
 *                      lambda memory length products                      
 */
void InitializeWeights( void );



/* Main loop and update functions */

/** Loop: loops (making moves, updating stats etc.) until the system is 
 *         considered frozen according to the stop criterion            
 */
void Loop( void );

/** UpdateS: update inverse temperature S at every Sskip step */
void UpdateS( void );

/** UpdateStats: updates mean, variance and acc_ratio after tau moves
 *                it needs i to do sanity check in parallel code     
 */
void UpdateStats( int i );

/** UpdateParameter: update parameters A, B, D and E and the estimators
 *                    for mean and standard deviation for the current S
 */
void UpdateParameter( void );

/**  Frozen: returns TRUE if frozen, FALSE otherwise; 'frozen' is defined  
 *           by 'freeze_count', 'stop_flag' and 'criterion' (see also sa.h 
 *           for a more extensive comment on this)                         
 */
int Frozen( void );



/* Functions that are needed for equilibration runs only */

/** FixTLoop: loops (making moves, updating stats etc.) keeping the temp-
 *             erature fixed for an equilibration run                    
*/
void FixTLoop( void );



/** GetEquil: returns the results of an equilibration run */
void GetEquil( double *equil_var );

/** SetEquilibrate: simply makes the equil_param struct static to lsa.c */
void SetEquilibrate( ChuParam ep );


/* functions which communicate with other source files, these are needed
 * for reading/writing .state files                                       
 */

/** GetLamstats: returns Lam statistics in an array of doubles; used to
 *                store Lam statistics in a state file                 
 */
double *GetLamstats( void );

/**  GetTimes: returns a two-element array with the current wallclock and  
 *             user time to be saved in the state file; for parallel code  
 *             we average the times for all processes                      
 */
double *GetTimes( void );

/** SetOutname: sets the output filename in lsa.c; this is necessary to 
 *               have the diverse .log and .ac and .mb etc files have the
 *               name of the output file if -w is chosen                 
 */
void SetOutname( char *outname );



/* functions which restore things in lsa.c upon restart from state file */

/** RestoreLamstats: restores static Lam statistics in lsa.c from an 
 *                    array of doubles; used to restore runs from a state 
 *                    file.
 */
void RestoreLamstats( double *stats );

/**  RestoreLog: restores the .log (and the .llog files) upon restart */
void RestoreLog( void );

/** RestoreTimes: restores the wallclock and user times if -t is used */
void RestoreTimes( double *delta );




/* functions to write the .log */

/**  WriteLandscape: write iterations, temperature, dS/S, energy,        
 *                   delta_energy,  mean, std deviation, estimate_mean,  
 *                   estimate_sd, acceptance ratio                       
 *              to look at the landscape of the problem Will be either   
 *              the entire landscape or only the accepted landscape      
 */
void WriteLandscape( char *landfile, int iteration, double delta_energy );

/**  InitLandscape: sets flag for printing landscape output and acceptance 
 *                 landscape and initializes the landscape file names      
 *                 called from xxx_sa.c to make filenames static and       
 *                 set landscape flag                                      
 */
void InitLandscape( int value, char *file );

/** WriteLog: writes things like mean and variation, Lam estimators, dS,   
 *             alpha and acceptance ratio to the log files and to stdout   
 *             (if -l is chosen).                                          
 */
void WriteLog( void );

/** PrintLog: actually prints the log to wherever it needs to be printed */
void PrintLog( FILE * outptr, int local_flag );



/* PROBLEM SPECIFIC FUNCTIONS THAT NEED TO BE DEFINED OUTSIDE LSA.C ********/

/* miscellaneous functions that usually live in <problem>_sa.c */

/**  ParseCommandLine: well, parses the command line and returns an index  
 *                     to the 1st argument after the command line options  
 */
int ParseCommandLine( int argc, char **argv );

/*  InitialMove: initializes the following stuff: 
 *                - reads in Lam and other annealing parameters (passed to *
 *                  lsa.c through the state_ptr; the first three arguments *
 *                  are used to open the right data file etc.)             *
 *                - initializes the cost function, establishes link be-    *
 *                  tween cost function and annealer and passes the init-  *
 *                  tial energy to lsa.c by p_chisq)                       *
 *                - initializes move generation in move(s).c               *
 *                - sets initial energy by evaluating cost function for    * 
 *                  the first time                                         *
 *                then it returns the initial temperature to the caller    *
 */
//double InitialMove(int argc, char **argv, int opt_index, NucStatePtr state_ptr, double *p_chisq);

/** MoveSA: This function actually does almost everything.
 * First it creates a static Input structure 'inp', where it puts all the 
 * information from the input file. This part is executed only once (when init == 1).
 * Then it creates a ScoreOutput structure 'out' where the score, penalty 
 * and residual vectors will be stored.
 * At the end it runs the score function, where all the calculation is done.
 */
//double MoveSA( NucStatePtr state_ptr, DistParms * distp, ScoreOutput * out, Files * files, int init, int jacobian );


void FinalMove( void );

/** WriteTimes: writes the timing information to wherever it needs to be 
 *               written to at the end of a run                            
 */
void WriteTimes( double *times );



/* move generation functions that are used in lsa.c (live in move(s).c) */

/** GenerateMove: evaluates the old energy, changes a parameter, then eval- 
 *               uates the new energy; returns the difference between old  
 *               and new energy to the caller                              
 */
//double GenerateMove( Files * files, DistParms * distp, ScoreOutput * out );

/** AcceptMove: sets new energy as the old energy for the next step and 
 *               keeps track of the number of successful moves          
 */
void AcceptMove( void );

/**  RejectMove: simply resets the tweaked parameter to the pretweak value */
void RejectMove( void );




/* a function that writes the .state file (should live in savestate.c) */

/**  StateWrite: collects Lam statistics, move state and the state of the 
 *               erand48 random number generator and writes all that into 
 *               the state file, which can then be used to restore the run
 *               in case it gets interrupted                              
 */
void StateWrite( char *infile );

#endif
