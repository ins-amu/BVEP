/**
 *                                                               
 *   @file MPI.h                                                 
 *                                                               
 *****************************************************************
 *                                                               
 *   written by John Reinitz                                     
 *   modified by King-Wai Chu and Johannes Jaeger                
 *                                                               
 *   last modified, 06/13/2002                                   
 *                                                               
 *****************************************************************
 *                                                               
 * IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL  
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND   
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!  
 *                                                               
 *****************************************************************
 *                                                               
 * MPI.h contains structs and constants that are specific to     
 * parallel annealing code using MPI.                            
 *                                                               
 * This includes prototypes for all the tuning functions below.  
 *                                                               
 * There are two problem-specific functions declared below that  
 * need to be defined in move(s).c.                              
 *                                                               
 *****************************************************************
 *                                                               
 * NOTE: this header only contains prototypes for functions used 
 *       for parallel annealing code only; all prototypes of     
 *       functions that include serial code need to go into sa.h 
 *                                                               
 */

#ifndef MPI_INCLUDED_SA
#define MPI_INCLUDED_SA



/*** CONSTANTS *************************************************************/

extern const int MAX_MIX;       /* max number of mixes during a tuning run */
                    /* set this to a lower number if you run out of memory */
extern const int GROUP_SIZE;    /* group size for calculating upper bound */

extern const int STOP_TUNE_CNT; /* stop tune count */
extern const int STOP_TUNE_CRIT;        /* tuning stop criterion */

extern const int LSTAT_LENGTH;  /* length of Lam msg array when annealing */
extern const int LSTAT_LENGTH_TUNE;     /* length of Lam msg array when tuning */



/*** PARALLEL GLOBALS ******************************************************/

int myid;                       /* id of local node (processor) */
int nnodes;                     /* number of nodes (processors), also known as P */

int tuning;                     /* flag for switching on tuning mode */
int covar_index;                /* covariance sample index for tuning (in 'tau' units) */
int write_tune_stat;            /* how often to write tuning statistics */
int auto_stop_tune;             /* auto stop tune flag to stop tuning runs early */
int write_llog;                 /* flag for writing local log files */



/*** FUNCTION PROTOTYPES ***************************************************/

/* lsa.c: parallel non-tuning funcs: update func for local Lam parameters */

/** UpdateLParameter: update local parameters l_A, l_B, l_D and l_E and 
 *                    the local estimators for mean and standard deviation 
 *                    for both upper and lower bounds of M_Opt for the     
 *                    current S                                            
 */
void UpdateLParameter( void );



/* lsa.c: parallel non-tuning funcs: mixing functions */

/** DoMix: does the mixing; sends move state and local Lam stats to the 
 *          dance partner(s)                                            
 */
void DoMix( void );

/** DoFixMix: for equilibration run, we only need to pass the move state 
 *             since Lam stats are not needed at constant temperature    
 */
void DoFixMix( void );

/** MakeLamMsg: packages local Lam stats into send buffer */
void MakeLamMsg( double **sendbuf );

/** AcceptLamMsg: receives new energy and Lam stats upon mixing */
void AcceptLamMsg( double *recvbuf );



/* lsa.c: tuning functions */

/** InitTuning: sets up/restores structs and variables for tuning runs */
void InitTuning( void );

/** DoTuning: calculates the cross-correlation (for lower bound) and the 
 *             variance of local means (for upper bound) for a sub_tune_   
 *             interval; the results are added to arrays, which span the   
 *             whole tuning_interval; when these values get written to the 
 *             'bound' files, they will be divided by the number of mixes  
 *             we have done already to average them out (see also King-Wai 
 *             Chu's thesis, p. 63)                                        
 */
void DoTuning( void );

/** WriteTuning: writes tuning stats to files ever tune_interval */
void WriteTuning( void );

/** StopTuning: is to tuning runs what Frozen() is to a normal annealing 
 *               run: it basically checks if the tuning stop criterion   
 *               applies and returns true if that's the case             
 */
int StopTuning( void );




/* move(s).c: functions for communicating move state for mixing */

/** MakeStateMsg: function to prepare a message which is then passed 
 *                 to other nodes via MPI telling the other nodes about    
 *                 move stats; since we don't know about the move state    
 *                 structs, but can assume that we'll only have to send    
 *                 longs and doubles, we split the message in two arrays,  
 *                 one for the longs and one for the doubles; then we re-  
 *                 turn the arrays and their sizes for lsa.c to send them  
 *                 to the dance partners                                   
 *                                                                         
 *                 note that all the arguments need to get passed by refe- 
 *                 rence since we need to allocate the arrays according to 
 *                 their problem-specific size in move(s).c                
 */
void MakeStateMsg( long **longbuf, int *lsize, double **doublebuf, int *dsize );

/** AcceptMsg: communicates a message about move stats received via MPI 
 *              to move(s).c; see the comment for MakeStateMsg for the ra- 
 *              tionale behind the two arrays that are passed              
 */
void AcceptStateMsg( long *longbuf, double *doublebuf );

#endif
