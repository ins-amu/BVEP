/**
 * @file score.h
 * @author JR, modified by Yoginho
 *                                                                  
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Stuff that is needed for reading and defining facts and limits,
 * and for scoring functions. 
 *
 * The functions declared here initialize or manipulate facts or     
 * data time tables, read and initialize limits and penalty (if  
 * needed and do the actual scoring. It should be included in code 
 * that does scoring (e.g. printscore or annealing code).                                           
 *                                                               
 * The following describes the search space to the scoring function.       
 * There are two ways to specify limits. Lambda, R, and d always get a     
 * range for each variable---an upper & lower limit. Elements of T, h      
 * and m that contribute to u in g(u) can be given a range.                
 * This is probably the way to go as it definitely results in an           
 * ergodic search space. However, as I write this code (10/95) we have     
 * only one set of runs using this method. The alternative, which has      
 * been mostly used, is to treat T, h and m with a penalty function on     
 * u of the form:                                                          
 *                                                                         
 *           | 0 if exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +               
 *           |                                     (m_i mmax)^2 - 1 < 0    
 * penalty = |                                                             
 *           | exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +(m_i mmax)^2        
 *           |                                                otherwise    
 *                                                                         
 * where vmax_i is the highest level of gene i in the data, similarly      
 * for mmax and bcd. This method can be non-ergodic if h and T or m are    
 * being altered and variables change one at a time. In any case, one      
 * scheme or the other must be used and the programmer must insure that    
 * one or another set of pointers (see below) are NULL.                    
 *                                                                         
 * NOTE: - Lambda is NOT the same stuff as lambda in equation params or    
 *         the lambda of the Lam schedule (don't get confused!)            
 */


#ifndef SCORE_INCLUDED
#define SCORE_INCLUDED

#include "maternal.h"
#include "fly_io.h"
#include "integrate.h"          /* for blastoderm and EPSILON and stuff */
#include "solvers.h"            /* for compare() */
#include "ioTools.h"
#include "global.h"

extern const int SLEEP_LGTH;
extern const int NPOINTS;

/** Yousong's GutInfo struct for writing square diff guts */
typedef struct GutInfo {
    /** for setting the gut flag in score.c */
    int flag; 
    /** gut output precision */
    int ndigits; 
} GutInfo;


/* FUNCTION PROTOTYPES *****************************************************/

/** InitScoring: intializes a) facts-related structs and TTimes and 
 *                           b) parameter ranges for the Score function.   
 */
Scoring InitScoring( FILE * fp, int method, Input * inp );

/** InitTweak: installs tweak as a static variable in translate.c; tweak 
 *              is read from the $tweak section in newstyle data files     
 */
Tweak InitTweak( FILE * fp, int *mask, TheProblem defs );

/** InitFacts: puts facts records into the appropriate DataTable. 
 *              Returns a pointer to a Facts struture, which contains a    
 *              Datatable, which in turn points to a sized array of        
 *              DataRecords, one for each time step.                       
 */
Facts InitFacts( FILE * fp, Input * inp );

/** InitWeights: puts facts records into the appropriate DataTable. 
 *              Returns a pointer to a DataTable, which in turn points to  
 *              a sized array of DataRecords                               
 */
Weights InitWeights( FILE * fp, Input * inp );

/** getFacts: gets facts records from the appropriate DataTable. 
 *              Returns a pointer to a Facts structure, which contains a   
 *              DataTable, which in turn points to a sized array of        
 *              DataRecords, one for each time step.                       
 */
DataTable getFact( int i, GenoType * facttype );

/** InitLimits: reads limits section from the data file into the struct 
 *               limits, which is static to score.c. Then, it initializes  
 *               the penalty function if necessary.                        
 *   NOTE:       lambda limits are stored as protein half lives in the     
 *               data file and therefore get converted upon reading        
 */
SearchSpace *InitLimits( FILE * fp, Input * inp );

/** InitPenalty: initializes vmax[] and mmax static to score.c; these 
 *                variables are used to calculate the penalty function     
 *                for scoring.                                             
 *         NOTE: Should only get called, if penalty function is used       
 */
void InitPenalty( FILE * fp, TheProblem defs, SearchSpace * limits );

/** InitTTs: Initializes the time points for which we need model output 
 *            i.e. the time points for which we have data.                 
 *      NOTE: we do not allow data at t=0 but it still gets included into  
 *            TTs                                                          
 */
GenoType *InitTTs( GenoType * facttype, int nalleles );

/** InitStepsize: the only thing this function does is putting stepsize 
 *                 and accuracy in a structure to be passed to the solver 
 *                 by the Score() function                             
 */
Step_Acc InitStepsize( double step, double accuracy, FILE * slog, char *infile );


/* Actual Scoring Functions */

double checkBound( TheProblem defs, SearchSpace limits );

/** Score: as the name says, score runs the simulation, gets a solution 
 *          and then compares it to the data using the Eval least squares  
 *          function                                                       
 *   NOTE:  both InitZygote and InitScoring have to be called first!       
 */
void Score( Input * inp, ScoreOutput * out, int jacobian );

double ScoreNoCheck( void );

/** Eval: scores the summed squared differences between equation solution 
 *         and data. Because the times for states written to the Solution  
 *         structure are read out of the data file itself, we do not check 
 *         for consistency of times in this function---all times with data 
 *         will be in the table, but the table may also contain additional 
 *         times.                                                          
 */
void Eval( ScoreEval * eval, NArrPtr * Solution, int gindex, Input * inp );

/*** Scoregut functions */

/** SetGuts: sets the gut info in score.c for printing out guts */
void SetGuts( int srsflag, int ndigits );

/** GutEval: this is the same as Eval, i.e it calculates the summed squa- 
 *            red differences between equation solution and data, with the 
 *            addition that individual squared differences between data-   
 *            points are written to STDOUT in the unfold output format     
 */
void GutEval( ScoreEval * eval, NArrPtr * Solution, int gindex, Input * inp );

/* Functions that return info about score.c-specific stuff */

/** GetTTimes: this function returns the times for which there's data 
 *              for a given genotype                                       
 *     CAUTION: InitTTs has to be called first!                            
 */
DArrPtr GetTTimes( char *genotype );

/** GetLimits:  returns a pointer to the static limits struct in score.c 
 *     CAUTION:  InitScoring must be called first!                         
 */
SearchSpace *GetLimits( Input * inp );

    /** GetPenalty: calculates penalty from static limits, vmax and mmax 
 *      CAUTION: InitPenalty must be called first!                         
 */
double GetPenalty( Input * inp, SearchSpace * limits );

double GetCurPenalty( void );

/* A function for converting penalty to explicit limits */

/** Penalty2Limits: uses the inverse function of g(u) to calculate upper 
 *                   and lower limits for T, m and h using the penalty     
 *                   lambda parameter; these limits can only be an appro-  
 *                   ximation, since all T, m and h are added up to yield  
 *                   u for g(u); we try to compensate for this summation   
 *                   dividing the limits by sqrt(n); this function then    
 *                   sets the penalty vector to NULL and supplies explicit 
 *                   limits for T, m and h, which can be used for scram-   
 *                   bling parameters and such                             
 *          CAUTION: this func DOES NOT RESET pen_vec, caller must do this 
 ***************************************************************************
 *                                                                         
 * A short comment on penalty limits (JJ, Aug 7, 2001):                    
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
 * See JJs lab notes for further detail on g(u)-inverse and such.          
 */
void Penalty2Limits( SearchSpace * limits, TheProblem defs );

/* Functions used to read stuff into structs */

/** Takes a Dlist and returns the corresponding DataTable structure we use 
 * for facts data.                          
 *                                                                         
 * An extensive comment about indices: (by JJ) 
 *                                                                         
 * All data points with -1 as a value are NOT read from the  
 * data file. The difference between such ignored and zero   
 * values is crucial: -1 data points WILL NOT BE COMPARED TO 
 * simulation data, whereas 0 means NO PROTEIN AT THAT TIME  
 * IN THAT NUCLEUS.                                          
 * Index numbers help maintain the integrity of the data.    
 * An index number is defined as the array index at which a  
 * protein concentration would be if the data was complete,  
 * i.e. available for all nuclei at all times. In this way   
 * a sparse set of data can be compared to a complete set of 
 * simulation output.                                        
 * Thus, indices are defined as starting from 1 for each     
 * DataRecord (each time step) and increase by one for each  
 * gene in each nucleus in the order predefined by JR.       
 */
DataTable *List2Facts( Dlist * inlist, int ngenes );

/** getMediansFromDataTable: Function to calculate medians per gene from data */
double *getMediansFromDataTable( DataTable * data, int ngenes );

/** Rescale: Function to rescale data - for every gene multiply data per previously calculated multiplier */
void Rescale( DataTable * data, double *multipliers, int ngenes );

/*** PrintFact:  ***************************************************************************/

int isResComp(  );
void setCostType( char *cost );

void PrintFacts( FILE * fp, int ndigits, int columns );
void initResComp( NArrPtr Solution );

/** FreeFacts: Function to fre tha DataTable object */
void FreeFacts( DataTable * D );

/** InitHistory: Initializing the full set of nuclei based on the lineages of the history,
       please make sure that all of the alleles' lineages are the same */
InterpObject *InitHistory( FILE * fp, Input * inp );
InterpObject *InitExternalInputs( FILE * fp, Input * inp );
void GetInterp( FILE * fp, char *title, Input * inp, int num_genes, DataTable ** interp_tables );
void FreeHistory( int nalleles, InterpObject * polations );
void FreeExternalInputs( int nalleles, InterpObject * extinp_polations );

#endif
