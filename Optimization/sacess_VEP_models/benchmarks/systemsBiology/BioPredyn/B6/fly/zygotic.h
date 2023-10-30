/**
 * @file zygotic.h                                             
 * @authors JR, modified by Yoginho, Manu, Damjan and Anton
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Declarations of functions that deal with the right hand side of the 
 * ODEs. 
 *  
 * We have an initializing function called \c InitZygote() that reads 
 * parameters, defs data and search space limits from a data file and 
 * installs a couple of things like  
 * the solver and temporary arrays for the dvdt function.        
 * Then we have dvdt_orig, which is the inner loop of the model. 
 * It propagates the equations from one time point to another    
 * and gets called whenever Blastoderm is in PROPAGATE mode.     
 * Last but not least, we have a few mutator functions in this   
 * file. They change the local lparm EqParm struct.
 */

#ifndef ZYGOTIC_INCLUDED
#define ZYGOTIC_INCLUDED

#include <stdio.h>
#include <time.h>
#include <sys/resource.h>

#include "maternal.h"


/*** CONSTANTS *************************************************************/

/* these are the propagation rules for dvdt_orig */
extern const int INTERPHASE;
extern const int MITOSIS;

/*** ENUM ******************************************************************/

/** This is the g(u)-function enum which describes the different types of   
 * g(u) functions we can use in derivative functions                       
 */
typedef enum GFunc {
    Sqrt,
    Tanh,
    Exp,
    Hvs,
    Kolja,
} GFunc;

/*** A GLOBAL **************************************************************/

/* The g(u) function we're using - by default gofu takes the 0th element of 
   the enum, in our case Sqrt */
extern GFunc gofu;   

/* Derivative and Jacobian */
extern void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * );
extern void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * );
extern void ( *dd ) ( double *, double **, double, double *, int, SolverInput *, Input * );

/*** FUNCTION PROTOTYPES ***************************************************/

/* Initialization Functions */

/** InitZygote: makes pm and pd visible to all functions in zygotic.c and 
 *               reads EqParms and TheProblem. It then initializes bicoid  
 *               and bias (including BTimes) in maternal.c. Lastly, it     
 *               allocates memory for structures used by the derivative    
 *               function DvdtOrig that are static to zygotic.c.           
 */
Zygote InitZygote( FILE * fp, void ( *pd ) (  ), void ( *pj ) (  ), Input * inp, char *section_title );

/* Cleanup functions */

/** FreeZygote: frees memory for D */
void FreeZygote( void );

/** FreeMutant: frees mutated parameter struct */
void FreeMutant( EqParms lparm );


/* Derivative Function(s) */

/* Derivative functions calculate the derivatives for the solver. **********/

/** DvdtOrig: the original derivative function; implements the equations 
 *             as published in Reinitz & Sharp (1995), Mech Dev 49, 133-58
 *             plus different g(u) functions as used by Yousong Wang in    
 *             spring 2002.                                                
 */
void DvdtOrig( double *v, double t, double *vdot, int n, SolverInput * si, Input * inp );

/**  DvdtDelay: the delay derivative function; implements the equations 
 *             as discussed in the delay report,                    
 *             plus different g(u) functions.					    
 */
void DvdtDelay( double *v, double **vd, double t, double *vdot, int n, SolverInput * si, Input * inp );


/** Dvdt_sqrt: reimplementation of part of DvdtOrig that should make the 
 *              maintenance easier and gave a small speed up (~12%). The   
 *              function is now subdivided into 3 subfunctions that allow  
 *              easy reuse with the krylov preconditioner.                 
 */
void Dvdt_sqrt( double *v, double t, double *vdot, int n, SolverInput * si, Input * inp );

/* the preconditioner version lacks diffusion */
/*void Dvdt_sqrt_precond( double *v, double t, double *vdot, int n, int allele, Input * inp );*/

/** subfunction for (regulated) production */
void Dvdt_production( double *v, double t, double *vdot, int n, double *v_ext, DArrPtr bcd, int m, Input * inp );

/** subfunction for the degradation */
void Dvdt_degradation( double *v, double t, double *vdot, int n, Input * inp );

/** subfunction doing the diffusion */
void Dvdt_diffusion( double *v, double t, double *vdot, int n, int m, double *D, Input * inp );


/* Jacobian Function(s) */

/* Calculate the Jacobian for a given model at a given time; these funcs ***
 * are used by certain implicit solvers                                    */

/**  JacobnOrig: Jacobian function for the DvdtOrig model; calculates the 
 *               Jacobian matrix (matrix of partial derivatives) for the   
 *               equations at a give time t; input concentrations come in  
 *               v (of size n), the Jacobian is returned in jac; note that 
 *               all dfdt's are zero in our case since our equations are   
 *               autonomous (i.e. have no explicit t in them)              
 ***************************************************************************
 *                                                                         
 * The Equation: dv/dx = Rg(u) + D() + lambda()                            
 *                                                                         
 * The Jacobian: (in this example: nnucs = 3, ngenes = 3                   
 *                                                                         
 *                        i = 1       i = 2        i = 3                   
 *         b =          1   2   3   1   2   3    1   2   3                 
 *                                                                         
 *         a = 1      X11 Y12 Y13  D1   0   0    0   0   0                 
 * i = 1   a = 2      Y21 X22 Y23   0  D2   0    0   0   0                 
 *         a = 3      Y31 Y32 X33   0   0  D3    0   0   0                 
 *                                                                         
 *         a = 1       D1   0   0 X11 Y12 Y13   D1   0   0                 
 * i = 2   a = 2        0  D2   0 Y21 X22 Y23    0  D2   0                 
 *         a = 3        0   0  D3 Y31 Y32 X33    0   0  D3                 
 *                                                                         
 *         a = 1        0   0   0  D1   0   0  X11 Y12 Y13                 
 * i = 3   a = 2        0   0   0   0  D2   0  Y21 X22 Y23                 
 *         a = 3        0   0   0   0   0   0  Y31 Y32 Y33                 
 *                                                                         
 * Where:  Xab = Rg'(u) - 2D{a} - lambda{a}                                
 *         Yab = Rg'(u)                                                    
 *         Da  = D{a}                                                      
 *                                                                         
 */
void JacobnOrig( double t, double *v, double *dfdt, double **jac, int n, SolverInput * si, Input * inp );


/*** GUTS FUNCTIONS ********************************************************/

/**  CalcGuts: calculates guts for genotpye 'gtype' using unfold output in 
 *             'table' and the guts string 'gutsdefs'; it returns the num- 
 *             ber of columns we'll need to print and the guts table in    
 *             'gtable'                                                    
 */
int CalcGuts( int gindex, char *gtype, InterpObject * hist_interrp, InterpObject * extinp_interrp,
              NArrPtr table, NArrPtr * gtable, char *gutsdefs, Input * inp );

/** CalcRhs: calculates the components of the right-hand-side (RHS) of 
 *            the equatiion which make up guts (e.g. regulatory contribu-  
 *            tions of specific genes, diffusion or protein decay); this   
 *            function makes use of the derivative function and calculates 
 *            everything outside g(u) in reverse, i.e. starting from the   
 *            derivative and calculating the desired gut properties back-  
 *            wards from there; everything within g(u) is simply recon-    
 *            structed from concentrations and parameters                  
 */
void CalcRhs( double *v, double t, double *guts, int n, int gn, int numguts, int which, unsigned long *gutcomps, SolverInput * si, Input * inp );

/** ParseString: parses a line of the $gutsdefs section (dataformatX.X 
 *                has the details); this function then returns two things: 
 *                1) the number of 'words' separated by whitespace in that 
 *                   string is the function's return value                 
 *                2) the 'words' themselves are returned in arginp         
 */
int ParseString( char *v, char **arginp );

/** GetGutsComps: takes a geneID string and an array of ID strings; then 
 *                 returns an array of long ints with bit flags set that   
 *                 tell CalcGuts() what to calculate for each column of    
 *                 guts output; this function also returns a pointer to    
 *                 the current gene in the gene ID string                  
 */
char *GetGutsComps( char *geneidstring, char *egeneidstring, char **specsin, unsigned long *specsout );


/* Mutator functions */

/** Mutate: calls mutator functions according to genotype string */
EqParms Mutate( char *g_type, EqParms parm, TheProblem * defs );

/** T_Mutate: mutates genes by setting all their T matrix entries 
 *             to zero. Used to simulate mutants that express a   
 *             non-functional protein.                            
 */
void T_Mutate( int gene, int ngenes, EqParms * lparm );

/** R_Mutate: mutates genes by setting their promotor strength R 
 *             to zero, so there will be no transcription at all           
 *             anymore. Used to simulate mutants that don't pro-           
 *             duce any protein anymore.                                   
 */
void R_Mutate( int gene, EqParms * lparm );

/** RT_Mutate: mutates gene by setting both promoter strength R 
 *              and T matrix entries to zero. Useful, if you want          
 *              to be really sure that there is NO protein pro-            
 *              duction and that maternal contribution also don't          
 *              contribute anything to gene interactions.                  
 */
void RT_Mutate( int gene, int ngenes, EqParms * lparm );

/** CopyParm: copies all the parameters into the lparm struct */
EqParms CopyParm( EqParms orig_parm, TheProblem * defs );


/* A function that sets static stuff in zygotic.c */

/** SetRule: sets the static variable rule to MITOSIS or INTERPHASE */
void SetRule( int r );


/* A function that return static stuff from zygotic.c */

/*** GetParameters: returns the parm struct to the caller; note that this **
 *                  function returns the ORIGINAL PARAMETERS as they are   *
 *                  in the data file and NOT THE MUTATED ONES; this is im- *
 *                  portant to prevent limit violations in Score()         *
 ***************************************************************************/

/*EqParms setParameter(EqParms origParm, TheProblem *defs);
EqParms *GetParameters(Input *inp);
 */

/** GetMutParameters: same as above but returns the mutated copy of the 
 *                     parameter struct; important for writing guts        
 */
EqParms *GetMutParameters( void );


/***** WriteDerivLog: write to solver log file *****************************/

/* void WriteDerivLog(char *deriv, int rule, int num_nuc); */

/** Function given the starting (fb) and ending (fa) time structures,
calculates the use cpu time in seconds */
double tvsub( struct rusage a, struct rusage b );

#endif
