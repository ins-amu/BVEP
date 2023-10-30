/**
 * @file integrate.h                                           
 * @author JR, modified by Yoginho
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Has problem-specific stuff needed for right hand side of the ODEs.
 * 
 * The function Blastoderm runs the fly model and   
 * calls the appropriate solver for propagating the equations.   
 * PrintBlastoderm formats and prints the output of the Blasto-  
 * derm function and FreeSolution frees the solution allocated   
 * by Blastoderm. integrate.h also comes with a few utility      
 * functions for TLists which are linked lists used to initia-   
 * lize the time/mode table for Blastoderm.                      
 *                                                               
 * @note all right-hand-of-ODE-specific stuff is in zygotic.h    
 */

#ifndef INTEGRATE_INCLUDED
#define INTEGRATE_INCLUDED

#include <stdio.h>
#include <float.h>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#include "error.h"
#include "solvers.h" 
#include "maternal.h"           /* for GetBTimes to fetch bias info */
#include "zygotic.h"            /* still needed for mutators */
#include "score.h"
#include "global.h"
#include "maternal.h"

/* CONSTANTS FOR DBL_EPSILON ***********************************************/
/* EPSILON is used as minimal time step and corresponds to the minimal     */
/* double (or float) that can still be handled by C                        */

/* maybe the following declarations 
 * can be restricted, but for now they go here */

/* DBL_EPSILON is about 2 x 10^-16. so BIG_EPSILON is ~ 2 x 10^-10 min. */

#define      EPSILON     (DBL_EPSILON * 1000.)
/* This assumes all times < 100. min. !! */

#ifdef       FLOAT_EVERYTHING
#undef       EPSILON
#define      EPSILON     (FLT_EPSILON * 100.)   /* see above */
#endif

extern const double BIG_EPSILON;
extern const double HALF_EPSILON;

/* MORE CONSTANTS: OPS FOR BLASTODERM -- WITH PRIORITIES *******************/

extern const int ADD_BIAS;      /* can happen any step */
extern const int NO_OP;         /* ops below: first found set executed      */
extern const int DIVIDE;        /* NO_OP does nothing, DIVIDE divides nucs  */
extern const int PROPAGATE;     /* and PROPAGATE propagates the equations   */
extern const int MITOTATE;      /* and MITOTATE carries forward the system by epsilon (to handle rounding errors from the solver   */

/* A STRUCT ****************************************************************/

typedef struct TList {          /* Tlist is a linked list we use to set up */
    double time;                /* the solution struct in integrate.c; the */
    int n;                      /* TList has an element for each time con- */
    int op;                     /* taining the number of elements for the  */
    struct TList *next;         /* solution struct at that time (n) and a  */
} TList;                        /* rule (op) to tell the solver what to do */




/*** NEW GLOBAL ************************************************************/

extern void ( *ps ) ( double *, double *, double, double, double, double, int, FILE *, SolverInput * si, Input * );
/* this is the solver */



/* FUNCTION PROTOTYPES *****************************************************/

/* Blastoderm Functions */

/**  Blastoderm: runs embryo model and returns an array of concentration 
 *               arrays for each requested time (given by TabTimes) using  
 *               stephint as a suggested stepsize and accuracy as the re-  
 *               quired numerical accuracy (in case of adaptive stepsize   
 *               solvers or global stepsize control) for the solver.       
 *         NOTE: TabTimes *must* start from 0 and have increasing times.   
 *               It includes times when bias is added, cell division times 
 *               and times for which we have data and ends with the time   
 *               of gastrulation.                                          
 */
NArrPtr Blastoderm( int genindex, char *genotype, Input * inp, FILE * slog );

//double *BlastodermJac( int genindex, char *genotype, DArrPtr tabtimes, double stephint, double accuracy, FILE * slog );

/**  ConvertAnswer: little function that gets rid of bias times, division 
 *                  times and such and only returns the times in the tab-  
 *                  times struct as its output; this is used to produce    
 *                  unfold output that contains only the requested times;  
 *                  it also makes sure that we return the right solution   
 *                  at cell division, i.e. the solution right *after* the  
 *                  cells have divided                                     
 */
NArrPtr ConvertAnswer( NArrPtr answer, DArrPtr tabtimes );

/** FreeSolution: frees memory of the solution structure created by 
 *                 Blastoderm() or gut functions                           
 */
void FreeSolution( NArrPtr * solution );





/* TList Utility Functions */

/** InitTList: initializes TList and adds first (t=0) and last 
 *              (t=gastrulation) element of Tlist.             
 */
TList *InitTList( Zygote * zyg, int *nnucs );

/** InsertTList: takes pointer to first element of TList plus time and 
 *                desired op for new TList element and inserts a new TList 
 *                element for time at the appropriate place within the     
 *                linked list. This function returns a pointer to the      
 *                first element of the TList if insertion worked out fine. 
 */
TList *InsertTList( Zygote * zyg, TList * first, double time, int op );

/** CountEntries: counts how many entries we have in a TList */
int CountEntries( TList * first );

/** FreeTList: frees the memory of a TList */
void FreeTList( TList * first );

/** Dat2NArrPtr: Converts a data table into NarrPtr and also returns the index when
the number of nuclei is maximum */
NArrPtr Dat2NArrPtr( DataTable * table, int *maxind );

/** DoInterp: Sets up the interpolation functions for the lineage that has most
nuclei, these functions will be used later in conjunction with
Go_Backward to return history for particular times */
void DoInterp( DataTable * interp_dat, InterpObject * interp_res, int num_genes, Zygote * zyg );

/** SetFactDiscons: make one object from History and ExternalInputs */
FactDiscons SetFactDiscons( InterpObject * hist_interp_object, InterpObject * extinp_interp_object );

/** FreeFactDiscons: free FactDiscons object. FactDiscons object is made from History and ExternalInputs */
void FreeFactDiscons( double *fact_discons );

void Go_Forward( double *output, double *input, int output_ind, int input_ind, Zygote * zyg, int num_genes );
void Go_Backward( double *output, double *input, int output_ind, int input_ind, Zygote * zyg, int num_genes );
void History( double t, double t_size, double *yd, int n, InterpObject hist_interp_object, int ngenes, Zygote * zyg );
double *GetFactDiscons( int *sss, FactDiscons fd );
void FreeInterpObject( InterpObject * interp_obj );
void ExternalInputs( double t, double t_size, double *yd, int n, InterpObject extinp_interp_object, int egenes, Zygote * zyg );
//void TestInterp( int num_genes, int type );
void FreeExternalInputTemp( void );
void FreeHistoryTemp( void );


#endif
