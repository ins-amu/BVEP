/**
 * @file integrate.c                                        
 * @author JR, modified by Yoginho
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Implementation of functions needed for the right hand side of ODEs. 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/times.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include "error.h"              /* for error handling functions */
//#include "maternal.h"           /* for GetBTimes to fetch bias info */
#include "integrate.h"          /* obviously */
//#include "solvers.h"            /* for .slog filename */
//#include "zygotic.h"            /* still needed for mutators */
//#include "score.h"


/**********************************************************************
 **********************************************************************
 * The routines below, and all executible programs which use them,    *
 * are dedicated to Jerry Garcia, who died on August 9, 1995.         *
 * If you are using any of these functions, please do not             *
 * alter the following external array, which should remain            *
 * readable in the executable, using the strings(1) or what(1) cmds.  *
 * You may remove the const qualifier for non-ANSI compilers, and     *
 * you are invited to add your own functions  and  copy  this         *
 * notice to your own code. JR, August 14, 1995. Jerry, we miss you.  *
 **********************************************************************
 **********************************************************************/


/* CONSTANTS */
const double BIG_EPSILON = ( EPSILON * 1000. );
const double HALF_EPSILON = ( EPSILON * .5 );

/* MORE CONSTANTS: OPS FOR BLASTODERM -- WITH PRIORITIES *******************/

const int ADD_BIAS = 1;         /* can happen at any step */
const int NO_OP = 2;            /* ops below: first found set executed      */
const int DIVIDE = 4;           /* NO_OP does nothing, DIVIDE divides nucs  */
const int PROPAGATE = 8;        /* and PROPAGATE propagates the equations   */
const int MITOTATE = 16;        /* and MITOTATE carries forward the system by epsilon (to handle rounding errors from the solver   */


void ( *ps ) ( double *, double *, double, double, double, double, int, FILE *, SolverInput * si, Input * );


const char Jerry[] = "@(#) In Memoriam Jerome John Garcia, 8/1/42-8/9/95";

// size of jacobian???
int jacSize;

/* variables used for timing */

/* these variables are used to evaluate real time using time() */
extern double start;            /* wallclock time before run */
extern double finish;           /* wallclock time after run */
/* these structs are used to evaluate user time using times() or MPI_Wtime */
extern struct tms *cpu_start;   /* user time before run */
extern struct tms *cpu_finish;

double *delta;

// For the delay solver
//static double *fact_discons, fact_discons_size;

/*** BLASTODERM FUNCTIONS **************************************************/

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
NArrPtr
Blastoderm( int genindex, char *genotype, Input * inp, FILE * slog ) {

    SolverInput si;
    const double epsilon = EPSILON;     /* epsilons: very small in- */
    const double big_epsilon = BIG_EPSILON;     /* creases used for division */

    NArrPtr solution;           /* solution will be an array of */
    /* concs for each requested time */

    int *what2do = NULL;        /* what to do at each time step */

    double *divtable = NULL;    /* cell div times in reverse order */
    double *transitions = NULL; /* this is when cell divs start */
    double *durations = NULL;   /* durations of cell divisions */

    DArrPtr biastimes;          /* times a which bias is added */
    //DArrPtr   oldtimes;                  /* statically saves tabtimes */ //we want to keep this in memory to avoid initialization in case
    //our times are repeated
    DArrPtr bias;               /* bias for given time & genotype */

    //static int       allocate;             /* flag: need to allocate or not? */
    int i, ii, j;               /* loop counters */
    int k;                      /* index of gene k in current nuc */
    int ap;                     /* nuc. position on AP axis */
    int lin;                    /* first lineage number at each ccycle */

    int rule;                   /* MITOSIS or INTERPHASE? */


    TList *entries = NULL;      /* temp linked list for times and */
    TList *current;             /* ops for the solver */
    jacSize = 0;
    /* allocate transitions array (skip if ndivs == 0) */
    if( inp->zyg.defs.ndivs > 0 )
        transitions = ( double * ) calloc( inp->zyg.defs.ndivs, sizeof( double ) );
    /* for each genotype, the 'genotype' variable has to be made static to zy- *
     * gotic.c so that the derivative functions know which genotype they're    *
     * dealing with (i.e. they need to get the appropriate bcd gradient)       */
    InitDelaySolver(  );
    si.genindex = genindex;
    si.all_fact_discons = SetFactDiscons( &( inp->his[genindex] ), &( inp->ext[genindex] ) );

    /* INITIALIZATION OF THE MODEL STRUCTS AND ARRAYS ************************* */

    solution.array = NULL;
    solution.size = 0;

    what2do = NULL;

    biastimes.array = NULL;
    biastimes.size = 0;


    /* get bias times and initialize information about cell divisions */

    biastimes = GetBTimes( genotype, &( inp->zyg ) );
    if( !( biastimes.array ) )
        error( "Blastoderm: error getting bias times" );
    if( inp->zyg.defs.ndivs > 0 ) {
        if( !( divtable = inp->zyg.times.div_times ) )
            error( "Blastoderm: error getting division table" );
        if( !( durations = inp->zyg.times.div_duration ) )
            error( "Blastoderm: error getting division durations" );
        for( i = 0; i < inp->zyg.defs.ndivs; i++ ) {
            transitions[i] = divtable[i] - durations[i];
        }
    }
    /* entries is a linked list, which we use to set up the solution struct    *
     * and the what2do array; each of these needs an entry for:                *
     * - start and end (gastrulation) time                                     *
     * - times for mitoses: - beginning of mitosis                             *
     *                      - cell division time (still belongs to previous    *
     *                        cleavage cycle)                                  *
     *                      - time right after cell division (+EPSILON), be-   *
     *                        longs to new cell cycle with doubled nnucs       *
     * - times at which we add bias                                            *
     * - tabulated times for which we have data or which we want to display    */

    /* add start and end (gastrulation) time */
    entries = InitTList( &( inp->zyg ), inp->zyg.nnucs );
    /* add all times required for mitoses (skip this for 0 div schedule) */
    for( i = 0; i < inp->zyg.defs.ndivs; i++ ) {

        entries = InsertTList( &( inp->zyg ), entries, divtable[i], DIVIDE );
        if( GetNNucs( &( inp->zyg.defs ), inp->zyg.nnucs, divtable[i], &( inp->zyg.times ) ) ==
            GetNNucs( &( inp->zyg.defs ), inp->zyg.nnucs, ( divtable[i] + epsilon ), &( inp->zyg.times ) ) )
            error( "Blastoderm: epsilon of %g too small! %g ", epsilon, divtable[i] );
        entries = InsertTList( &( inp->zyg ), entries, divtable[i] + epsilon, PROPAGATE );
        entries = InsertTList( &( inp->zyg ), entries, transitions[i], MITOTATE );
        if( GetNNucs( &( inp->zyg.defs ), inp->zyg.nnucs, transitions[i], &( inp->zyg.times ) ) !=
            GetNNucs( &( inp->zyg.defs ), inp->zyg.nnucs, ( transitions[i] + epsilon ), &( inp->zyg.times ) ) )
            error( "Blastoderm: division within epsilon of %g! %g %g ", epsilon, transitions[i], durations[i] );
        entries = InsertTList( &( inp->zyg ), entries, transitions[i] + epsilon, PROPAGATE );
    }
    /* add bias times */
    
    for( i = 0; i < biastimes.size; i++ ) {
        entries = InsertTList( &( inp->zyg ), entries, biastimes.array[i], ADD_BIAS | PROPAGATE );
    }

    /* tabulated times */
    for( i = 0; i < inp->sco.facts.tt[genindex].ptr.times.size; i++ ) {
        entries = InsertTList( &( inp->zyg ), entries, inp->sco.facts.tt[genindex].ptr.times.array[i], PROPAGATE );
    }
    /* now we know the number of solutions we have to calculate, so we can     *
     * allocate and initialize the solution struct and the what2do array       */
    solution.size = CountEntries( entries );
    solution.array = ( NucState * ) calloc( solution.size, sizeof( NucState ) );
    what2do = ( int * ) calloc( solution.size, sizeof( int ) );
    current = entries;
    for( i = 0; i < solution.size; i++ ) {
        solution.array[i].time = current->time;
        solution.array[i].state.size = current->n;
        solution.array[i].state.array = ( double * ) calloc( current->n, sizeof( double ) );
        for( j = 0; j < current->n; j++ )       //check if this zeroing is really needed or it was just for testing purposes
            solution.array[i].state.array[j] = 0;
        what2do[i] = current->op;
        current = current->next;
    }
    /* free the linked list and save new tab times in oldtimes array */
    FreeTList( entries );

    /* RUNNING THE MODEL ****************************************************** */
    /* Before running the model, mutate zygotic params appropriately */
    FreeMutant( inp->lparm );
    inp->lparm = Mutate( genotype, inp->zyg.parm, &( inp->zyg.defs ) );
    /*if (debug) {
       fprintf(slog, "\n--------------------------------------------------");
       fprintf(slog, "--------------------------------------------------\n");
       fprintf(slog, "Blastoderm: mutated genotype to %s.\n", genotype);
       } */
    /* Below is the loop that evaluates the solution and saves it in the solu- *
     * tion struct      
     *                                                        */
    //printf("SOLUTION BEFORE %lg %lg %lg %lg\n", solution.array[0].state.array[0], solution.array[0].state.array[1], solution.array[0].state.array[2], solution.array[0].state.array[3]);
    for( i = 0; i < solution.size; i++ ) {
        //printf("%d) %lg %lg %lg %lg\n", i, solution.array[i].state.array[0], solution.array[i].state.array[1], solution.array[i].state.array[2], solution.array[i].state.array[3]);
        /* First we have to set the propagation rule in zygotic.c (either MITOSIS  *
         * or INTERPHASE) to make sure that the correct differential equations are *
         * propagated (rules are defined in zygotic.h)                             */
        rule = GetRule( solution.array[i].time, &( inp->zyg ) );
        si.time = solution.array[i].time;
        //printf("WHAT2DO = %d at step %d\n", what2do[i], i);
        /*if (debug) {
           if (rule == 0) { //rule = INTERPHASE
           fprintf(stderr, "Blastoderm: rule is INTERPHASE.\n");
           } else if (rule == 1 && !(what2do[i] & DIVIDE)) { //rule = MITOSIS and !(op = DIVIDE)
           fprintf(stderr, "Blastoderm: rule is MITOSIS.\n");
           } else if (what2do[i] & DIVIDE) { //rule = MITOSIS and op = DIVIDE
           fprintf(stderr, "Blastoderm: rule is DIVISION.\n");
           }
           } */

        /* ADD_BIAS is a special op in that it can be combined with any other op   *
         * (see also next comment); we can add (or subtract) protein conentrations *
         * at any time; this is for adding initial conditions for genes that are   *
         * both maternal and zygotic (e.g. cad and hb) or to simulated perturba-   *
         * tions like heat shocks or induced overexpression of certain genes;      *
         * note that the bias lives in maternal.c and has to be fetched from there */
        if( what2do[i] & ADD_BIAS ) {
            //printf("%d-%d ADD_BIAS\n", i, what2do[i]);  
            for( j = 0; j < biastimes.size; j++ ) {
                if( fabs( solution.array[i].time - biastimes.array[j] ) < big_epsilon ) {
                    bias = GetBias( biastimes.array[j], genindex, &( inp->zyg ) );      //bias is an array with concentrations
                    //Here we choose if we want to "add" or "set" bias concentrations;
                    //just uncomment the line you need (and comment the other one)
                    for( ii = 0; ii < bias.size; ii++ )
                        //solution.array[i].state.array[ii] += bias.array[ii]; //adding bias concentrations to the system
                        solution.array[i].state.array[ii] = bias.array[ii];    //setting bias concentrations to the system
                }
            }
            if( debug )
                fprintf( slog, "Blastoderm: added bias at time %f.\n", solution.array[i].time );
        }
        /* The ops below can be executed in addition to ADD_BIAS but they cannot   *
         * be combined between themselves; if more than one op is set, the prio-   * 
         * rities are as follows: NO_OP > DIVIDE > PROPAGATE; please make sure     *
         * that you do not set double ops, which will only obfuscate the code;     * 
         *                                                                         *
         * NO_OP simply does nothing (used for gastrulation time)                  */

        if( what2do[i] & NO_OP ) {
            //printf("%d-%d NO_OP\n", i, what2do[i]);
            ;
        }

        /* This is the DIVIDE rule: only symmetrical division is supported yet,    *
         * i.e. both daughter nuclei inherit all concentrations from their mother; *
         * we have to take special care of the fact that we might not need the     *
         * most anterior and most posterior daughter nuclei, depending on the ran- *
         * ge of nuclei in cycle 14 we've chosen to anneal on; e.g. if the most    *
         * anterior nucleus in cycle 14 has an odd lineage number, we need to push *
         * its sibling off the anterior limit when we divide after cycle 13; or in *
         * other words: our most anterior cell in cycle 14 is the posterior daugh- *
         * ter of our most anterior cell in cycle 13; therefore, we need to 'loose'*
         * its anterior sibling; the same applies for the posterior end, where we  *
         * want to loose the most posterior daughter cell if we don't need it any  *
         * more at the later cycle                                                 *
         *                                                                         *
         * This is implemented as follows below: lin is the lineage number of the  *
         * most anterior cell of the next cell cycle; if it's odd numbered, we'll  *
         * push off the most anterior daughter by subtracting the number of genes  * 
         * from the daughter indices (i.e. we shift the solution array for the     *
         * next cycle posteriorly by one nucleus); if on the other hand, the most  *
         * posterior nucleus lies outside our new array, we just forget about it   */
        else if( what2do[i] & DIVIDE ) {
            //printf("%d-%d DIVIDE\n", i, what2do[i]);
            lin = GetStartLin( solution.array[i + 1].time, inp->zyg.defs, inp->zyg.lin_start, &( inp->zyg.times ) );
            for( j = 0; j < solution.array[i].state.size; j++ ) {
                k = j % inp->zyg.defs.ngenes;   /* k: index of gene k in current nucleus */
                ap = j / inp->zyg.defs.ngenes;  /* ap: rel. nucleus position on AP axis */

                /* evaluate ii: index of anterior daughter nucleus */

                if( lin % 2 )
                    ii = 2 * ap * inp->zyg.defs.ngenes + k - inp->zyg.defs.ngenes;
                else
                    ii = 2 * ap * inp->zyg.defs.ngenes + k;

                /* skip the first most anterior daughter nucleus in case lin is odd */
                if( ii >= 0 )
                    solution.array[i + 1].state.array[ii] = solution.array[i].state.array[j];
                /* the second daughter only exists if it is still within the region */
                if( ii + inp->zyg.defs.ngenes < solution.array[i + 1].state.size )
                    solution.array[i + 1].state.array[ii + inp->zyg.defs.ngenes] = solution.array[i].state.array[j];
            }
            /* Divide the history of the delay solver */
            DivideHistory( solution.array[i].time, solution.array[i + 1].time, &( inp->zyg ) );

            if( debug )
                fprintf( slog, "Blastoderm: nuclear division %d -> %d nuclei.\n", solution.array[i].state.size / inp->zyg.defs.ngenes,
                         solution.array[i + 1].state.size / inp->zyg.defs.ngenes );
        }
        /* This is the MITOTATE rule, which advances the blastoderm from
           transitions[i] to transitions[i] + epsilon, from where we propagate
           again. This is there to ensure that rounding errors from the solver do
           not lead to incorrect rules being used at the end-points */

        else if( what2do[i] & MITOTATE ) {
            //printf("%d-%d MITOTATE\n", i, what2do[i]);
            for( j = 0; j < solution.array[i].state.size; j++ )
                solution.array[i + 1].state.array[j] = solution.array[i].state.array[j];
        }
        /* In case we have to PROPAGATE the differential equeations, we call the   *
         * solver; you have to make sure that the appropriate rule has been set in *
         * zygotic.c (use SetRule(), rules are MITOSIS or INTERPHASE), otherwise   *
         * you will propagate the wrong equations; the solver needs pointers to    *
         * arrays of concentration for start and end time as well as those start   *
         * and end times themselves; all solvers need a suggested stepsize, adap-  *
         * tive stepsize solvers need the accuracy argument; we also need the accu-*
         * racy argument for global stepsize control (if ever implemented); lastly *
         * we need to tell the solver how big input and output arrays are          */

        else if( what2do[i] & PROPAGATE ) {
            //printf("%d-%d PROPAGATE\n", i, what2do[i]);
            /*
               if (debug) {
               printf("From %d, %lg %lg %lg %lg\n", i, solution.array[i].state.array[0], solution.array[i].state.array[1], solution.array[i].state.array[2], solution.array[i].state.array[3]);
               }
             */
            ( *ps ) ( solution.array[i].state.array, solution.array[i + 1].state.array, solution.array[i].time, solution.array[i + 1].time, inp->ste.stepsize,
                      inp->ste.accuracy, solution.array[i].state.size, slog, &si, inp );
            /*
               if (debug) {
               printf("To %d, %lg %lg %lg %lg\n", i+1, solution.array[i+1].state.array[0], solution.array[i+1].state.array[1], solution.array[i+1].state.array[2], solution.array[i+1].state.array[3]);
               }
             */
            jacSize += solution.array[i].state.size;

        } else {                /* unknown op? -> error! */
            error( "op was %d!?", what2do[i] );
        }

    }                           //end for
    //printf("AFTER %lg %lg %lg %lg\n", solution.array[solution.size-1].state.array[0], solution.array[solution.size-1].state.array[1], solution.array[solution.size-1].state.array[2], solution.array[solution.size-1].state.array[3]);
    /* After having calculated the solution, free the mutant parameter structs * 
     * and return the result                                                   */

    FreeDelaySolver(  );
    FreeFactDiscons( si.all_fact_discons.fact_discons );
    free( what2do );
    free( transitions );
    return solution;
}

// don't really understand how or why this works like this
/*
 double*  BlastodermJac(int genindex, char *genotype, DArrPtr tabtimes,
           double stephint, double accuracy, TheProblem defs, FILE *slog)
{
  const  double    epsilon     = EPSILON;      // epsilons: very small in-
  const  double    big_epsilon = BIG_EPSILON; // creases used for division

  static NArrPtr   solution;               // solution will be an array of
                                          // concs for each requested time
  static int       *what2do;               // what to do at each time step
  int cc;
  double *dfdt;          
  double **jac;                                     // the Jacobian matrix
  double *finalJac;
  int jindex=0;
  static double    *divtable;           // cell div times in reverse order
  static double    *transitions;           // this is when cell divs start
  static double    *durations;              // durations of cell divisions

  static DArrPtr   biastimes;               // times a which bias is added
  static DArrPtr   oldtimes;                  // statically saves tabtimes

  static int       allocate;             // flag: need to allocate or not?

  int              i,ii,j,ij,ijj;                                // loop counters
  int              k;                    // index of gene k in current nuc
  int              ap;                         // nuc. position on AP axis
  int              lin;             // first lineage number at each ccycle

  int              rule;                         // MITOSIS or INTERPHASE?

  DArrPtr          bias;                 // bias for given time & genotype

  TList            *entries    = NULL;   // temp linked list for times and
  TList            *current;                         // ops for the solver


  if ( !(finalJac = (double *)calloc(jacSize, sizeof(double))))
    error("BaDe: error allocating jac.\n");
  
 allocate transitions array (skip if ndivs == 0) 
 
  if ( (transitions == NULL) && (defs.ndivs > 0) )
    transitions = (double *)calloc(defs.ndivs, sizeof(double));
 */

/* for each genotype, the 'genotype' variable has to be made static to zy- *
 * gotic.c so that the derivative functions know which genotype they're    *
 * dealing with (i.e. they need to get the appropriate bcd gradient)       */


//InitGenotype(genindex);
//update: we don't use this static variable anymore

/* Blastoderm() first checks if tabtimes has the same size as in the pre-  *
 * vious run (remembered by the static oldtimes array) and if all tab ti-  *
 * mes are the same; if this is the case, the initialization is skipped    */

/*
   if ( tabtimes.size == oldtimes.size ) {
   for ( i=0; i<tabtimes.size; i++ ) 
     if ( tabtimes.array[i] != oldtimes.array[i] ) 
   allocate = 1;
 } else
   allocate = 1;
 */
/* INITIALIZATION OF THE MODEL STRUCTS AND ARRAYS **************************/
/*
// free structures that are already present

  if ( allocate ) {
    if ( solution.array != NULL ) { 
      FreeSolution(solution);   
      solution.array = NULL;      
      solution.size  = 0;
    }
    if (what2do != NULL) { 
      free(what2do);     
      what2do = NULL;   
    }
    if (biastimes.array != NULL) {  
      free(biastimes.array);      
      biastimes.array = NULL;   
      biastimes.size  = 0;
    }
    if (oldtimes.array != NULL) { 
      free(oldtimes.array);      
      oldtimes.array = NULL;     
      oldtimes.size  = 0;
    }
 
// get bias times and initialize information about cell divisions

    biastimes = GetBTimes(genotype);
    if ( !(biastimes.array) )
      error("Blastoderm: error getting bias times");
    
    if ( defs.ndivs > 0 ) {
      if ( !(divtable =inp->zyg.times.div_times) )
    error("Blastoderm: error getting division table");
      if ( !(durations =inp->zyg.times.div_duration) )
    error("Blastoderm: error getting division durations");
      for(i=0; i<defs.ndivs; i++) 
    transitions[i] = divtable[i] - durations[i];
    }
 */
/* entries is a linked list, which we use to set up the solution struct    *
 * and the what2do array; each of these needs an entry for:                *
 * - start and end (gastrulation) time                                     *
 * - times for mitoses: - beginning of mitosis                             *
 *                      - cell division time (still belongs to previous    *
 *                        cleavage cycle)                                  *
 *                      - time right after cell division (+EPSILON), be-   *
 *                        longs to new cell cycle with doubled nnucs       *
 * - times at which we add bias                                            *
 * - tabulated times for which we have data or which we want to display    */
/*
// add start and end (gastrulation) time

    entries = InitTList();

// add all times required for mitoses (skip this for 0 div schedule)

    if ( defs.ndivs > 0 )
      for (i=0; i<defs.ndivs; i++) { 
    entries = InsertTList(entries, divtable[i], DIVIDE);
    if( (GetNNucs(divtable[i])) == (GetNNucs(divtable[i] + epsilon)) )
      error("Blastoderm: epsilon of %g too small!", epsilon);
    entries = InsertTList(entries, divtable[i] + epsilon, PROPAGATE);
    entries = InsertTList(entries, transitions[i], PROPAGATE);
      }

// add bias times

    for(i=0; i < biastimes.size; i++)
      entries = 
    InsertTList(entries, biastimes.array[i], ADD_BIAS | PROPAGATE);
    
// tabulated times

    for(i=0; i < tabtimes.size; i++)
      entries = InsertTList(entries, tabtimes.array[i], PROPAGATE);
 */
/* now we know the number of solutions we have to calculate, so we can     *
 * allocate and initialize the solution struct and the what2do array       */

/*
 solution.size  = CountEntries(entries);
 solution.array = (NucState *)calloc(solution.size, sizeof(NucState));
    
 what2do = (int *)calloc(solution.size, sizeof(int));
    
 current = entries;       
 for(i=0; i<solution.size; i++) {          
   solution.array[i].time = current->time;   
   solution.array[i].state.size = current->n; 
   solution.array[i].state.array =            
 (double *)calloc(current->n, sizeof(double));
   for(j=0; j<current->n; j++)
 solution.array[i].state.array[j] = 0;
   what2do[i] = current->op;            
   current = current->next;          
 }
   
// free the linked list and save new tab times in oldtimes array

 FreeTList(entries);
 entries = NULL;
    
 oldtimes.size  = tabtimes.size;         
 oldtimes.array = (double *)calloc(tabtimes.size, sizeof(double));
 for(i=0; i< tabtimes.size; i++)
   oldtimes.array[i] = tabtimes.array[i];

 allocate = 0; 
   
// if we're not allocating: reset initial conditions; we'll *add* them again below

} else
    
 for(i=0; i < solution.array[0].state.size; i++)
   solution.array[0].state.array[i] = 0;
  
 */
/* RUNNING THE MODEL *******************************************************/

/*
// Before running the model, mutate zygotic params appropriately

  Mutate(genotype);

  if ( debug ) {
    fprintf(slog, "\n--------------------------------------------------");
    fprintf(slog, "--------------------------------------------------\n");
    fprintf(slog, "Blastoderm: mutated genotype to %s.\n", genotype);
  }

// Below is the loop that evaluates the solution and saves it in the solution struct

  for(i=0; i<solution.size; i++) {
 */

/* First we have to set the propagation rule in zygotic.c (either MITOSIS  *
 * or INTERPHASE) to make sure that the correct differential equations are *
 * propagated (rules are defined in zygotic.h)                             */

/*  rule = GetRule(solution.array[i].time);
  SetRule(rule); 

  if ( debug ) {
      if ( rule == 0 ) {
          fprintf(slog, "Blastoderm: rule is INTERPHASE.\n");
      } else if ( rule == 1 && !(what2do[i] & DIVIDE) ) {
          fprintf(slog, "Blastoderm: rule is MITOSIS.\n");
      } else if ( what2do[i] & DIVIDE ) {
          fprintf(slog, "Blastoderm: rule is DIVISION.\n");
      }
  }
 */

/* ADD_BIAS is a special op in that it can be combined with any other op   *
 * (see also next comment); we can add (or subtract) protein conentrations *
 * at any time; this is for adding initial conditions for genes that are   *
 * both maternal and zygotic (e.g. cad and hb) or to simulated perturba-   *
 * tions like heat shocks or induced overexpression of certain genes;      *
 * note that the bias lives in maternal.c and has to be fetched from there */

/*
 if (what2do[i] & ADD_BIAS) {
   for (j=0; j<biastimes.size; j++) {
 if ( fabs(solution.array[i].time - biastimes.array[j]) 
      < big_epsilon ) {
   bias = GetBias(biastimes.array[j], genindex);
   for (ii=0; ii < bias.size; ii++)
     solution.array[i].state.array[ii] += bias.array[ii];
 }
   }

   if ( debug )
 fprintf(slog, "Blastoderm: added bias at time %f.\n",
     solution.array[i].time);
 }
 */

/* The ops below can be executed in addition to ADD_BIAS but they cannot   *
 * be combined between themselves; if more than one op is set, the prio-   * 
 * rities are as follows: NO_OP > DIVIDE > PROPAGATE; please make sure     *
 * that you do not set double ops, which will only obfuscate the code;     * 
 *                                                                         *
 * NO_OP simply does nothing (used for gastrulation time)                  */

/*
    if (what2do[i] & NO_OP)
   ;
 */

/* This is the DIVIDE rule: only symmetrical division is supported yet,    *
 * i.e. both daughter nuclei inherit all concentrations from their mother; *
 * we have to take special care of the fact that we might not need the     *
 * most anterior and most posterior daughter nuclei, depending on the ran- *
 * ge of nuclei in cycle 14 we've chosen to anneal on; e.g. if the most    *
 * anterior nucleus in cycle 14 has an odd lineage number, we need to push *
 * its sibling off the anterior limit when we divide after cycle 13; or in *
 * other words: our most anterior cell in cycle 14 is the posterior daugh- *
 * ter of our most anterior cell in cycle 13; therefore, we need to 'loose'*
 * its anterior sibling; the same applies for the posterior end, where we  *
 * want to loose the most posterior daughter cell if we don't need it any  *
 * more at the later cycle                                                 *
 *                                                                         *
 * This is implemented as follows below: lin is the lineage number of the  *
 * most anterior cell of the next cell cycle; if it's odd numbered, we'll  *
 * push off the most anterior daughter by subtracting the number of genes  * 
 * from the daughter indices (i.e. we shift the solution array for the     *
 * next cycle posteriorly by one nucleus); if on the other hand, the most  *
 * posterior nucleus lies outside our new array, we just forget about it   */

/*
    else if ( what2do[i] & DIVIDE ) {

      lin = GetStartLin(solution.array[i+1].time, &(inp->zyg.times)); 
      for (j=0; j < solution.array[i].state.size; j++) {

    k  = j % defs.ngenes;     // k: index of gene k in current nucleus
    ap = j / defs.ngenes;      // ap: rel. nucleus position on AP axis

// evaluate ii: index of anterior daughter nucleus

    if ( lin % 2 )
      ii = 2 * ap * defs.ngenes + k - defs.ngenes;
    else 
      ii = 2 * ap * defs.ngenes + k;

// skip the first most anterior daughter nucleus in case lin is odd

        if ( ii >= 0 )
      solution.array[i+1].state.array[ii] = 
        solution.array[i].state.array[j]; 
	
// the second daughter only exists if it is still within the region

    if ( ii + defs.ngenes < solution.array[i+1].state.size )
      solution.array[i+1].state.array[ii+defs.ngenes] =  
        solution.array[i].state.array[j]; 
      }

      if ( debug )
    fprintf(slog, "Blastoderm: nuclear division %d -> %d nuclei.\n",
        solution.array[i].state.size / defs.ngenes,
        solution.array[i+1].state.size / defs.ngenes);

    }
 */

/* In case we have to PROPAGATE the differential equeations, we call the   *
 * solver; you have to make sure that the appropriate rule has been set in *
 * zygotic.c (use SetRule(), rules are MITOSIS or INTERPHASE), otherwise   *
 * you will propagate the wrong equations; the solver needs pointers to    *
 * arrays of concentration for start and end time as well as those start   *
 * and end times themselves; all solvers need a suggested stepsize, adap-  *
 * tive stepsize solvers need the accuracy argument; we also need the accu-*
 * racy argument for global stepsize control (if ever implemented); lastly *
 * we need to tell the solver how big input and output arrays are          */

/*
    else if ( what2do[i] & PROPAGATE ){
      
      (*ps)(solution.array[i].state.array,
        solution.array[i+1].state.array,
        solution.array[i].time, 
        solution.array[i+1].time, 
        stephint, accuracy,
        solution.array[i].state.size,
        slog);
      
      int ns=solution.array[i].state.size;
      if ( !(jac = (double **)calloc(ns, sizeof(double *))))
    error("BaDe: error allocating jac.\n");
      
      for (ij=0; ij<ns; ij++)
    if ( !(jac[ij] = (double *)calloc(ns, sizeof(double))))
      error("BaDe: error allocating jac's second dimension.\n");
      
      
      JacobnOrig(solution.array[i+1].time,solution.array[i+1].state.array,dfdt,jac,ns);
      for (ij=0;ij<ns;ij++){
    for (ijj=0;ijj<ns;ijj++){
      cc=jindex+(ns*ij)+ijj;
      finalJac[cc]=jac[ij][ijj];
    }
      }
    }
    
// unknown op? -> error!

    else
      error("op was %d!?", what2do[i]);
  }
  
// After having calculated the solution, free the mutant parameter structs and return the result 
 
  FreeMutant();           
  return finalJac;        

}

 */

/**  FreeSolution: frees memory of the solution structure created by 
 *                 Blastoderm() or gut functions                           
 */
void
FreeSolution( NArrPtr * solution ) {
    int i;                      /* loop counter */
    if( solution != NULL ) {
        for( i = 0; i < solution->size; i++ ) {
            free( solution->array[i].state.array );
        }
        free( solution->array );
    }
}

/**  ConvertAnswer: little function that gets rid of bias times, division 
 *                  times and such and only returns the times in the tab-  
 *                  times struct as its output; this is used to produce    
 *                  unfold output that contains only the requested times;  
 *                  it also makes sure that we return the right solution   
 *                  at cell division, i.e. the solution right *after* the  
 *                  cells have divided                                     
 */
NArrPtr
ConvertAnswer( NArrPtr answer, DArrPtr tabtimes ) {
    int i;                      /* loop counter for answer */
    int outindex = 0;           /* loop counter for tabtimes */
    double nexttime;            /* used below to avoid memory errors */
    NArrPtr outtab;             /* answer to be returned */

    /* allocate the outtab array */

    outtab.size = tabtimes.size;
    outtab.array = ( NucState * ) calloc( outtab.size, sizeof( NucState ) );

    /* the following loop goes through all answer elements and copies only     */
    /* those which correspond to a tabtime into the outtab struct              */

    for( i = 0; ( i < answer.size ) && ( outindex < tabtimes.size ); i++ ) {

        /* this kludge makes sure that we don't bomb below when there is no next   *
         * time anymore                                                            */
        
        if( i == ( answer.size - 1 ) )
            nexttime = 999999999;
        else                    /* if statement */
            nexttime = answer.array[i + 1].time;

        /* this if makes sure that we only return the requested times and that     *
         * we return the solution *after* the cells have divided for times that    *
         * correspond exactly to the cell division time (i.e. the end of mitosis); *
         * note that times before and after cell divisions are separated by EPIS-  *
         * LON and that BIG_EPSILON is always (much) bigger than EPSILON (but      *
         * still small enough for all practical purposes...)                       */
        //printf("answertime = %lg, tabtimestime = %lg, nextime = %lg, diff<bige %lg, diff>bige %lg\n",answer.array[i].time, tabtimes.array[outindex], nexttime, (fabs(answer.array[i].time - tabtimes.array[outindex])), (fabs(answer.array[i].time - nexttime)));
        if( ( fabs( answer.array[i].time - tabtimes.array[outindex] )
              < BIG_EPSILON )
            && ( fabs( answer.array[i].time - nexttime )
                 > BIG_EPSILON ) ) {
            outtab.array[outindex].time = answer.array[i].time;
            outtab.array[outindex].state.size = answer.array[i].state.size;
            outtab.array[outindex].state.array = ( double * ) calloc( answer.array[i].state.size, sizeof( double ) );
            memcpy( outtab.array[outindex].state.array, answer.array[i].state.array, answer.array[i].state.size * sizeof( double ) );
            outindex++;
            //printf("outtab.array[%d].time=%lg size=%d, state0=%lg\n", outindex, outtab.array[outindex].time, outtab.array[outindex].state.size, outtab.array[outindex].state.array[0]);
        }
    }
    return outtab;
}



/*** THE FOLLOWING FUNCTIONS ARE FOR HANDLING TLIST, a linked list used to *
 *   initialize the structure that tells the solver for which time there's *
 *   data, how many nuclei there are and what to do.                       *
 ***************************************************************************/

/**  InitTList: initializes TList and adds first (t=0) and last 
 *              (t=gastrulation) element of Tlist.                         
 */
TList *
InitTList( Zygote * zyg, int *nnucs ) {
    TList *first, *last;        /* first and last element of TList */

    /* initialize the first element (t=0) */


    first = ( TList * ) malloc( sizeof( TList ) );
    first->time = 0. + EPSILON;
    first->n = zyg->defs.ngenes * GetNNucs( &( zyg->defs ), nnucs, 0, &( zyg->times ) );
    first->op = PROPAGATE;

    /* initialize the last element (t=gast_time) */

    last = ( TList * ) malloc( sizeof( TList ) );
    if( !( last->time = zyg->times.gast_time ) )
        error( "InitTList: error getting gastrulation time" );
    last->n = zyg->defs.ngenes * GetNNucs( &( zyg->defs ), nnucs, last->time, &( zyg->times ) );
    last->op = NO_OP;

    /* link the two elements and return the list */

    first->next = last;
    last->next = NULL;

    return first;
}

/**  InsertTList: takes pointer to first element of TList plus time and 
 *                desired op for new TList element and inserts a new TList 
 *                element for time at the appropriate place within the     
 *                linked list. This function returns a pointer to the      
 *                first element of the TList if insertion worked out fine. 
 */
TList *
InsertTList( Zygote * zyg, TList * first, double time, int op ) {
    TList *current;             /* used to step through TList */
    TList *new;                 /* used to allocate memory for new element */

    int n;                      /* how many nuclei at 'time' */


    n = zyg->defs.ngenes * GetNNucs( &( zyg->defs ), zyg->nnucs, time, &( zyg->times ) );       /* get number of nuclei for 'time' */

    if( first == NULL )
        error( "InsertTList: TList pointer is NULL at time %f!", time );

    /* the following loop steps through the linked list and places the new     *
     * element at the appropriate position depending on its time               */

    current = first;
    do {

        if( fabs( time - current->time ) < BIG_EPSILON ) {

            if( current->n == n ) {     /* if new time really close to an existing time point and if no cell division happened... */
                if( current->op == MITOTATE ) { /* if we are skipping epsilon after MITOSIS start */
                    new = ( TList * ) malloc( sizeof( TList ) );
                    new->time = time;   /* -> allocate new element for right after */
                    new->n = n; /* cell division has occurred */
                    new->op = op;
                    new->next = current->next;
                    current->next = new;
                    return first;
                }
                current->op |= op;      /* -> just add new op to existing one - from the resulting number we can obtain */
                return first;   /* uniquely which operations belongs to that timepoint */
            } else if( current->n < n ) {       /* if, on the other hand cell div HAS happened */
                new = ( TList * ) malloc( sizeof( TList ) );
                new->time = time;       /* -> allocate new element for right after */
                new->n = n;     /* cell division has occurred */
                new->op = op;
                new->next = current->next;
                current->next = new;
                return first;
            } else {            /* a sudden reduction of nuclei will be hard to explain */
                error( "InsertTList: sudden reduction of nuclei at time %g!", current->time );
            }

        } else if( time > current->time ) {     /* new time > than current time */
            if( fabs( time - current->next->time ) < BIG_EPSILON ) {    /* is it */
                current = current->next;        /* really close to the next time point? *//* -> go there! */
            } else if( time > current->next->time ) {   /* or if time is >> */
                current = current->next;        /* than the next time point *//* -> go there too! */
            } else if( time < current->next->time ) {   /* but if time is < */
                new = ( TList * ) malloc( sizeof( TList ) );    /* the next time point */
                new->time = time;       /* -> allocate new element */
                new->n = n;
                new->op = op;
                new->next = current->next;
                current->next = new;
                return first;
            } else {            /* if all the above don't apply, there's something wrong */
                error( "InsertTList: impossible error at time %g!", current->time );
            }
        } else {                /* if time < current time, there's something wrong too */
            error( "InsterTList: we missed our exit at time %g %g!", current->time, time );
        }
    } while( current != NULL );
    return current;
}

/**  CountEntries: counts how many entries we have in a TList */
int
CountEntries( TList * first ) {
    int n = 0;
    while( first != NULL ) {
        n++;
        first = first->next;
    }
    return n;
}

/**  FreeTList: frees the memory of a TList */
void
FreeTList( TList * first ) {
    if( first->next != NULL )
        FreeTList( first->next );

    free( first );
}

double *
GetFactDiscons( int *sss, FactDiscons fd ) {
    *sss = fd.fact_discons_size;
    return fd.fact_discons;

}

/** Dat2NArrPtr: Converts a data table into NarrPtr and also returns the index when
the number of nuclei is maximum */
NArrPtr
Dat2NArrPtr( DataTable * table, int *maxind ) {

    int i, j, max;
    NArrPtr input;

    max = 0;
    *maxind = 0;
    input.size = table->size;
    input.array = ( NucState * ) calloc( input.size, sizeof( NucState ) );

    for( i = 0; i < table->size; i++ ) {

        // printf("Heloo: %d %f\n", table->record[i].size, table->record[i].time);

        input.array[i].time = table->record[i].time + EPSILON;
        input.array[i].state.size = table->record[i].size;

        if( max < input.array[i].state.size ) {
            max = input.array[i].state.size;
            *maxind = i;
        }

        input.array[i].state.array = ( double * )
            calloc( input.array[i].state.size, sizeof( double ) );
        for( j = 0; j < table->record[i].size; j++ )
            input.array[i].state.array[j] = table->record[i].array[j].conc;

    }

    return input;
}

/** DoInterp: Sets up the interpolation functions for the lineage that has most
 * nuclei, these functions will be used later in conjunction with Go_Backward to 
 * return history for particular times 
 */
void
DoInterp( DataTable * interp_dat, InterpObject * interp_res, int num_genes, Zygote * zyg ) {
    double *z;
    int i, j, kk;
    NArrPtr Nptrfacts;
    int maxind;
    double *x, *y;
    int accepted, currsize;
    gsl_spline *temp_spline;
    gsl_interp_accel *temp_acc;
    Nptrfacts = Dat2NArrPtr( interp_dat, &maxind );
    //printf("MAXIND %d\n", maxind);

    interp_res->maxsize = Nptrfacts.array[maxind].state.size;
    interp_res->maxtime = Nptrfacts.array[maxind].time;

    //printf("maxsize:%d, maxtime:%f\n",interp_res->maxsize,interp_res->maxtime);
    /*  PrintBlastoderm(stdout, Nptrfacts, "Wloo", 2, num_genes); */
    interp_res->func.size = Nptrfacts.size;
    interp_res->func.array = ( NucState * ) calloc( interp_res->func.size, sizeof( NucState ) );
    interp_res->slope.size = Nptrfacts.size;
    interp_res->slope.array = ( NucState * ) calloc( interp_res->slope.size, sizeof( NucState ) );
    interp_res->fact_discons_size = Nptrfacts.size;
    interp_res->fact_discons = ( double * ) calloc( interp_res->fact_discons_size, sizeof( double ) );

    //test
    /*for (i = 0; i < interp_res->func.size; i++) {
       printf("times[%d]: intres_time=%lg, nptrfacts_time=%lg\n", i, interp_res->func.array[i].time, Nptrfacts.array[i].time);
       } */
    //test

    for( i = 0; i < interp_res->func.size; i++ ) {

        interp_res->func.array[i].time = Nptrfacts.array[maxind].time;
        interp_res->slope.array[i].time = Nptrfacts.array[i].time;

        interp_res->fact_discons[i] = Nptrfacts.array[i].time;

        interp_res->func.array[i].state.size = Nptrfacts.array[maxind].state.size;
        interp_res->slope.array[i].state.size = Nptrfacts.array[maxind].state.size;

        interp_res->func.array[i].state.array = ( double * ) calloc( interp_res->func.array[i].state.size, sizeof( double ) );
        interp_res->slope.array[i].state.array = ( double * ) calloc( interp_res->slope.array[i].state.size, sizeof( double ) );

        //printf("DoInterp(): Going from: %d to %d\n",Nptrfacts.array[i].state.size, interp_res->func.array[i].state.size);



        if( interp_res->func.array[i].state.size >= Nptrfacts.array[i].state.size ) {
            //printf("DoInterp(): intres_time=%lg, nptrfacts_time=%lg\n", interp_res->func.array[i].time, Nptrfacts.array[i].time);
            Go_Forward( interp_res->func.array[i].state.array, Nptrfacts.array[i].state.array,
                        GetStartLinIndex( interp_res->func.array[i].time, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( Nptrfacts.array[i].time,
                                                                                                                               &( zyg->defs ),
                                                                                                                               &( zyg->times ) ), zyg,
                        num_genes );
        } else
            Go_Backward( interp_res->func.array[i].state.array, Nptrfacts.array[i].state.array,
                         GetStartLinIndex( interp_res->func.array[i].time, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( Nptrfacts.array[i].time,
                                                                                                                                &( zyg->defs ),
                                                                                                                                &( zyg->times ) ), zyg,
                         num_genes );

        interp_res->func.array[i].time = Nptrfacts.array[i].time;

    }
    for( i = 0; i < interp_res->maxsize; i++ ) {
        /* If the last data point is zero, make it equal to the last
         * non-zero value, so that the concentration at a time 
         * after the last non -1 is maintained at that value */

        if( interp_res->func.array[Nptrfacts.size - 1].state.array[i] == -1. ) {

            kk = Nptrfacts.size - 2;
            while( ( interp_res->func.array[kk].state.array[i] == -1. ) && ( kk >= 0 ) )
                kk--;
            if( kk == -1 ) {

                printf( "All data is -1!\n" );
                exit( 1 );

            }

            interp_res->func.array[Nptrfacts.size - 1].state.array[i] = interp_res->func.array[kk].state.array[i];

        }
        /* If the first value in a data set is -1, set it to 0. */

        if( interp_res->func.array[0].state.array[i] == -1. )
            interp_res->func.array[0].state.array[i] = 0.;

        /* Now weed-out the -1s so you can interpolate linearly
         * between time points for which you have data */

        x = ( double * ) calloc( Nptrfacts.size, sizeof( double ) );
        y = ( double * ) calloc( Nptrfacts.size, sizeof( double ) );
        z = ( double * ) calloc( 1, sizeof( double ) );

        accepted = 0;
        currsize = Nptrfacts.size;

        for( j = 0; j < Nptrfacts.size; j++ ) {

            if( interp_res->func.array[j].state.array[i] != -1. ) {

                x[accepted] = Nptrfacts.array[j].time;
                y[accepted] = interp_res->func.array[j].state.array[i];
                accepted++;

            } else {

                currsize--;
                x = realloc( x, currsize * sizeof( double ) );
                y = realloc( y, currsize * sizeof( double ) );

            }
        }
        temp_acc = gsl_interp_accel_alloc(  );
        temp_spline = gsl_spline_alloc( gsl_interp_linear, currsize );
        gsl_spline_init( temp_spline, x, y, currsize );
        free( x );
        free( y );
        for( j = 0; j < Nptrfacts.size; j++ ) {
            gsl_spline_eval_e( temp_spline, interp_res->slope.array[j].time, temp_acc, z );     //doesn't work without _e... mistery
            //interp_res->func.array[j].state.array[i] = gsl_spline_eval(temp_spline, interp_res->slope.array[j].time, temp_acc);
            interp_res->func.array[j].state.array[i] = z[0];
        }

        free( z );
        for( j = 0; j < Nptrfacts.size; j++ ) {
            if( j < Nptrfacts.size - 1 )
                interp_res->slope.array[j].state.array[i] = ( interp_res->func.array[j + 1].state.array[i] - interp_res->func.array[j].state.array[i] ) /
                    ( interp_res->slope.array[j + 1].time - interp_res->slope.array[j].time );
            else
                interp_res->slope.array[j].state.array[i] = 0.;
        }
        gsl_interp_accel_free( temp_acc );
        gsl_spline_free( temp_spline );
    }
    FreeSolution( &Nptrfacts );
    return;
}

/*void TestInterp( int num_genes, int type ) {

    NArrPtr interp_out;
    double t;
	
    double ttee[ 10 ] = { 0.0, 10.550, 24.225, 30.475, 36.725, 42.975, 49.225,
        55.475, 61.725, 67.975};
    double ttee1[ 19 ] = { 0.0, 5.0, 10.550, 20.0, 24.225, 27.0, 30.475, 32.0,
        36.725, 40.0, 42.975, 49.0, 49.225, 52.225, 55.475, 58.0, 61.725, 
        65.33, 67.975 };
    
    double ttee2[3]={-30.0, -20.0, -5.0};
    int i;

    interp_out.size = 3;
    interp_out.array = (NucState *) calloc(interp_out.size, sizeof(NucState));
    t = 0.0+EPSILON;
    for (i=0; i < interp_out.size; i++)
    {
        t = ttee2[i]+EPSILON;
        t += 5.0;
        interp_out.array[i].time=0.0;
        interp_out.array[i].state.size=Index2NNuc(GetStartLinIndex(t, &(zyg->defs), &(zyg->times)))*num_genes;
        interp_out.array[i].state.array=(double *) calloc(interp_out.array[i].state.size, sizeof(double));

        if (!type)
                    History(t, t, interp_out.array[i].state.array, interp_out.array[i].state.size);
        else 
                    ExternalInputs(t, t, interp_out.array[i].state.array, interp_out.array[i].state.size);
    }	

	
    PrintBlastoderm(stdout, interp_out, "gloo", 2, num_genes);
    FreeSolution(interp_out);

}*/

/** SetFactDiscons: make one object from History and ExternalInputs */
FactDiscons
SetFactDiscons( InterpObject * hist_interp_object, InterpObject * extinp_interp_object ) {

    int i;

    FactDiscons all_fact_discons;

    all_fact_discons.fact_discons_size = hist_interp_object->fact_discons_size + extinp_interp_object->fact_discons_size;

    all_fact_discons.fact_discons = ( double * ) calloc( all_fact_discons.fact_discons_size, sizeof( double ) );

    for( i = 0; i < hist_interp_object->fact_discons_size; i++ )
        all_fact_discons.fact_discons[i] = hist_interp_object->fact_discons[i];

    for( i = 0; i < extinp_interp_object->fact_discons_size; i++ )
        all_fact_discons.fact_discons[i + hist_interp_object->fact_discons_size] = extinp_interp_object->fact_discons[i];

    qsort( ( void * ) all_fact_discons.fact_discons, all_fact_discons.fact_discons_size, sizeof( double ),
           ( int ( * )( const void *, const void * ) ) compare );
    return all_fact_discons;
}

/** FreeFactDiscons: free FactDiscons object. FactDiscons object is made from History and ExternalInputs */
void
FreeFactDiscons( double *fact_discons ) {
    free( fact_discons );
}

void
History( double t, double t_size, double *yd, int n, InterpObject hist_interp_object, int ngenes, Zygote * zyg ) {

    int j, k;
    double *blug;
    double t_interp, t_diff;

    /*  printf("Going from %d to %d, time:%f time for size:%f\n",maxsize,n,t,t_size); */

    blug = ( double * ) calloc( hist_interp_object.maxsize, sizeof( double ) );

    k = -1;
    do {
        k++;

    } while( ( k < hist_interp_object.slope.size ) && ( t > hist_interp_object.slope.array[k].time ) );
    if( k == 0 )
        t_interp = hist_interp_object.slope.array[k].time;
    else {

        k--;
        t_interp = t;

    }

    t_diff = t_interp - hist_interp_object.slope.array[k].time;

    for( j = 0; j < hist_interp_object.maxsize; j++ ) {

        blug[j] = hist_interp_object.func.array[k].state.array[j] + hist_interp_object.slope.array[k].state.array[j] * t_diff;

    }

    if( n >= hist_interp_object.maxsize ) {
        Go_Forward( yd, blug, GetStartLinIndex( t_size, &( zyg->defs ), &( zyg->times ) ),
                    GetStartLinIndex( hist_interp_object.maxtime, &( zyg->defs ), &( zyg->times ) ), zyg, ngenes );
    } else
        Go_Backward( yd, blug, GetStartLinIndex( t_size, &( zyg->defs ), &( zyg->times ) ),
                     GetStartLinIndex( hist_interp_object.maxtime, &( zyg->defs ), &( zyg->times ) ), zyg, ngenes );
    free( blug );

    return;
}

void
ExternalInputs( double t, double t_size, double *yd, int n, InterpObject extinp_interp_object, int egenes, Zygote * zyg ) {

    int j, k;
    double *blug;
    double t_interp, t_diff;

    /*  printf("Going from %d to %d, time:%f time for size:%f\n",maxsize,n,t,t_size); */

    blug = ( double * ) calloc( extinp_interp_object.maxsize, sizeof( double ) );

    k = -1;

    do {
        k++;
    } while( ( k < extinp_interp_object.slope.size ) && ( t > extinp_interp_object.slope.array[k].time ) );

    if( k == 0 )
        t_interp = extinp_interp_object.slope.array[k].time;
    else {
        k--;
        t_interp = t;
    }

    t_diff = t_interp - extinp_interp_object.slope.array[k].time;

    for( j = 0; j < extinp_interp_object.maxsize; j++ ) {
        blug[j] = extinp_interp_object.func.array[k].state.array[j] + extinp_interp_object.slope.array[k].state.array[j] * t_diff;
    }
    if( n >= extinp_interp_object.maxsize ) {
        Go_Forward( yd, blug, GetStartLinIndex( t_size, &( zyg->defs ), &( zyg->times ) ),
                    GetStartLinIndex( extinp_interp_object.maxtime, &( zyg->defs ), &( zyg->times ) ), zyg, egenes );
    } else
        Go_Backward( yd, blug, GetStartLinIndex( t_size, &( zyg->defs ), &( zyg->times ) ),
                     GetStartLinIndex( extinp_interp_object.maxtime, &( zyg->defs ), &( zyg->times ) ), zyg, egenes );
    free( blug );

    return;
}

void
FreeInterpObject( InterpObject * interp_obj ) {
    FreeSolution( ( &( interp_obj->func ) ) );
    FreeSolution( ( &( interp_obj->slope ) ) );

    if( interp_obj->fact_discons )
        free( interp_obj->fact_discons );
}

void
Go_Forward( double *output, double *input, int output_ind, int input_ind, Zygote * zyg, int num_genes ) {

    double *y;
    int output_lin, size, newsize;
    int k, ap, j, ii;

    y = NULL;
    size = 0;

    output_lin = Index2StartLin( output_ind, zyg->full_lin_start );
    newsize = Index2NNuc( output_ind, zyg->full_nnucs ) * num_genes;

    //printf("output lineage start, indices:%d, %d, %d\n", output_lin, input_ind, output_ind);

    if( output_ind < input_ind - 1 ) {
        size = Index2NNuc( output_ind + 1, zyg->full_nnucs ) * num_genes;
        y = ( double * ) calloc( size, sizeof( double ) );
        /*              printf("Passing on to another fwd with targets %d %d %d\n",size,output_ind+1,input_ind); */
        //printf("Go_FW\n");
        Go_Forward( y, input, output_ind + 1, input_ind, zyg, num_genes );
    } else if( output_ind == input_ind - 1 ) {
        size = Index2NNuc( input_ind, zyg->full_nnucs ) * num_genes;
        y = ( double * ) calloc( size, sizeof( double ) );
        /*              printf("Goin' to do the tranfer:%d %d\n",size,newsize); */
        y = memcpy( y, input, size * sizeof( double ) );
    } else if( output_ind == input_ind ) {
        output = memcpy( output, input, newsize * sizeof( double ) );
        return;
    } else {
        //printf("====\n");
        //printf("input: input_ind=%d, full_nnucs=%d\n", input_ind, zyg->full_nnucs[input_ind]);
        //printf("output: output_ind=%d, full_nnucs=%d\n", output_ind, zyg->full_nnucs[output_ind]);

        error( "Go_Forward(): You are trying to go from nnucs %d to %d!", Index2NNuc( input_ind, zyg->full_nnucs ), Index2NNuc( output_ind, zyg->full_nnucs ) );
    }

    for( j = 0; j < size; j++ ) {

        k = j % num_genes;      /* k: index of gene k in current nucleus */
        ap = j / num_genes;     /* ap: rel. nucleus position on AP axis */

        /* evaluate ii: index of anterior daughter nucleus */

        if( output_lin % 2 )
            ii = 2 * ap * num_genes + k - num_genes;
        else
            ii = 2 * ap * num_genes + k;

        /* skip the first most anterior daughter nucleus in case output_lin is odd */

        if( ii >= 0 )
            output[ii] = y[j];

        /* the second daughter only exists if it is still within the region */

        if( ii + num_genes < newsize )
            output[ii + num_genes] = y[j];
    }

    free( y );

    return;
}

void
Go_Backward( double *output, double *input, int output_ind, int input_ind, Zygote * zyg, int num_genes ) {

    double *y;
    int output_lin, input_lin, size, newsize;
    int k, ap, j, ii;

    y = NULL;
    size = 0;
    
    output_lin = Index2StartLin( output_ind, zyg->full_lin_start );
    input_lin = Index2StartLin( input_ind, zyg->full_lin_start );
    newsize = Index2NNuc( output_ind, zyg->full_nnucs ) * num_genes;

    /*  printf("output lineage start, indices:%d, %d, %d\n",output_lin,output_ind,input_ind); */

    if( output_ind > input_ind + 1 ) {
        size = Index2NNuc( output_ind - 1, zyg->full_nnucs ) * num_genes;
        y = ( double * ) calloc( size, sizeof( double ) );
        /*              printf("Passing on to another bkd with targets %d %d %d\n",size,output_ind-1,input_ind); */
        Go_Backward( y, input, output_ind - 1, input_ind, zyg, num_genes );
        input_lin = Index2StartLin( output_ind - 1, zyg->full_lin_start );
    } else if( output_ind == input_ind + 1 ) {
        size = Index2NNuc( input_ind, zyg->full_nnucs ) * num_genes;
        y = ( double * ) calloc( size, sizeof( double ) );
        /*              printf("Goin' to do the tranfer\n"); */
        memcpy( y, input, size * sizeof( double ) );
    } else if( output_ind == input_ind ) {
        memcpy( output, input, newsize * sizeof( double ) );
        return;
    } else
        error( "Go_Backward(): You are trying to go from nnucs %d to %d!", Index2NNuc( input_ind, zyg->full_nnucs ),
               Index2NNuc( output_ind, zyg->full_nnucs ) );

    for( j = 0; j < newsize; j++ ) {

        k = j % num_genes;      /* k: index of gene k in current nucleus */
        ap = j / num_genes;     /* ap: rel. nucleus position on AP axis */

        /* evaluate ii: index of anterior daughter nucleus */

        if( input_lin % 2 )
            ii = 2 * ap * num_genes + k - num_genes;
        else
            ii = 2 * ap * num_genes + k;

        /* skip the first most anterior daughter nucleus in case output_lin is odd */

        if( ii < 0 )
            output[j] = y[ii + num_genes];

        if( ( ii >= 0 ) && ( ii + num_genes < size ) )
            output[j] = .5 * ( y[ii] + y[ii + num_genes] );

        if( ii + num_genes >= size )
            output[j] = y[ii];

    }

    free( y );

    return;
}
