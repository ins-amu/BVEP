/**
 * @file zygotic.c                                                   
 * @authors JR, modified by Yoginho, 
 *          additional g(u)'s by Yousong Wang in Feb 2002,    
 *          guts code by Manu in June/July 2002
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Implementation of functions that deal with the right hand side 
 * of the ODEs.
 * 
 * We have an initializing function called          
 * InitZygote that reads parameters, defs data and search space  
 * limits from a data file and installs a couple of things like  
 * the solver and temporary arrays for the dvdt function.        
 * Then we have the DvdtOrig's, which represent the inner loop   
 * of the model. Each of these derivative functions calculates   
 * the derivatives of a particular right-hand-side at a given    
 * time point. They are called by the solver and need to know if 
 * we're in mitosis or interphase, since the equations differ    
 * in each case.                                                 
 * Last but not least, we have a few mutator functions in this   
 * file. They change the local lparm EqParm struct.              
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "error.h"              /* for error handling */
#include "maternal.h"           /* need this for BArrPtr (for the Bicoid structure) */
#include "solvers.h"            /* for p_deriv */
#include "integrate.h"          /* for ExternalInputs */
#include "zygotic.h"            /* obviously */
#include "fly_io.h"             /* i/o of parameters and data */
#include "ioTools.h"


//test
#include <sys/time.h>           /* for time calculation */
struct timeval start, end;      /* time_t is defined on <time.h> and <sys/types.h> as long */
//endtest

// following is a superfast vector square root function available for DEC
/*#ifdef ALPHA_DU
extern void vsqrt_();
extern void vexp_();
#endif*/


/*** CONSTANTS *************************************************************/

/* these are the propagation rules for dvdt_orig */
const int INTERPHASE = 0;
const int MITOSIS = 1;

/* diffusion coefficients */
static double *D;

/* following arrays are used by DvdtOrig */
double *vinput;                 /* vinput, bot2 and bot are used for */
double *bot2, *bot;             /* storing intermediate stuff for vector */


GFunc gofu;
void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * );
void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * );
void ( *dd ) ( double *, double **, double, double *, int, SolverInput *, Input * );

/*** INITIALIZATION FUNCTIONS **********************************************/

/** InitZygote: makes pm and pd visible to all functions in zygotic.c and 
 *               reads EqParms and TheProblem; it then initializes bicoid  
 *               and bias (including BTimes) in maternal.c; lastly, it     
 *               allocates memory for structures used by the derivative    
 *               functions                                                 
 */
Zygote
InitZygote( FILE * fp, void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * ),
            void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * ), Input * inp, char *section_title ) {

    /***************************************************************************
     * - defs:     is the global problem struct; must be visible to scoring    *
     *             functions, to Blastoderm in integrate.c, to stuff in mater- *
     *             nal.c and to functions in zygotic.c (need any more justifi- *
     *             cation for this to be global?)                     *
     ***************************************************************************/

    Zygote zyg;
    /* install the dvdt function: makes pd global (solvers need to access it) */

    p_deriv = pd;
    p_jacobn = pj;
    d_deriv = dd;               //delayed derivative

    zyg.ndp = 0;
    zyg.nalleles = 0;
    zyg.full_lin_start = NULL;
    zyg.full_nnucs = NULL;


    /* read equation parameters and the problem */
    zyg.defs = ReadTheProblem( fp );
    D = ( double * ) calloc( zyg.defs.ngenes, sizeof( double ) );       /* contains info about diffusion sched. */

    /* install bicoid and bias and nnucs in maternal.c */
    zyg.bcdtype = InitBicoid( fp, &zyg );

    zyg.bias = InitBias( fp, &zyg );
    zyg.nnucs = InitNNucs( &zyg );      //nnucs for every cellcycle
    zyg.times = ReadDivTimes( fp, zyg.defs );   //Look for times section in the input file, if not there take ones from maternal.c //Under Construction

    // read initial state
    zyg.parm = ReadParameters( fp, zyg.defs, section_title );

    return zyg;
}

/*** CLEANUP FUNCTIONS *****************************************************/

/** FreeZygote: frees memory for D */
void
FreeZygote( void ) {
    free( D );
}

/** FreeMutant: frees mutated parameter struct */
void
FreeMutant( EqParms lparm ) {
    free( lparm.R );
    free( lparm.T );
    free( lparm.E );
    free( lparm.m );
    free( lparm.h );
    free( lparm.d );
    free( lparm.lambda );
    free( lparm.tau );
}

/*** DERIVATIVE FUNCTIONS **************************************************/

/*** Derivative functions: *************************************************
 *                                                                         *
 *   These functions calculate derivatives for the solver function based   *
 *   on the state variables in array v for the time t; the output is       *
 *   written into array vdot; n describes the size of both v and vdot      *
 *                                                                         *
 *   There are different g(u) functions that can be used in each deriva-   *
 *   tive function (see comments in the code below). They are set by the   *
 *   command line option -g.                                               * 
 *                                                                         *
 *   There are two rules for the derivative function, one for INTERPHASE   *
 *   and one for MITOSIS. The difference between the two is that only de-  *
 *   cay and diffusion happen during mitosis and the regulation term is    *
 *   only included during INTERPHASE. The rule gets set in Blastoderm      *
 *   (using SetRule) to ensure that the derivative function is in the      *
 *   right rule.                                                           *
 *                                                                         *
 *   JR: We get rid of 3 of 4 divisions in the inner loop by nesting. We   *
 *   get rid of an if by adding one iteration of the loop before and one   *
 *   after for the diffusion (wall) boundary conditions.                   *
 *                                                                         *
 ***************************************************************************/

/** DvdtOrig: the original derivative function; implements the equations 
 *             as published in Reinitz & Sharp (1995), Mech Dev 49, 133-58 
 *             plus different g(u) functions as used by Yousong Wang in    
 *             spring 2002.                                                
 */
void
DvdtOrig( double *v, double t, double *vdot, int n, SolverInput * si, Input * inp ) {

    double vinput1 = 0;
    int m;                      /* number of nuclei */
    int ap;                     /* nuclear index on AP axis [0,1,...,m-1] */
    int i, j;                   /* local loop counters */
    int k;                      /* index of gene k in a specific nucleus */
    int base, base1;            /* index of first gene in a specific nucleus */
    int *l_rule;                /* for autonomous implementation */
#ifdef ALPHA_DU
    int incx = 1;               /* increment step size for vsqrt input array */
    int incy = 1;               /* increment step size for vsqrt output array */
#endif

    static int num_nucs = 0;    /* store the number of nucs for next step */
    static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
    static DArrPtr bcd;         /* pointer to appropriate bicoid struct */
    double *v_ext;              /* array to hold the external input
                                   concentrations at time t */
    int allele = si->genindex;

    vinput = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot2 = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );


    /* get D parameters and bicoid gradient according to cleavage cycle */
    /* n is the total number of genes */
    m = n / inp->zyg.defs.ngenes;       /* m is the number of nuclei */
    /* inp->zyg.defs.ngenes is the number of
     * genes per nucleus */
    if( m != num_nucs ) {       /* time-varying quantities only vary by ccycle */

        GetD( t, inp->lparm.d, D, &( inp->zyg ) );
        /* get diff coefficients, according to diff schedule */
        if( num_nucs > m )      /* started a new iteration in score */
            bcd_index = 0;      /* -> start all over again! */
        num_nucs = m;           /* store # of nucs for next step */
        bcd = GetBicoid( t, allele, inp->zyg.bcdtype, &( inp->zyg ) );  /* get bicoid gradient */
        if( bcd.size != num_nucs )
            error( "DvdtOrig: %d nuclei don't match Bicoid!", num_nucs );
        bcd_index++;            /* store index for next bicoid gradient */
    }
    l_rule = ( int * ) calloc( inp->zyg.defs.ngenes, sizeof( int ) );
    for( i = 0; i < inp->zyg.defs.ngenes; i++ )
        l_rule[i] = !( Theta( t, &( inp->zyg ) ) );     // Theta(u) = false while interphase
    /* l_rule is zero during mitosis, in order
     * to put to zero the regulation part of the
     * equation. Remember, no regulation during
     * mitosis */
    // Here we retrieve the external input concentrations into v_ext
    v_ext = ( double * ) calloc( m * inp->zyg.defs.egenes, sizeof( double ) );
    ExternalInputs( t, t, v_ext, m * inp->zyg.defs.egenes, inp->ext[allele], inp->zyg.defs.egenes, &( inp->zyg ) );     //here we assume that all the genotypes use the same external inps, so we take the first one -- ask Yogi 2
    /* This is how it works (by JR): 

       ap      nucleus position on ap axis
       base    index of first gene (0) in current nucleus
       k       index of gene k in current nucleus 

       Protein synthesis terms are calculated according to g(u)

       First we do loop for vinput contributions; vinput contains
       the u that goes into g(u)

       Then we do a separate loop or vector func for sqrt or exp

       Then we do vdot contributions (R and lambda part)

       Then we do the spatial part (Ds); we do a special case for 
       each end 

       Note, however, that before the real loop, we have to do a 
       special case for when rhere is one nuc, hence no diffusion

       These loops look a little funky 'cause we don't want any 
       divides'                                                       */


    /***************************************************************************
     *                                                                         *
     *        g(u) = 1/2 * ( u / sqrt(1 + u^2) + 1)                            *
     *                                                                         *
     ***************************************************************************/

    if( gofu == Sqrt ) {

        gettimeofday( &start, NULL );   //start the timer

        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;
                //printf("HK=%lg\n", inp->lparm.h[k]);
                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */
                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[base1 + j];
                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];
                bot2[i] = 1 + vinput1 * vinput1;
                vinput[i] = vinput1;
            }
        }

        gettimeofday( &end, NULL );     //start the timer
        //printf("beforehalfG %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

#ifdef ALPHA_DU
        vsqrt_( bot2, &incx, bot, &incy, &n );  /* superfast DEC vector function */
#else
        for( i = 0; i < n; i++ )        /* slow traditional style sqrt */
            bot[i] = sqrt( bot2[i] );
#endif
        gettimeofday( &end, NULL );     //start the timer
        //printf("halfG %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */
        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */
            register double vdot1, g1;
            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */
            register double vdot1, g1;
            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;          /* v[i-inp->zyg.defs.ngenes]-v[i] == 0 */
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }

            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {     /* last: posterior-most nucleus */
                k = i - base;   /* v[i+inp->zyg.defs.ngenes]-v[i] == 0 */
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

        gettimeofday( &end, NULL );     //start the timer
        //printf("endG %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));


        /***************************************************************************
         *                                                                         *
         *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
         *                                                                         *
         ***************************************************************************/

    } else if( gofu == Tanh ) {
        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];

                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = tanh( vinput[i] ) + 1;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = tanh( vinput[i] ) + 1;
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = tanh( vinput[i] ) + 1;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = tanh( vinput[i] ) + 1;
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }
        /***************************************************************************
         *                                                                         *
         *        g(u) = 1 / (1 + exp(-2u))                                        *
         *                                                                         *
         ***************************************************************************/

    } else if( gofu == Exp ) {
        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];

                vinput[i] = -2.0 * vinput1;
            }
        }

        /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
        vexp_( vinput, &incx, bot, &incy, &n ); /* superfast DEC vector function */
#else
        for( i = 0; i < n; i++ )        /* slow traditional style exp */
            bot[i] = exp( vinput[i] );
#endif

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 / ( 1 + bot[i] );
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 / ( 1 + bot[i] );
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 / ( 1 + bot[i] );
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 / ( 1 + bot[i] );
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 0 if u<0, 1 if u>=0 (Heaviside function)                  *
         *                                                                         *
         *        this makes the model quasi boolean and the equations locally     *
         *        linear for both u>0 and u<0                                      *
         ***************************************************************************/

    } else if( gofu == Hvs ) {
        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];
                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    if( vinput[i] >= 0. )
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                if( vinput[i] >= 0. )
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    if( vinput[i] >= 0. )
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                if( vinput[i] >= 0. )
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

    } else if( gofu == Kolja ) {
        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];

                vinput[i] = vinput1;    //u
            }
        }

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = vinput[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;

                vdot1 = -inp->lparm.lambda[k] * v[i];
                //printf("afterdecay vdot0=%lg vi=%lg\n", vdot1, v[i]);
                g1 = vinput[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }

            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = vinput[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }

            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = vinput[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

    } else
        error( "DvdtOrig: unknown g(u)" );

    /* during mitosis only diffusion and decay happen */
    free( l_rule );
    free( v_ext );
    free( vinput );
    free( bot2 );
    free( bot );
    return;
}



/*** JACOBIAN FUNCTION(S) **************************************************/

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
void
JacobnOrig( double t, double *v, double *dfdt, double **jac, int n, SolverInput * si, Input * inp ) {


    int m;                      /* number of nuclei */
    int ap;                     /* nuclear index on AP axis [0,1,...,m-1] */
    int i, j;                   /* local loop counters */
    int k, kk;                  /* index of gene k in a specific nucleus */
    int base;                   /* indices of 1st gene in a specific nucleus */
    int rule;                   /* propagation rule for JacobnOrig */
#ifdef ALPHA_DU
    int incx = 1;               /* increment step size for vsqrt input array */
    int incy = 1;               /* increment step size for vsqrt output array */
#endif

    static int num_nucs = 0;    /* store the number of nucs for next step */
    static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
    DArrPtr bcd = ( const struct DArrPtr ){ 0 }; 
                                /* pointer to appropriate bicoid struct */

    int allele = si->genindex;

    vinput = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot2 = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );

    /* get D parameters and bicoid gradient according to cleavage cycle */

    m = n / inp->zyg.defs.ngenes;       /* m is the number of nuclei */
    if( m != num_nucs ) {       /* time-varying quantities only vary by ccycle */
        GetD( t, inp->lparm.d, D, &( inp->zyg ) );
        /* get diff coefficients, according to diff schedule */
        if( num_nucs > m )      /* started a new iteration in score */
            bcd_index = 0;      /* -> start all over again! */
        num_nucs = m;           /* store # of nucs for next step */
        bcd = GetBicoid( t, allele, inp->zyg.bcdtype, &( inp->zyg ) );  /* get bicoid gradient */
        if( bcd.size != num_nucs )
            error( "JacobnOrig: %d nuclei don't match Bicoid!", num_nucs );
        bcd_index++;            /* store index for next bicoid gradient */
    }

    rule = GetRule( t, &( inp->zyg ) );

    /*** INTERPHASE rule *******************************************************/

    if( rule == INTERPHASE ) {

        register double vinput1 = 0;    /* used to calculate u */
        register double gdot1, vdot1;   /* used to calculate Xab's and Yab's */


        /***************************************************************************
         *                                                                         *
         *  g(u)  = 1/2 * ( u / sqrt(1 + u^2) + 1)                                 *
         *  g'(u) = 1/2 * ( 1 / ((1 + u^2)^3/2)) * T{ab}                           *
         *                                                                         *
         ***************************************************************************/

        if( gofu == Sqrt ) {

            /* this looks confusing, but it's actually quite simple: the two loops be- *
             * low are in reality just one that loops over the first dimension of the  *
             * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
             * long introductory comment above); ap keeps track of the nucleus number, * 
             * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

            for( base = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, ap++ ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;

                    vinput1 = inp->lparm.h[k];
                    vinput1 += inp->lparm.m[k] * bcd.array[ap];

                    for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                        vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];

                    bot2[i] = 1 + vinput1 * vinput1;
                }
            }

            /* now calculate sqrt(1+u^2); store it in bot[] */

#ifdef ALPHA_DU
            vsqrt_( bot2, &incx, bot, &incy, &n );      /* superfast DEC vector function */
#else
            for( i = 0; i < n; i++ )    /* slow traditional style sqrt */
                bot[i] = sqrt( bot2[i] );
#endif

            /* resume loop after vector sqrt above; we finish calculating g'(u) and    *
             * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;

                    gdot1 = 1 / ( bot[i] * bot2[i] );
                    gdot1 *= inp->lparm.R[k] * 0.5;

                    for( j = base; j < base + inp->zyg.defs.ngenes; j++ ) {

                        kk = j - base;

                        vdot1 = inp->lparm.T[( k * inp->zyg.defs.ngenes ) + kk] * gdot1;
                        if( k == kk ) {
                            if( n > inp->zyg.defs.ngenes ) {
                                if( base > 0 && base < n - inp->zyg.defs.ngenes ) {
                                    vdot1 -= 2. * D[k];
                                } else {
                                    vdot1 -= D[k];
                                }
                                vdot1 -= inp->lparm.lambda[k];
                            }
                        }
                        jac[i][j] = vdot1;

                    }

                    if( base > 0 )
                        jac[i][i - inp->zyg.defs.ngenes] = D[k];

                    if( base < n - inp->zyg.defs.ngenes )
                        jac[i][i + inp->zyg.defs.ngenes] = D[k];

                }
            }

            /***************************************************************************
             *                                                                         *
             * g(u)  = 1/2 * (tanh(u) + 1) )  or  g(u)  = 1 / ( 1 + e^(-2u))           *
             * g'(u) = Sech(u)^2 /2           or  g'(u) = 2e^(-2u) / (1 + e^(-2u))^2   *
             *                                                                         *
             ***************************************************************************
             *                                                                         *
             * These are actually the same function in different guises; we implement  *
             * Exp below since it's faster                                             *
             *                                                                         *
             ***************************************************************************/

        } else if( ( gofu == Tanh ) || ( gofu == Exp ) ) {

            /* this looks confusing, but it's actually quite simple: the two loops be- *
             * low are in reality just one that loops over the first dimension of the  *
             * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
             * long introductory comment above); ap keeps track of the nucleus number, * 
             * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

            for( base = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, ap++ ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;

                    vinput1 = inp->lparm.h[k];
                    vinput1 += inp->lparm.m[k] * bcd.array[ap];

                    for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                        vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * v[base + j];

                    bot[i] = -2.0 * vinput1;
                }
            }

            /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
            vexp_( bot, &incx, bot2, &incy, &n );       /* superfast DEC vector function */
#else
            for( i = 0; i < n; i++ )    /* slow traditional style exp */
                bot2[i] = exp( bot[i] );
#endif

            /* resume loop after vector exp above; we finish calculating g'(u) and     *
             * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;

                    gdot1 = 2. * bot2[i];
                    gdot1 /= ( 1. + bot2[i] ) * ( 1. + bot2[i] );
                    gdot1 *= inp->lparm.R[k];

                    for( j = base; j < base + inp->zyg.defs.ngenes; j++ ) {

                        kk = j - base;

                        vdot1 = inp->lparm.T[( k * inp->zyg.defs.ngenes ) + kk] * gdot1;
                        if( k == kk ) {
                            if( n > inp->zyg.defs.ngenes ) {
                                if( base > 0 && base < n - inp->zyg.defs.ngenes ) {
                                    vdot1 -= 2. * D[k];
                                } else {
                                    vdot1 -= D[k];
                                }
                                vdot1 -= inp->lparm.lambda[k];
                            }
                        }
                        jac[i][j] = vdot1;

                    }

                    if( base > 0 )
                        jac[i][i - inp->zyg.defs.ngenes] = D[k];

                    if( base < n - inp->zyg.defs.ngenes )
                        jac[i][i + inp->zyg.defs.ngenes] = D[k];

                }
            }

            /*** semi-implicit solvers are NOT allowed with heaviside g(u) *************/

        } else if( gofu == Hvs ) {
            error( "JacobnOrig: can't use semi-implicit solver on heaviside g(u)!" );

        } else {
            error( "JacobnOrig: unknown g(u) function!\n" );
        }

        /* during mitosis only diffusion and decay happen */

    } else if( rule == MITOSIS ) {

        register double vdot1;

        for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;

                for( j = base; j < base + inp->zyg.defs.ngenes; j++ ) {
                    kk = j - base;

                    if( k == kk ) {
                        vdot1 = -inp->lparm.lambda[k];
                        if( n < inp->zyg.defs.ngenes ) {
                            if( base > 0 && base < n - inp->zyg.defs.ngenes ) {
                                vdot1 -= 2. * D[k];
                            } else {
                                vdot1 -= D[k];
                            }
                        }
                    } else {
                        vdot1 = 0.;
                    }

                    jac[i][j] = vdot1;

                }

                if( base > 0 )
                    jac[i][i - inp->zyg.defs.ngenes] = D[k];

                if( base < n - inp->zyg.defs.ngenes )
                    jac[i][i + inp->zyg.defs.ngenes] = D[k];

            }
        }

    } else
        error( "JacobnOrig: Bad rule %i sent to JacobnOrig", rule );
    free( vinput );
    free( bot2 );
    free( bot );
    return;
}



/*** GUTS FUNCTIONS ********************************************************/

/** CalcGuts: calculates guts for genotpye 'gtype' using unfold output in 
 *             'table' and the guts string 'gutsdefs'; it returns the num- 
 *             ber of columns we'll need to print and the guts table in    
 *             'gtable'                                                    
 */
int
CalcGuts( int gindex, char *gtype, InterpObject * hist_interrp, InterpObject * extinp_interrp, NArrPtr table, NArrPtr * gtable, char *gutsdefs, Input * inp ) {
    int m;                      // number of nuclei at each time 
    int i;                      // just for the loops 
    int numguts;                // number of guts columns requested 
    int which;                  // for which gene to calculate guts 
    int rule;                   // propagation rule 
    unsigned long *gutcomps = NULL;     // bit flags for each word of ID string 

    NArrPtr gutsy;              // temporary place for guts, to be returned 

    char **w = NULL;            // w(ord) string array after parsing 
    char **dummy = NULL;        // beginning of w(ord) string array 
    char *targetptr;            // for parsing the gutsdefs, and a counter 

    SolverInput si;


    // allocate memory for parsed gut strings, then parse the string

    dummy = ( char ** ) calloc( MAX_RECORD, sizeof( char * ) );
    w = dummy;
    numguts = ( ParseString( gutsdefs, w ) - 1 );       // w contains array of words

    // if string is empty -> return zero

    if( numguts == -1 ) {
        while( *++w )
            free( *w );
        free( *dummy );
        free( dummy );
        return numguts + 1;
    }
    // if no guts specified but gene name in string -> error
    if( numguts == 0 )
        error( "CalcGuts: no guts specified for target gene %s", w[0] );

    /* allocate the gutcomps array and fill it with the bit flags correspon-   *
     * ding to the things that need to get calculated below; each word in an   *
     * ID string corresponds to one long int in the gutcomps array; targetptr  *
     * points to the name of the gene in the geneID string for which we cal-   *
     * culate guts                                                             */

    gutcomps = ( unsigned long * ) calloc( numguts, sizeof( unsigned long ) );
    if( !( targetptr = GetGutsComps( inp->zyg.defs.gene_ids, inp->zyg.defs.egene_ids, w, gutcomps ) ) )
        error( "CalcGuts: target gene (%s) is not defined in the geneID string", w[0] );

    // allocate gutsy array and set the right genotype in score.c & zygotic.c

    gutsy.size = table.size;
    gutsy.array = ( NucState * ) calloc( table.size, sizeof( NucState ) );
    /*InitDelaySolver(); */
    inp->his = hist_interrp;
    inp->ext = extinp_interrp;
    si.all_fact_discons = SetFactDiscons( inp->his, inp->ext );
    si.genindex = gindex;
    inp->lparm = Mutate( gtype, inp->zyg.parm, &( inp->zyg.defs ) );

    // which tells us which gene we calculate guts for

    which = targetptr - inp->zyg.defs.gene_ids;

    // free gut ID strings, we don't need them anymore since we have gutcomps

    while( *++w )
        free( *w );
    free( *dummy );
    free( dummy );

    // then start calculating guts

    for( i = 0; i < table.size; i++ ) {
        si.time = gutsy.array[i].time;

        gutsy.array[i].time = table.array[i].time;

        // calculate size of 2D array that holds guts and allocate memory

        m = table.array[i].state.size / inp->zyg.defs.ngenes;
        gutsy.array[i].state.size = numguts * m;
        gutsy.array[i].state.array = ( double * ) calloc( gutsy.array[i].state.size, sizeof( double ) );

        // set whether it is MITOSIS or INTERPHASE when calculating guts

        rule = GetRule( gutsy.array[i].time, &( inp->zyg ) );

        // call the function that calculates the internals of the RHS

        CalcRhs( table.array[i].state.array, gutsy.array[i].time, gutsy.array[i].state.array, table.array[i].state.size,
                 gutsy.array[i].state.size, numguts, which, gutcomps, &si, inp );

    }

    // clean up

    FreeMutant( inp->lparm );
    /*FreeDelaySolver(); */
    FreeFactDiscons( si.all_fact_discons.fact_discons );
    free( gutcomps );
    // return the guts and number of columns for PrintGuts()

    *gtable = gutsy;
    return numguts;
}

/** CalcRhs: calculates the components of the right-hand-side (RHS) of 
 *            the equatiion which make up guts (e.g. regulatory contribu-  
 *            tions of specific genes, diffusion or protein decay); this   
 *            function makes use of the derivative function and calculates 
 *            everything outside g(u) in reverse, i.e. starting from the   
 *            derivative and calculating the desired gut properties back-  
 *            wards from there; everything within g(u) is simply recon-    
 *            structed from concentrations and parameters                  
 */
void
CalcRhs( double *v, double t, double *guts, int n, int gn, int numguts, int which, unsigned long *gutcomps, SolverInput * si, Input * inp ) {

    int m;                      // number of nuclei
    int ap;                     // nuclear index on AP axis [0,1,...,m-1]

    int i, j;                   // loop counters
    int base, base1, base2;     // index counters for data, guts, and external inputs respectively
    unsigned long *tempctr;     // temporary counter for traversing gutcomps

    double *deriv;              // array for the derivatives at a specific time

    static int num_nucs = 0;    // store the number of nucs for next step
    static int bcd_index = 0;   // the *next* array in bicoid struct for bcd

    static DArrPtr bcd;         // pointer to appropriate bicoid struct
    double *v_ext;              // array to hold the external input concentrations at time t

    int allele = si->genindex;

    /* the following calcululates appropriate D's according to the diffusion   *
     * schedule, gets the bicoid gradient and complains if anything is wrong;  *
     * this all only happens after a cell division or at the beginning         */


    m = n / inp->zyg.defs.ngenes;       // m is the current number of nuclei

    if( m != num_nucs ) {
        GetD( t, inp->lparm.d, D, &( inp->zyg ) );
        num_nucs = m;
        bcd = GetBicoid( t, allele, inp->zyg.bcdtype, &( inp->zyg ) );
        if( bcd.size != num_nucs )
            error( "CalcRhs: %d nuclei don't match Bicoid!", num_nucs );
        bcd_index++;
    }
    // call the derivative function to get the total derivative

    deriv = ( double * ) calloc( n, sizeof( double ) );
    p_deriv( v, t, deriv, n, si, inp );

    // Here we retrieve the external input concentrations into v_ext

    v_ext = ( double * ) calloc( m * inp->zyg.defs.egenes, sizeof( double ) );
    ExternalInputs( t, t, v_ext, m * inp->zyg.defs.egenes, inp->ext[allele], inp->zyg.defs.egenes, &( inp->zyg ) );

    /* all the code below calculates the requested guts (by checking the bits  *
     * set in the gutcomps array); it does so by forward reconstructing all    *
     * the guts that are part of u and then reverse calculating (starting from *
     * the derivative) the guts outside g(u) such that guts work for any g(u); *
     * most of this is directly derived from the derivative function(s) above  */

    /* the next block of code calculates guts which are within u, i.e. all the *
     * regulatory contributions or the threshold of a certain gene             */

    for( base = 0, base1 = 0, base2 = 0, ap = 0; /*base < n, */ base1 < gn;
         base += inp->zyg.defs.ngenes, base1 += numguts, base2 += inp->zyg.defs.egenes, ap++ ) {

        tempctr = gutcomps;
        for( i = base1; i < base1 + numguts; i++ ) {
            for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                if( ( *tempctr >> j ) % 2 )
                    guts[i] += inp->lparm.T[( which * inp->zyg.defs.ngenes ) + j] * v[base + j];

            for( j = 0; j < inp->zyg.defs.egenes; j++ )
                if( ( *tempctr >> ( inp->zyg.defs.ngenes + j ) ) % 2 )
                    guts[i] += inp->lparm.E[( which * inp->zyg.defs.egenes ) + j] * v_ext[base2 + j];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes ) ) % 2 )
                guts[i] += inp->lparm.m[which] * bcd.array[ap];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 1 ) ) % 2 )
                guts[i] += inp->lparm.h[which];
            tempctr++;
        }
    }

    // the code below 'reverse calculates' everything outside g(u)

    if( n == inp->zyg.defs.ngenes ) {   // first part: one nuc, no diffusion

        // I don't trust this for loop!! What is the end condition?
        for( base = 0, base1 = 0;
             /*base < n, */ base1 < gn;
             base += inp->zyg.defs.ngenes, base1 += numguts ) {

            tempctr = gutcomps;
            for( i = base1; i < base1 + numguts; i++ ) {

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 2 ) ) % 2 )
                    guts[i] += -inp->lparm.lambda[which] * v[base + which];

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 5 ) ) % 2 )
                    guts[i] += deriv[base + which];

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 6 ) ) % 2 )
                    guts[i] += deriv[base + which]
                        + inp->lparm.lambda[which] * v[base + which];

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 7 ) ) % 2 )
                    guts[i] += ( deriv[base + which]
                                 + inp->lparm.lambda[which] * v[base + which] )
                        / inp->lparm.R[which];

                tempctr++;
            }
        }

    } else {                    // then for multiple nuclei -> diffusion

        tempctr = gutcomps;     // first anterior-most nucleus
        for( i = 0; i < numguts; i++ ) {

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 2 ) ) % 2 )
                guts[i] += -inp->lparm.lambda[which] * v[which];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 4 ) ) % 2 )
                guts[i] += D[which] * ( v[which + inp->zyg.defs.ngenes] - v[which] );

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 5 ) ) % 2 )
                guts[i] += deriv[which];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 6 ) ) % 2 )
                guts[i] += deriv[which]
                    + inp->lparm.lambda[which] * v[which]
                    - D[which] * ( v[which + inp->zyg.defs.ngenes] - v[which] );

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 7 ) ) % 2 )
                guts[i] += ( deriv[which]
                             + inp->lparm.lambda[which] * v[which]
                             - D[which] * ( v[which + inp->zyg.defs.ngenes] - v[which] ) )
                    / inp->lparm.R[which];

            tempctr++;
        }

        // then middle nuclei
        for( base = inp->zyg.defs.ngenes, base1 = numguts;
             /*base < n - inp->zyg.defs.ngenes, */ base1 < gn - numguts;
             base += inp->zyg.defs.ngenes, base1 += numguts ) {

            tempctr = gutcomps;
            for( i = base1; i < base1 + numguts; i++ ) {

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 2 ) ) % 2 )
                    guts[i] += -inp->lparm.lambda[which] * v[base + which];

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 3 ) ) % 2 )
                    guts[i] += D[which] * ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] );

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 4 ) ) % 2 )
                    guts[i] += D[which] * ( v[base + which + inp->zyg.defs.ngenes] - v[base + which] );

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 5 ) ) % 2 )
                    guts[i] += deriv[base + which];

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 6 ) ) % 2 )
                    guts[i] += deriv[base + which]
                        + inp->lparm.lambda[which] * v[base + which]
                        - D[which] * ( ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] )
                                       + ( v[base + which + inp->zyg.defs.ngenes] - v[base + which] ) );

                if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 7 ) ) % 2 )
                    guts[i] += ( deriv[base + which]
                                 + inp->lparm.lambda[which] * v[base + which]
                                 - D[which] * ( ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] )
                                                + ( v[base + which + inp->zyg.defs.ngenes] - v[base + which] ) ) )
                        / inp->lparm.R[which];

                tempctr++;
            }
        }
        // last: posterior-most nucleus
        tempctr = gutcomps;
        for( i = base1; i < base1 + numguts; i++ ) {

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 2 ) ) % 2 )
                guts[i] += -inp->lparm.lambda[which] * v[base + which];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 3 ) ) % 2 )
                guts[i] += D[which] * ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] );

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 5 ) ) % 2 )
                guts[i] += deriv[base + which];

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 6 ) ) % 2 )
                guts[i] += deriv[base + which]
                    + inp->lparm.lambda[which] * v[base + which]
                    - D[which] * ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] );

            if( ( *tempctr >> ( inp->zyg.defs.ngenes + inp->zyg.defs.egenes + 7 ) ) % 2 )
                guts[i] += ( deriv[base + which]
                             + inp->lparm.lambda[which] * v[base + which]
                             - D[which] * ( v[base + which - inp->zyg.defs.ngenes] - v[base + which] ) )
                    / inp->lparm.R[which];

            tempctr++;
        }
    }

    free( deriv );
    free( v_ext );
    return;
}

/** ParseString: parses a line of the $gutsdefs section (dataformatX.X 
 *                has the details); this function then returns two things: 
 *                1) the number of 'words' separated by whitespace in that 
 *                   string is the function's return value                 
 *                2) the 'words' themselves are returned in arginp         
 */
int
ParseString( char *v, char **arginp ) {
    const int IN = 1;           /* are we inside or outside a word? */
    const int OUT = 0;

    int state = OUT;            /* state can be either 'IN' or 'OUT' above */
    int n = 0;                  /* counter for words */
    char *mover = v;            /* pointer-counter for pos within whole string */
    char **myargs = arginp;     /* pointer-counter for different words */
    char *tempo = NULL;         /* pointer-counter for pos within current word */

    /* parsing starts here */

    n = 0;
    state = OUT;
    while( ( *mover != '\0' ) && ( *mover != '\n' ) ) { /* loop thru whole line */
        if( *mover == ' ' || *mover == '\t' ) { /* whitespace -> done with word */
            if( state == IN ) { /* word just finished: */
                *tempo = '\0';  /* finish string by adding \0 */
                myargs++;       /* go to next word */
            }
            state = OUT;        /* whitespace -> outside word */
        } else {
            if( state == OUT ) {        /* a new word */
                state = IN;     /* we're inside a word now */
                *myargs = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
                tempo = *myargs;        /* allocate and assign tempo to new string */
                ++n;
            }
            *tempo = *mover;    /* characters copied here */
            tempo++;            /* tempo keeps track of pos within word */
        }
        mover++;                /* mover keeps track of pos within whole string */
    }

    if( state == IN ) {         /* if we finished within a word, add a final \0 */
        *tempo = '\0';
        myargs++;
    }

    return n;                   /* return the number of words in the string */
}

/** GetGutsComps: takes a geneID string and an array of ID strings; then 
 *                 returns an array of long ints with bit flags set that   
 *                 tell CalcGuts() what to calculate for each column of    
 *                 guts output; this function also returns a pointer to    
 *                 the current gene in the gene ID string                  
 */
char *
GetGutsComps( char *geneidstring, char *egeneidstring, char **specsin, unsigned long *specsout ) {
    const char *gutsids = "BALXYDJZU";  /* legal guts IDs */

    char *strcntr = NULL;       /* ptr-counter within ID string */
    char *ids = NULL;           /* string with all legal IDs */
    char **strarrcntr = NULL;   /* ptr-counter for gut ID strings */
    unsigned long *arrcntr = NULL;      /* ptr-counter for gutcomps */
    char *i = NULL;             /* ptr for finding IDs in ids string */


    /* ids = geneIDs plus gutsIDs */

    ids = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    ids = strcat( strcat( strcat( ids, geneidstring ), egeneidstring ), gutsids );

    /* set pointer-counters to both ID strings and gutcomps */

    arrcntr = specsout;
    strarrcntr = specsin;
    strarrcntr++;               /* skip first entry (gene name) */

    /* the following loop parses the words of the ID string and sets according *
     * bits in the gutscomp array                                              */

    // why the assignment?? Does it terminate??
    while( ( strcntr = *strarrcntr ) ) {
        while( *strcntr != '\0' ) {
            i = strchr( ids, *strcntr );
            if( i == NULL )
                error( "GetGutsComps: unknown gene in gutsdefs string" );
            else if( *strcntr == 'U' ) {
                if( ( *( strcntr + 1 ) != '\0' ) || ( strcntr != *strarrcntr ) ) {
                    error( "GetGutsComps: U not alone in gutsdefs string" );
                }
                *arrcntr = ( 1 << ( strlen( geneidstring ) + strlen( egeneidstring ) + 2 ) ) - 1;
            } else {
                *arrcntr |= ( 1 << ( i - ids ) );
            }
            strcntr++;
        }
        arrcntr++;
        strarrcntr++;
    }

    /* clean up and return a pointer to the gut gene in the gene ID string */
    free( ids );
    return strchr( geneidstring, **specsin );
}




/*** MUTATOR FUNCTIONS *****************************************************/

/** Mutate: calls mutator functions according to genotype string */
EqParms
Mutate( char *g_type, EqParms parm, TheProblem * defs ) {
    int i;
    int c;

    char *record;

    EqParms lparm;

    lparm = CopyParm( parm, defs );     /* make local copy of parameters to be mutated */

    record = g_type;
    c = ( int ) *record;

    for( i = 0; c != '\0'; i++, c = ( int ) *( ++record ) ) {
        if( c == 'W' )
            continue;
        else if( c == 'R' )
            R_Mutate( i, &lparm );
        else if( c == 'S' )
            RT_Mutate( i, defs->ngenes, &lparm );
        else if( c == 'T' )
            T_Mutate( i, defs->ngenes, &lparm );
        else
            error( "Mutate: unrecognized letter in genotype string!" );
    }
    return lparm;
}

/** T_Mutate: mutates genes by setting all their T matrix entries 
 *             to zero. Used to simulate mutants that express a            
 *             non-functional protein.                                     
 */
void
T_Mutate( int gene, int ngenes, EqParms * lparm ) {

    int i;

    for( i = 0; i < ngenes; i++ )
        lparm->T[( i * ngenes ) + gene] = 0;
}

/** R_Mutate: mutates genes by setting their promotor strength R 
 *             to zero, so there will be no transcription at all           
 *             anymore. Used to simulate mutants that don't pro-           
 *             duce any protein anymore.                                   
 */
void
R_Mutate( int gene, EqParms * lparm ) {
    lparm->R[gene] = 0;
}

/** RT_Mutate: mutates gene by setting both promoter strength R 
 *              and T matrix entries to zero. Useful, if you want          
 *              to be really sure that there is NO protein pro-            
 *              duction and that maternal contribution also don't          
 *              contribute anything to gene interactions.                  
 */
void
RT_Mutate( int gene, int ngenes, EqParms * lparm ) {
    int i;

    for( i = 0; i < ngenes; i++ )
        lparm->T[( i * ngenes ) + gene] = 0;
    lparm->R[gene] = 0;
    /* Don't need to zero param.thresh in (trans acting) mutants */
}

/** CopyParm: copies all the parameters into the lparm struct */
EqParms
CopyParm( EqParms orig_parm, TheProblem * defs ) {
    int i, j;                   /* local loop counters */

    EqParms l_parm;             /* copy of parm struct to be returned */

    l_parm.R = ( double * ) calloc( defs->ngenes, sizeof( double ) );
    l_parm.T = ( double * ) calloc( defs->ngenes * defs->ngenes, sizeof( double ) );
    l_parm.E = ( double * ) calloc( defs->ngenes * defs->egenes, sizeof( double ) );
    l_parm.m = ( double * ) calloc( defs->ngenes, sizeof( double ) );
    l_parm.h = ( double * ) calloc( defs->ngenes, sizeof( double ) );
    if( ( defs->diff_schedule == 'A' ) || ( defs->diff_schedule == 'C' ) ) {
        l_parm.d = ( double * ) malloc( sizeof( double ) );
    } else {
        l_parm.d = ( double * ) calloc( defs->ngenes, sizeof( double ) );
    }
    l_parm.lambda = ( double * ) calloc( defs->ngenes, sizeof( double ) );
    l_parm.tau = ( double * ) calloc( defs->ngenes, sizeof( double ) );

    for( i = 0; i < defs->ngenes; i++ ) {
        l_parm.R[i] = orig_parm.R[i];
        for( j = 0; j < defs->ngenes; j++ )
            l_parm.T[( i * defs->ngenes ) + j] = orig_parm.T[( i * defs->ngenes ) + j];
        for( j = 0; j < defs->egenes; j++ )
            l_parm.E[( i * defs->egenes ) + j] = orig_parm.E[( i * defs->egenes ) + j];
        l_parm.m[i] = orig_parm.m[i];
        l_parm.h[i] = orig_parm.h[i];
        l_parm.lambda[i] = orig_parm.lambda[i];
        l_parm.tau[i] = orig_parm.tau[i];
    }

    if( ( defs->diff_schedule == 'A' ) || ( defs->diff_schedule == 'C' ) ) {
        l_parm.d[0] = orig_parm.d[0];
    } else {
        for( i = 0; i < defs->ngenes; i++ )
            l_parm.d[i] = orig_parm.d[i];
    }

    return l_parm;
}


/*** GetMutParameters: same as above but returns the mutated copy of the ***
 *                     parameter struct; important for writing guts        *
 ***************************************************************************/

/*EqParms *GetMutParameters(void)
{
  return &lparm;
}*/

/** Function given the starting (fb) and ending (fa) time structures,
calculates the use cpu time in seconds */
double
tvsub( struct rusage a, struct rusage b ) {
    double fa, fb;

    fa = a.ru_utime.tv_sec + a.ru_stime.tv_sec + ( a.ru_utime.tv_usec + a.ru_stime.tv_usec ) / 1e6;
    fb = b.ru_utime.tv_sec + b.ru_stime.tv_sec + ( b.ru_utime.tv_usec + b.ru_stime.tv_usec ) / 1e6;

    return fa - fb;
}

/** DvdtDelay: the delay derivative function; implements the equations 
 *             as discussed in the delay report,                    
 *             plus different g(u) functions.					    
 */
void
DvdtDelay( double *v, double **vd, double t, double *vdot, int n, SolverInput * si, Input * inp ) {
    double vinput1 = 0;
    int m;                      /* number of nuclei */
    int ap;                     /* nuclear index on AP axis [0,1,...,m-1] */
    int i, j;                   /* local loop counters */
    int k;                      /* index of gene k in a specific nucleus */
    int base, base1;            /* index of first gene in a specific nucleus */
    int *l_rule;                /* for autonomous implementation */
#ifdef ALPHA_DU
    int incx = 1;               /* increment step size for vsqrt input array */
    int incy = 1;               /* increment step size for vsqrt output array */
#endif

    static int num_nucs = 0;    /* store the number of nucs for next step */
    static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
    static DArrPtr bcd;         /* pointer to appropriate bicoid struct */
    double **v_ext;             /* array to hold the external input
                                   concentrations at time t */
    int allele = si->genindex;

    vinput = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot2 = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );
    bot = ( double * ) calloc( inp->zyg.defs.ngenes * inp->zyg.defs.nnucs, sizeof( double ) );

    /* get D parameters and bicoid gradient according to cleavage cycle */
    m = n / inp->zyg.defs.ngenes;       /* m is the number of nuclei */
    if( m != num_nucs ) {       /* time-varying quantities only vary by ccycle */
        GetD( t, inp->lparm.d, D, &( inp->zyg ) );
        /* get diff coefficients, according to diff schedule */
        if( num_nucs > m )      /* started a new iteration in score */
            bcd_index = 0;      /* -> start all over again! */
        num_nucs = m;           /* store # of nucs for next step */
        bcd = GetBicoid( t, allele, inp->zyg.bcdtype, &( inp->zyg ) );  /* get bicoid gradient */
        if( bcd.size != num_nucs )
            error( "DvdtDelay: %d nuclei don't match Bicoid!", num_nucs );
        bcd_index++;            /* store index for next bicoid gradient */
    }

    l_rule = ( int * ) calloc( inp->zyg.defs.ngenes, sizeof( int ) );

    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {
        l_rule[i] = !Theta( t - inp->lparm.tau[i], &( inp->zyg ) );     /*for autonomous equations */
        if( debug ) {
            if( l_rule[i] == 0 )
                printf( "L_RULE = %d for t=%lg, tau=%lg\n", l_rule[i], t, inp->lparm.tau[i] );
        }

    }

    //Theta(u) = false while interphase

    /* l_rule is zero during mitosis, in order
     * to put to zero the regulation part of the
     * equation. Remember, no regulation during
     * mitosis */

    /* Here we retrieve the external input concentrations into v_ext */

    v_ext = ( double ** ) calloc( inp->zyg.defs.ngenes, sizeof( double * ) );
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {

        v_ext[i] = ( double * ) calloc( m * inp->zyg.defs.egenes, sizeof( double ) );
        ExternalInputs( t - inp->lparm.tau[i], t, v_ext[i], m * inp->zyg.defs.egenes, inp->ext[allele], inp->zyg.defs.egenes, &( inp->zyg ) );

    }
    /* This is how it works (by JR): 

       ap      nucleus position on ap axis
       base    index of first gene (0) in current nucleus
       k       index of gene k in current nucleus 

       Protein synthesis terms are calculated according to g(u)

       First we do loop for vinput contributions; vinput contains
       the u that goes into g(u)

       Then we do a separate loop or vector func for sqrt or exp

       Then we do vdot contributions (R and lambda part)

       Then we do the spatial part (Ds); we do a special case for 
       each end 

       Note, however, that before the real loop, we have to do a 
       special case for when there is one nuc, hence no diffusion

       These loops look a little funky 'cause we don't want any 
       divides'                                                       */


    /***************************************************************************
     *                                                                         *
     *        g(u) = 1/2 * ( u / sqrt(1 + u^2) + 1)                            *
     *                                                                         *
     ***************************************************************************/
    if( gofu == Sqrt ) {

        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[k][base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * vd[k][base + j];

                bot2[i] = 1 + vinput1 * vinput1;
                vinput[i] = vinput1;
            }
        }

        /* now calculate sqrt(1+u^2); store it in bot[] */
#ifdef ALPHA_DU
        vsqrt_( bot2, &incx, bot, &incy, &n );  /* superfast DEC vector function */
#else
        for( i = 0; i < n; i++ )        /* slow traditional style sqrt */
            bot[i] = sqrt( bot2[i] );
#endif
        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
         *                                                                         *
         ***************************************************************************/

    } else if( gofu == Tanh ) {

        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[k][base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * vd[k][base + j];

                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = tanh( vinput[i] ) + 1;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = tanh( vinput[i] ) + 1;
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = tanh( vinput[i] ) + 1;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = tanh( vinput[i] ) + 1;
                vdot1 += l_rule[k] * inp->lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1 / (1 + exp(-2u))                                        *
         *                                                                         *
         ***************************************************************************/

    } else if( gofu == Exp ) {

        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[k][base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * vd[k][base + j];

                vinput[i] = -2.0 * vinput1;
            }
        }

        /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
        vexp_( vinput, &incx, bot, &incy, &n ); /* superfast DEC vector function */
#else
        for( i = 0; i < n; i++ )        /* slow traditional style exp */
            bot[i] = exp( vinput[i] );
#endif

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 / ( 1 + bot[i] );
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 / ( 1 + bot[i] );
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = 1 / ( 1 + bot[i] );
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = 1 / ( 1 + bot[i] );
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 0 if u<0, 1 if u>=0 (Heaviside function)                  *
         *                                                                         *
         *        this makes the model quasi boolean and the equations locally     *
         *        linear for both u>0 and u<0                                      *
         ***************************************************************************/

    } else if( gofu == Hvs ) {

        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[k][base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * vd[k][base + j];
                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    if( vinput[i] >= 0. )
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                if( vinput[i] >= 0. )
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    if( vinput[i] >= 0. )
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                if( vinput[i] >= 0. )
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

    /***************************************************************************
     *                                                                         *
     *        g(u) = u                                                         *
     *                                                                         *
     *        we don't have a sigmoid funcion                                  *
     *                                                                         *
     ***************************************************************************/

    } else if( gofu == Kolja ) {
        for( base = 0, base1 = 0, ap = 0; base < n; base += inp->zyg.defs.ngenes, base1 += inp->zyg.defs.egenes, ap++ ) {
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {

                k = i - base;

                vinput1 = inp->lparm.h[k];
                vinput1 += inp->lparm.m[k] * bcd.array[ap];     /* ap is nuclear index */

                for( j = 0; j < inp->zyg.defs.egenes; j++ )
                    vinput1 += inp->lparm.E[( k * inp->zyg.defs.egenes ) + j] * v_ext[k][base1 + j];

                for( j = 0; j < inp->zyg.defs.ngenes; j++ )
                    vinput1 += inp->lparm.T[( k * inp->zyg.defs.ngenes ) + j] * vd[k][base + j];

                vinput[i] = vinput1;    //u
            }
        }
        if( debug ) {
            printf( "vinput after u = %lg\n", vinput[0] );
        }
        if( n == inp->zyg.defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base = 0; base < n; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = vinput[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* first anterior-most nucleus */
                k = i;

                vdot1 = -inp->lparm.lambda[k] * v[i];
                //printf("afterdecay vdot0=%lg vi=%lg\n", vdot1, v[i]);
                g1 = vinput[i];
                if( debug ) {
                    printf( "i=%d, vdot1=%lg, g1=%lg, l_rule[k]=%d, lambda=%lg R=%lg\n", i, vdot1, g1, l_rule[k], -inp->lparm.lambda[k], inp->lparm.R[k] );
                }
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                if( debug ) {
                    printf( "beforediffusion vdot0=%lg\n", vdot1 );
                }
                vdot1 += D[i] * ( v[i + inp->zyg.defs.ngenes] - v[i] );
                if( debug ) {
                    printf( "afterdiffusion vdot0=%lg D=%lg\n", vdot1, D[i] );
                }
                vdot[i] = vdot1;
            }
            if( debug ) {
                printf( "1nuc vdot0=%lg\n", vdot[0] );
            }

            /* then middle nuclei */
            for( base = inp->zyg.defs.ngenes; base < n - inp->zyg.defs.ngenes; base += inp->zyg.defs.ngenes ) {
                for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                    k = i - base;
                    vdot1 = -inp->lparm.lambda[k] * v[i];
                    g1 = vinput[i];
                    vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                    vdot1 += D[k] * ( ( v[i - inp->zyg.defs.ngenes] - v[i] ) + ( v[i + inp->zyg.defs.ngenes] - v[i] ) );
                    vdot[i] = vdot1;
                }
            }

            /* last: posterior-most nucleus */
            for( i = base; i < base + inp->zyg.defs.ngenes; i++ ) {
                k = i - base;
                vdot1 = -inp->lparm.lambda[k] * v[i];
                g1 = vinput[i];
                vdot1 += l_rule[k] * inp->lparm.R[k] * g1;
                vdot1 += D[k] * ( v[i - inp->zyg.defs.ngenes] - v[i] );
                vdot[i] = vdot1;
            }
        }

    } else {
        error( "DvdtDelay: unknown g(u)" );
    }

    /* during mitosis only diffusion and decay happen */
    if( debug ) {
        printf( "vdot0=%lg, vdot1=%lg, vdot2=%lg, vdot3=%lg\n", vdot[0], vdot[1], vdot[2], vdot[3] );
    }

    free( l_rule );

    for( i = 0; i < inp->zyg.defs.ngenes; i++ )
        free( v_ext[i] );

    free( v_ext );
    free( vinput );
    free( bot2 );
    free( bot );
    return;
}

/////////////////////////////////////////////////////////////////////////////

/** Dvdt_sqrt: reimplementation of part of DvdtOrig that should make the 
 *              maintenance easier and gave a small speed up (~12%). The   
 *              function is now subdivided into 3 subfunctions that allow  
 *              easy reuse with the krylov preconditioner.                 
 */
void
Dvdt_sqrt( double *v, double t, double *vdot, int n, SolverInput * si, Input * inp ) {

    //printf("+++++V = %f %f %f %f %f %f\n", v[0], v[1], v[2], v[3], v[4], v[5]);

    // number of nuclei
    const int MM = n / inp->zyg.defs.ngenes;
    // counter
    int i, rule;

    int allele = si->genindex;

    // array to hold the external input concentrations at time t
    double *v_ext;
    // store the number of nucs for next step
    static int num_nucs = 0;
    // the next array in bicoid struct for Bcd
    static int bcd_index = 0;
    // pointer to appropriate bicoid struct
    static DArrPtr bcd;

    if( num_nucs != MM ) {
        // get D parameters according to cleavage cycle
        // NOTE: time-varying quantities only vary by ccycle
        GetD( t, inp->lparm.d, D, &( inp->zyg ) );
        if( num_nucs > MM )     /* started a new iteration in score */
            bcd_index = 0;      /* -> start all over again! */
        num_nucs = MM;          /* store # of nucs for next step */
        bcd = GetBicoid( t, allele, inp->zyg.bcdtype, &( inp->zyg ) );  /* get bicoid gradient */
        if( bcd.size != num_nucs )
            error( "DvdtDelay: %d nuclei don't match Bicoid!", num_nucs );
        bcd_index++;            /* store index for next bicoid gradient */
    }

    // here we retrieve the external input concentrations into v_ext
    v_ext = ( double * ) calloc( MM * inp->zyg.defs.egenes, sizeof( double ) );
    ExternalInputs( t, t, v_ext, MM * inp->zyg.defs.egenes, inp->ext[allele], inp->zyg.defs.egenes, &( inp->zyg ) );

    //printf("TIME = %lg\n", si->time);
    rule = GetRule( si->time, &( inp->zyg ) );

    // in interphase there is gene product synthesis (rna or protein)
    if( rule == INTERPHASE ) {
        Dvdt_production( v, t, vdot, n, v_ext, bcd, MM, inp );
    } else {
        // set vdot to zero: no production
        for( i = 0; i < n; ++i )
            vdot[i] = 0.0;
    }
    Dvdt_degradation( v, t, vdot, n, inp );
    Dvdt_diffusion( v, t, vdot, n, MM, D, inp );
    free( v_ext );
}

/*
void Dvdt_sqrt_precond(double *v, double t, double *vdot, int n, int allele, Input *inp) {

    // number of nuclei
    const int MM = n / inp->zyg.defs.ngenes;
    // counter
    int i, rule;
    // array to hold the external input concentrations at time t
    double *v_ext;

    rule = GetRule(t, &(inp->zyg));
    // Here we retrieve the external input concentrations into v_ext
    v_ext = (double *) calloc(MM * inp->zyg.defs.egenes, sizeof (double));
    ExternalInputs(t, t, v_ext, MM * inp->zyg.defs.egenes, inp->ext[allele], inp->zyg.defs.egenes, &(inp->zyg));

    if (rule == INTERPHASE) {
        Dvdt_production(v, t, vdot, n, v_ext, MM, inp);
    } else {
        // set vdot to zero: no production
        for (i = 0; i < n; ++i) vdot[ i ] = 0.0;
    }
    Dvdt_degradation(v, t, vdot, n, inp);
    // and no diffusion
    free(v_ext);
}
*/

/** subfunction for (regulated) production */
void
Dvdt_production( double *v, double t, double *vdot, int n, double *v_ext, DArrPtr bcd, int MM, Input * inp ) {
    // nuclear index on AP axis [0,1,...,m-1]
    int ap;
    // local loop counters
    int i, j, aux1, aux2;
    // index of first gene in a specific nucleus
    int base, base_e;
    // auxilary var
    double aux;
    // forall nuclei do production/regulation
    for( ap = 0; ap < MM; ++ap ) {
        // external genes
        base_e = ap * inp->zyg.defs.egenes;
        // gap genes
        base = ap * inp->zyg.defs.ngenes;
        // per nucleus
        for( i = 0; i < inp->zyg.defs.ngenes; ++i ) {

            aux = inp->lparm.h[i];
            aux += inp->lparm.m[i] * bcd.array[ap];
            aux1 = ( i * inp->zyg.defs.egenes );
            aux2 = ( i * inp->zyg.defs.ngenes );
            for( j = 0; j < inp->zyg.defs.egenes; ++j ) {
                aux += inp->lparm.E[aux1 + j] * v_ext[base_e + j];
            }
            for( j = 0; j < inp->zyg.defs.ngenes; ++j ) {
                aux += inp->lparm.T[aux2 + j] * v[base + j];
                //if ((ap == 0) && (i == 0))
                //      printf("T%d=%f v=%f\n", j, inp->lparm.T[ aux2 + j ], v[ base + j ]);
            }
            //printf("aux=%lg interm=%lg rez=%lg\n", aux, (sqrt(1 + aux * aux)), (aux / sqrt(1 + aux * aux)));
            vdot[base + i] = inp->lparm.R[i] * 0.5 * ( 1 + aux / sqrt( 1 + aux * aux ) );
        }
    }
}

/** subfunction for the degradation */
void
Dvdt_degradation( double *v, double t, double *vdot, int n, Input * inp ) {

    // forall nuclei*ngenes do degradation
    int base = 0, i;

    for( i = 0; i < n; ++i ) {
        // cycle through gap genes
        if( base == inp->zyg.defs.ngenes )
            base = 0;
        // degrade
        vdot[i] -= inp->lparm.lambda[base] * v[i];
        // increment
        ++base;
    }
}

/** subfunction doing the diffusion */
void
Dvdt_diffusion( double *v, double t, double *vdot, int n, int MM, double *D, Input * inp ) {
    const int GG = inp->zyg.defs.ngenes;
    // counters and auxilary vars
    int i, ap, base, aux1;
    // and forall nuclei diffusion (one nucleus special case)
    if( n != GG ) {

        // anterior most nucleus
        for( i = 0; i < GG; ++i ) {
            vdot[i] += D[i] * ( v[i + GG] - v[i] );
        }

        // then middle nuclei
        for( ap = 1; ap < MM - 1; ++ap ) {
            // gap genes
            base = ap * GG;
            for( i = 0; i < GG; ++i ) {
                vdot[base + i] += D[i] * ( v[base + i - GG] + v[base + i + GG] - 2 * v[base + i] );
            }
        }

        // last: posterior-most nucleus
        base = ( MM - 1 ) * GG;
        for( i = 0; i < GG; ++i ) {
            aux1 = base + i;
            vdot[aux1] += D[i] * ( v[aux1 - GG] - v[aux1] );
        }

    }                           // else we skip diffusion
}



