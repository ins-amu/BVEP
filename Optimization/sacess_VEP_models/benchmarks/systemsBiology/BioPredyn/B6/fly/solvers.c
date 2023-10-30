/**
 * @file solvers.c                                             
 * @author JR, modified by Yoginho,               
 * Additional solvers by: Joel Linton (Rk4),                    
 *                        Johannes Jaeger (Rk2, Meuler, Heun),  
 *                        Konstantin Kozlov (Milne, Adams), 
 *                        Marcel Wolf (Rkck, Rkf),    
 *                        Manu (Adams, BuSt),
 *                        Anton Crombach (Krylov)
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Contains the solver function implementations (can be toggled by -s 
 * option) that propagate the equations.
 *                                                               
 * @note \em ONLY general solvers allowed here; they \em MUST comply to 
 * the generic solver interface                         
 *                                                               
 * @warning The functions in this file return vin and vout unchanged iff  
 * stepsize = 0. It is the responsibility of the calling function to be sure 
 * that this is not a problem.                 
 *                                                                
 * \par Usage:
 *                                                               
 *  All solvers adhere to the same interface:                    
 *                                                               
 *  Solver(double *vin, double *vout, double tin, double tout,   
 *         double stephint, double accuracy, int n);             
 *                                                               
 *  Arguments:                                                   
 *                                                               
 *  - vin       array of dependent variables at 'tin' (start)    
 *  - vout      array of dependent variables at 'tout' (end);    
 *              this is what's returned by the solver            
 *  - tin       time at start                                    
 *  - tout      time at end                                      
 *  - stephint  suggested stepsize for fixed stepsize solvers;    
 *              see extensive comment below; the embedded Rk     
 *              solvers (Rkck, Rkf) take stephint as their ini-  
 *              tial stepsize                                    
 *  - accuracy  accuracy for adaptive stepsize solvers; accuracy 
 *              is always relative, i.e. 0.1 means 10% of the    
 *              actual v we've evaluated                         
 *  - n         size of vin and vout arrays; the user is respon- 
 *              sible for checking that the number of elements   
 *              does not change in the model between vin/vout    
 *                                                               
 *  Note that if tin = tout the solvers will just returns with-  
 *  out doing anything. The caller is responsible for handling   
 *  such a situation.                                            
 *                                                               
 * \par Solver Naming Conventions: 
 *                                                               
 *  There are various contradictory solver nomenclatures in the  
 *  literature. If not stated otherwise, we use the ones from    
 *  the following two books, on which most of our solvers are    
 *  based:                                                       
 *                                                               
 *  Press, Teukolsky, Vetterling & Flannery (1992). Numerical    
 *  Recipes in C. Cambridge University Press, Cambridge, U.K.    
 *                                                               
 *  Burden & Faires (1993). Numerical Analysis, 5th Edition.     
 *  PWS Publishing Co., Boston MA, U.S.A.                        
 *                                                               
 *  
 * \par Important Comment: 
 *                                                               
 * stephint and the stepsize of fixed-stepsize solvers:          
 *                                                               
 * The stepsize which is passed to a fixed-stepsize solvers is   
 * actually a stephint, i.e. only an approximation to the actual  
 * stepsize which depends on the length of the interval over     
 * which the equations are propagated:                            
 *                                                               
 * if the remainder of (deltat div stephint) is:                 
 *                                                               
 *   0                                   stepsize = stephint     
 *   > 0 && < stephint/2                 stepsize > stephint     
 *   >= stephint/2 && < stephint         stepsize < stephint     
 *                                                               
 * The exact amount of the difference between stepsize and step- 
 * hint depends on the size of deltat and stephint:              
 *  - the larger deltat, the smaller the maximum possible diffe- 
 *    rence for a given stephint                                 
 *  - the larger stephint, the larger the maximum possible di-   
 *    fference for a given stephint                              
 *                                                               
 * Since mitoses only last about 3-5 minutes, these differences  
 * are quite significant for large stepsizes ( > 0.1 ). For the  
 * 4div schedule (newstyle datafiles) we have maximum relative   
 * differences between stephint and actual stepsize of about:     
 *                                                               
 *   stephint                     max. rel. difference           
 *                                                               
 *        1.0                                    10.0%           
 *        0.5                                     5.8%           
 *        0.3                                     4.3%           
 *        0.1                                     1.0%           
 *                                                               
 * See ../doc/step4div.ps for more details.                      
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

//#include <error.h>
#include <solvers.h>
#include <maternal.h>
#include <zygotic.h>
#include <integrate.h>

//#include "fly_opt.h"

/*** STATIC VARIABLES AND MACROS *******************************************/

double *d;                      /* D's used for Neville extrapolation in BuSt() */
double *hpoints;                /* stepsizes h (=H/n) which we try in BuSt() */
double maxdel, mindel;
int numdel;                     /* delay parameters used by DCERk32, y_delayed */
int gridstart;                  /* for the heuristic */
double *delay;                  /* static array set in SoDe, used by DCERk32 */
int gridpos;                    /* where you are in the grid */
double *tdone;                  /* the grid */
double **derivv1;               /* intermediate derivatives for the Cash-Karp formula */
double **derivv2;
double **derivv3;
double **derivv4;
double **vdonne;

void ( *d_deriv ) ( double *, double **, double, double *, int, SolverInput *, Input * );
void ( *p_deriv ) ( double *, double, double *, int, SolverInput *, Input * );
void ( *p_jacobn ) ( double, double *, double *, double **, int, SolverInput *, Input * );

/* three macros used in various solvers below */

double dqrarg;

#define DSQR(a) ((dqrarg=(a)) == 0.0 ? 0.0 : dqrarg*dqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1 = (a), dmaxarg2 = (b), (dmaxarg1) > \
(dmaxarg2) ?  (dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1 = (a), dminarg2 = (b), (dminarg1) < \
(dminarg2) ?  (dminarg1) : (dminarg2))

static Input *inp;

static SolverInput *si;

/*** Krylov solver variables added by Anton Crombach, October 2010 *********/

/* extra data needed to be shared among different functions */
//static ExtraData edata = NULL;

/* memory for the solver to use */
static void *cvode_mem = NULL;
/* memory that holds the current state of the system */
static N_Vector vars = NULL;
/* number of equations (is also length of `vars') */
static int neq = -1;



/*** SOLVERS ***************************************************************/

/** Euler: propagates vin (of size n) from tin to tout by the Euler 
 *          method; the result is returned by vout                  
 */
void
Euler( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *deriv;              /* the derivatives at the beginning of a step */

    double deltat;              /* tout - tin */
    double t;                   /* current time */

    double m;                   /* tmp var to calculate precise stepsize from stephint */
    double stepsize;            /* real stepsize */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* total number of steps we have to take */

    int nd = 0;                 /* number of deriv evaluations */

    /* steps too small: just return */

    if( tin == tout )
        return;

    /* if steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv = ( double * ) calloc( n, sizeof( double ) );

    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */
    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Euler: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    for( step = 0; step < nsteps; step++ ) {    /* loop for steps */

        p_deriv( vnow, t, deriv, n, si, inp );  /* call derivative func to evaluate deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )
            vnext[i] = vnow[i] + stepsize * deriv[i];   /* Euler formula */

        t += stepsize;          /* go to next step */

        if( step < nsteps - 2 ) {       /* CASE 1: many steps to go */
            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
            vnow = v[toggle];   /* set vout */
            vnext = vout;

        } else if( step > nsteps - 2 ) {        /* CASE 3: just did final iteration */
            free( v[0] );       /* clean up and go home! */
            free( v[1] );
            free( v );
            free( deriv );
        }
    }

    if( debug )
        WriteSolvLog( "Euler", tin, tout, stepsize, nsteps, nd, slog );

    return;
}

/** Meuler: propagates vin (of size n) from tin to tout by the Modified 
 *           Euler method (this is NOT the midpoint method, see Rk2()); 
 *           the result is returned by vout                             
 */
void
Meuler( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's for midpoint */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* the derivatives at midpoint of a step */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double th;                  /* time at end of a step (t+stepsize) */

    double m;                   /* tmp var to calculate precise stepsize from stephint */
    double stepsize;            /* real stepsize */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* total number of steps we have to take */

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );


    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */
    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Meuler: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    for( step = 0; step < nsteps; step++ ) {    /* loop for steps */

        th = t + stepsize;      /* time of step endpoint */

        p_deriv( vnow, t, deriv1, n, si, inp ); /* first call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )        /* evaluate v's at endpoint */
            vtemp[i] = vnow[i] + stepsize * deriv1[i];

        p_deriv( vtemp, th, deriv2, n, si, inp );       /* second call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )        /* Modified Euler formula */
            vnext[i] = vnow[i] + ( stepsize / 2. ) * ( deriv1[i] + deriv2[i] );

        t += stepsize;          /* go to next step */

        if( step < nsteps - 2 ) {       /* CASE 1: many steps to go */
            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
            vnow = v[toggle];   /* set vout */
            vnext = vout;

        } else if( step > nsteps - 2 ) {        /* CASE 3: just did final iteration */
            free( v[0] );       /* clean up and go home! */
            free( v[1] );
            free( v );
            free( vtemp );
            free( deriv1 );
            free( deriv2 );
        }
    }

    if( debug )
        WriteSolvLog( "Meuler", tin, tout, stepsize, nsteps, nd, slog );

    return;
}

/** Heun: propagates vin (of size n) from tin to tout by Heun's method 
 *         the result is returned by vout                              
 */
void
Heun( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's for midpoint */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* the derivatives at midpoint of a step */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double thh;                 /* time at 2/3 of the step */

    double m;                   /* tmp var to calculate precise stepsize from stephint */
    double stepsize;            /* real stepsize */
    double hh;                  /* 2/3 the stepsize */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* total number of steps we have to take */

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );


    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */
    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Heun: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 2. / 3.;    /* evaluate 2/3 of stepsize */

    for( step = 0; step < nsteps; step++ ) {    /* loop for steps */

        thh = t + hh;           /* time at 2/3 of the step */

        p_deriv( vnow, t, deriv1, n, si, inp ); /* first call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )        /* evaluate v's at 2/3 of the step */
            vtemp[i] = vnow[i] + hh * deriv1[i];

        p_deriv( vtemp, thh, deriv2, n, si, inp );      /* second call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )        /* Heun's formula */
            vnext[i] = vnow[i] + ( stepsize / 4. ) * ( deriv1[i] + 3. * deriv2[i] );

        t += stepsize;          /* go to next step */

        if( step < nsteps - 2 ) {       /* CASE 1: many steps to go */
            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
            vnow = v[toggle];   /* set vout */
            vnext = vout;

        } else if( step > nsteps - 2 ) {        /* CASE 3: just did final iteration */
            free( v[0] );       /* clean up and go home! */
            free( v[1] );
            free( v );
            free( vtemp );
            free( deriv1 );
            free( deriv2 );
        }
    }

    if( debug )
        WriteSolvLog( "Heun", tin, tout, stepsize, nsteps, nd, slog );

    return;
}

/** Rk2: propagates vin (of size n) from tin to tout by the Midpoint or 
 *        Second-Order Runge-Kutta method; the result is returned by vout  
 */
void
Rk2( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's for midpoint */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* the derivatives at midpoint of a step */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double thh;                 /* time at half of the step */

    double m;                   /* tmp var to calculate precise stepsize from stephint */
    double stepsize;            /* real stepsize */
    double hh;                  /* half the stepsize */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* total number of steps we have to take */

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );


    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */
    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Rk2: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;        /* evaluate half of stepsize */

    for( step = 0; step < nsteps; step++ ) {    /* loop for steps */

        thh = t + hh;           /* time of interval midpoints */

        p_deriv( vnow, t, deriv1, n, si, inp ); /* first call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )        /* evaluate v's at midpoint */
            vtemp[i] = vnow[i] + hh * deriv1[i];

        p_deriv( vtemp, thh, deriv2, n, si, inp );      /* second call to deriv */

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )
            vnext[i] = vnow[i] + stepsize * deriv2[i];  /* Midpoint formula */

        t += stepsize;          /* go to next step */

        if( step < nsteps - 2 ) {       /* CASE 1: many steps to go */
            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
            vnow = v[toggle];   /* set vout */
            vnext = vout;

        } else if( step > nsteps - 2 ) {        /* CASE 3: just did final iteration */
            free( v[0] );       /* clean up and go home! */
            free( v[1] );
            free( v );
            free( vtemp );
            free( deriv1 );
            free( deriv2 );
        }
    }

    if( debug )
        WriteSolvLog( "Rk2", tin, tout, stepsize, nsteps, nd, slog );

    return;
}

/** Rk4: propagates vin (of size n) from tin to tout by the Fourth-Order 
 *        Runge-Kutta method; the result is returned by vout               
 ***************************************************************************
 *                                                                         
 * written by Joel Linton (somewhere around 1998)                          
 * fixed and modified by Yoginho (somewhere around 2001)                   
 *                                                                         
 */
void
Rk4( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* the derivatives at test-midpoint 1 */
    double *deriv3;             /* the derivatives at test-midpoint 2 */
    double *deriv4;             /* the derivatives at test-endpoint */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double th;                  /* time at the end of the step */
    double thh;                 /* time for interval midpoints */

    double m;                   /* used to calculate number of steps */
    double stepsize;            /* real stepsize */
    double hh;                  /* half the stepsize */
    double h6;                  /* one sixth of the stepsize (for Rk4 formula) */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* number of steps we have to take */

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );
    deriv3 = ( double * ) calloc( n, sizeof( double ) );
    deriv4 = ( double * ) calloc( n, sizeof( double ) );

    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */
    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Rk4: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    for( step = 0; step < nsteps; step++ ) {    /* loop for steps */

        thh = t + hh;           /* time of interval midpoints */
        th = t + stepsize;      /* time at end of the interval */

        /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

        p_deriv( vnow, t, deriv1, n, si, inp );

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )
            vtemp[i] = vnow[i] + hh * deriv1[i];
        p_deriv( vtemp, thh, deriv2, n, si, inp );

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )
            vtemp[i] = vnow[i] + hh * deriv2[i];
        p_deriv( vtemp, thh, deriv3, n, si, inp );

        if( debug )
            nd++;

        for( i = 0; i < n; i++ )
            vtemp[i] = vnow[i] + stepsize * deriv3[i];
        p_deriv( vtemp, th, deriv4, n, si, inp );

        if( debug )
            nd++;

        /* ... then feed them to the Fourth-Order Runge-Kutta formula */

        for( i = 0; i < n; i++ )
            vnext[i] = vnow[i]
                + h6 * ( deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i] );

        /* next step */

        t += stepsize;


        if( step < nsteps - 2 ) {       /* CASE 1: many steps to go */
            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
            vnow = v[toggle];   /* set vout */
            vnext = vout;

        } else if( step > nsteps - 2 ) {        /* CASE 3: just did final iteration */
            free( v[0] );       /* clean up and go home! */
            free( v[1] );
            free( v );
            free( vtemp );
            free( deriv1 );
            free( deriv2 );
            free( deriv3 );
            free( deriv4 );
        }
    }

    if( debug )
        WriteSolvLog( "Rk4", tin, tout, stepsize, nsteps, nd, slog );

    return;
}

/** Rkck: propagates vin (of size n) from tin to tout by the Runge-Kutta 
 *         Cash-Karp method, which is an adaptive-stepsize Rk method; it   
 *         uses a fifth-order Rk formula with an embedded forth-oder for-  
 *         mula for calucalting the error; its result is returned by vout  
 ***************************************************************************
 *                                                                         
 * This solver was written by Marcel Wolf, Spring 2002.                    
 *                                                                         
 */
void
Rkck( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {

    int i;                      /* local loop counter */
    double **v; /** used for storing intermediate steps */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *verror;             /* error estimate */
    double verror_max;          /* the maximum error in verror */

    double *deriv1;             /* intermediate derivatives for the Cash-Karp formula */
    double *deriv2;
    double *deriv3;
    double *deriv4;
    double *deriv5;
    double *deriv6;

    double t;                   /* the current time */

    double h = stephint;        /* initial stepsize */
    double hnext;               /* used to calculate next stepsize */
    double SAFETY = 0.9;        /* safety margin for decreasing stepsize */



    /* declare and initialize Cash-Karp parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;

    static double
        b21 = 0.2,
        b31 = 3.0 / 40.0,
        b32 = 9.0 / 40.0,
        b41 = 0.3,
        b42 = -0.9,
        b43 = 1.2,
        b51 = -11.0 / 54.0,
        b52 = 2.5,
        b53 = -70.0 / 27.0, b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0;

    static double c1 = 37.0 / 378.0,
        /*  c2  =     0.0 */
        c3 = 250.0 / 621.0, c4 = 125.0 / 594.0,
        /*  c5  =     0.0 */
        c6 = 512.0 / 1771.0;

    double dc1 = c1 - 2825.0 / 27648.0,
        /*  dc2 =          0.0 */
        dc3 = c3 - 18575.0 / 48384.0, dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.0 / 14336.0, dc6 = c6 - 0.25;



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    verror = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );
    deriv3 = ( double * ) calloc( n, sizeof( double ) );
    deriv4 = ( double * ) calloc( n, sizeof( double ) );
    deriv5 = ( double * ) calloc( n, sizeof( double ) );
    deriv6 = ( double * ) calloc( n, sizeof( double ) );

    t = tin;
    vnow = vin;
    vnext = v[0];

    /* initial stepsize cannot be bigger than total time */

    if( tin + h >= tout )
        h = tout - tin;
    while( t < tout ) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */
        while( 1 ) {

            /* do the Runge-Kutta thing here: calulate intermediate derivatives */
            p_deriv( vnow, t, deriv1, n, si, inp );
            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b21 * deriv1[i] );
            p_deriv( vtemp, t + a2 * h, deriv2, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b31 * deriv1[i] + b32 * deriv2[i] );
            p_deriv( vtemp, t + a3 * h, deriv3, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b41 * deriv1[i] + b42 * deriv2[i] + b43 * deriv3[i] );
            p_deriv( vtemp, t + a4 * h, deriv4, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b51 * deriv1[i] + b52 * deriv2[i] + b53 * deriv3[i]
                                           + b54 * deriv4[i] );
            p_deriv( vtemp, t + a5 * h, deriv5, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b61 * deriv1[i] + b62 * deriv2[i] + b63 * deriv3[i]
                                           + b64 * deriv4[i] + b65 * deriv5[i] );
            p_deriv( vtemp, t + a6 * h, deriv6, n, si, inp );

            /* ... then feed them to the Cash-Karp formula */
            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h * ( c1 * deriv1[i] + c3 * deriv3[i] + c4 * deriv4[i]
                            + c6 * deriv6[i] );

            /* calculate the error estimate using the embedded formula */

            for( i = 0; i < n; i++ )
                verror[i] = h * ( dc1 * deriv1[i] + dc3 * deriv3[i]
                                  + dc4 * deriv4[i] + dc5 * deriv5[i]
                                  + dc6 * deriv6[i] );

            /* find the maximum error */

            verror_max = 0.;
            for( i = 0; i < n; i++ ) {
                //printf("calculation was wade for %d: %lg / %lg = %lg abs=%lg\n", i, verror[i], vnext[i], (verror[i]/vnext[i]), (fabs(verror[i]/vnext[i])));
                if( vnext[i] != 0. ) {
                    verror_max = DMAX( fabs( verror[i] / vnext[i] ), verror_max );
                    
                } else {
                    verror_max = DMAX( verror[i] / DBL_EPSILON, verror_max );
                }
                //printf("verrormax(%d)=%lg\n", i, verror_max);
            }

            /* scale error according to desired accuracy */
            //printf("verrormax_b = %lg\n", verror_max);
            verror_max /= accuracy;
            //printf("verrormax_a = %lg\n", verror_max);

            /* compare maximum error to the desired accuracy; if error < accuracy, we  *
             * are done with this step; otherwise, the stepsize has to be reduced and  *
             * the step repeated; for detailed comments on the approximation involving *
             * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

            if( verror_max <= 1.0 )
                break;

            hnext = SAFETY * h * pow( verror_max, -0.25 );
            /* decrease stepsize by no more than a factor of 10; check for underflows */
            h = ( hnext > 0.1 * h ) ? hnext : 0.1 * h;
            //printf("h=%lg\n", h);
            if( h < DBL_EPSILON )
                error( "Rkck: stepsize underflow" );
        }
        /* advance the current time by last stepsize */

        t += h;

        if( t >= tout )
            break;              /* that was the last iteration */

        /* increase stepsize according to error (5th order) for next iteration */

        h = h * pow( verror_max, -0.20 );

        /* make sure t does not overstep tout */

        if( t + h >= tout )
            h = tout - t;

        /* toggle vnow and vnext between v[0] and v[1] */

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];
    }
    /* copy the last result to vout after the final iteration */

    memcpy( vout, vnext, sizeof( *vnext ) * n );

    free( v[0] );
    free( v[1] );
    free( v );
    free( vtemp );
    free( verror );
    free( deriv1 );
    free( deriv2 );
    free( deriv3 );
    free( deriv4 );
    free( deriv5 );
    free( deriv6 );

    //exit(1);
}

/** Rkf: propagates vin (of size n) from tin to tout by the Runge-Kutta 
 *        Fehlberg method, which is a the original adaptive-stepsize Rk    
 *        method (Cash-Karp is an improved version of this); it uses a     
 *        fifth-order Rk formula with an embedded forth-oder formula for   
 *        calucalting the error; its result is returned by vout            
 ***************************************************************************
 *                                                                         
 * This solver was written by Marcel Wolf, Spring 2002.                    
 *                                                                         
 */
void
Rkf( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* used for storing intermediate steps */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */

    double *verror;             /* error estimate */
    double verror_max;          /* the maximum error in verror */

    double *deriv1;             /* intermediate derivatives for the Fehlberg formula */
    double *deriv2;
    double *deriv3;
    double *deriv4;
    double *deriv5;
    double *deriv6;

    double t;                   /* the current time */

    double h = stephint;        /* initial stepsize */
    double hnext;               /* used to calculate next stepsize */
    const double SAFETY = 0.9;  /* safety margin for decreasing stepsize */



    /* declare and initialize Fehlberg parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double a2 = 0.25, a3 = 0.375, a4 = 12.0 / 13.0, a5 = 1.0, a6 = 0.5;

    static double
        b21 = 0.25,
        b31 = 3.0 / 32.0,
        b32 = 9.0 / 32.0,
        b41 = 1932.0 / 2197.0,
        b42 = -7200.0 / 2197.0,
        b43 = 7296.0 / 2197.0,
        b51 = 439.0 / 216.0,
        b52 = -8.0,
        b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0, b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

    static double c1 = 25.0 / 216.0,
        /*  c2  =     0.0 */
        c3 = 1408.0 / 2565.0, c4 = 2197.0 / 4104.0, c5 = -0.2;
    /*  c6  =     0.0 */

    static double dc1 = 1.0 / 360.0,
        /*  dc2 =     0.0 */
        dc3 = -128.0 / 4275.0, dc4 = -2197.0 / 75240.0, dc5 = 1.0 / 50.0, dc6 = 2.0 / 55.0;



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    verror = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );
    deriv3 = ( double * ) calloc( n, sizeof( double ) );
    deriv4 = ( double * ) calloc( n, sizeof( double ) );
    deriv5 = ( double * ) calloc( n, sizeof( double ) );
    deriv6 = ( double * ) calloc( n, sizeof( double ) );

    t = tin;
    vnow = vin;
    vnext = v[0];

    /* initial stepsize cannot be bigger than total time */

    if( tin + h >= tout )
        h = tout - tin;

    while( t < tout ) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */

        while( 1 ) {

            /* do the Runge-Kutta thing here: calulate intermediate derivatives */

            p_deriv( vnow, t, deriv1, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b21 * deriv1[i] );
            p_deriv( vtemp, t + a2 * h, deriv2, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b31 * deriv1[i] + b32 * deriv2[i] );
            p_deriv( vtemp, t + a3 * h, deriv3, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b41 * deriv1[i] + b42 * deriv2[i] + b43 * deriv3[i] );
            p_deriv( vtemp, t + a4 * h, deriv4, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b51 * deriv1[i] + b52 * deriv2[i] + b53 * deriv3[i]
                                           + b54 * deriv4[i] );
            p_deriv( vtemp, t + a5 * h, deriv5, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + h * ( b61 * deriv1[i] + b62 * deriv2[i] + b63 * deriv3[i]
                                           + b64 * deriv4[i] + b65 * deriv5[i] );
            p_deriv( vtemp, t + a6 * h, deriv6, n, si, inp );

            /* ... then feed them to the Fehlberg formula */

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h * ( c1 * deriv1[i] + c3 * deriv3[i] + c4 * deriv4[i]
                            + c5 * deriv5[i] );

            /* calculate the error estimate using the embedded formula */

            for( i = 0; i < n; i++ )
                verror[i] = fabs( h * ( dc1 * deriv1[i] + dc3 * deriv3[i]
                                        + dc4 * deriv4[i] + dc5 * deriv5[i]
                                        + dc6 * deriv6[i] ) );

            /* find the maximum error */

            verror_max = 0.;
            for( i = 0; i < n; i++ )
                if( vnext[i] != 0. )
                    verror_max = DMAX( verror[i] / vnext[i], verror_max );
                else
                    verror_max = DMAX( verror[i] / DBL_EPSILON, verror_max );

            /* scale error according to desired accuracy */

            verror_max /= accuracy;

            /* compare maximum error to the desired accuracy; if error < accuracy, we  *
             * are done with this step; otherwise, the stepsize has to be reduced and  *
             * the step repeated; for detailed comments on the approximation involving *
             * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

            if( verror_max <= 1.0 )
                break;

            hnext = SAFETY * h * pow( verror_max, -0.25 );

            /* decrease stepsize by no more than a factor of 10; check for underflows */

            h = ( hnext > 0.1 * h ) ? hnext : 0.1 * h;
            if( h < DBL_EPSILON )
                error( "Rkf: stepsize underflow" );

        }

        /* advance the current time by last stepsize */

        t += h;

        if( t >= tout )
            break;              /* that was the last iteration */

        /* increase stepsize according to error for next iteration */

        h = h * pow( verror_max, -0.20 );

        /* make sure t does not overstep tout */

        if( t + h >= tout )
            h = tout - t;

        /* toggle vnow and vnext between v[0] and v[1] */

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

    }

    /* copy the last result to vout after the final iteration */

    memcpy( vout, vnext, sizeof( *vnext ) * n );

    free( v[0] );
    free( v[1] );
    free( v );
    free( vtemp );
    free( verror );
    free( deriv1 );
    free( deriv2 );
    free( deriv3 );
    free( deriv4 );
    free( deriv5 );
    free( deriv6 );

}

/**    Milne: propagates vin (of size n) from tin to tout by Milne-Simpson 
 *            which is a predictor-corrector method; the result is retur-  
 *            ned by vout                                                  
 ***************************************************************************
 *                                                                         
 * This solver was implemented by Konstantin Koslov, Dec 2001/Jan 2002     
 *                                                                         
 ***************************************************************************
 *                                                                         
 * THIS SOLVER SEEMS TO BE BUGGY FOR SOME REASON, DO NOT USE IT!!!!!       
 *                                                                         
 */
void
Milne( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* array to store intermediate results */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* current protein concs */
    double *vnext;              /* next protein concs */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* derivatives at test-midpoint 1 */
    double *deriv3;             /* derivatives at test-midpoint 2 */
    double *deriv4;             /* derivatives at test-endpoint */

    double **history_dv;        /* history of the derivatives */
    double **history_v;         /* history of the v: v_{i-3} v_{i-1} */
    double *dblank;             /* pointer for history_dv rotation */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double th;                  /* time at end of interval */
    double thh;                 /* time for interval midpoints */

    double m;                   /* used to calculate number of steps */
    double stepsize;            /* real stepsize */
    double hh;                  /* half the stepsize */
    double h6;                  /* one sixth of the stepsize (for Rk4 formula) */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* number of steps we have to take */

    double hp;                  /* multiplier for the predictor derivs */
    double hc;                  /* multiplier for the corrector derivs */

    double mistake;             /* error estimator */
    double corr;                /* temp variable for the corrector */
    double ee;                  /* temp variable for the error estimator */
    double cc29;                /* 1/29: used for calculating ee */

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout ) {
        return;
    }

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );
    deriv3 = ( double * ) calloc( n, sizeof( double ) );
    deriv4 = ( double * ) calloc( n, sizeof( double ) );

    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */

    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Milne: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    /* do we need Milne? no, not for less than 4 steps! */

    if( nsteps < 4 ) {

        for( step = 0; step < nsteps; step++ ) {        /* loop for steps */

            thh = t + hh;       /* time of interval midpoints */
            th = t + stepsize;  /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            p_deriv( vnow, t, deriv1, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv1[i];
            p_deriv( vtemp, thh, deriv2, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv2[i];
            p_deriv( vtemp, thh, deriv3, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            p_deriv( vtemp, th, deriv4, n, si, inp );

            if( debug )
                nd++;


            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h6 * ( deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i] );

            /* next step */

            t += stepsize;

            if( step < nsteps - 2 ) {   /* CASE 1: many steps to go */
                vnow = v[toggle];       /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];

            } else if( step == nsteps - 2 ) {   /* CASE 2: next iteration = final */
                vnow = v[toggle];       /* set vout */
                vnext = vout;

            } else if( step > nsteps - 2 ) {    /* CASE 3: just did final iteration */
                free( v[0] );   /* clean up and go home! */
                free( v[1] );
                free( v );
                free( deriv1 );
                free( deriv2 );
                free( deriv3 );
                free( deriv4 );
                free( vtemp );
            }
        }

        /* more than 3 steps: yes, we need Milne */

    } else {

        history_dv = ( double ** ) calloc( 3, sizeof( double * ) );
        history_v = ( double ** ) calloc( 4, sizeof( double * ) );

        for( step = 0; step < 3; step++ ) {
            history_dv[step] = ( double * ) calloc( n, sizeof( double ) );
        }

        for( step = 0; step < 4; step++ ) {
            history_v[step] = ( double * ) calloc( n, sizeof( double ) );
        }

        mistake = 0.;
        cc29 = 1. / 29.;

        /* loop for initial steps using Rk4 */

        for( step = 0; step < 3; step++ ) {

            thh = t + hh;       /* time of interval midpoints */
            th = t + stepsize;  /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */
            p_deriv( vnow, t, deriv1, n, si, inp );

            if( debug )
                nd++;

            if( step > 0 )
                for( i = 0; i < n; i++ )
                    history_dv[step - 1][i] = deriv1[i];

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv1[i];
            p_deriv( vtemp, thh, deriv2, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )    /* evaluate deriv at second midpoint */
                vtemp[i] = vnow[i] + hh * deriv2[i];
            p_deriv( vtemp, thh, deriv3, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )    /* evaluate deriv at end point */
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            p_deriv( vtemp, th, deriv4, n, si, inp );

            if( debug )
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h6 * ( deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i] );

            /* next step */

            t += stepsize;

            for( i = 0; i < n; i++ )    /* save old v values in history_v */
                history_v[step][i] = vnow[i];

            if( step == 1 ) {
                vnow = vnext;
                vnext = vout;
            } else {
                vnow = v[toggle];       /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];
            }
        }

        /* we have calculated 4 initial points with rk4: start multistepping */

        hc = stepsize / 3.0;    /* stepsize multipliers: for the corrector */
        hp = 4. * hc;           /* ... and for the predictor */

        for( i = 0; i < n; i++ )        /* save current values of v in history_v */
            history_v[3][i] = vout[i];

        for( step = 3; step < nsteps; step++ ) {        /* loop for Milne steps */

            th = t + stepsize;  /* time at end of the interval */

            /* do the Milne thing here */

            p_deriv( vout, t, deriv1, n, si, inp );     /* evaluate deriv at start of interval */

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )    /* save current derivs in history_dv */
                history_dv[2][i] = deriv1[i];

            /* Predictor step (P): extrapolate derivatives */

            for( i = 0; i < n; i++ )
                v[0][i] = history_v[0][i]
                    + hp * ( 2 * history_dv[0][i] - history_dv[1][i]
                             + 2 * history_dv[2][i] );

            /* Evaluate the extrapolated derivatives (E) */

            p_deriv( v[0], th, deriv1, n, si, inp );

            if( debug )
                nd++;

            /* Corrector step (C): interpolate current v (+ calculate error estimator) */

            for( i = 0; i < n; i++ ) {

                corr = history_v[2][i]
                    + hc * ( history_dv[1][i] + 4 * history_dv[2][i] + deriv1[i] );
                ee = fabs( corr - v[0][i] ) * cc29;     /* error estimator */

                if( ee > mistake )
                    mistake = ee;
                vout[i] = corr;

            }

            t += stepsize;      /* go to next step */

            /* update the arrays */

            if( step <= nsteps - 2 ) {  /* CASE 1: many steps to go */

                dblank = history_dv[0];
                history_dv[0] = history_dv[1];
                history_dv[1] = history_dv[2];
                history_dv[2] = dblank;

                dblank = history_v[0];
                history_v[0] = history_v[1];
                history_v[1] = history_v[2];
                history_v[2] = history_v[3];
                history_v[3] = dblank;

                for( i = 0; i < n; i++ )
                    history_v[3][i] = vout[i];

            } else {            /* CASE 2: just did final iteration */

                free( v[0] );   /* clean up and go home! */
                free( v[1] );
                free( v );
                free( history_v[0] );
                free( history_v[1] );
                free( history_v[2] );
                free( history_v[3] );
                free( history_v );
                free( history_dv[0] );
                free( history_dv[1] );
                free( history_dv[2] );
                free( history_dv );
                free( deriv1 );
                free( deriv2 );
                free( deriv3 );
                free( deriv4 );
                free( vtemp );
            }
        }

        /* use below for sanity check if Milne behaves unstably */

        /* 
           fprintf(stdout,"mistake=%.15f\n", mistake);
           fflush(stdout);
           getchar();
         */

    }

    if( debug )
        WriteSolvLog( "Milne", tin, tout, stepsize, nsteps, nd, slog );

    return;

}

/** Adams: propagates vin (of size n) from tin to tout by Adams-Moulton 
 *            which is an implicit predictor-corrector method of second    
 *            order; the result is returned by vout                        
 ***************************************************************************
 *                                                                         
 * This solver was implemented by Konstantin Koslov, Spring 2002           
 * Slightly modified by Manu, July 2002                                    
 *                                                                         
 */
void
Adams( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i;                      /* local loop counter */

    double **v;                 /* array to store intermediate results */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* current protein concs */
    double *vnext;              /* next protein concs */

    double *deriv1;             /* the derivatives at the beginning of a step */
    double *deriv2;             /* derivatives at test-midpoint 1 */
    double *deriv3;             /* derivatives at test-midpoint 2 */
    double *deriv4;             /* derivatives at test-endpoint */

    double **history_dv;        /* history of the derivatives */
    double *dblank;             /* pointer for history_dv rotation */

    double deltat;              /* tout - tin */
    double t;                   /* current time */
    double th;                  /* time at end of interval */
    double thh;                 /* time for interval midpoints */

    double m;                   /* used to calculate number of steps */
    double stepsize;            /* real stepsize */
    double hh;                  /* half the stepsize */
    double h6;                  /* one sixth of the stepsize (for Rk4 formula) */
    int step;                   /* loop counter for steps */
    int nsteps;                 /* number of steps we have to take */

    double mistake;             /* error estimators */
    //  double eps1;
    //  double eps2;

    int nd = 0;                 /* number of deriv evaluations */



    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );
    deriv1 = ( double * ) calloc( n, sizeof( double ) );
    deriv2 = ( double * ) calloc( n, sizeof( double ) );
    deriv3 = ( double * ) calloc( n, sizeof( double ) );
    deriv4 = ( double * ) calloc( n, sizeof( double ) );

    deltat = tout - tin;        /* how far do we have to propagate? */
    m = floor( deltat / stephint + 0.5 );       /* see comment on stephint above */

    if( m < 1. )
        m = 1.;                 /* we'll have to do at least one step */

    stepsize = deltat / m;      /* real stepsize calculated here */
    nsteps = ( int ) m;         /* int number of steps */
    t = tin;                    /* set current time */

    if( t == t + stepsize )
        error( "Adams: stephint of %g too small!", stephint );

    vnow = vin;

    if( nsteps == 1 )           /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    /* do we need Adams? no, not for less than 4 steps! */

    if( nsteps < 4 ) {

        for( step = 0; step < nsteps; step++ ) {        /* loop for steps */

            thh = t + hh;       /* time of interval midpoints */
            th = t + stepsize;  /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            p_deriv( vnow, t, deriv1, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv1[i];
            p_deriv( vtemp, thh, deriv2, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv2[i];
            p_deriv( vtemp, thh, deriv3, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            p_deriv( vtemp, th, deriv4, n, si, inp );

            if( debug )
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h6 * ( deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i] );

            /* next step */

            t += stepsize;

            if( step < nsteps - 2 ) {   /* CASE 1: many steps to go */
                vnow = v[toggle];       /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];

            } else if( step == nsteps - 2 ) {   /* CASE 2: next iteration = final */
                vnow = v[toggle];       /* set vout */
                vnext = vout;

            } else if( step > nsteps - 2 ) {    /* CASE 3: just did final iteration */
                free( v[0] );   /* clean up and go home! */
                free( v[1] );
                free( v );
                free( vtemp );
                free( deriv1 );
                free( deriv2 );
                free( deriv3 );
                free( deriv4 );
            }
        }

    } else {

        /* more than 3 steps: yes, we need Adams */

        history_dv = ( double ** ) calloc( 4, sizeof( double * ) );
        for( step = 0; step < 4; step++ ) {
            history_dv[step] = ( double * ) calloc( n, sizeof( double ) );
        }

        mistake = 0.;

        /* loop for initial steps using Rk4 */

        for( step = 0; step < 3; step++ ) {

            thh = t + hh;       /* time of interval midpoints */
            th = t + stepsize;  /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            p_deriv( vnow, t, deriv1, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                history_dv[step][i] = deriv1[i];

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv1[i];
            p_deriv( vtemp, thh, deriv2, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + hh * deriv2[i];
            p_deriv( vtemp, thh, deriv3, n, si, inp );

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            p_deriv( vtemp, th, deriv4, n, si, inp );

            if( debug )
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + h6 * ( deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i] );

            /* next step */

            t += stepsize;

            vnow = v[toggle];   /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        }

        /* we have calculated 4 initial points with rk4: start multistepping */

        hh = stepsize / 24.0;   /* reset hh to 1/24 of stepsize */

        if( nsteps == 4 )
            vnext = vout;

        for( step = 3; step < nsteps; step++ ) {        /* loop for steps */

            th = t + stepsize;  /* time at end of the interval */

            /* do the Adams thing here */

            p_deriv( vnow, t, deriv1, n, si, inp );     /* evaluate deriv at start of interval */

            if( debug )
                nd++;

            for( i = 0; i < n; i++ )
                history_dv[3][i] = deriv1[i];

            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + hh * ( 55 * history_dv[3][i] - 59 * history_dv[2][i]
                             + 37 * history_dv[1][i] - 9 * history_dv[0][i] );

            p_deriv( vnext, th, deriv1, n, si, inp );

            if( debug )
                nd++;


            /* the following evaluates the error, but then nothing is done with it */

            /*    for(i=0; i<n; i++) {              
               eps1 = hh * (           deriv1[i] - 3 * history_dv[3][i]
               + 3 * history_dv[2][i] -     history_dv[1][i]);
               eps2 = fabs(eps1);
               if (eps2 > mistake)
               mistake=eps2;
               } 
             */
            for( i = 0; i < n; i++ )
                vnext[i] = vnow[i]
                    + hh * ( 9 * deriv1[i] + 19 * history_dv[3][i]
                             - 5 * history_dv[2][i] + history_dv[1][i] );

            t += stepsize;      /* go to next step */

            if( step <= nsteps - 2 ) {  /* CASE 1: many steps to go */
                dblank = history_dv[0];
                history_dv[0] = history_dv[1];
                history_dv[1] = history_dv[2];
                history_dv[2] = history_dv[3];
                history_dv[3] = dblank;

                if( step < nsteps - 2 ) {
                    vnow = v[toggle];   /* toggle results from v[0] and v[1] */
                    toggle++;
                    toggle %= 2;
                    vnext = v[toggle];

                } else if( step == nsteps - 2 ) {       /* CASE 2: next iteration = final */
                    vnow = v[toggle];   /* set vout */
                    vnext = vout;
                }

            } else {            /* CASE 2: just did final iteration */

                free( v[0] );   /* clean up and go home! */
                free( v[1] );
                free( v );
                free( history_dv[0] );
                free( history_dv[1] );
                free( history_dv[2] );
                free( history_dv[3] );
                free( history_dv );
                free( deriv1 );
                free( deriv2 );
                free( deriv3 );
                free( deriv4 );
                free( vtemp );
            }
        }
    }

    if( debug )
        WriteSolvLog( "Adams", tin, tout, stepsize, nsteps, nd, slog );

    return;

}

/** BuSt: propagates v(t) from t1 to t2 by Bulirsch-Stoer; this method 
 *           uses Richardson extrapolation to estimate v's at a hypothe-   
 *           tical stepsize of 0; the extrapolation also yields an error   
 *           estimate, which is used to adapt stepsize and change the or-  
 *           der of the method as required; the result is returned by vout 
 ***************************************************************************
 *                                                                         
 * This solver was implemented by Manu, July 2002                          
 *                                                                         
 */
void
BuSt( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {

    /* bsstep is the BuSt stepper function */
    void bsstep( double *v, double *deriv, int n, double *t, double htry, double accuracy, double *hdid, double *hnext, SolverInput * si, Input * inp );

    int i;                      /* local loop counter */

    double t;                   /* current time */

    double ht;                  /* stepsize of next step */
    double hd;                  /* stepsize of step we just did */

    double *initderiv;          /* array for initial derivatives */


    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough -> initialize variables */

    initderiv = ( double * ) calloc( n, sizeof( double ) );

    t = tin;

    /* initial stepsize is either stephint or tout-tin */

    ht = DMIN( stephint, ( tout - tin ) );

    /* we're saving old v's in vin and propagate vout below */

    for( i = 0; i < n; i++ )
        vout[i] = vin[i];

    /* this loop steps through the whole interval tout-tin */

    do {

        p_deriv( vout, t, initderiv, n, si, inp );
        bsstep( vout, initderiv, n, &t, ht, accuracy, &hd, &ht, si, inp );

        if( t + ht > tout )
            ht = tout - t;

    } while( t < tout );

    /* clean up and go home */

    free( initderiv );

    return;
}

/** bsstep: this is the Bulirsch-Stoer stepper function; it propagates 
 *           the equations at a given overall stepsize and order, then     
 *           evaluates the error and repeats the step with increased order 
 *           if necessary, until the required accuracy is achieved; its    
 *           arguments are:                                                
 *           - v:        input and output v's over a give total stepsize   
 *           - deriv:    array of initial derivatives (at tin)             
 *           - n:        size of v and deriv                               
 *           - t:        pointer to the current time                       
 *           - htry:     initial stepsize to be tried                      
 *           - accuracy: relative accuracy to be achieved                  
 *           - hdid:     returns the stepsize of the step we just did      
 *           - hnext:    returns the stepsize of the next step to be done  
 */
void
bsstep( double *v, double *deriv, int n, double *t, double htry, double accuracy, double *hdid, double *hnext, SolverInput * si, Input * inp ) {

    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *            mmid: implements the modified midpoint method                *
     *          pzextr: implements the Richardson extrapolation (polynomial)   */

    void mmid( double *vin, double *vout, double *deriv, double tin, double htot, int nstep, int n, SolverInput * si, Input * inp );

    void pzextr( int iest, double hest, double *yest, double *yz, double *dy, int nv, const int KMAXX );

    /*** constants *************************************************************/

    const int KMAXX = 8;        /* this is how many successive stepsize */
    const int IMAXX = KMAXX + 1;        /* refinements we're going to try */

    const double SAFE1 = 0.25;  /* safety margins for errors and h reduction */
    const double SAFE2 = 0.7;

    const double REDMAX = 1.0e-5;       /* boundaries for stepsize reduction */
    const double REDMIN = 0.7;

    const double SCALMX = 0.1;  /* maximum scaling value for stepsize */

    /*** variables *************************************************************/

    int i, k, km;               /* loop counters */

    double *vsav;               /* v's at beginning of interval */
    double *vseq;               /* v's at end of interval for a given stepsize */

    double *verror;             /* error estimates */
    double verror_max;          /* the maximum error in verror */
    double *err;                /* normalized error estimates */

    double h;                   /* stepsize for solver method used below */
    double hest;                /* square of h passed to extrapolation function */
    double red = 0.0;           /* reduction factor for reducing stepsize */
    double scale = 0.0;         /* scaling factor for next stepsize */

    double work;                /* temp vars for calculating row for convergence */
    double wrkmin;
    double fact;

    static double old_accuracy = -1.0;  /* used to save old accuracy */
    static double tnew;         /* used to save old start time */

    /* the following two arrays are used for Deuflhard's error estimation; a   *
     * contains the work coefficients and alf (alpha) the correction factors   */

    static double *a = NULL;
    static double **alf = NULL;

    double accuracy1;           /* error (< accuracy) used to calculate alphas */

    static int kmax;            /* used for finding kopt */
    static int kopt;            /* optimal row number for convergence */

    /* sequence of separate attempts to cross interval htot with increasing    *
     * values of nsteps as suggested by Deuflhard (Num Rec, p. 726)            */

    static int nseq[] = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18 };

    /* miscellaneous flags */

    static int first = 1;       /* is this the first try for a given step? */
    int reduct;                 /* flag indicating if we have reduced stepsize yet */
    int exitflag = 0;           /* exitflag: when set, we exit (!) */

    /* static global arrays */

    extern double *d;           /* D's used for extrapolation in pzextr */
    extern double *hpoints;     /* stepsizes h (=H/n) which we have tried */

    /* allocate arrays */

    if( !( err = ( double * ) calloc( KMAXX, sizeof( double ) ) ) )
        error( "BuSt: error allocating err.\n" );

    if( !( verror = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BuSt: error allocating verror.\n" );

    if( !( vsav = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BuSt: error allocating vsav.\n" );

    if( !( vseq = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BuSt: error allocating vseq.\n" );

    if( !( d = ( double * ) calloc( KMAXX * n, sizeof( double ) ) ) )
        error( "BuSt: error allocating d.\n" );

    if( !( hpoints = ( double * ) calloc( KMAXX, sizeof( double ) ) ) )
        error( "BuSt: error allocating hpoints.\n" );

    /* set initial stepsize to try to htry */

    h = htry;

    /* new accuracy? -> initialize (here, this only applies at start) */

    if( accuracy != old_accuracy ) {

        *hnext = tnew = -1.0e29;        /* init these to impossible value */

        /* allocate memory if necessary */

        if( !a ) {
            a = ( double * ) calloc( IMAXX + 1, sizeof( double ) );
            alf = ( double ** ) calloc( KMAXX + 1, sizeof( double * ) );
            for( i = 0; i < KMAXX + 1; i++ )
                alf[i] = ( double * ) calloc( KMAXX + 1, sizeof( double ) );
        }

        /* initialize the work coefficients */

        a[1] = nseq[1] + 1;
        for( k = 1; k <= KMAXX; k++ )
            a[k + 1] = a[k] + nseq[k + 1];

        /* initialize the correction factors (alpha) */

        accuracy1 = SAFE1 * accuracy;   /* accuracy1 used to calculate alphas */
        for( i = 2; i <= KMAXX; i++ )
            for( k = 1; k < i; k++ )
                alf[k][i] = pow( accuracy1, ( ( a[k + 1] - a[i + 1] ) / ( ( a[i + 1] - a[1] + 1.0 ) * ( 2 * k + 1 ) ) ) );
        old_accuracy = accuracy;

        /* determine optimal row number for convergence */

        for( kopt = 2; kopt < KMAXX; kopt++ )
            if( a[kopt + 1] > a[kopt] * alf[kopt - 1][kopt] )
                break;
        kmax = kopt;

    }

    /* actual stepping starts here: first save starting values of v[] in vsav  */

    for( i = 0; i < n; i++ )
        vsav[i] = v[i];

    /* new integration or stepsize? -> re-establish the order window */

    if( *t != tnew || h != ( *hnext ) ) {
        first = 1;
        kopt = kmax;
    }

    reduct = 0;                 /* we have not reduced stepsize yet */

    for( ;; ) {

        /* try increasing numbers of steps over the interval h */

        for( k = 1; k <= kmax; k++ ) {

            tnew = ( *t ) + h;
            if( tnew == ( *t ) )
                error( "BuSt: stepsize underflow in bsstep\n" );

            /* propagate the equations from t to t+h with nseq[k] steps using vsav and *
             * deriv as input; mmid returns vseq                                       */

            mmid( vsav, vseq, deriv, *t, h, nseq[k], n, si, inp );

            /* extrapolate v's for h->0 based on the different stepsizes (hest's) we   *
             * have tried already; the result is returned in v (errors in verror);     *
             * hest is squared since the error series is even                          */

            hest = DSQR( h / nseq[k] );
            pzextr( k - 1, hest, vseq, v, verror, n, KMAXX );

            if( k > 1 ) {

                /* find the maximum error */

                verror_max = 0.;
                for( i = 0; i < n; i++ )
                    if( v[i] != 0. )
                        verror_max = DMAX( fabs( verror[i] / v[i] ), verror_max );
                    else
                        verror_max = DMAX( fabs( verror[i] / DBL_EPSILON ), verror_max );

                /* scale error according to desired accuracy */

                verror_max /= accuracy;

                /* compute normalized error estimates (epsilons) */

                km = k - 1;
                err[km - 1] = pow( verror_max / SAFE1, 1.0 / ( 2 * km + 1 ) );

            }

            /* are we in the order window? -> converged */

            if( k > 1 && ( k >= kopt - 1 || first ) ) {

                if( verror_max < 1.0 ) {        /* exit if accuracy good enough */
                    exitflag = 1;
                    break;
                }

                if( k == kmax || k == kopt + 1 ) {      /* stepsize reduction possible? */
                    red = SAFE2 / err[km - 1];
                    break;
                } else if( k == kopt && alf[kopt - 1][kopt] < err[km - 1] ) {
                    red = 1.0 / err[km - 1];
                    break;
                } else if( kopt == kmax && alf[km][kmax - 1] < err[km - 1] ) {
                    red = alf[km][kmax - 1] * SAFE2 / err[km - 1];
                    break;
                } else if( alf[km][kopt] < err[km - 1] ) {
                    red = alf[km][kopt - 1] / err[km - 1];
                    break;
                }
            }
        }

        if( exitflag )
            break;

        /* reduce stepsize by at least REDMIN and at most REDMAX, then try again */

        red = DMIN( red, REDMIN );
        red = DMAX( red, REDMAX );

        h *= red;
        reduct = 1;
    }

    /* we've taken a successful step */

    *t = tnew;
    *hdid = h;
    first = 0;
    wrkmin = 1.0e35;

    /* compute optimal row for convergence and corresponding stepsize */

    for( i = 1; i <= km; i++ ) {

        fact = DMAX( err[i - 1], SCALMX );
        work = fact * a[i + 1];

        if( work < wrkmin ) {

            scale = fact;
            wrkmin = work;
            kopt = i + 1;

        }
    }

    *hnext = h / scale;

    /* check for possible order increase, but not if stepsize was just reduced */

    if( kopt >= k && kopt != kmax && !reduct ) {

        fact = DMAX( scale / alf[kopt - 1][kopt], SCALMX );

        if( a[kopt + 1] * fact <= wrkmin ) {
            *hnext = h / fact;
            kopt++;
        }
    }

    /* clean up */

    free( d );
    free( hpoints );
    free( vseq );
    free( verror );
    free( err );
    free( vsav );
}

/** pzextr: uses polynomial extrapolation (Neville's algorithm) to evalu- 
 *           ate v's at a the hypothetical stepsize 0; this is called Ri-  
 *           chardson extrapolation; the arguments are:                    
 *           - iest: how many steps have we done already before?           
 *           - hest: current stepsize h                                    
 *           - vest: v's obtained using current stepsize h                  
 *           - vout: extrapolated v's that will be returned                
 *           - dv:   array of error estimates to be returned               
 *           - n:    size of verst, vout and dv arrays                     
 *           Neville's algorithm uses a recursive approach to determine a  
 *           suitable Lagrange polynomial for extrapolation by traversing  
 *           a tableau of differences between Lagrange polynomials of in-  
 *           creasing order; these differences are called C and D below    
 *           (see Numerical Recipes in C, Chapter 3.1 for details)         
 */
void
pzextr( int iest, double hest, double *vest, double *vout, double *dv, int n, const int KMAXX ) {
    int i, j;                   /* loop counters */

    double q;                   /* temp vars for calculating C's and D's */
    double f1, f2;
    double delta;
    double *c;                  /* C's used for extrapolation */

    extern double *d;           /* D's used for extrapolation */
    extern double *hpoints;     /* stepsizes h (=H/n) which we have tried */

    if( !( c = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "pzextr: error allocating c.\n" );

    hpoints[iest] = hest;       /* stepsizes h (=H/n) which we have tried */

    for( i = 0; i < n; i++ )
        dv[i] = vout[i] = vest[i];

    /* first time around? -> store first estimate in first column of d */

    if( iest == 0 )
        for( i = 0; i < n; i++ )
            d[i * KMAXX] = vest[i];

    /* more than one point? -> do polynomial extrapolation to h->0 */

    else {

        for( i = 0; i < n; i++ )
            c[i] = vest[i];

        for( i = 1; i <= iest; i++ ) {

            /* calculate new values of C's and D's to traverse Neville's tableau */

            delta = 1.0 / ( hpoints[iest - i] - hest );
            f1 = hest * delta;
            f2 = hpoints[iest - i] * delta;

            /* propagate tableau one diagonal more; the error is calculated using the  *
             * values of the D's */

            for( j = 0; j < n; j++ ) {
                q = d[j * KMAXX + ( i - 1 )];
                d[j * KMAXX + ( i - 1 )] = dv[j];
                delta = c[j] - q;
                dv[j] = f1 * delta;
                c[j] = f2 * delta;
                vout[j] += dv[j];


            }
        }

        /* save current D's for future calls to this function in the appropriate   *
         * column of the D tableau                                                 */

        for( i = 0; i < n; i++ )
            d[i * KMAXX + iest] = dv[i];

    }

    free( c );
}

/** mmid: implements the modified midpoint method used by BuSt; it subdi- 
 *         vides a large step (htot) into nstep intervals and uses a mid-  
 *         point method over the whole of htot, with a stepsize of 2*h     
 *         except for the first and the last step, hence *modified* mid-   
 *         point method); the nice thing about this method is that the     
 *         error follows a power law depending on the stepsize, i.e. its   
 *         error converges to zero really fast as h is diminished; this    
 *         allows for good extrapolation to h=0 (see pzextr() above)       
 */
void
mmid( double *vin, double *vout, double *deriv, double tin, double htot, int nstep, int n, SolverInput * si, Input * inp ) {

    int i, j;                   /* loop counters */

    double *vm;                 /* v's at beginning of 2*h */
    double *vn;                 /* v's at midpoint of 2*h */
    double swap;                /* tmp var used for swapping vm and vn */

    double t;                   /* current time */

    double h;                   /* small (h) stepsize (equals htot / nstep) */
    double h2;                  /* double the small (h) stepsize */


    vm = calloc( n, sizeof( double ) );
    vn = calloc( n, sizeof( double ) );

    /* calculate h from H and n */

    h = htot / nstep;

    /* calculate the first step (according to Euler) */

    for( i = 0; i < n; i++ ) {
        vm[i] = vin[i];
        vn[i] = vin[i] + h * deriv[i];
    }

    /* do the intermediate steps */

    h2 = 2.0 * h;               /* mmid uses 2 * h as stepsize for interm. steps */
    t = tin + h;
    p_deriv( vn, t, vout, n, si, inp );

    for( i = 2; i <= nstep; i++ ) {
        for( j = 0; j < n; j++ ) {
            swap = vm[j] + h2 * vout[j];
            vm[j] = vn[j];
            vn[j] = swap;
        }

        t += h;
        p_deriv( vn, t, vout, n, si, inp );
    }

    /* calculate the last step */

    for( i = 0; i < n; i++ )
        vout[i] = 0.5 * ( vm[i] + vn[i] + h * vout[i] );

    /* clean up */

    free( vm );
    free( vn );
}

/**    BaDe: propagates v(t) from t1 to t2 by Bader-Deuflhard; this method 
 *           uses Richardson extrapolation to estimate v's at a hypothe-   
 *           tical stepsize of 0; the extrapolation also yields an error   
 *           estimate, which is used to adapt stepsize and change the or-  
 *           der of the method as required; the result is returned by vout 
 ***************************************************************************
 *                                                                         
 * This solver was implemented by Yogi, based on BuSt, Aug 2002            
 *                                                                         
 */
void
BaDe( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {

    /* stifbs is the BaDe stepper function */

    void stifbs( double *v, double *deriv, int n, double *t, double htry, double accuracy, double *hdid, double *hnext, SolverInput * si, Input * inp );

    int i;                      /* local loop counter */

    double t;                   /* current time */

    double ht;                  /* stepsize of next step */
    double hd;                  /* stepsize of step we just did */

    double *initderiv;          /* array for initial derivatives */


    /* the do-nothing case; too small steps dealt with under usual */

    if( tin == tout )
        return;

    /* the usual case: steps big enough -> initialize variables */

    initderiv = ( double * ) calloc( n, sizeof( double ) );

    t = tin;

    /* initial stepsize is either stephint or tout-tin */

    ht = DMIN( stephint, ( tout - tin ) );

    /* we're saving old v's in vin and propagate vout below */

    for( i = 0; i < n; i++ )
        vout[i] = vin[i];

    /* this loop steps through the whole interval tout-tin */

    do {

        p_deriv( vout, t, initderiv, n, si, inp );
        stifbs( vout, initderiv, n, &t, ht, accuracy, &hd, &ht, si, inp );

        if( t + ht > tout )
            ht = tout - t;

    } while( t < tout );

    /* clean up and go home */

    free( initderiv );

    return;
}

/**  stifbs: this is the Bader-Deuflhard stepper function; it propagates 
 *           the equations at a given overall stepsize and order, then     
 *           evaluates the error and repeats the step with increased order 
 *           if necessary, until the required accuracy is achieved; its    
 *           arguments are:                                                
 *           - v:        input and output v's over a give total stepsize   
 *           - deriv:    array of initial derivatives (at tin)             
 *           - n:        size of v and deriv                               
 *           - t:        pointer to the current time                       
 *           - htry:     initial stepsize to be tried                      
 *           - accuracy: relative accuracy to be achieved                  
 *           - hdid:     returns the stepsize of the step we just did      
 *           - hnext:    returns the stepsize of the next step to be done  
 */
void
stifbs( double *v, double *deriv, int n, double *t, double htry, double accuracy, double *hdid, double *hnext, SolverInput * si, Input * inp ) {

    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *           simpr: implements the semi-implicit mid-point rule            *
     *          pzextr: implements the Richardson extrapolation (polynomial)   */

    void simpr( double *vin, double *vout, double *deriv, double *dfdt,
                double **jac, double tin, double htot, int nstep, int n, SolverInput * si, Input * inp );

    void pzextr( int iest, double hest, double *yest, double *yz, double *dy, int nv, const int KMAXX );

    /*** constants *************************************************************/

    const int KMAXX = 7;        /* this is how many successive stepsize */
    const int IMAXX = KMAXX + 1;        /* refinements we're going to try */

    const double SAFE1 = 0.25;  /* safety margins for errors and h reduction */
    const double SAFE2 = 0.7;

    const double REDMAX = 1.0e-5;       /* boundaries for stepsize reduction */
    const double REDMIN = 0.7;

    const double SCALMX = 0.1;  /* maximum scaling value for stepsize */

    /*** variables *************************************************************/

    int i, k, km;               /* loop counters */

    double *vsav;               /* v's at beginning of interval */
    double *vseq;               /* v's at end of interval for a given stepsize */

    double *verror;             /* error estimates */
    double verror_max;          /* the maximum error in verror */
    double *err;                /* normalized error estimates */

    double *dfdt;
    double **jac;               /* the Jacobian matrix */

    double h;                   /* stepsize for solver method used below */
    double hest;                /* square of h passed to extrapolation function */
    double red = 0.0;           /* reduction factor for reducing stepsize */
    double scale = 0.0;         /* scaling factor for next stepsize */

    double work;                /* temp vars for calculating row for convergence */
    double wrkmin;
    double fact;

    static double old_accuracy = -1.0;  /* used to save old accuracy */
    static double tnew;         /* used to save old start time */
    static int nold = -1;       /* for saving old value of n */

    /* the following two arrays are used for Deuflhard's error estimation; a   *
     * contains the work coefficients and alf (alpha) the correction factors   */

    static double *a = NULL;
    static double **alf = NULL;

    double accuracy1;           /* error (< accuracy) used to calculate alphas */

    static int kmax;            /* used for finding kopt */
    static int kopt;            /* optimal row number for convergence */

    /* sequence of separate attempts to cross interval htot with increasing    *
     * values of nsteps as suggested by Deuflhard (Num Rec, p. 726)            */

    static int nseq[] = { 0, 2, 6, 10, 14, 22, 34, 50, 70 };

    /* miscellaneous flags */

    static int first = 1;       /* is this the first try for a given step? */
    int reduct;                 /* flag indicating if we have reduced stepsize yet */
    int exitflag = 0;           /* exitflag: when set, we exit (!) */

    /* static global arrays */

    extern double *d;           /* D's used for extrapolation in pzextr */
    extern double *hpoints;     /* stepsizes h (=H/n) which we have tried */

    /* allocate arrays */

    if( !( d = ( double * ) calloc( KMAXX * n, sizeof( double ) ) ) )
        error( "BaDe: error allocating d.\n" );

    if( !( dfdt = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BaDe: error allocating dfdt.\n" );

    if( !( jac = ( double ** ) calloc( n, sizeof( double * ) ) ) )
        error( "BaDe: error allocating jac.\n" );

    for( i = 0; i < n; i++ )
        if( !( jac[i] = ( double * ) calloc( n, sizeof( double ) ) ) )
            error( "BaDe: error allocating jac's second dimension.\n" );

    if( !( err = ( double * ) calloc( KMAXX, sizeof( double ) ) ) )
        error( "BaDe: error allocating err.\n" );

    if( !( verror = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BaDe: error allocating verror.\n" );

    if( !( vsav = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BaDe: error allocating vsav.\n" );

    if( !( vseq = ( double * ) calloc( n, sizeof( double ) ) ) )
        error( "BaDe: error allocating vseq.\n" );

    if( !( hpoints = ( double * ) calloc( KMAXX, sizeof( double ) ) ) )
        error( "BaDe: error allocating hpoints.\n" );

    /* set initial stepsize to try to htry */

    h = htry;

    /* allocate memory if necessary */

    if( !a ) {
        a = ( double * ) calloc( IMAXX + 1, sizeof( double ) );
        alf = ( double ** ) calloc( KMAXX + 1, sizeof( double * ) );
        for( i = 0; i < KMAXX + 1; i++ )
            alf[i] = ( double * ) calloc( KMAXX + 1, sizeof( double ) );
    }

    /* new accuracy? -> initialize (here, this only applies at start) */

    if( ( accuracy != old_accuracy ) || ( n != nold ) ) {

        *hnext = tnew = -1.0e29;        /* init these to impossible value */

        /* initialize the work coefficients */

        a[1] = nseq[1] + 1;
        for( k = 1; k <= KMAXX; k++ )
            a[k + 1] = a[k] + nseq[k + 1];

        /* initialize the correction factors (alpha) */

        accuracy1 = SAFE1 * accuracy;   /* accuracy1 used to calculate alphas */
        for( i = 2; i <= KMAXX; i++ )
            for( k = 1; k < i; k++ )
                alf[k][i] = pow( accuracy1, ( ( a[k + 1] - a[i + 1] ) / ( ( a[i + 1] - a[1] + 1.0 ) * ( 2 * k + 1 ) ) ) );
        old_accuracy = accuracy;
        nold = n;

        /* add cost of Jacobian evaluations to work coefficients */

        a[1] += n;
        for( k = 1; k <= KMAXX; k++ )
            a[k + 1] = a[k] + nseq[k + 1];

        /* determine optimal row number for convergence */

        for( kopt = 2; kopt < KMAXX; kopt++ )
            if( a[kopt + 1] > a[kopt] * alf[kopt - 1][kopt] )
                break;
        kmax = kopt;

    }

    /* actual stepping starts here: first save starting values of v[] in vsav  */

    for( i = 0; i < n; i++ )
        vsav[i] = v[i];

    /* evaluate jacobian */

    p_jacobn( *t, v, dfdt, jac, n, si, inp );

    /* new integration or stepsize? -> re-establish the order window */

    if( *t != tnew || h != ( *hnext ) ) {
        first = 1;
        kopt = kmax;
    }

    reduct = 0;                 /* we have not reduced stepsize yet */

    for( ;; ) {

        /* try increasing numbers of steps over the interval h */

        for( k = 1; k <= kmax; k++ ) {

            tnew = ( *t ) + h;
            if( tnew == ( *t ) )
                error( "BaDe: stepsize underflow in stifbs.\n" );

            /* propagate the equations from t to t+h with nseq[k] steps using vsav and *
             * deriv as input; dfdt are the function derivs and jac is the Jacobian    *
             * trix; simpr returns vseq                                                */

            simpr( vsav, vseq, deriv, dfdt, jac, *t, h, nseq[k], n, si, inp );

            /* extrapolate v's for h->0 based on the different stepsizes (hest's) we   *
             * have tried already; the result is returned in v (errors in verror);     *
             * hest is squared since the error series is even                          */

            hest = DSQR( h / nseq[k] );
            pzextr( k - 1, hest, vseq, v, verror, n, KMAXX );

            if( k > 1 ) {

                /* find the maximum error */

                verror_max = 0.;
                for( i = 0; i < n; i++ )
                    if( v[i] != 0. )
                        verror_max = DMAX( fabs( verror[i] / v[i] ), verror_max );
                    else
                        verror_max = DMAX( fabs( verror[i] / DBL_EPSILON ), verror_max );

                /* scale error according to desired accuracy */

                verror_max /= accuracy;

                /* compute normalized error estimates (epsilons) */

                km = k - 1;
                err[km - 1] = pow( verror_max / SAFE1, 1.0 / ( 2 * km + 1 ) );

            }

            /* are we in the order window? -> converged */

            if( k > 1 && ( k >= kopt - 1 || first ) ) {

                if( verror_max < 1.0 ) {        /* exit if accuracy good enough */
                    exitflag = 1;
                    break;
                }

                if( k == kmax || k == kopt + 1 ) {      /* stepsize reduction possible? */
                    red = SAFE2 / err[km - 1];
                    break;
                } else if( k == kopt && alf[kopt - 1][kopt] < err[km - 1] ) {
                    red = 1.0 / err[km - 1];
                    break;
                } else if( kopt == kmax && alf[km][kmax - 1] < err[km - 1] ) {
                    red = alf[km][kmax - 1] * SAFE2 / err[km - 1];
                    break;
                } else if( alf[km][kopt] < err[km - 1] ) {
                    red = alf[km][kopt - 1] / err[km - 1];
                    break;
                }
            }
        }

        if( exitflag )
            break;

        /* reduce stepsize by at least REDMIN and at most REDMAX, then try again */

        red = DMIN( red, REDMIN );
        red = DMAX( red, REDMAX );

        h *= red;
        reduct = 1;
    }

    /* we've taken a successful step */

    *t = tnew;
    *hdid = h;
    first = 0;
    wrkmin = 1.0e35;

    /* compute optimal row for convergence and corresponding stepsize */

    for( i = 1; i <= km; i++ ) {

        fact = DMAX( err[i - 1], SCALMX );
        work = fact * a[i + 1];

        if( work < wrkmin ) {

            scale = fact;
            wrkmin = work;
            kopt = i + 1;

        }
    }

    *hnext = h / scale;

    /* check for possible order increase, but not if stepsize was just reduced */

    if( kopt >= k && kopt != kmax && !reduct ) {

        fact = DMAX( scale / alf[kopt - 1][kopt], SCALMX );

        if( a[kopt + 1] * fact <= wrkmin ) {
            *hnext = h / fact;
            kopt++;
        }
    }

    /* clean up */

    for( i = 0; i < n; i++ )
        free( jac[i] );
    free( jac );

    free( d );
    free( dfdt );
    free( hpoints );
    free( vseq );
    free( verror );
    free( err );
    free( vsav );
}

/**  simpr: implements the semi-implicit midpoint rule used by BaDe; it 
 *          subdivides a large step (htot) into nstep intervals and uses a 
 *          semi-implicit midpoint rule over the whole of htot, with a     
 *          stepsize of 2*h except for the first and the last step; the    
 *          nice thing about this method is that the error follows a power 
 *          law depending on the stepsize, i.e. its error converges to 0   
 *          really fast as h is diminished; this allows for good extrapo-  
 *          lation to h=0 (see pzextr() above)                             
 */
void
simpr( double *vin, double *vout, double *deriv, double *dfdt, double **jac, double tin, double htot, int nstep, int n, SolverInput * si, Input * inp ) {
    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *          lubcmp: does DU decomposition of a matrix                      *
     *          lubksb: solves linear system a * b                             */

    void ludcmp( double **a, int n, int *indx, double *d );
    void lubksb( double **a, int n, int *indx, double *b );

    /*** variables *************************************************************/

    int i, j;                   /* loop counters */

    double t;                   /* current time */
    double h;                   /* small (h) stepsize (equals htot / nstep) */
    double d;                   /* d indicates if eve/odd rows were switched in ludcmp */

    int *indx;                  /* array for permutated row indices of matrix a */

    double *del;                /* delta: difference in v between steps */
    double *vtemp;              /* array for temp derivs and v's */

    double **a;                 /* matrix [1 - hf'] */

    /* allocate arrays */

    indx = ( int * ) calloc( n, sizeof( int ) );

    del = ( double * ) calloc( n, sizeof( double ) );
    vtemp = ( double * ) calloc( n, sizeof( double ) );

    a = ( double ** ) calloc( n, sizeof( double * ) );
    for( i = 0; i < n; i++ )
        a[i] = ( double * ) calloc( n, sizeof( double ) );

    /* calculate h from H and n */

    h = htot / nstep;

    /* set up the matrix [1 - hf'] */

    for( i = 0; i < n; i++ ) {
        for( j = 0; j < n; j++ )
            a[i][j] = -h * jac[i][j];
        ++a[i][i];
    }

    /* the following is slightly bizarre */

    /* do LU decomposition of matrix [1 - hf']; this is needed for all steps   *
     * below, which use linearization of the equations to get estimates of the *
     * derivatives at the endpoint of a step (which is done in lubksb)         */

    ludcmp( a, n, indx, &d );

    /* do the first step */

    /* set up the right-hand side for the first step; use vout for temp sto-   *
     * rage; note that since our equations are autonomous, all dfdt's are 0.0; *
     * lubksb solves the linear system given by a and vout and then returns    *
     * the solution vector in vout; this vector contains the difference be-    *
     * tween this step's v's and the next steps' v's which are stored in del   */

    for( i = 0; i < n; i++ )
        vout[i] = h * ( deriv[i] + h * dfdt[i] );
    lubksb( a, n, indx, vout );

    for( i = 0; i < n; i++ )
        vtemp[i] = vin[i] + ( del[i] = vout[i] );

    /* do the intermediate steps */

    t = tin + h;
    p_deriv( vtemp, t, vout, n, si, inp );

    for( i = 2; i <= nstep; i++ ) {

        /* set up the right hand side */

        for( j = 0; j < n; j++ )
            vout[j] = h * vout[j] - del[j];
        lubksb( a, n, indx, vout );

        /* take the step */

        for( j = 0; j < n; j++ )
            vtemp[j] += ( del[j] += 2.0 * vout[j] );

        /* go to next step */

        t += h;
        p_deriv( vtemp, t, vout, n, si, inp );
    }

    /* do the last step */

    /* set up the right-hand side for the last step */

    for( i = 0; i < n; i++ )
        vout[i] = h * vout[i] - del[i];
    lubksb( a, n, indx, vout );

    /* take the last step */

    for( i = 0; i < n; i++ )
        vout[i] += vtemp[i];

    /* clean up */

    for( i = 0; i < n; i++ )
        free( a[i] );
    free( a );

    free( vtemp );
    free( del );
    free( indx );

}

/**  ludcmp: does an LU decomposition of matrix a, which is needed for in- 
 *           verting the matrix; LU decomposition dissects the matrix into 
 *           a lower and an upper triangular matrix that when multiplied,  
 *           equal matrix a again; this function uses Crout's algorithm    
 *           for the decomposition; the result is returned in a and con-   
 *           tains the LU decomposition of a rowwise permutation of a;     
 *           indx returns the order of the permutated rows; n is the size  
 *           of matrix a (as in n x n); d returns +1 if number of row in-  
 *           terchanges was even, -1 if odd.                               
 */
void
ludcmp( double **a, int n, int *indx, double *d ) {
    int i, j, k;                /* loop counters */
    int imax;                   /* index of row of largest element in a column */

    double big;                 /* largest value, used for pivoting */
    double dum;                 /* dummy used for pivoting */
    double sum;                 /* temp var for summing */
    double temp;                /* temp var for absolute values */
    double *vv;                 /* stores the implicit scaling information */

    const double TINY = 1.0e-20;        /* this is a tiny number */

    imax = 0;
    
    vv = ( double * ) calloc( n, sizeof( double ) );

    *d = 1.0;                   /* no interchanges yet */

    /* loop over rows to get the implicit scaling information */

    for( i = 0; i < n; i++ ) {

        big = 0.0;
        for( j = 0; j < n; j++ )
            if( ( temp = fabs( a[i][j] ) ) > big )
                big = temp;
        /* whole row == 0 -> singular matrix */
        if( big == 0.0 )
            error( "ludcmp: cannot decompose singular matrix" );

        vv[i] = 1.0 / big;      /* save the scaling */
    }

    /* loop over columns of Crout's method; this is done in two loops which    *
     * reflect the two different formulas of Crout's method, except that the   *
     * first formula is used for the diagonal elements in the second loop; in  *
     * the second loop we're also finding the largest element and its row num- *
     * ber for pivoting                                                        */

    for( j = 0; j < n; j++ ) {

        for( i = 0; i < j; i++ ) {
            sum = a[i][j];
            for( k = 0; k < i; k++ )
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }

        big = 0.0;              /* used for search of largest pivot element */

        for( i = j; i < n; i++ ) {

            sum = a[i][j];      /* summing happens here */
            for( k = 0; k < j; k++ )
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;

            if( ( dum = vv[i] * fabs( sum ) ) >= big ) {        /* pivoting stuff here */
                big = dum;
                imax = i;
            }

        }

        /* pivoting: do we need to interchange rows? if yes -> do so! */

        if( j != imax ) {
            for( k = 0; k < n; k++ ) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }

            *d = -( *d );       /* change parity of d */
            vv[imax] = vv[j];   /* rearrange scale factors */
        }

        indx[j] = imax;         /* store index of permutated row */

        /* if the pivot element is zero, the matrix is singular; TINY avoids div/0 *
         * below (some applications can deal with singular matrices, we don't)     */

        if( a[j][j] == 0.0 )
            a[j][j] = TINY;

        /* divide by the pivot element */

        if( j != n - 1 ) {
            dum = 1.0 / ( a[j][j] );
            for( i = j + 1; i < n; i++ )
                a[i][j] *= dum;
        }
        /* end of loop: next column */
    }

    free( vv );
}

/**  lubksb: does forward and backsubstitution of matrix a; in fact, a is 
 *           not passed in its original form but as the LU decomposition   
 *           of a rowwise permutation of a as returned by the function     
 *           ludcmp(); the right hand side vector is called b, which also  
 *           returns the solution vector to the calling function; indx is  
 *           a vector that indicates the order of rows in the permutated   
 *           matrix; n is the dimension of the matrix a (as in n x n).     
 */
void
lubksb( double **a, int n, int *indx, double *b ) {
    int i, j;                   /* counter variables */
    int ii = -1;                /* index of first non-zero element of b */
    int ip;                     /* index of permutated matrix a */

    double sum;                 /* temp var for summing */

    /* first loop does the forward substitution; we don't worry about the dia- *
     * gonal elements of the a matrix since here they're all 1 by definition   *
     * (see description of Crout's algorithm in Num Rec, chapter 2.3)          */

    for( i = 0; i < n; i++ ) {

        ip = indx[i];           /* get the index of the (real) first row */
        sum = b[ip];            /* get the according b */
        b[ip] = b[i];           /* un-permutate the b-vector */
        /* start summing at first non-zero element of b */
        if( ii >= 0 )
            for( j = ii; j <= i - 1; j++ )
                sum -= a[i][j] * b[j];
        else if( sum )
            ii = i;

        b[i] = sum;             /* b's are now in the right (un-permutated) order */

    }

    /* second loop does the backsubstitution */

    for( i = n - 1; i >= 0; i-- ) {
        sum = b[i];
        for( j = i + 1; j <= n - 1; j++ )
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

/*** DEBUGGING AND SOLVER LOG FUNCS ***************************************/

void
WriteSolvLog( char *solver, double tin, double tout, double h, int n, int nderivs, FILE * slog ) {
    double nds;                 /* Number of Derivative evaluations per Step */
    double ttot;                /* total time of propagation interval */

    nds = ( double ) nderivs / ( double ) n;
    ttot = tout - tin;

    fprintf( slog, "%s:   tin = %7.3f   tout = %7.3f   ttot = %7.3f   ", solver, tin, tout, ttot );
    fprintf( slog, "h = %5.3f   nsteps = %4d    nderivs/step = %4.1f\n", h, n, nds );
}

/**
 * Delay solver - by Manu?
 */
void
SoDe( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp ) {
    int i, j;
    EqParms *lp;
    double *shift_discs, *fact_discons;
    double *discon_array;
    int num_discons, fact_discons_size;

    double **vees, *tees;
    //printf("running delay solver\n");
    /*lp = GetMutParameters(); */
    lp = &( inp->lparm );
    fact_discons = GetFactDiscons( &fact_discons_size, si->all_fact_discons );

    //printf("%d %f %f %f\n",fact_discons_size,fact_discons[0],fact_discons[1], fact_discons[2]);

    maxdel = 0.;
    mindel = 1000.;

    for( j = 0; j < inp->zyg.defs.ngenes; j++ ) {

        if( lp->tau[j] > maxdel )
            maxdel = lp->tau[j];
        if( lp->tau[j] < mindel )
            mindel = lp->tau[j];

    }
    numdel = inp->zyg.defs.ngenes;
    delay = lp->tau;

    //printf("maxdel:%f, mindel:%f and numdel:%d\n",maxdel,mindel,numdel);

    j = -1;
    do {
        j++;

    } while( ( j < fact_discons_size ) && ( tout > fact_discons[j] ) );

    shift_discs = ( double * ) calloc( j + 1, sizeof( double ) );
    shift_discs[0] = 0.;

    for( i = 0; i < j; i++ ) {
        shift_discs[i + 1] = fact_discons[i] - tin;
    }
    discon_array = Construct_Discont_Array( tout - tin, lp->tau, inp->zyg.defs.ngenes, shift_discs, j + 1, &num_discons );
    for( i = 0; i < num_discons; i++ ) {
        discon_array[i] += tin;
    }



    vees = ( double ** ) calloc( 2, sizeof( double * ) );
    tees = ( double * ) calloc( 2, sizeof( double ) );

    tees[0] = tin;
    tees[1] = tout;
    vees[0] = vin;
    vees[1] = vout;

    DCERk32( vees, n, tees, 2, discon_array, num_discons, stephint, accuracy, si, inp );

    free( shift_discs );
    if( discon_array )
        free( discon_array );
    free( vees );
    free( tees );

    return;

}

int
compare( double *x, double *y ) {

    if( *x > *y )
        return 1;
    else if( *x < *y )
        return -1;
    else
        return 0;

}

double *
Construct_Discont_Array( double range, double *taus, int n, double *starts, int sn, int *disc_size ) {

    int M = 0;
    long m;
    int i, j, k;
    double *delay_array = NULL;
    for( i = 0; i < n; i++ ) {
        if( taus[i] != 0. ) {
            m = ( ( long ) floor( range / taus[i] ) ) + 1;
        } else
            m = 1;

        if( m > 5 )
            m = 5;

        //printf("%f %f m=%d\n",range/taus[i],floor(range/taus[i]),m);

        for( k = 0; k < sn; k++ ) {
            if( starts[k] > -taus[i] ) {

                delay_array = ( double * ) realloc( delay_array, ( M + m ) * sizeof( double ) );
                for( j = M; j < M + m; ++j ) {
                    delay_array[j] = starts[k] + ( j - M + 1 ) * taus[i];
                }
                M += m;
            }
        }
    }
    /*  for (i=0;i<M;i++) printf("%1.14f\n",delay_array[i]); */

    qsort( ( void * ) delay_array, M, sizeof( double ), ( int ( * )( const void *, const void * ) ) compare );

    /*  for (i=0;i<M;i++) printf("%1.14f\n",delay_array[i]); */
    qsort( ( void * ) taus, n, sizeof( double ), ( int ( * )( const void *, const void * ) ) compare );
    /* Now lets remove duplicates */

    i = 0;

    while( i < M ) {
        if( i < M - 1 ) {
            if( ( fabs( delay_array[i] - delay_array[i + 1] ) < 1E-10 )
                || ( delay_array[i] > range ) ) {
                memmove( ( delay_array + i ), ( delay_array + i + 1 ), ( M - i - 1 ) * sizeof( double ) );
                /*                              printf("Shifted %d elements to %d\n",M-i-1,i); */
                M--;
                i--;
            }
        } else if( ( i == M - 1 ) && ( delay_array[i] > range ) ) {
            M--;
            i--;
        }

        i++;
        delay_array = ( double * ) realloc( delay_array, M * sizeof( double ) );
    }
    /*  for (i=0;i<M;i++) printf("Filtered:%1.14f\n",delay_array[i]); */

    *disc_size = M;
    return delay_array;

}

void
CE( double t, double *vans, double tbegin, double *v_at_tbegin, double ech, double *d1, double *d2, double *d3, double *d4, int n ) {
    static double e1 = -4.0 / 3.0, e2 = 5.0 / 9.0, e3 = -2.0 / 3.0, e4 = -8.0 / 9.0;

    int i;
    double sigma;               /* for the continuous ext. v(n+sigma) */
    double sigma_sq;            /* sigma's square */
    double dt, f1, f2, f3, f4;

    sigma = ( t - tbegin ) / ech;
    sigma_sq = sigma * sigma;
    dt = t - tbegin;

    f1 = dt * ( 1 + e1 * sigma + e2 * sigma_sq );
    f2 = dt * ( sigma + e3 * sigma_sq );
    f3 = dt * ( e1 * sigma - e4 * sigma_sq );
    f4 = dt * ( sigma_sq - sigma );

    for( i = 0; i < n; i++ )
        vans[i] = v_at_tbegin[i] + f1 * d1[i];

    for( i = 0; i < n; i++ )
        vans[i] = vans[i] + f2 * d2[i];

    for( i = 0; i < n; i++ )
        vans[i] = vans[i] - f3 * d3[i];

    for( i = 0; i < n; i++ )
        vans[i] = vans[i] + f4 * d4[i];

    return;

}

int
y_delayed( double ***vd, int n, double *rktimes, double *tau,
           double *grid, double **vdone, double **deriv1, double **deriv2, double
           **deriv3, double **deriv4, int gridsize, double accu, SolverInput * si, Input * inp ) {

    int i, j, vc, dc, it;
    double t, ech;
    double *vtemp, *vnext, *vprev, *dummy;
    double *drv1, *drv2, *drv3, *drv4;
    double verror_max;

    //  static double a1 = 0.0, a2 = 0.5, a3 = 0.75, a4 = 1.0;

    static double b21 = 0.5,
        /*  b31 =     0.0, */
        b32 = 0.75;

    static double c1 = 2.0 / 9.0, c2 = 1.0 / 3.0, c3 = 4.0 / 9.0;

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    vnext = ( double * ) calloc( n, sizeof( double ) );
    vprev = ( double * ) calloc( n, sizeof( double ) );
    drv1 = ( double * ) calloc( n, sizeof( double ) );
    drv2 = ( double * ) calloc( n, sizeof( double ) );
    drv3 = ( double * ) calloc( n, sizeof( double ) );
    drv4 = ( double * ) calloc( n, sizeof( double ) );

    ech = rktimes[3] - rktimes[0];

    /* First lets initialiaze all the the vdelays based on what is
       available to us without iteration. We will carry out iteration
       subsequently, only if t-mindel lies within the integration interval */

    for( vc = 0; vc < 4; vc++ ) {

        t = rktimes[vc];

        for( dc = 0; dc < numdel; dc++ )
            if( tau[dc] == 0. )
                vd[vc][dc] = memcpy( vd[vc][dc], vdone[gridsize - 1], sizeof( double ) * n );
            else if( t - tau[dc] <= grid[0] )
                History( t - tau[dc], t, vd[vc][dc], n, inp->his[si->genindex], inp->zyg.defs.ngenes, &( inp->zyg ) );
            else if( t - tau[dc] <= grid[gridsize - 1] ) {

                j = 0;
                while( ( j < gridsize ) && ( t - tau[dc] > grid[j] ) )
                    j++;

                if( j >= gridsize ) {
                    printf( "y_past:time requested not in grid, bailing!\n" );
                    exit( 1 );
                }

                if( ( j == gridsize - 1 ) && ( t - tau[dc] == grid[gridsize - 1] ) ) {
                    vd[vc][dc] = memcpy( vd[vc][dc], vdone[gridsize - 1], sizeof( double ) * n );
                    /*                          for (i=0; i<n; i++)
                       printf("t=%f, t-del=%f, vdone[%d]=%f,"
                       "vpast[%d]=%f\n", 
                       t,t-tau[dc],j,vdone[j][i],i,vd[vc][dc][i]); */
                } else {

                    CE( t - tau[dc], vd[vc][dc], grid[j - 1], vdone[j - 1], grid[j] - grid[j - 1], deriv1[j - 1], deriv2[j - 1], deriv3[j - 1], deriv4[j - 1],
                        n );
                    /*                                  for (i=0; i<n; i++)
                       printf("t=%f, t-del=%f, vdone[%d]=%f," 
                       "vpast[%d]=%f\n",
                       t,t-tau[dc],j,vdone[j-1][i],i,vd[vc][dc][i]); */
                }
            } else if( t - tau[dc] > grid[gridsize - 1] ) {

                CE( t - tau[dc], vd[vc][dc], grid[gridsize - 2], vdone[gridsize - 2], grid[gridsize - 1] - grid[gridsize - 2], deriv1[gridsize - 2],
                    deriv2[gridsize - 2], deriv3[gridsize - 2], deriv4[gridsize - 2], n );
                /*                              for (i=0; i<n; i++)
                   printf("IterInit [%.6f,%.6f], t-del=%.6f,"
                   "vd[%d][%d]=%.6f\n", 
                   grid[gridsize-2],grid[gridsize-1],
                   t-tau[dc],vc,dc,vd[vc][dc][i]); */
            } else {

                printf( "Requested point not to be found\n" );
                exit( 1 );

            }
    }

    /* Now lets do the promised iteration, if required ofcourse! */
    if( rktimes[3] - mindel > grid[gridsize - 1] ) {
        d_deriv( vdone[gridsize - 1], vd[0], rktimes[0], drv1, n, si, inp );

        it = 0;
        verror_max = 100.;

        while( ( it < 5 ) && ( verror_max > 0.01 * accu ) ) {

            /*                  printf("Iteration No.%d, error: %f\n",it,verror_max); */

            for( i = 0; i < n; i++ )
                vtemp[i] = vdone[gridsize - 1][i] + ech * ( b21 * drv1[i] );
            d_deriv( vtemp, vd[1], rktimes[1], drv2, n, si, inp );

            for( i = 0; i < n; i++ )
                vtemp[i] = vdone[gridsize - 1][i] + ech * ( b32 * drv2[i] );
            d_deriv( vtemp, vd[2], rktimes[2], drv3, n, si, inp );

            for( i = 0; i < n; i++ )
                vnext[i] = vdone[gridsize - 1][i] + ech * ( c1 * drv1[i] + c2 * drv2[i] + c3 * drv3[i] );
            d_deriv( vnext, vd[3], rktimes[3], drv4, n, si, inp );

            for( vc = 0; vc < 4; vc++ ) {

                t = rktimes[vc];

                for( dc = 0; dc < numdel; dc++ )
                    if( t - tau[dc] > grid[gridsize - 1] ) {
                        CE( t - tau[dc], vd[vc][dc], grid[gridsize - 1], vdone[gridsize - 1], ech, drv1, drv2, drv3, drv4, n );

                        /*                                      for (i=0; i<n; i++)
                           printf("[%.6f,%.6f], t-del=%.6f,"
                           "vd1[%d]=%.6f,vd2[%d]=%.6f, vd3[%d]=%.6f," 
                           "vd4[%d]=%.6f, vnext[%d]=%.6f, t=%.6f\n",
                           rktimes[0],rktimes[3],
                           t-tau[dc],dc,vd[0][dc][0], 
                           dc,vd[1][dc][0], dc,vd[2][dc][0],  
                           dc,vd[3][dc][0],i,vnext[i],
                           grid[gridsize-1]); */
                    }
            }

            if( it ) {
                verror_max = 0.;

                for( i = 0; i < n; i++ ) {
                    if( fabs( vprev[i] ) > 10000. * DBL_EPSILON )
                        verror_max = DMAX( fabs( ( vnext[i] - vprev[i] ) / vprev[i] ), verror_max );
                    else
                        verror_max = DMAX( fabs( ( vnext[i] - vprev[i] ) / ( 10000. * DBL_EPSILON ) ), verror_max );
                }
            }

            dummy = vnext;
            vnext = vprev;
            vprev = dummy;

            it++;
        }

        if( it == 5 )
            return 1;
    }

    free( vtemp );
    free( vnext );
    free( vprev );
    free( drv1 );
    free( drv2 );
    free( drv3 );
    free( drv4 );

    return 0;
}

/**  DCERk3(2): propagates v[0] (of size n) according to tarray  by  
the Runge-Kutta 3(2) pair with continuous extension storing the     
result in vatt. Initial conditions are specified in vatt[0],        
corresponding to tarray[0]											
*/
void
DCERk32( double **vatt, int n, double *tarray, int tpoints, double *darray, int dpoints, double stephint, double accuracy, SolverInput * si, Input * inp ) {
    int i, dc, vc;              /* local loop counters */
    int tpos = 1;               /* where in the t array we are, starting
                                   at the value right after the starting time */

    int dpos = 0;               /* where is the discontinuities array we
                                   are (darray) */
    int DISCON = 0;             /* if we just hit a discontinuity */

    double **v; /** used for storing intermediate steps */
    int toggle = 0;             /* used to toggle between v[0] and v[1] */

    double *vtemp;              /* guessed intermediate v's */

    double *vnow;               /* ptr to v at current time */
    double *vnext;              /* ptr to v at current time + stepsize */
    double **v_delayed[4];      /* ptr to v at t-delay */

    double *verror;             /* error estimate */
    double verror_max;          /* the maximum error in verror */

    //  double *tempswap;

    double t;                   /* the current time */
    double tms[4];              /* array for the rk step times */
    double h = stephint;        /* initial stepsize */
    double hnext;               /* used to calculate next stepsize */
    const double SAFETY = 0.9;  /* safety margin for decrsg stepsize */
    const double ERRCON = 5.832e-3;     /* to prevent huge jumps, see
                                           num recipes pg. 719 */



    /* declare and initialize parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double step_exponent = -1.0 / 3.0;

    static double
        /*  a1  =      0.0, */
      a2 = 0.5, a3 = 0.75;
    /*  a4  =     1.0, */

    static double b21 = 0.5,
        /*  b31 =     0.0, */
        b32 = 0.75;

    static double c1 = 2.0 / 9.0, c2 = 1.0 / 3.0, c3 = 4.0 / 9.0;


    double dc1 = c1 - 7.0 / 24.0, dc2 = c2 - 0.25, dc3 = c3 - c2, dc4 = -0.125;


    /* the do-nothing case; too small steps dealt with under usual */

    if( tarray[0] >= tarray[tpoints - 1] )
        return;

    /* the usual case: steps big enough */

    v = ( double ** ) calloc( 2, sizeof( double * ) );

    vtemp = ( double * ) calloc( n, sizeof( double ) );
    verror = ( double * ) calloc( n, sizeof( double ) );
    v[0] = ( double * ) calloc( n, sizeof( double ) );
    v[1] = ( double * ) calloc( n, sizeof( double ) );

    for( vc = 0; vc < 4; vc++ ) {
        v_delayed[vc] = ( double ** ) calloc( numdel, sizeof( double * ) );
        for( dc = 0; dc < numdel; dc++ )
            v_delayed[vc][dc] = ( double * ) calloc( n, sizeof( double ) );
    }

    t = tarray[0];
    vnow = vatt[0];
    vnext = v[0];

    /* allocate more grid and derivs */
    gridpos++;
    if( !( tdone = ( double * ) realloc( tdone, ( gridpos + 1 ) * sizeof( double ) ) ) ) {
        printf( "Could not allocate memory for tdone\n" );
        exit( 1 );
    }

    if( !( vdonne = ( double ** ) realloc( vdonne, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
        printf( "Could not allocate memory for vdonne\n" );
        exit( 1 );
    }

    if( !( derivv1 = ( double ** ) realloc( derivv1, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
        printf( "Could not allocate memory for derivv1\n" );
        exit( 1 );
    }

    if( !( derivv2 = ( double ** ) realloc( derivv2, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
        printf( "Could not allocate memory for derivv2\n" );
        exit( 1 );
    }

    if( !( derivv3 = ( double ** ) realloc( derivv3, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
        printf( "Could not allocate memory for derivv3 \n" );
        exit( 1 );
    }

    if( !( derivv4 = ( double ** ) realloc( derivv4, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
        printf( "Could not allocate memory for derivv4\n" );
        exit( 1 );
    }
    derivv1[gridpos] = ( double * ) calloc( n, sizeof( double ) );
    derivv2[gridpos] = ( double * ) calloc( n, sizeof( double ) );
    derivv3[gridpos] = ( double * ) calloc( n, sizeof( double ) );
    derivv4[gridpos] = ( double * ) calloc( n, sizeof( double ) );
    vdonne[gridpos] = ( double * ) calloc( n, sizeof( double ) );
    tdone[gridpos] = t;
    memcpy( vdonne[gridpos], vnow, sizeof( *vnow ) * n );


    /* initial stepsize cannot be bigger than total time */
    if( tarray[0] + h >= tarray[tpoints - 1] )
        h = tarray[tpoints - 1] - tarray[0];

    /*  printf("%f\n",step_exponent); */

    /*  printf("disconny:%d %d\n",dpoints,dpos); */

    if( ( dpoints ) && ( dpos < dpoints ) )
        if( ( darray[dpos] > t ) && ( darray[dpos] <= t + h ) ) {

            /*                  printf("\n***Hit discont at %f ***\n",darray[dpos]); */

            h = darray[dpos] - t;
            dpos++;
            /*                  printf("dpos:%d\n",dpos); */
            DISCON = 1;
        }

    if( ( h > mindel ) && ( h <= 2. * mindel ) ) {

        if( DISCON ) {

            dpos--;
            /*                  printf("dpos:%d\n",dpos); */
            DISCON = 0;

        }
        h = .5 * h;
    }

    tms[0] = t;
    tms[1] = t + a2 * h;
    tms[2] = t + a3 * h;
    tms[3] = t + h;
    /* we need to calculate derivv1 only the first time, since if the
       previous step was a success, we can use the last derivv4, and if it is
       was a failure, we don't have to recalculate it */
    while( y_delayed( v_delayed, n, tms, delay, tdone, vdonne, derivv1, derivv2, derivv3, derivv4, gridpos + 1, accuracy, si, inp ) ) {
        printf( "Rejected Iteration for [%f,%f]!\n", t, t + h );
        h = 0.5 * h;
        tms[0] = t;
        tms[1] = t + a2 * h;
        tms[2] = t + a3 * h;
        tms[3] = t + h;

        if( DISCON ) {

            dpos--;
            /*                          printf("dpos:%d\n",dpos); */
            DISCON = 0;

        }
    }

    d_deriv( vnow, v_delayed[0], t, derivv1[gridpos], n, si, inp );

    while( t < tarray[tpoints - 1] ) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */

        while( 1 ) {

            tms[0] = t;
            tms[1] = t + a2 * h;
            tms[2] = t + a3 * h;
            tms[3] = t + h;

            if( !y_delayed( v_delayed, n, tms, delay, tdone, vdonne, derivv1, derivv2, derivv3, derivv4, gridpos + 1, accuracy, si, inp ) ) {

                /* do the Runge-Kutta thing here: calulate intermediate 
                   derivatives */
                for( i = 0; i < n; i++ )
                    vtemp[i] = vnow[i] + h * ( b21 * derivv1[gridpos][i] );

                d_deriv( vtemp, v_delayed[1], t + a2 * h, derivv2[gridpos], n, si, inp );

                for( i = 0; i < n; i++ )
                    vtemp[i] = vnow[i] + h * ( b32 * derivv2[gridpos][i] );

                d_deriv( vtemp, v_delayed[2], t + a3 * h, derivv3[gridpos], n, si, inp );

                /* ... then feed them to the Rk32 formula */

                for( i = 0; i < n; i++ )
                    vnext[i] = vnow[i] + h * ( c1 * derivv1[gridpos][i] + c2 * derivv2[gridpos][i] + c3 * derivv3[gridpos][i] );

                /* Now calculate k4 for the embedded 4-stage formula, if this step is
                   succesful, it will get used as derivv1 (k1) in the next step */
                //              for (i=0; i<10; i++) {
                //                  printf("%d gridpos=%d; derivv1=%lg; derivv2=%lg; derivv3=%lg; h=%lg\n", i, gridpos, derivv1[gridpos][i], derivv2[gridpos][i], derivv3[gridpos][i], h);
                //              }
                d_deriv( vnext, v_delayed[3], t + h, derivv4[gridpos], n, si, inp );
                /* calculate the error estimate using the embedded formula */

                for( i = 0; i < n; i++ ) {
                    verror[i] = h * ( dc1 * derivv1[gridpos][i] + dc2 * derivv2[gridpos][i] + dc3 * derivv3[gridpos][i] + dc4 * derivv4[gridpos][i] );
                }
                //printf("verror=%lg, h=%lg, dc1=%lg, deriv1=%lg, dc2=%lg, deriv2=%lg, dc3=%lg, deriv3=%lg, dc4=%lg, deriv4=%lg\n", verror[0], h, dc1, derivv1[gridpos][0], dc2, derivv2[gridpos][0], dc3, derivv3[gridpos][0], dc4, derivv4[gridpos][0]);
                /* find the maximum error */
                verror_max = 0.;
                for( i = 0; i < n; i++ ) {
                    if( fabs( vnext[i] ) > 10000. * DBL_EPSILON ) {
                        verror_max = DMAX( fabs( verror[i] / vnext[i] ), verror_max );
                        //printf("VERROR THEN i=%d verror=%lg vnext=%lg\n", i, verror[i], vnext[i]);
                    } else {
                        verror_max = DMAX( fabs( verror[i] / ( 10000. * DBL_EPSILON ) ), verror_max );
                        //printf("VERROR ELSE i=%d verror=%lg vnext=%lg\n", i, verror[i], vnext[i]);
                    }
                }

                /*      for (i=0; i<n; i++) 
                   printf("The deriv1[%d]=%.10f deriv2[%d]=%.10f deriv3[%d]=%.10f "
                   "deriv4[%d]=%.10f verror[%d]=%.10f verror_max=%.10f \n",i,
                   derivv1[gridpos][i],i, derivv2[gridpos][i],i,derivv3[gridpos][i],i,derivv4[gridpos][i],i,
                   verror[i],verror_max);
                 */
                /* scale error according to desired accuracy */

                verror_max /= accuracy;
                /* compare maximum error to the desired accuracy; if error < accuracy, we  *
                 * are done with this step; otherwise, the stepsize has to be reduced and  *
                 * the step repeated; for detailed comments on the approximation involving *
                 * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

                if( verror_max <= 1.0 ) {
                    break;
                }

                /*        printf("Step size %f rejected!\n",h); */

                /* kludge for going back to the original position in the discontinuity
                   array if the step is rejected */
                if( DISCON ) {

                    dpos--;
                    /*          printf("dpos:%d\n",dpos); */
                    DISCON = 0;

                }

                hnext = SAFETY * h * pow( verror_max, step_exponent );

                /* decrease stepsize by no more than a factor of 10; check for underflows */

                h = ( hnext > 0.1 * h ) ? hnext : 0.1 * h;

                if( ( h > mindel ) && ( h <= 2. * mindel ) ) {

                    if( DISCON ) {

                        dpos--;
                        /*                      printf("dpos:%d\n",dpos); */
                        DISCON = 0;

                    }
                    h = .5 * h;

                }

                /*        printf("Step size %f suggested!\n",h); */

                if( h < DBL_EPSILON )
                    error( "SoDe: stepsize underflow" );

            } else {
                /*              printf("Rejected Iteration for [%f,%f]!\n",t,t+h); */
                h = 0.5 * h;

                if( ( h > mindel ) && ( h <= 2. * mindel ) )
                    h = .5 * h;

                if( DISCON ) {

                    dpos--;
                    /*          printf("dpos:%d\n",dpos); */
                    DISCON = 0;

                }
            }

        }


        /* advance the current time by last stepsize */

        /*      printf("The stepsize at time %f is %f error is %f\n",
           t,h,verror_max*accuracy); */




        /* Now we are going to do the continuous extension business, we will
           look for t < t* <= t+h, and for all such t*, we will calculate the CE
           and add the entried to the the vatt array. If t* = t+h, we just
           exchange the pointers */
        while( ( tarray[tpos] < t + h ) && ( tpos < tpoints ) ) {

            CE( tarray[tpos], vatt[tpos], t, vnow, h, derivv1[gridpos], derivv2[gridpos], derivv3[gridpos], derivv4[gridpos], n );
            /*          printf("Vatt: %d %f %f %f %f %f\n",
               tpos,tarray[tpos],vatt[tpos][0],vnow[0],t,t+h); */

            tpos += 1;
        }

        if( ( tarray[tpos] == t + h ) && ( tpos < tpoints ) ) {

            memcpy( vatt[tpos], vnext, sizeof( double ) * n );
            tpos += 1;
        }

        /* toggle vnow and vnext between v[0] and v[1], vnow is just 
           vnow_temp right now, as we need the old vnow in the Continuous 
           Extension stuff */

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

        t += h;


        /*    if (t >= tarray[tpoints - 1])
           break;                                that was the last iteration */

        /* increase stepsize according to error (5th order) for next iteration */

        /*      printf("ERRCON:%f\n",ERRCON); */

        if( verror_max > ERRCON ) {

            h = SAFETY * h * pow( verror_max, step_exponent );

        } else {
            h = 5.0 * h;
        }

        /* make sure t does not overstep last time point */

        if( t + h >= tarray[tpoints - 1] ) {
            h = tarray[tpoints - 1] - t;
        }

        DISCON = 0;

        if( ( dpoints ) && ( dpos < dpoints ) )
            if( ( darray[dpos] > t ) && ( darray[dpos] <= t + h ) ) {

                /*                      printf("\n***Hit discont at %f ***\n",darray[dpos]); */

                h = darray[dpos] - t;
                dpos++;
                /*                      printf("dpos:%d\n",dpos); */
                DISCON = 1;
            }


        if( ( h > mindel ) && ( h <= 2. * mindel ) ) {

            if( DISCON ) {

                dpos--;
                /*                      printf("dpos:%d\n",dpos); */
                DISCON = 0;

            }
            h = .5 * h;
        }


        /* allocate more grid and derivs */

        gridpos++;
        if( !( tdone = ( double * ) realloc( tdone, ( gridpos + 1 ) * sizeof( double ) ) ) ) {
            printf( "Could not allocate memory for tdone\n" );
            exit( 1 );
        }

        if( !( vdonne = ( double ** ) realloc( vdonne, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
            printf( "Could not allocate memory for vdonne\n" );
            exit( 1 );
        }

        if( !( derivv1 = ( double ** ) realloc( derivv1, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
            printf( "Could not allocate memory for derivv1\n" );
            exit( 1 );
        }

        if( !( derivv2 = ( double ** ) realloc( derivv2, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
            printf( "Could not allocate memory for derivv2\n" );
            exit( 1 );
        }

        if( !( derivv3 = ( double ** ) realloc( derivv3, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
            printf( "Could not allocate memory for derivv3 \n" );
            exit( 1 );
        }

        if( !( derivv4 = ( double ** ) realloc( derivv4, ( gridpos + 1 ) * sizeof( double * ) ) ) ) {
            printf( "Could not allocate memory for derivv4\n" );
            exit( 1 );
        }
        derivv1[gridpos] = ( double * ) calloc( n, sizeof( double ) );
        derivv2[gridpos] = ( double * ) calloc( n, sizeof( double ) );
        derivv3[gridpos] = ( double * ) calloc( n, sizeof( double ) );
        derivv4[gridpos] = ( double * ) calloc( n, sizeof( double ) );
        vdonne[gridpos] = ( double * ) calloc( n, sizeof( double ) );
        tdone[gridpos] = t;

        if( t - maxdel > tdone[0] ) {
            while( t - maxdel >= tdone[gridstart] )
                gridstart++;

            gridstart--;
        }

        memcpy( vdonne[gridpos], vnow, sizeof( *vnow ) * n );


        /* put present derivv4 into future derivv1 */

        memcpy( derivv1[gridpos], derivv4[gridpos - 1], sizeof( **derivv4 ) * n );
    }

    /*       for (j=0; j<=gridpos;j++)
       {
       for (i=0; i<n; i++) 
       printf("At time %f, The deriv1[%d][%d]=%.10f deriv2[%d][%d]"
       "=%.10f deriv3[%d][%d]=%.10f deriv4[%d][%d]=%.10f\n",tdone[j],j,i,
       derivv1[j][i],j,i, derivv2[j][i],j,i,derivv3[j][i],j,i,
       derivv4[j][i]);
       }
     */

    /*  printf("vs on grid:\n"); */

    /*           for (j=0; j<=gridpos;j++)
       {
       printf("%f %f\n",tdone[j],vdonne[j][0]);

       }
     */


    free( v[0] );
    free( v[1] );

    for( vc = 0; vc < 4; vc++ ) {
        for( dc = 0; dc < numdel; dc++ )
            free( v_delayed[vc][dc] );
        free( v_delayed[vc] );
    }

    /* Put zeroes in derivv1[gridpos], all the rest are zero anyways */

    free( derivv1[gridpos] );
    derivv1[gridpos] = ( double * ) calloc( n, sizeof( double ) );


    free( v );
    free( vtemp );
    free( verror );

}

void
FreeDelaySolver( void ) {

    int j;

    for( j = 0; j <= gridpos; j++ ) {
        free( derivv1[j] );
        free( derivv2[j] );
        free( derivv3[j] );
        free( derivv4[j] );
        free( vdonne[j] );
    }

    free( derivv1 );
    free( derivv2 );
    free( derivv3 );
    free( derivv4 );
    free( vdonne );
    free( tdone );

}

void
InitDelaySolver( void ) {

    gridpos = -1;
    gridstart = 0;
    tdone = NULL;
    vdonne = NULL;
    derivv1 = NULL;
    derivv2 = NULL;
    derivv3 = NULL;
    derivv4 = NULL;

}

void
DivideHistory( double t1, double t2, Zygote * zyg ) {
    double *blug;
    int i, size;

    if( ( size =
          GetNNucs( &( zyg->defs ), zyg->nnucs, t2, &( zyg->times ) ) * zyg->defs.ngenes ) > GetNNucs( &( zyg->defs ), zyg->nnucs, t1,
                                                                                                       &( zyg->times ) ) * zyg->defs.ngenes )
        for( i = 0; i <= gridpos; i++ ) {
            blug = ( double * ) calloc( size, sizeof( double ) );
            Go_Forward( blug, vdonne[i], GetStartLinIndex( t2, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( t1, &( zyg->defs ), &( zyg->times ) ), zyg,
                        zyg->defs.ngenes );

            free( vdonne[i] );
            vdonne[i] = blug;

            blug = ( double * ) calloc( size, sizeof( double ) );
            Go_Forward( blug, derivv1[i], GetStartLinIndex( t2, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( t1, &( zyg->defs ), &( zyg->times ) ), zyg,
                        zyg->defs.ngenes );

            free( derivv1[i] );
            derivv1[i] = blug;

            blug = ( double * ) calloc( size, sizeof( double ) );
            Go_Forward( blug, derivv2[i], GetStartLinIndex( t2, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( t1, &( zyg->defs ), &( zyg->times ) ), zyg,
                        zyg->defs.ngenes );

            free( derivv2[i] );
            derivv2[i] = blug;

            blug = ( double * ) calloc( size, sizeof( double ) );
            Go_Forward( blug, derivv3[i], GetStartLinIndex( t2, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( t1, &( zyg->defs ), &( zyg->times ) ), zyg,
                        zyg->defs.ngenes );

            free( derivv3[i] );
            derivv3[i] = blug;

            blug = ( double * ) calloc( size, sizeof( double ) );
            Go_Forward( blug, derivv4[i], GetStartLinIndex( t2, &( zyg->defs ), &( zyg->times ) ), GetStartLinIndex( t1, &( zyg->defs ), &( zyg->times ) ), zyg,
                        zyg->defs.ngenes );

            free( derivv4[i] );
            derivv4[i] = blug;

        }

}

/**  Krylov: propagates vin (of size n) from tin to tout by BDF (Backward  *
 *           Differential Formulas and use of a Newton-Krylov method with  *
 *           preconditioning to avoid the costly computation of the        *
 *           jacobian.                                                     *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Anton Crombach, October 2010                 *
 *                                                                         *
 ***************************************************************************/
/*
 * I realised I cannot evaluate only part of the derivative (f.i. only the
 * production/decay part)... so let's see what happens if I evaluate the full
 * derivative for preconditioning. That actually works better.
 *
 * There is one problem hidden here. If tout coincides with the change from
 * INTERPHASE to MITOSIS, the krylov solver somehow sets all gene product
 * concentrations to zero. Luckily we normally don't want to get output at 
 * these times (f.i. t = 16.000)... so all should be fine. It is most likely 
 * a problem with the discontinuity caused by the sudden stop of gene product
 * regulation -- only decay and diffusion remain in MITOSIS.
 * Damjan: Is all this true for the Band solver too? If so, it should be added to the Krylov (Band) solver function description
 */
/*
 * 
 //the "old" Krylov solver not used anymore. Now we use the Krylov Band solver
 * 
void Krylov(double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE *slog, SolverInput *sinput, Input *input) {
    int flag, i, j;
    realtype t, tstop;
    double *divtimes, *divdurations;

    inp = input;
    si = sinput;
    // If nothing to do, return 
    if (fabs(tin - tout) < 1e-6) return;
    if (cvode_mem != 0) {
        FreeKrylovSolver();
    }
    InitKrylovVariables(vin, n);
    InitKrylovSolver(tin, stephint, accuracy, accuracy);
    // Krylov solver looks ahead and then gets confused by the change
    //   in number of equations, so we need to set a stop time beyond which it
    //   is not allowed to look 
    i = inp->zyg.defs.ndivs - 1;
    j = -1;
    divtimes = inp->zyg.times.div_times;
    divdurations = inp->zyg.times.div_duration;
    while (i != j)
        if (tout > divtimes[ i ])
            --i;
        else
            j = i;

    // there is still a division to happen
    if (i > -1) {
        // perhaps mitosis even before
        tstop = divtimes[ i ];
        if (tout < divtimes[ i ] - divdurations[ i ])
            tstop = divtimes[ i ] - divdurations[ i ];
    } else {
        // otherwise don't look beyond gastrulation
        tstop = inp->zyg.times.gast_time;
    }
    CVodeSetStopTime(cvode_mem, tstop);

    // old code, works if networks would always be well-behaved 
    flag = CVode(cvode_mem, tout, vars, &t, CV_NORMAL);
    if (CheckFlag(&flag, "CVode", 1)) return;

    // copy vars into vout 
    for (i = 0; i < n; ++i) {
        vout[ i ] = NV_Ith_S(vars, i);
    }
}

int my_f(realtype t, N_Vector y, N_Vector ydot, void *extra_data) { //to call the derivative         //needed by the old Krylov solver
    // wrapper function 
    int i, n = NV_LENGTH_S(y);
    ExtraData edata = (ExtraData) extra_data;

    p_deriv(NV_DATA_S(y), t, NV_DATA_S(ydot), n, si, inp);

    // caching interaction data 
    for (i = 0; i < n; ++i) {
        edata->fsave[ i ] = NV_Ith_S(ydot, i);
    }

    return 0;
}

int Precond(realtype tn, N_Vector c, N_Vector fc, booleantype jok,
        booleantype *jcurPtr, realtype gamma, void *extra_data, N_Vector vtemp1,
        N_Vector vtemp2, N_Vector vtemp3) {                                      //needed by the old Krylov solver

    const int NRNUC = GetNNucs(&(inp->zyg.defs), inp->zyg.nnucs, tn, &(inp->zyg.times));

    realtype ***P;
    int **pivot;
    int n, i, j, ier, flag, ii, jj;

    realtype fac, r, r0, save, srur;
    realtype *f1, *fsave, *cdata, *err_data;

    ExtraData edata = (ExtraData) extra_data;
    void *cvode_mem = edata->cvode_mem;
    N_Vector err_weights = edata->err_weights;

    flag = CVodeGetErrWeights(cvode_mem, err_weights);
    if (CheckFlag(&flag, "CVodeGetErrWeights", 1)) return 1;
    err_data = NV_DATA_S(err_weights);

    cdata = NV_DATA_S(c);

    P = edata->P;
    pivot = edata->pivot;
    srur = SQRT(UNIT_ROUNDOFF);
    fsave = edata->fsave;

    //Approximate each diagonal block of Jacobian.
    //Here, fsave contains the base value of the rate vector and
    //r0 is a minimum increment factor for the difference quotient.

    f1 = NV_DATA_S(vtemp1);
    fac = N_VWrmsNorm(fc, err_weights);
    r0 = RCONST(1000.0) * ABS(gamma) * UNIT_ROUNDOFF * neq * fac;
    if (fabs(r0 - RCONST(0.0)) < 1e-6) r0 = RCONST(1.0);

    // get all interactions
    p_deriv(cdata, tn, f1, neq, si, inp);

    for (n = 0; n < NRNUC; ++n) {
        // Generate i-th diagonal block
        for (i = 0; i < inp->zyg.defs.ngenes; i++) {
            // Generate the j-th column as a difference quotient
            ii = n * inp->zyg.defs.ngenes + i;
            save = cdata[ ii ];

            r = MAX(srur * ABS(save), r0 / err_data[ ii ]);
            cdata[ ii ] += r;
            fac = -gamma / r;
            for (j = 0; j < inp->zyg.defs.ngenes; ++j) {
                jj = n * inp->zyg.defs.ngenes + j;
                P[ n ][ i ][ j ] = fac * (f1[ jj ] - fsave[ jj ]);
            }
            cdata[ ii ] = save;
        }
    }

    // Add identity matrix and do LU decompositions on blocks.
    for (i = 0; i < NRNUC; ++i) {
        denseAddIdentity(P[ i ], inp->zyg.defs.ngenes);
        ier = denseGETRF(P[ i ], inp->zyg.defs.ngenes, inp->zyg.defs.ngenes, pivot[ i ]);
        if (ier != 0) return 1;
    }

    *jcurPtr = TRUE;
    return 0;
}
                                                                        
int PSolve(realtype tn, N_Vector c, N_Vector fc, N_Vector r,             
        N_Vector z, realtype gamma, realtype delta, int lr,
        void *extra_data, N_Vector vtemp) {                              //needed by the old Krylov solver

    const int NRNUC = GetNNucs(&(inp->zyg.defs), inp->zyg.nnucs, tn, &(inp->zyg.times));

    realtype ***P;
    int **pivot;
    int i;
    ExtraData edata = (ExtraData) extra_data;

    N_VScale(RCONST(1.0), r, z);

    // call GSIter for Gauss-Seidel iterations
    //GSIter(gamma, z, vtemp, wdata);

    // Do backsolves for inverse of block-diagonal preconditioner factor
    P = edata->P;
    pivot = edata->pivot;
    for (i = 0; i < NRNUC; ++i) {
        denseGETRS(P[ i ], inp->zyg.defs.ngenes, pivot[ i ], &(NV_DATA_S(z)[ i ]));
    }

    return 0;
}

// This is the old slow Krylov solver
int InitKrylovSolver(realtype tzero, double stephint, double rel_tol, double abs_tol) {  //needed by the old Krylov solver
    int flag;
    neq = inp->zyg.defs.ngenes * GetNNucs(&(inp->zyg.defs), inp->zyg.nnucs, tzero, &(inp->zyg.times));
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (CheckFlag((void*) cvode_mem, "CVodeCreate", 0)) {
        printf("Error creating ODE solver\n");
        return 1;
    }

    // Call CVodeSetUserData to register the user data that will be passed
    //  around.
     
    edata = NewExtraData(GetNNucs(&(inp->zyg.defs), inp->zyg.nnucs, tzero, &(inp->zyg.times)));
    flag = CVodeSetUserData(cvode_mem, edata);
    if (CheckFlag(&flag, "CVodeSetUserData", 1)) return 1;
    edata->cvode_mem = cvode_mem;

    // Call CVodeInit to initialize the integrator memory and specify the
    // user's right hand side function, the inital time tzero, and
    // the initial dependent variable vector vars_
    // 
    flag = CVodeInit(cvode_mem, my_f, tzero, vars);
    if (CheckFlag(&flag, "CVodeInit", 1)) {
        printf("Error setting up ODE solver\n");
        return 1;
    }

    // Call CVodeSStolerances to specify the scalar relative tolerance
    // and scalar absolute tolerance
    //
    flag = CVodeSStolerances(cvode_mem, rel_tol, abs_tol);
    if (CheckFlag(&flag, "CVodeSStolerances", 1)) {
        printf("Error setting up tolerances\n");
        return 1;
    }

    // Call CVSpgmr to specify the linear solver CVSPGMR
    // with left preconditioning and the maximum Krylov dimension maxl
    //
    flag = CVSpgmr(cvode_mem, PREC_LEFT, 0);
    if (CheckFlag(&flag, "CVSPGMR", 1)) {
        printf("Error setting up linear solver CVSPGMR\n");
        return 1;
    }
    // Set the preconditioner solve and setup functions 
    flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
    if (CheckFlag(&flag, "CVSpilsSetPreconditioner", 1)) {
        printf("Error setting preconditioner\n");
        return 1;
    }

    CVodeSetInitStep(cvode_mem, stephint);
    //printf( "Krylov Solver successfully (re)initialized\n" );
    return 0;
}

void FreeKrylovSolver(void) { //before initializing again    //needed by the old Krylov solver
    if (cvode_mem != NULL) {
        CVodeFree(&cvode_mem);
        cvode_mem = NULL;
    }
    if (vars != NULL) {
        N_VDestroy_Serial(vars);
        vars = NULL;
    }
    if (edata != NULL) {     
        FreeExtraData(edata);
        edata = NULL;
    }
}

ExtraData NewExtraData(int len) {    //needed by the old Krylov solver

    int i;
    ExtraData edata;

    edata = (ExtraData) malloc(sizeof *edata);
    // set length of P and pivot
    edata->size = len;
    // allocating memory
    edata->P = (realtype ***) malloc(len * sizeof ( realtype **));
    edata->pivot = (int **) malloc(len * sizeof ( int *));
    for (i = 0; i < len; ++i) {
        (edata->P)[ i ] = newDenseMat(inp->zyg.defs.ngenes, inp->zyg.defs.ngenes);
        (edata->pivot)[ i ] = newIntArray(inp->zyg.defs.ngenes);
    }
    edata->fsave = (realtype*) malloc(neq * sizeof ( realtype));
    edata->err_weights = N_VNew_Serial(neq);

    return edata;
}

void FreeExtraData(ExtraData edata) {   //needed by the old Krylov solver
    
    int i;
    for (i = 0; i < edata->size; ++i) {
        destroyMat((edata->P)[ i ]);
        destroyArray((edata->pivot)[ i ]);
    }
    free(edata->P);
    free(edata->pivot);
    free(edata->fsave);
    N_VDestroy_Serial(edata->err_weights);
    free(edata);
}
*/

int
InitKrylovVariables( double *vin, int n ) {

    int i;

    vars = N_VNew_Serial( n );
    if( CheckFlag( ( void * ) vars, "N_VNew_Serial", 0 ) )
        return 1;
    for( i = 0; i < n; ++i ) {
        NV_Ith_S( vars, i ) = vin[i];
    }
    return 0;
}

int
CheckFlag( void *flagvalue, char *funcname, int opt ) {
    // function taken from examples supplied by Sundial
    int *errflag;

    // Check if SUNDIALS function returned NULL pointer - no memory allocated
    if( opt == 0 && flagvalue == 0 ) {
        fprintf( stderr, "SUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname );
        return 1;
    } else if( opt == 1 ) {
        // Check if flag < 0
        errflag = ( int * ) flagvalue;
        if( *errflag < 0 ) {
            fprintf( stderr, "SUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag );
            return 1;
        }
    } else if( opt == 2 && flagvalue == 0 ) {
        // Check if function returned NULL pointer - no memory allocated
        fprintf( stderr, "MEMORY_ERROR: %s() failed - returned 0 pointer\n\n", funcname );
        return 1;
    }
    return 0;
}

/**  Krylov: propagates vin (of size n) from tin to tout by BDF (Backward  
 *           Differential Formulas and use of a Newton-Krylov method with  
 *           preconditioning to avoid the costly computation of the        
 *           jacobian.                                                     
 ***************************************************************************
 *                                                                         
 * This solver was written by Anton Crombach, October 2010                 
 *                                                                          
 * I realised I cannot evaluate only part of the derivative (f.i. only the
 * production/decay part)... so let's see what happens if I evaluate the full
 * derivative for preconditioning. That actually works better.
 * 
 * This is actually the Krylov Band solver. */
void
Krylov( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * sinput, Input * input ) {
    int flag, i, j;
    realtype t, tstop;
    double *divtimes, *divdurations;
    inp = input;
    si = sinput;

    /* If nothing to do, return */
    if( fabs( tin - tout ) < 1e-6 )
        return;
    if( cvode_mem != 0 ) {
        FreeBandSolver(  );
    }
    InitKrylovVariables( vin, n );
    InitBandSolver( tin, stephint, accuracy, accuracy );

    /* Krylov solver looks ahead and then gets confused by the change
       in number of equations, so we need to set a stop time beyond which it
       is not allowed to look */
    i = inp->zyg.defs.ndivs - 1;
    j = -1;

    divtimes = inp->zyg.times.div_times;
    divdurations = inp->zyg.times.div_duration;
    // bounded linear search for any division that still has to happen
    while( i != j )
        if( tout > divtimes[i] )
            --i;
        else
            j = i;
    // there is still a division to happen
    if( i > -1 ) {
        // perhaps mitosis even before
        tstop = divtimes[i];
        if( tout < divtimes[i] - divdurations[i] ) {
            tstop = divtimes[i] - divdurations[i];
        }
    } else {
        // otherwise don't look beyond gastrulation
        tstop = inp->zyg.times.gast_time;
    }
    CVodeSetMaxStep(cvode_mem,(realtype)1000000);
    CVodeSetMaxStep(cvode_mem,(realtype)50);
    CVodeSetMaxErrTestFails(cvode_mem, 1000);

    CVodeSetStopTime( cvode_mem, tstop );
    /* old code, works if networks would always be well-behaved */
    flag = CVode( cvode_mem, tout, vars, &t, CV_NORMAL );
    if( CheckFlag( &flag, "CVode", 1 ) )
        return;
    /* copy vars into vout */
    for( i = 0; i < n; ++i ) {
        vout[i] = NV_Ith_S( vars, i );
    }
}

/** wrapper function - to call the derivative */
int
my_f_band( realtype t, N_Vector y, N_Vector ydot, void *extra_data ) {  
    
    int n = NV_LENGTH_S( y );
    p_deriv( NV_DATA_S( y ), t, NV_DATA_S( ydot ), n, si, inp );
    return 0;
}

int
InitBandSolver( realtype tzero, double stephint, double rel_tol, double abs_tol ) {
    int flag;
    neq = inp->zyg.defs.ngenes * GetNNucs( &( inp->zyg.defs ), inp->zyg.nnucs, tzero, &( inp->zyg.times ) );
    cvode_mem = CVodeCreate( CV_BDF, CV_NEWTON );
    if( CheckFlag( ( void * ) cvode_mem, "CVodeCreate", 0 ) ) {
        printf( "Error creating ODE solver\n" );
        return 1;
    }


    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function, the inital time tzero, and
     * the initial dependent variable vector vars_
     */
    flag = CVodeInit( cvode_mem, my_f_band, tzero, vars );
    if( CheckFlag( &flag, "CVodeInit", 1 ) ) {
        printf( "Error setting up ODE solver\n" );
        return 1;
    }

    /* Call CVodeSStolerances to specify the scalar relative tolerance
     * and scalar absolute tolerance
     */
    flag = CVodeSStolerances( cvode_mem, rel_tol, abs_tol );
    if( CheckFlag( &flag, "CVodeSStolerances", 1 ) ) {
        printf( "Error setting up tolerances\n" );
        return 1;
    }

    flag = CVBand( cvode_mem, neq, inp->zyg.defs.ngenes + 1, inp->zyg.defs.ngenes + 1 );
    if( CheckFlag( &flag, "CVBand", 1 ) ) {
        printf( "Error setting band linear solver\n" );
        return 1;
    }

    /* Set step size hint, pass 0.0 to use internal estimate */
    //stephint = 0.0;
    CVodeSetInitStep( cvode_mem, stephint );

    //printf( "Krylov Solver successfully (re)initialized\n" );
    return 0;
}

void
FreeBandSolver( void ) {        //before initializing again

    if( cvode_mem != NULL ) {
        CVodeFree( &cvode_mem );
        cvode_mem = NULL;
    }
    if( vars != NULL ) {
        N_VDestroy_Serial( vars );
        vars = NULL;
    }
}

/*
void gaussSeidel( realtype gamma, N_Vector z, N_Vector aux, int nrnuc ) {
    // perform max GS_ITER_MAX=5 Gauss-Seidel iterations to compute an
    // approximation to P^inv * z, where P = I - gamma *J_diff, and J_diff
    // represents the diffusion part of the jacobian.
    //
    // The solution/answer is stored in z, aux is a temporary vector.
    //
    // Note: we use some vector methods defined in reactions.hh
    
    const double D[ 4 ] = { 0.237, 0.3, 0.115, 0.3 };
    
    const int GS_ITER_MAX = 5;
    const int NN = NV_LENGTH_S( z );
    const int GG = defs.ngenes;
    
    int i, iter, ap, base;    
    realtype *xd = NV_DATA_S( aux );
    realtype *zd = NV_DATA_S( z );
    realtype beta[ GG ], beta2[ GG ], coef[ GG ];
    realtype aux1;
    
    // write matrix as P = D - L - U
    // and load local arrays
    // coef is inverse of diagonal item
    for( i = 0; i < GG; ++i ) {
        realtype temp = 1.0 / ( 1.0 + 2.0 * gamma * D[ i ] );
        
        beta[ i ] = gamma * D[ i ] * temp;
        beta2[ i ] = 2 * gamma * D[ i ] * temp;
        coef[ i ] = temp;
    }

    // initialize loop; set in each nucleus the diffusion terms
    for( i = 0; i < nrnuc; ++i ) {
        base = i * GG;
        vec_prod( xd+base, coef, zd+base, GG );
    }
    
    // and do the loop
    for( iter = 0; iter < GS_ITER_MAX; ++iter ) {
        
        // calculate (D-inverse)*U*x if not first iteration
        if( iter > 0 ) {
            // most anterior nucleus
            vec_prod( xd+base, beta2, xd+base+GG, GG );
            
            // middle
            for( ap = 1; ap < nrnuc-1; ++ap ) {
                // gap genes
                base = ap * GG;
                vec_prod( xd+base, beta, xd+base+GG, GG );
            }
            
            // most posterior nucleus
            base = (nrnuc-1) * GG;
            vec_zero( xd+base, GG );
        }
        
        // overwrite x with [ ( I - (D-inverse)*L )-inverse ]*x
        
        // skip anterior nucleus
        for( ap = 1; ap < nrnuc-1; ++ap ) {
            // gap genes
            base = ap * GG;
            vec_inc_by_prod( xd+base, beta, xd+base-GG, GG );
        }
        // most posterior nucleus is different (why?)
        base = (nrnuc-1) * GG;
        vec_inc_by_prod( xd+base, beta2, xd+base-GG, GG );
        
        // and increment x to z
        N_VLinearSum( 1.0, z, 1.0, aux, z );
    }
}

void vec_prod( realtype u[], realtype v[], realtype w[], int n ) {
    int i;
    for( i = 0; i < n; ++i ) u[ i ] = v[ i ] * w[ i ];
}

void vec_inc_by_prod( realtype u[], realtype v[], realtype w[], int n ) {
    int i;
    for( i = 0; i < n; ++i ) u[ i ] += v[ i ] * w[ i ];
}

void vec_zero( realtype u[], int n ) {
    int i;
    for( i = 0; i < n; ++i ) u[ i ] = 0;
}
 */

/** get some info and write it out */
void
writeInfo(  ) {
    long int lenrw, leniw;
    long int lenrwLS, leniwLS;
    long int nst, nfe, nsetups, nni, ncfn, netf;
    long int nli, npe, nps, ncfl, nfeLS;
    int flag;

    flag = CVodeGetWorkSpace( cvode_mem, &lenrw, &leniw );
    CheckFlag( &flag, "CVodeGetWorkSpace", 1 );
    flag = CVodeGetNumSteps( cvode_mem, &nst );
    CheckFlag( &flag, "CVodeGetNumSteps", 1 );
    flag = CVodeGetNumRhsEvals( cvode_mem, &nfe );
    CheckFlag( &flag, "CVodeGetNumRhsEvals", 1 );
    flag = CVodeGetNumLinSolvSetups( cvode_mem, &nsetups );
    CheckFlag( &flag, "CVodeGetNumLinSolvSetups", 1 );
    flag = CVodeGetNumErrTestFails( cvode_mem, &netf );
    CheckFlag( &flag, "CVodeGetNumErrTestFails", 1 );
    flag = CVodeGetNumNonlinSolvIters( cvode_mem, &nni );
    CheckFlag( &flag, "CVodeGetNumNonlinSolvIters", 1 );
    flag = CVodeGetNumNonlinSolvConvFails( cvode_mem, &ncfn );
    CheckFlag( &flag, "CVodeGetNumNonlinSolvConvFails", 1 );

    flag = CVSpilsGetWorkSpace( cvode_mem, &lenrwLS, &leniwLS );
    CheckFlag( &flag, "CVSpilsGetWorkSpace", 1 );
    flag = CVSpilsGetNumLinIters( cvode_mem, &nli );
    CheckFlag( &flag, "CVSpilsGetNumLinIters", 1 );
    flag = CVSpilsGetNumPrecEvals( cvode_mem, &npe );
    CheckFlag( &flag, "CVSpilsGetNumPrecEvals", 1 );
    flag = CVSpilsGetNumPrecSolves( cvode_mem, &nps );
    CheckFlag( &flag, "CVSpilsGetNumPrecSolves", 1 );
    flag = CVSpilsGetNumConvFails( cvode_mem, &ncfl );
    CheckFlag( &flag, "CVSpilsGetNumConvFails", 1 );
    flag = CVSpilsGetNumRhsEvals( cvode_mem, &nfeLS );
    CheckFlag( &flag, "CVSpilsGetNumRhsEvals", 1 );

    fprintf( stderr, "! Statistics..\n" );
    fprintf( stderr, "lenrw   = %5ld     leniw   = %5ld     ", lenrw, leniw );
    fprintf( stderr, "lenrwls = %5ld     leniwls = %5ld\n", lenrwLS, leniwLS );
    fprintf( stderr, "nst     = %5ld\n", nst );
    fprintf( stderr, "nfe     = %5ld     nfeLS   = %5ld     nfetot  = %5ld\n", nfe, nfeLS, ( nfe + nfeLS ) );
    fprintf( stderr, "nni     = %5ld     nli     = %5ld\n", nni, nli );
    fprintf( stderr, "nsetups = %5ld     netf    = %5ld\n", nsetups, netf );
    fprintf( stderr, "npe     = %5ld     nps     = %5ld\n", npe, nps );
    fprintf( stderr, "ncfn    = %5ld     ncfl    = %5ld\n", ncfn, ncfl );
}
