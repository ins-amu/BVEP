/**
 * @file solvers.h
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Contains the interface to the solver functions.
 * 
 * In other words, solver function prototypes and the p_deriv global which 
 * serves as a definition of the interface to the derivative func.  
 *                                                               
 * @note \em ONLY general solvers allowed here; they \em MUST comply to the 
 * generic solver interface (see solvers.c for details)
 */

#ifndef SOLVERS_INCLUDED
#define SOLVERS_INCLUDED

/*
 * added by Anton Crombach, October 2010
 */
#include <cvode/cvode.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_spgmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include "maternal.h"

/*** GLOBAL VARIABLES ******************************************************/

extern void ( *p_deriv ) ( double *, double, double *, int, SolverInput *, Input * );
extern void ( *p_deriv_lite ) ( double *, double, double *, int, SolverInput *, Input * );
extern void ( *d_deriv ) ( double *, double **, double, double *, int, SolverInput *, Input * );
extern void ( *p_jacobn ) ( double, double *, double *, double **, int, SolverInput *, Input * );




/** FUNCTION PROTOTYPES ***************************************************/

/*** Euler: propagates vin (of size n) from tin to tout by the Euler 
 *          method; the result is returned by vout                         
 */
void Euler( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



/** Meuler: propagates vin (of size n) from tin to tout by the Modified 
 *           Euler method (this is NOT the midpoint method, see Rk2()); 
 *           the result is returned by vout                             
 */
void Meuler( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



/** Heun: propagates vin (of size n) from tin to tout by Heun's method 
 *         the result is returned by vout                                  
 */
void Heun( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



/** Rk2: propagates vin (of size n) from tin to tout by the Midpoint or 
 *        Second-Order Runge-Kutta method; the result is returned by vout  
 */
void Rk2( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



/** Rk4: propagates vin (of size n) from tin to tout by the Fourth-Order 
 *        Runge-Kutta method; the result is returned by vout             
 ***************************************************************************
 *                                                                         
 * written by Joel Linton (somewhere around 1998)                          
 * fixed and modified by Yoginho (somewhere around 2001)                   
 *                                                                         
 */
void Rk4( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );


/** Rkck: propagates vin (of size n) from tin to tout by the Runge-Kutta 
 *         Cash-Karp method, which is an adaptive-stepsize Rk method; it   
 *         uses a fifth-order Rk formula with an embedded forth-oder for-  
 *         mula for calucalting the error; its result is returned by vout  
 ***************************************************************************
 *                                                                         
 * This solver was written by Marcel Wolf, Spring 2002.                    
 *                                                                         
 */
void Rkck( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );


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
void Rkf( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



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
void Milne( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



/**    Adams: propagates vin (of size n) from tin to tout by Adams-Moulton 
 *            which is an implicit predictor-corrector method of second    
 *            order; the result is returned by vout                        
 ***************************************************************************
 *                                                                         
 * This solver was implemented by Konstantin Koslov, Spring 2002           
 * Slightly modified by Manu, July 2002                                    
 *                                                                         
 */
void Adams( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



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
void BuSt( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );



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
void BaDe( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );


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
void Krylov( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );

int InitKrylovVariables( double *vin, int n );

int InitBandSolver( realtype tzero, double stephint, double rel_tol, double abs_tol );
void FreeBandSolver( void );

int CheckFlag( void *flagvalue, char *funcname, int opt );

/** wrapper function - to call the derivative */
int my_f_band( realtype t, N_Vector y, N_Vector ydot, void *extra_data );

/** WriteSolvLog: write to solver log file */
void WriteSolvLog( char *solver, double tin, double tout, double h, int n, int nderivs, FILE * slog );

/**
 * Delay solver - by Manu?
 */
void SoDe( double *vin, double *vout, double tin, double tout, double stephint, double accuracy, int n, FILE * slog, SolverInput * si, Input * inp );

double *Construct_Discont_Array( double range, double *taus, int n, double *starts, int sn, int *disc_size );

int compare( double *x, double *y );

int y_delayed( double ***vd, int n, double *rktimes, double *tau, double *grid, double **vdone, double **deriv1, double **deriv2, double **deriv3,
               double **deriv4, int gridsize, double accu, SolverInput * si, Input * inp );

void DCERk32( double **vatt, int n, double *tarray, int tpoints, double *darray, int dpoints, double stephint, double accuracy, SolverInput * si, Input * inp );

void CE( double t, double *vans, double tbegin, double *v_at_tbegin, double ech, double *d1, double *d2, double *d3, double *d4, int n );

void DivideHistory( double t1, double t2, Zygote * zyg );
void FreeDelaySolver( void );
void InitDelaySolver( void );

#endif
