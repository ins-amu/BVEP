/* 
 * File:   serSolvers.h
 * Author: david
 *
 * Created on 6 de junio de 2013, 12:09
 */

#ifndef SERSOLVERS_H
#define	SERSOLVERS_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef MPI2
    #include <structure_paralleltestbed.h>
    
int  execute_Solver(experiment_total *, result_solver *, double, void* function( double *,  void*));

int execute_parallel_solver (experiment_total *, result_solver *, long , double, void* function( double *,  void*));

int execute_serial_solver(experiment_total *, result_solver *, long, double, void* function( double *,  void*));

#endif

#ifdef	__cplusplus
}
#endif

#endif	/* SERSOLVERS_H */

