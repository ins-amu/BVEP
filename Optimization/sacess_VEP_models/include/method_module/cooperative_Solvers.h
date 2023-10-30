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
    
    
int cooperative_asynchronous_DE(experiment_total *, double(*fitnessfunction)(double*,void *), result_solver *, double , double , int);

int cooperative_synchronous_DE (experiment_total *, double(*fitnessfunction)(double*,void *), result_solver *, double , double , int);


#endif

#ifdef	__cplusplus
}
#endif

#endif	/* SERSOLVERS_H */

