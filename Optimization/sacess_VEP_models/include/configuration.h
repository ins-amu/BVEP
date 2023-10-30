/* 
 * File:   configuration.h
 * Author: david
 *
 * Created on 11 de junio de 2013, 19:53
 */

#ifndef CONFIGURATION_H
#define	CONFIGURATION_H

#ifdef	__cplusplus
extern "C" {
#endif


    
#include <math.h>
#include <time.h>   
#include <sys/time.h>     
    

#define URAND_BBOB   ((double)rand()/((double)RAND_MAX + 1.0))
#define INITRAND_BBOB  srand(time(NULL))
#define INITRAND_PAR_BBOB(id) srand(time(NULL)*(id+1))
#define ISIZE_BBOB 1
#define ISIZE_SYSTBIO   1
#define EPSILON_LOG 1e-6
    
#define e_definition expl(1.0) 

#ifdef	__cplusplus
}
#endif

#endif	/* CONFIGURATION_H */

