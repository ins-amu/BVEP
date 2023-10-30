/** 
 * @file:   ggn.h
 * @author: Damjan Cicin-Sain
 * @email: damjan.cicin@crg.es
 *
 * Created on June 30, 2010, 12:10 PM
 */

#ifndef _GGN_H
#define	_GGN_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "global.h"
#include "mathLib.h"            /* Trunc function is in there */    
#include "fly_sa.h"             /* MoveX function is in here */


void returnbounds( double *, double * );
void returnupperbound(double *);
void returnlowerbound(double *);

double ggn(double x[], int );

double fitnessfunctionB6(double *);

#ifdef	__cplusplus
}
#endif

#endif	/* _GGN_H */

