/** 
 * File:   fly_sa.h
 * Author: Damjan Cicin-Sain
 * Contact: damjan.cicin@crg.es
 *
 * Created on June 17, 2013, 5:26 PM
 */

//#include "global.h"

#ifndef _FLY_SA_H
#define	_FLY_SA_H

#ifdef	__cplusplus
extern "C" {
#endif

    
/** MoveX: This function actually does almost everything.
 * First it creates a static Input structure 'inp', where it puts all the 
 * information from the input file. This part is executed only once (when init == 1).
 * Then it creates a ScoreOutput structure 'out' where the score, penalty 
 * and residual vectors will be stored.
 * At the end it runs the score function, where all the calculation is done.
 * This function is called from the mex file and it is used to connect with the matlab
 * ssm code. Once that the code is translated to c, we plan to call MoveSA instead of 
 * MoveX. 
 */
void
MoveX( double *x, int *mask, ScoreOutput * out, Files * files, int init, int jacobian, int solver );


#ifdef	__cplusplus
}
#endif

#endif	/* _FLY_SA_H */

