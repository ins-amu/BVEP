/**
 * @file translate.c                                           
 * @author JR, modified by Yoginho
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                                 
 * @brief Makes the list of pointers that the annealer uses to make moves.
 *
 * SENDING PARAMATERS TO THE ANNEALER                            
 *                                                               
 * The idea is to conceal as much   
 * problem specific details as possible from the annealer. Some  
 * of the more general ways of doing this must wait for future   
 * code. e.g., when there are multiple possible right hand-side- 
 * of-the-ODE derivative functions, each of these must have a    
 * translation function, and possibly a small family of other    
 * functions as well. These would need to be handled as arrays   
 * of pointers to functions for each choice, etc. Also, the      
 * ParamList structure now has ranges, because it may be important 
 * in some contexts for the annealer to know the limits of  
 * the search space. If penalty funcs are being used (pen_vec != 
 * NULL), a little more work is required, and the coder might,   
 * say, write a function that goes in this file which calculates 
 * the effective range (say for penalty <= 10000) of one param   
 * while asuming the values of the others to be fixed.           
 *                                                               
 * SPECIFYING WHICH PARAMETERS ARE ANNEALED, WHICH HELD CONSTANT 
 *                                                               
 * Parameters are specified as fixed one by one in the parameter 
 * file, by setting their entries in the $tweak section to zero. 
 * this entry is then used to construct the array of pointers to 
 * the parameters that are to be tweaked.                        
 *                                                               
 * The full GPL copyright notice can be found in lsa.c           
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "translate.h"


/* Translation: installing an array of pointer to params-to-be-tweaked */

/** @brief Creates an array of pointers that point to the parameters 
 * to be tweaked.
 *
 * The PArrPtr that is returned also includes pointers to the 
 * corresponding parameter ranges for each parameter, although 
 * this feature is not yet used anywhere in the annealing code
 *
 * CAUTION: InitZygote and InitScoring have to be called first!        
 */
PArrPtr
Translate( Input * inp ) {
    int i, j;
    int index = 0;

    PArrPtr plist;              /* local copy of the PArrPtr */
    ParamList *p;               /* local copy of parameter list */

    EqParms *parm;              /* local copy of the EqParms */
    SearchSpace *limits;        /* local copy of the SearchSpace */

    int max_size = 0;           /* maximum size of parameter array */
    int max_elem = 0;           /* maximum # of elements of parameter array */

    /* initially we allocate max_elem members of the PArrPtr array in max_size */
    /* bytes. The actual space needed will be less than or equal to this.      */

    max_elem = inp->zyg.defs.ngenes * inp->zyg.defs.ngenes + inp->zyg.defs.ngenes * inp->zyg.defs.egenes + 6 * inp->zyg.defs.ngenes;
    max_size = max_elem * sizeof( ParamList );

    /* initialize the ParamList */
    p = ( ParamList * ) malloc( max_size );

    for( i = 0; i < max_elem; ++i ) {
        p[i].param = NULL;
        p[i].param_range = NULL;
        //p[ i ].scale = NULL;
    }
    plist.size = 0;

    /* Get limits and parameters */
    parm = &( inp->zyg.parm );  /* get pointer to EqParm struct in zygotic.c */
    limits = GetLimits( inp );  /* get pointer to SearchSpace in score.c */

    /* if we are using a penalty, not all of Limits will have been allocated   */
    /* and we must take care not to take indices of unallocated arrays!!!      */

    if( limits->pen_vec != NULL )       /* we are using a penalty, hence pen_vec */
        plist.pen_vec = limits->pen_vec;
    else                        /* or: explicit ranges for everything */
        plist.pen_vec = NULL;
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {       /* pointers to R stuff */
        if( inp->twe.Rtweak[i] == 1 ) {
            p[index].param = &parm->R[i];
            p[index].param_range = limits->Rlim[i];
            index++;
        }
    }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to T stuff */
        for( j = 0; j < inp->zyg.defs.ngenes; j++ )
            if( inp->twe.Ttweak[i * inp->zyg.defs.ngenes + j] == 1 ) {
                p[index].param = &( parm->T[j + i * inp->zyg.defs.ngenes] );
                if( plist.pen_vec == NULL )
                    p[index].param_range = limits->Tlim[j + i * inp->zyg.defs.ngenes];
                index++;
            }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to E stuff */
        for( j = 0; j < inp->zyg.defs.egenes; j++ )
            if( inp->twe.Etweak[i * inp->zyg.defs.egenes + j] == 1 ) {
                p[index].param = &( parm->E[j + i * inp->zyg.defs.egenes] );
                if( plist.pen_vec == NULL )
                    p[index].param_range = limits->Elim[j + i * inp->zyg.defs.egenes];
                index++;
            }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to m stuff */
        if( inp->twe.mtweak[i] == 1 ) {
            p[index].param = &parm->m[i];
            if( plist.pen_vec == NULL )
                p[index].param_range = limits->mlim[i];
            index++;
        }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to h stuff */
        if( inp->twe.htweak[i] == 1 ) {
            p[index].param = &parm->h[i];
            if( plist.pen_vec == NULL )
                p[index].param_range = limits->hlim[i];
            index++;
        }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to d stuff */
        if( inp->twe.dtweak[i] == 1 ) {
            p[index].param = &parm->d[i];
            p[index].param_range = limits->dlim[i];
            index++;
            if( ( inp->zyg.defs.diff_schedule == 'A' ) || ( inp->zyg.defs.diff_schedule == 'C' ) )
                break;
        }
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to lambda stuff */
        if( inp->twe.lambdatweak[i] == 1 ) {
            p[index].param = &parm->lambda[i];
            p[index].param_range = limits->lambdalim[i];
            index++;
        }

    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) /* pointers to tau stuff */
        if( inp->twe.tautweak[i] == 1 ) {
            p[index].param = &parm->tau[i];
            p[index].param_range = limits->taulim[i];
            index++;
        }
    // reallocate memory for smaller-than-or-equal-to-maxsize array
    p = ( ParamList * ) realloc( ( void * ) p, ( index * sizeof( ParamList ) ) );

    // finish up initializing PArrPtr
    plist.size = index;
    plist.array = p;
    return plist;
}
