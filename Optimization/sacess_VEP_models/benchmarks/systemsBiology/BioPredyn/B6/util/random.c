/**
 *                                                               
 *   @file random.c                                              
 *                                                               
 *****************************************************************
 *                                                               
 *   written by Yoginho                                          
 *                                                               
 *****************************************************************
 *                                                               
 *   functions for initializing and running dSFMT() random       
 *   number generator                                            
 *                                                               
 *****************************************************************
 *                                                               
 * Copyright (C) 1989-2003 John Reinitz                          
 * the full GPL copyright notice can be found in lsa.c           
 *                                                               
 */

#include <stdlib.h>
#include <dSFMT.h>
#include <dSFMT_str_state.h>
#include <random.h>



/* STATIC VARIABLES ********************************************************/
/* an array needed by dSFMT */

static dsfmt_t dsfmt;

/*** RANDOM NUMBER FUNCTIONS ***********************************************/

/** InitRand: initializes dSFMT random number generator by making seed 
 *             static to random.c                                       
 */
void
InitRand( int seed ) {
   
    dsfmt_init_gen_rand(&dsfmt, seed);
}

/** RestoreRand: restores dSFMT random number generator state by making seed 
 *             static to random.c 
 */
void
RestoreRand(char* restored) {
   dsfmt_str_to_state(&dsfmt, restored, NULL);
}


/** RandomReal: returns a random real number between 0 and 1 using the 
 *               random number generator of choice                      
 */
double
RandomReal( void ) {
    double i;

    i = dsfmt_genrand_close_open(&dsfmt);
    
    return ( i );
}

/** RandomInt: returns a random integer between 0 and max using the 
 *              random number generator of choice
 */
int
RandomInt( int max ) {
    register int i;

    i = ( ( int ) ( dsfmt_genrand_close_open(&dsfmt) * max ) ) % max;

    return ( i );
}

/** GetDSFMTState: returns the dSFMT state as a string, whish is used to 
 *                  initialize dSFMT(); used for saving the dSFMT state    
 *                  in a state file                                        
 */
char*
GetDSFMTState( void ) {
    char* p;
    char* prefix = NULL;
    
    p = dsfmt_state_to_str(&dsfmt, prefix);

    return p;
}