/**
 *                                                               
 *   @file random.h                                              
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
 *****************************************************************/

#ifndef RANDOM_INCLUDED
#define RANDOM_INCLUDED


/*** FUNCTION PROTOTYPES ***************************************************/

/** InitRand: initializes dSFMT random number generator by making seed 
 *             static to random.c                                       
 */
void InitRand( int seed );

/** RestoreRand: restores dSFMT random number generator state by making seed
 *             static to random.c 
 */
void RestoreRand( char* restored );

/** RandomReal: returns a random real number between 0 and 1 using the 
 *               random number generator of choice                      
 */
double RandomReal( void );

/** RandomInt: returns a random integer between 0 and max using the 
 *              random number generator of choice                          
 */
int RandomInt( int max );

/** GetDSFMTState: returns the dSFMT state as a string, whish is used to 
 *                  initialize dSFMT(); used for saving the dSFMT state    
 *                  in a state file                                        
 */
char *GetDSFMTState( void );

#endif
