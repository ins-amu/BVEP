/**
 *                                                               
 *   @file error.h                                               
 *                                                               
 *****************************************************************
 *                                                               
 *   written by JR, modified by Yoginho                          
 *                                                               
 *****************************************************************
 *                                                               
 *   contains error and warning functions                        
 *                                                               
 */
#pragma once

#ifndef ERROR_INCLUDED
#define ERROR_INCLUDED

#include <global.h>


/*** FUNCTION PROTOTYPES ***************************************************/

/** ERROR HANDLING FUNCTIONS 
 *   The two routines 'error' and 'warning' print error messages with value.
 *   'error'then exits, while 'warning' returns to the calling function.   
 *                                                                         
 *   Both functions are called analogous to the way you would call printf. 
 *   The only conversion specs that are legal are %c, %d, %f, %g and %s.   
 *   No modifiers (h,l or float precision or size specifiers are allowed). 
 *                                                                          
 *   No % sign means just print msg and exit.                              
 *                                                                         
 */
void error( const char *format, ... );

/** ERROR HANDLING FUNCTIONS 
 *   The two routines 'error' and 'warning' print error messages with value.
 *   'error'then exits, while 'warning' returns to the calling function.   
 *                                                                         
 *   Both functions are called analogous to the way you would call printf. 
 *   The only conversion specs that are legal are %c, %d, %f, %g and %s.   
 *   No modifiers (h,l or float precision or size specifiers are allowed). 
 *                                                                                                                                                  
 *   No % sign means just print msg and exit.                              
 *                                                                         
 */
void warning( const char *format, ... );

/** file_error prints an file handling error using perror(); it only 
 *   needs the name of the calling function as an argument                 
 */
void file_error( const char *call_name );

/** PrintMsg: prints a string to stderr, then quits with exit status 0;
 *             useful for printing help and usage messages             
 */
void PrintMsg( const char *msg, int exit_status );

#endif
