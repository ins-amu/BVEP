/**
 *   @file ioTools.h
 * 
 *   Copyright (C) 2006 by Yves Fomekong Nanfack, F.219,+31 20 525 7530,   
 *   yvesf@science.uva.nl   *
 *                                                                         
 *   This program is free software; you can redistribute it and/or modify  
 *   it under the terms of the GNU General Public License as published by  
 *   the Free Software Foundation; either version 2 of the License, or     
 *   (at your option) any later version.                                   
 *                                                                         
 *   This program is distributed in the hope that it will be useful,       
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 *   GNU General Public License for more details.                          
 *                                                                         
 *   You should have received a copy of the GNU General Public License     
 *   along with this program; if not, write to the                         
 *   Free Software Foundation, Inc.,                                       
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
 */


#ifndef IOTOOLS_INCLUDED
#define IOTOOLS_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "error.h"
#include <global.h>


/* A FUNCTION WHICH IS NEEDED BY ALL OTHER READING FUNCS *******************/

/** FindSection: This function finds a given section of the input file & 
 *                returns a pointer positioned to the first record of that
 *                section. Section titles should be passed without the 
 *                preceding '$'. If it can't find the right section, the 
 *                function returns NULL.                                      
 */
FILE *FindSection( FILE * fp, char *input_section );

/** KillSection: erases the section with 'title' from the file 'fp' */
void KillSection( char *filename, char *title );

#endif
