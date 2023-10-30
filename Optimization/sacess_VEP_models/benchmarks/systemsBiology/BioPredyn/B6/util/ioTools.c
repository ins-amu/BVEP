/**
 *   @file ioToold.c
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



#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ioTools.h"


/** FindSection: This function finds a given section of the input file & 
 *                returns a pointer positioned to the first record of that
 *                section. Section titles should be passed without the 
 *                preceding '$'. If it can't find the right section, the 
 *                function returns NULL.                                      
 */
FILE *
FindSection( FILE * fp, char *input_section ) {
    int c;                      /* input happens character by character */
    int nsought;                /* holds length of section title */
    char *base;                 /* string for section title */
    long looksite;              /* file position indicator */
    rewind( fp );               /* start looking at the beginning of the file */
    nsought = strlen( input_section );
    base = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    /*** while loop goes through the file character by character... ************/
    while( ( c = getc( fp ) ) != EOF ) {        /* ...until if finds a '$' */
        if( c == '$' ) {        /* found a sectioning control string */
            looksite = ftell( fp );     /* where are we? */
            base = fgets( base, MAX_RECORD, fp );       /* get sect title (without $) */
            if( !( strncmp( base, input_section, nsought ) ) ) {
                fseek( fp, looksite, 0 );       /* found the sought string: reposi- */
                fscanf( fp, "%*s\n" );  /* tion the pointer to '$', then */
                free( base );
                return ( fp );  /* advance the pointer to after the */
            } /* section title and return it */
            else {              /* didn't find it: skip this control */
                fseek( fp, looksite, 0 );       /* record, keep looking */
                fscanf( fp, "%*s" );    /* NOTE: "%*s" advances pointer */
            }                   /* without assignment */
        }
    }

    free( base );
    return ( NULL );            /* couldn't find the right section */
}

/** KillSection: erases the section with 'title' from the file 'fp' */
void
KillSection( char *filename, char *title ) {
    size_t length;              /* length of title string */

    char *fulltitle;            /* title incl. $ */
    char *temp;                 /* temporary file name */
    char *record;               /* record to be read and written */
    char *record_ptr;           /* pointer used to remember record for 'free' */

    char *shell_cmd;            /* used by system below */

    FILE *fp;                   /* name of file where section needs to be killed */
    FILE *tmpfile;              /* name of temporary file */


    fulltitle = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    temp = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    record = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    shell_cmd = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    record_ptr = record;        /* this is to remember record for 'free' */

    fp = fopen( filename, "r" );        /* open file for reading */
    if( !fp )
        error( "KillSection: error opening file %s", filename );

    temp = strcpy( temp, "killXXXXXX" );        /* required by mkstemp() */
    if( mkstemp( temp ) == -1 ) {       /* get unique name for temp file */
        error( "KillSection: error creating temporary file" );
    }

    tmpfile = fopen( temp, "w" );       /* ... and open it for writing */
    if( !tmpfile )
        error( "KillSection: error opening temporary file" );

    if( !FindSection( fp, title ) )
        error( "KillSection: section to be killed not found" );
    rewind( fp );

    /* remove section by simply ignoring it while copying the file */

    fulltitle = strcpy( fulltitle, "$" );
    fulltitle = strcat( fulltitle, title );
    length = strlen( fulltitle );

    while( strncmp( ( record = fgets( record, MAX_RECORD, fp ) ), fulltitle, length ) )
        fputs( record, tmpfile );

    while( strncmp( ( record = fgets( record, MAX_RECORD, fp ) ), "$$", 2 ) );

    do {
        record = fgets( record, MAX_RECORD, fp );
        if( !record )
            break;
    } while( strncmp( record, "$", 1 ) );

    if( record )
        fputs( record, tmpfile );

    while( ( record = fgets( record, MAX_RECORD, fp ) ) )
        fputs( record, tmpfile );

    fclose( fp );
    fclose( tmpfile );

    /* rename tmpfile into new file */

    sprintf( shell_cmd, "cp -f %s %s", temp, filename );

    if( -1 == system( shell_cmd ) )
        error( "KillSection: error renaming temp file %s", temp );

    if( remove( temp ) )
        error( "KillSection: temp file %s could not be deleted", temp );


    /* clean up */

    free( record_ptr );
    free( temp );
    free( shell_cmd );
    free( fulltitle );
}
