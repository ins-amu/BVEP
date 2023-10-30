
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "benchmarkshelper.h"

/* Checks if sDir exists, 
   creates it if not
   checks if is writable thereafter
   Fatal ERROR if anything fails
*/
void dirOK(char *sDir)
{
    char sLoc[256], fileNameLoc[1024];
    if (strlen(sDir) > 1024)
        ERROR("pathName too long %s\n",sDir); 
    sprintf(sLoc, "test -d %s", sDir); 
    if ( system(sLoc) )  /* does  NOT exist */
    {
        sprintf(sLoc, "mkdir %s", sDir);
        system(sLoc);
    }
    /* could we make the dir ? */
    sprintf(sLoc, "test -d %s", sDir);
    if ( system(sLoc) )  /* does  NOT exist */
        ERROR("Failed to create dir %s", sDir);
    /* is dir writable */
    createFullFileName(fileNameLoc, sDir, (char*)"toto");
    sprintf(sLoc, "touch %s", fileNameLoc);
    if ( system(sLoc) )  /* failed */
        ERROR("Problem writing in dir %s", sDir);
    else
        unlink(fileNameLoc);
    return;
}
