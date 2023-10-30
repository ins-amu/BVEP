/*
 * findStates.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */
 
/*$Id: findStates.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <stdio.h>
#include <stdlib.h>

int *findStates(int **adjMatrix, int n)
{
    int *stateVec=malloc(n*sizeof(int));
    int i,j;
    for (j = 0; j <n; j++)
    {
        stateVec[j]=0;
        for(i=0;i<n;i++)
        {
            if(adjMatrix[i][j])
            {
                stateVec[j]=1;
            }
        }
     }

   return(stateVec);
}
