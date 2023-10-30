/*$Id: getStateIndex.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <stdio.h>
#include <stdlib.h>

int *getStateIndex(int **adjMatrix, int n)
{
    int* indexVec=(int*)malloc(n*sizeof(int));
    int i,j;
    int stateNumber=0;
    int first;
    for (j = 0; j <n; j++)
    {
    	indexVec[j]=-1;
    	first=1;
        for(i=0;i<n;i++)
        {
            if(adjMatrix[i][j] && first)
            {
            	indexVec[j]=stateNumber++;
            	first=0;
            }
        }
     }
    /*
    printf("State Index\n");
    for (j = 0; j <n; j++)
    {
    	Rprintf("species %d \t\t index %d\n",j,indexVec[j]);
    }
    */

   return(indexVec);
}
