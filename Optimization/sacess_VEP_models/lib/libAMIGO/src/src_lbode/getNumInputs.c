/*$Id: getNumInputs.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <stdio.h>
#include <stdlib.h>

int* getNumInputs(int **adjMatrix,int n)
{
    int* numInputs=(int*)malloc(n*sizeof(int));
    int i,j,count;
    for (j = 0; j <n; j++)
    {
        count=0;
        for(i=0;i<n;i++)
        {
            if(adjMatrix[i][j])
            {
                count++;
            }
        }
        numInputs[j]=count;
     }
    return(numInputs);
}
