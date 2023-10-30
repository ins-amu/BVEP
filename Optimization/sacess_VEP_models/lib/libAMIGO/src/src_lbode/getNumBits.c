/*$Id: getNumBits.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int* getNumBits(int* numInputs,int n)
{
    int *numBits=(int*)malloc(n*sizeof(int));
    int i;
    for (i = 0; i <n; i++)
    {
        if(numInputs[i]>0)
        {
            numBits[i]=(int)pow(2,(double)numInputs[i]);
        }
        else
        {
            numBits[i]=1;
        }
    }
    return(numBits);
}
