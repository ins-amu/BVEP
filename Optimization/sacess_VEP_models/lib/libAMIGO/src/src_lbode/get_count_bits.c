
/*$Id: get_count_bits.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int* get_count_bits(int n,int** truth_tables, int* numBits)
{
    int *count_bits=(int*)malloc(n*sizeof(int));
    int i,j;
    for (i = 0; i < n; ++i)
    {
    	count_bits[i]=0;
    	for (j = 0; j < numBits[i]; ++j)
    	{
    		if(truth_tables[i][j])count_bits[i]++;
    	}
	}
    return(count_bits);
}
