/*$Id: get_truth_tables_index.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int** get_truth_tables_index(int n,int** truth_tables, int* numBits,int* count_bits)
{
    int **truth_tables_index=(int**)malloc(n*sizeof(int*));
    int i,j,count;
    for (i = 0; i < n; ++i)
    {
    	count=0;
    	truth_tables_index[i]=(int*)malloc(count_bits[i]*sizeof(int));
    	for (j = 0; j < numBits[i]; ++j)
    	{
    		if(truth_tables[i][j])truth_tables_index[i][count++]=j;
    	}
	}
    return(truth_tables_index);
}
