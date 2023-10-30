
/*$Id: get_support_truth_tables.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int* decimal2binary(int decimal_value,int nBits);

int*** get_support_truth_tables(int n,int *nInputs)
{
	int i,j,k;
	int*** support_truth_tables=(int***)malloc(n*sizeof(int**));

	for (i = 0; i < n; ++i)
	{
		support_truth_tables[i]=(int**)malloc(pow(2,nInputs[i])*sizeof(int*));
		for (j = 0; j < pow(2,nInputs[i]); ++j)
		{
			support_truth_tables[i][j]= (int*)decimal2binary(j,nInputs[i]);
		}
	}
	return(support_truth_tables);
}

