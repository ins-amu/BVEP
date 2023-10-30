/*$Id: printNminiTerms.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>

void printNminiterms(int*** miniTerms,int* nInputs, int* nMiniterms, int nRows)
{
	int i,j,k;
	printf("\n");
	for (i = 0; i < nRows; i++)
	{
		printf("Number of miniterms:%d\n",nMiniterms[i]);
		printf("Number of n inputs:%d\n",nInputs[i]);
		printf("Species %d\n",i);
		for (j = 0; j <nInputs[i]; j++)
		{
			for (k = 0; k < nMiniterms[i]; ++k)
			{
				printf("%d\t",miniTerms[i][j][k]);
			}
			printf("\n");
		}

	}
}
