/*$Id: printTruthTables.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>

void printTruthTables(int** truthTables,int* nBits, int nRows)
{
	int i,j;
	printf("-----------------------------\n");
		for (i = 0; i < nRows; i++)
		{
			for (j = 0; j < nBits[i]; ++j)
			{
				printf("%d \n",truthTables[i][j]);
			}
			printf("------------------------\n");
		}

}
