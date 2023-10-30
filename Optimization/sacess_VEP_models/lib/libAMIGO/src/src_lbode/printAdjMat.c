
 /*$Id: printAdjMat.c 1372 2012-01-24 11:16:50Z cokelaer $ */
 
#include <stdio.h>
#include <stdlib.h>

void printAdjMat(int** adjMat,int nNodes)
{
	int i,j;

	for (i = 0; i < nNodes; i++)
	{
		for (j = 0; j < nNodes; j++)
		{
			printf("%d\t",adjMat[i][j]);
		}
		printf(";\n");
	}
}
