/*$Id: printInterMat.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>

void printInterMat(int** interMat,int nRows,int nCols)
{
	int i,j;
	for (i = 0; i < nRows; ++i)
	{
		for (j = 0; j < nCols; ++j)
		{
			printf("%d\t",interMat[i][j]);
		}
		printf("\n");
	}

}
