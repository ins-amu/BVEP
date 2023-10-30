
/*$Id: getAdjacencyMatrix.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <stdio.h>
#include <stdlib.h>

int** getAdjacencyMatrix(int** interMat,int nRows, int nCols)
{
	int i,j,k;

	int** adjMat= (int**) malloc(nRows * sizeof(int*));
	for (i = 0; i < nRows; i++)
	{
		adjMat[i] = (int*) malloc(nCols * sizeof(int));
		for (j = 0; j < nRows; j++)
		{
			adjMat[i][j]=0;
		}
	}

	for (i = 0; i < nRows; i++)
	{
		  for (j = 0; j < nCols; j++)
		  {
			 if(interMat[i][j]==1)
			 {
				 for (k = 0; k < nRows; k++)
				 {
					 if(interMat[k][j]==-1)
						 adjMat[k][i]=1;
				 }
			 }
		  }
	}

	/*printAdjMat(adjMat,nRows);*/

	return(adjMat);
}

