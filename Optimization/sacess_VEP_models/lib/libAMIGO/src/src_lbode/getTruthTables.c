/*$Id: getTruthTables.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <stdio.h>
#include <stdlib.h>


int* decimal2binary(int decimal_value,int nBits);

int** getTruthTables(int** adjMat,int** interMat,int** notMat,
		int* isState,int* nInputs,int *nBits,int nRows,int nCols)
{

	int i,j,k,m,n,counter1,input,miniterm;
	/*printInterMat(interMat,nRows,nCols);*/
	/*printAdjMat(adjMat,nRows);*/
	/*Tridimensional matrix of Size: states x nMiniTerms(state) x nInputs(state)*/
	int*** miniTerms;
	int* binary_value;
	int flag;
	int** truthTables=(int**)malloc(nRows*sizeof(int*));
	int** miniTermPositions=(int**)malloc(nRows*sizeof(int*));
	int** miniTermInputs=(int**)malloc(nRows*sizeof(int*));
	/*The number of miniTerms of each Row*/
	int* nMiniTerms=(int*)malloc(nRows*sizeof(int));

	/*Find the Inputs to a miniterm*/

	/*Find how many miniterms there per state*/
	for (i = 0; i < nRows; i++)
	{
		/*Find how many miniterms there are*/
		counter1=0;
		for (j = 0; j < nCols; j++){
			if(interMat[i][j]==1)counter1++;
		}
		nMiniTerms[i]=counter1;
		/*Now find their positions and store them in a matrix*/
		miniTermPositions[i]=(int*)malloc(counter1*sizeof(int));
		counter1=0;
		for (j = 0; j < nCols; j++){
			if(interMat[i][j]==1)
				miniTermPositions[i][counter1++]=j;
		}
	}
	/*Find the inputs to these miniterms. This info is
	the adjacency matrix*/
	for (j = 0; j < nRows; j++){

		counter1=0;
		miniTermInputs[j]=(int*)malloc(nInputs[j]*sizeof(int));
		for (i = 0; i < nRows; i++){
			if(adjMat[i][j]==1)
				miniTermInputs[j][counter1++]=i;
		}
	}

	/*Allocate memory for miniTerms*/
	miniTerms=(int***)malloc(nRows*sizeof(int**));
	for (i = 0; i < nRows; i++)
	{
		miniTerms[i]=(int**)malloc(nInputs[i]*sizeof(int*));

		for (j = 0; j < nInputs[i]; j++)
		{
			miniTerms[i][j]=(int*)malloc(nMiniTerms[i]*sizeof(int));
			for (k = 0; k < nMiniTerms[i]; k++)
			{
				miniTerms[i][j][k]=0;
			}

		}
	}


	/*These should be ordered according to the adjacency matrix*/
	/*We iterate through columns first(targets). See definition of adjacency*/
	/*Matrix in case of doubt*/
	for (i = 0; i < nRows; i++)
	{
		if(isState[i])
		{
			for (j = 0; j < nInputs[i];j++)
			{
				input=miniTermInputs[i][j];
				for (k = 0; k < nMiniTerms[i];k++)
				{
					miniterm=miniTermPositions[i][k];
					if(interMat[input][miniterm]==-1)
					{
						if(notMat[input][miniterm]==1)	miniTerms[i][j][k]=-1;
						else							 miniTerms[i][j][k]=1;
					}
					else miniTerms[i][j][k]=0;
				}
			}
		}
	}

	/*printNminiterms(miniTerms,nInputs,nMiniTerms,nRows);*/

	/*Note that we have the miniterms it is easier to parse to*/
	/*Truth tables*/
	for (i = 0; i < nRows; i++)
	{
		truthTables[i]=(int*)malloc(nBits[i]*sizeof(int));
		for (j = 0; j < nBits[i]; ++j) truthTables[i][j]=0;

		/*Check active inputs in the miniterm*/

		for (k = 0; k <nBits[i];k++)
		{
			for (j = 0; j < nMiniTerms[i]; ++j)
			{
				binary_value =  decimal2binary(k, nInputs[i]);
				flag=1;
				for (m = 0; m < nInputs[i]; m++)
				{
					if(miniTerms[i][m][j]==1)
					{
						if(binary_value[m]!=1) flag=0;
					}
					if(miniTerms[i][m][j]==-1)
					{
						if(binary_value[m]!=0) flag=0;
					}
				}
				free(binary_value);
				if(flag) truthTables[i][k]=1;
			}
		}
	}

	/*printTruthTables(truthTables,nBits,nRows);*/

	/*Release memory*/

	/*miniTerms*/

	for (i = 0; i < nRows; i++){
		for (j = 0; j < nInputs[i]; j++)
			free(miniTerms[i][j]);
		free(miniTerms[i]);}
	free(miniTerms);

	/*miniTermPositions*/
	for (i = 0; i < nRows; i++)
		free(miniTermPositions[i]);
	free(miniTermPositions);

	for (i = 0; i < nRows; i++)
		free(miniTermInputs[i]);
	free(miniTermInputs);

	free(nMiniTerms);

	return(truthTables);
}
