/*$Id: rhsODE.c 1372 2012-01-24 11:16:50Z cokelaer $*/


#include <AMIGO_model.h>
#include "CNOStructure.h"
#include "mex.h"

#define Ith(v,i) ( NV_DATA_S(v)[i] )

int lbode_AMIGO_rhs(realtype t, N_Vector y, N_Vector ydot, void *data)
{
	int i,j,k;
	AMIGO_model* amigo_model=(AMIGO_model*) data;
	CNOStructure* cno=(CNOStructure*) amigo_model->data;

	int countPar=0;
	double tempProd;
	double kHill,nHill;
	int countState=0;

	/*Loop through every column j in the Graph adjacency matrix*/
	for (j = 0; j <cno->nRows; j++){

		if(cno->isState[j]){

			Ith(ydot,countState)=0;

			for (i = 0; i < cno->numInputs[j]; ++i){

				nHill=amigo_model->pars[countPar++];
				kHill=amigo_model->pars[countPar++];

				if(cno->isState[cno->input_index[j][i]]){
					cno->hillFuncValues[i]=
						cno->transfer_function(
						Ith(y,cno->state_index[cno->input_index[j][i]]),
						nHill,
						kHill
						);
					
				}
				else{
					cno->hillFuncValues[i]=
						cno->transfer_function(
						cno->state_array[amigo_model->exp_num][cno->input_index[j][i]],
						nHill,
						kHill
						);
				}
			}

			/*For every bit in the truth table*/
			for (i = 0; i < cno->count_bits[j]; ++i){

				tempProd=1;
				for (k = 0; k < cno->numInputs[j]; k++){

					if(!cno->support_truth_tables[j][cno->truth_tables_index[j][i]][k]){

						tempProd*=(1-cno->hillFuncValues[k]);
					}
					else{
						tempProd*=cno->hillFuncValues[k];
					}
				}

				Ith(ydot,countState)+=tempProd;
			}
			Ith(ydot,countState)=
				(Ith(ydot,countState)-Ith(y,countState))
											*amigo_model->pars[countPar++]
												*(1-cno->inhibitor_array[amigo_model->exp_num][j]);

			countState++;
		}
	}

	return(0);
}

