/*$Id: sim_logic_ode.c 1433 2012-04-04 13:24:36Z davidh $*/

#include <stdio.h>
#include <stdlib.h>
#include "CNOStructure.h"
#include <AMIGO_problem.h>

#include <math.h>
#include <string.h>
#include "mex.h"
#include "mat.h"


int lbode_AMIGO_rhs(realtype t, N_Vector y, N_Vector ydot, void *data);
 
AMIGO_model* mxAllocateAMIGOmodel(mxArray* privstruct_ptr,mxArray* inputs_ptr, int exp_num);

void mxSolveAMIGOivp(AMIGO_problem* amigo_problem, int save2workspace, char* save2File);

void mxSolveAMIGO_FSA(AMIGO_problem* amigo_problem, int save2Workspace, char* save2File);

void mxSolveAMIGO_MKLSENS(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File);

int mxAMIGOsave2File(mxArray *outputs_ptr,char* save2File);

void mxSolveAMIGOnl2sol(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File);

void mxEvalAMIGOlsq(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File);

void mxEvalAMIGOllk(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	int countProtected=0;
	int i,n_exp,iexp, j, k, b, counter;
	int *indexStim;
	int *indexInh;
	double* timeSig;
	int **interMAT;
	int **notMAT;
	int** adjMatrix;
	double **valueINHIBITORS;
	double **valueSTIMULI;
	CNOStructure* cno;
	double** state_array;

	double** inhibitor_array;
	int maxNumInputs=-1;
	int sensi=1;
	double**** sensResults;
	double unknown_ICs;

	int dims[4];

	double *outputData;
	double *sensData;
	mxArray *pa;
	mxArray* privstruct_ptr;
	mxArray*inputs_ptr;
	char *buf;
	mwSize buflen;
	int status;

	AMIGO_problem* amigo_problem;
	AMIGO_model** amigo_models;

	int nRows =(int)mxGetPr(prhs[0])[0];

	int nCols=(int)mxGetPr(prhs[1])[0];

	int nStimuli=(int)mxGetPr(prhs[2])[0];

	int nInhibitors=(int)mxGetPr(prhs[3])[0];

	int nExperiments=(int)mxGetPr(prhs[4])[0];

	int transfer_function=(int)mxGetPr(prhs[5])[0];

	cno=(CNOStructure*)malloc(sizeof(CNOStructure));

	counter=0;
	indexInh=(int*)malloc(nInhibitors*sizeof(int));
	for (i = 0; i < nInhibitors; i++) {
		indexInh[i] =(int) mxGetPr(prhs[6])[counter++]-1;
	}

	counter=0;
	indexStim=(int*)malloc(nStimuli*sizeof(int));
	for (i = 0; i < nStimuli; i++) {
		indexStim[i] =(int) mxGetPr(prhs[7])[counter++]-1;
	}

	interMAT = (int**) malloc(nRows * sizeof(int*));
	for (i = 0; i < nRows; i++)interMAT[i] = (int*) malloc(nCols * sizeof(int));

	counter=0;
	for (i = 0; i < nCols; i++){
		for (j = 0; j < nRows; j++)interMAT[j][i] =(int) mxGetPr(prhs[8])[counter++];
	}

	notMAT = (int**)malloc(nRows * sizeof(int*));
	for (i = 0; i < nRows; i++)notMAT[i] = (int*)malloc(nCols*sizeof(int));

	counter=0;
	for (i = 0; i < nCols; i++){
		for (j = 0; j < nRows; j++)notMAT[j][i]=(int)mxGetPr(prhs[9])[counter++];
	}

	adjMatrix = (int**)malloc(nRows * sizeof(int*));
	for (i = 0; i < nRows; i++)adjMatrix[i] = (int*)malloc(nRows*sizeof(int));

	counter=0;
	for (i = 0; i < nRows; i++){
		for (j = 0; j < nRows; j++)adjMatrix[j][i]=(int)mxGetPr(prhs[10])[counter++];
	}

	counter=0;
	valueINHIBITORS = (double**)malloc(nExperiments * sizeof(double*));
	for (i = 0; i < nExperiments; i++)
		valueINHIBITORS[i] = (double*)malloc(nInhibitors*sizeof(double));

	counter=0;
	for (i = 0; i < nInhibitors; i++){
		for (j = 0; j < nExperiments; j++)valueINHIBITORS[j][i]=mxGetPr(prhs[11])[counter++];
	}

	counter=0;
	valueSTIMULI = (double**)malloc(nExperiments * sizeof(double*));
	for (i = 0; i < nExperiments; i++)
		valueSTIMULI[i] = (double*)malloc(nStimuli*sizeof(double));

	counter=0;
	for (i = 0; i < nStimuli; i++){
		for (j = 0; j < nExperiments; j++)valueSTIMULI[j][i]=mxGetPr(prhs[12])[counter++];
	}

	buflen = mxGetN(prhs[13])*sizeof(mxChar)+1;
	buf = mxMalloc(buflen);
	status = mxGetString(prhs[13], buf, buflen);

	cno->interMat=interMAT;
	cno->notMat=notMAT;
	cno->valueInhibitors=valueINHIBITORS;
	cno->valueStimuli=valueSTIMULI;    
	cno->indexStimuli=indexStim;
	cno->indexInhibitors=indexInh;

	cno->nRows=nRows;
	cno->nCols=nCols;
	cno->nStimuli=nStimuli;
	cno->nInhibitors=nInhibitors;

	cno->adjacencyMatrix=adjMatrix;
	cno->numInputs =(int*) getNumInputs(cno->adjacencyMatrix, cno->nRows);

	for (i = 0; i < nRows; ++i) {
		if(cno->numInputs[i]>maxNumInputs)maxNumInputs=cno->numInputs[i];

	}
	cno->maxNumInputs=maxNumInputs;

	cno->numBits =(int*) getNumBits(cno->numInputs, cno->nRows);
	cno->isState =(int*) findStates(cno->adjacencyMatrix, cno->nRows);

	cno->truthTables =(int**) getTruthTables(cno->adjacencyMatrix, cno->interMat,
		cno->notMat, cno->isState, cno->numInputs, cno->numBits, cno->nRows, cno->nCols);

	state_array = (double**)malloc(nExperiments*sizeof(double*));
	inhibitor_array = (double**)malloc(nExperiments*sizeof(double*));
	for (i = 0; i < nExperiments; ++i){
		state_array[i] = (double*)malloc(cno->nRows*sizeof(double));
		inhibitor_array[i]=(double*)malloc((cno->nRows)*sizeof(double));
	}

	cno->state_index=(int*)getStateIndex(cno->adjacencyMatrix, cno->nRows);
	cno->inhibitor_array=inhibitor_array;
	cno->state_array=state_array;

	cno->count_bits=(int*)get_count_bits(nRows, cno->truthTables, cno->numBits);

	cno->truth_tables_index=(int**)get_truth_tables_index(nRows,
		cno->truthTables, cno->numBits, cno->count_bits);

	cno->input_index=(int**)get_input_index(cno->adjacencyMatrix, nRows, cno->numInputs);

	counter=0;
	for (i = 0; i < nRows; ++i)
		if(cno->isState[i]) counter++;

	cno->nStates = counter;

	cno->support_truth_tables=(int***)get_support_truth_tables(nRows, cno->numInputs);

	if(transfer_function==1)
		cno->transfer_function = &linear_transfer_function;
	else if(transfer_function==2)
		cno->transfer_function = &hill_function;
	else
		cno->transfer_function = &normHill;

	for(iexp=0; iexp<nExperiments; iexp++){

		for(i=0; i<cno->nRows; i++){

			cno->state_array[iexp][i] = 1;
			cno->inhibitor_array[iexp][i] = 0;
		}

		for(i=0; i<cno->nRows; ++i)
		{
			if(cno->isState[i]){
				for (j = 0; j < cno->nInhibitors; j++){

					if(cno->indexInhibitors[j]==i){

						cno->inhibitor_array[iexp][i] = cno->valueInhibitors[iexp][j];
					}
				}
			}
			else
			{
				for (j = 0; j < cno->nStimuli; j++)
				{
					if(cno->indexStimuli[j]==i){
						cno->state_array[iexp][i] = cno->valueStimuli[iexp][j];
					}
				}
			}
		}
	}

	cno->hillFuncValues=(double*)malloc(sizeof(double)*cno->maxNumInputs);

	privstruct_ptr=mexGetVariable("caller", "privstruct");
	inputs_ptr=mexGetVariable("caller", "inputs");

	pa=mxGetField(privstruct_ptr,0, "n_exp");
	n_exp=(int)mxGetPr(pa)[0];

	amigo_models=(AMIGO_model**)malloc(sizeof(AMIGO_model*)*n_exp);

	for (i = 0; i < n_exp; i++){
		amigo_models[i]=mxAllocateAMIGOmodel(privstruct_ptr,inputs_ptr,i);
		amigo_models[i]->rhs=lbode_AMIGO_rhs;
		amigo_models[i]->data=cno;
	}

	amigo_problem=allocate_AMIGO_problem(n_exp,amigo_models);
	
	if (strcmp(buf, "IVP")==0){
		mxSolveAMIGOivp(amigo_problem,1,"");
	}else if(strcmp(buf, "FSA")==0){
		//arg1 amigo_prob, 2 use cvodes, 3 use mkl, arg4 save2Workspace arg5 filename
		mxSolveAMIGO_FSA(amigo_problem,1,"");
	}else if(strcmp(buf, "MKLSENS")==0){
#ifdef MKL
		//arg1 amigo_prob, 2 use cvodes, 3 use mkl, arg4 save2Workspace arg5 filename
		mxSolveAMIGO_MKLSENS(amigo_problem,1,"");
#else
		mexErrMsgTxt("Compile using MKL library and try again.");
#endif
	}else if(strcmp(buf, "NL2SOL")==0){

		mxSolveAMIGOnl2sol(amigo_problem,1,"");

	}else if(strcmp(buf, "LSQ")==0){

		mxEvalAMIGOlsq(amigo_problem,1,"");

	}else if(strcmp(buf, "LLK")==0){

		mxEvalAMIGOllk(amigo_problem,1,"");

	}

	free_AMIGO_problem(amigo_problem);

	for (i = 0; i < nRows; ++i) {
		for (j = 0; j < pow((int)2, cno->numInputs[i]); ++j) {
			free(cno->support_truth_tables[i][j]);
		}
		free(cno->support_truth_tables[i]);
	}
	free(cno->support_truth_tables);

	free(indexStim);
	free(indexInh);

	for (i = 0; i < nRows; i++) {
		free(cno->adjacencyMatrix[i]);
		free(cno->truthTables[i]);
		free(cno->truth_tables_index[i]);
		free(cno->input_index[i]);
	}
	free(cno->truthTables);
	free(cno->adjacencyMatrix);
	free(cno->truth_tables_index);
	free(cno->input_index);

	for (i = 0; i < nExperiments; i++)
		free(valueSTIMULI[i]);
	free(valueSTIMULI);

	for (i = 0; i < nExperiments; i++)
		free(valueINHIBITORS[i]);
	free(valueINHIBITORS);

	for (i = 0; i < nRows; i++)
		free(notMAT[i]);
	free(notMAT);

	for (i = 0; i < nRows; i++)
		free(interMAT[i]);
	free(interMAT);

	for (i = 0; i < nExperiments; ++i){
		free(state_array[i]);
		free(inhibitor_array[i]);
	}
	free(state_array);
	free(inhibitor_array);

	free(cno->state_index);
	free(cno->numBits);
	free(cno->numInputs);
	free(cno->count_bits);
	free(cno->hillFuncValues);
	mxFree(buf);

}
AMIGO_model* mxAllocateAMIGOmodel(mxArray* privstruct_ptr,mxArray* inputs_ptr, int exp_num){

	mxArray *pa;
	mxArray *cell_ptr;
	const mwSize  *dim;
	AMIGO_model* amigo_model;

	int n_states,n_observables,n_pars,n_opt_pars,n_times,n_opt_ics,n_controls,n_controls_t,i,j,counter;

	//inputs.model.n_st
	pa=mxGetField(inputs_ptr,0, "model");
	pa=mxGetField(pa,0, "n_st");
	n_states=mxGetPr(pa)[0];

	//N_OBS
	pa=mxGetField(privstruct_ptr,0, "n_obs");
	cell_ptr = mxGetCell(pa, exp_num);
	n_observables=(int)mxGetPr(cell_ptr)[0];

	//n pars
	pa=mxGetField(inputs_ptr,0, "model");
	pa=mxGetField(pa,0, "n_par");
	n_pars=(int)mxGetPr(pa)[0];

	//inputs.PEsol.n_theta
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "n_theta");
	n_opt_pars=(int)mxGetPr(pa)[0];

	// privstruct.n_s
	pa=mxGetField(privstruct_ptr,0, "n_s");
	cell_ptr = mxGetCell(pa, exp_num);
	n_times=(int)mxGetPr(cell_ptr)[0];
	//printf("n_times=%d\n",n_times);

	//inputs.PEsol.n_theta_y0
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "n_local_theta_y0");
	cell_ptr = mxGetCell(pa, exp_num);
	n_opt_ics=(int)mxGetPr(cell_ptr)[0];
	//printf("n_opt_ics=%d\n",n_opt_ics);

	//inputs.model.n_stimulus
	pa=mxGetField(inputs_ptr,0, "model");
	pa=mxGetField(pa,0, "n_stimulus");
	n_controls=(int)mxGetPr(pa)[0];
	//printf("n_controls=%d\n",n_controls);

	//privstruct.n_s
	pa=mxGetField(privstruct_ptr,0, "t_con");
	cell_ptr = mxGetCell(pa, exp_num);
	dim=mxGetDimensions(cell_ptr);
	n_controls_t=(int)dim[1];
	//printf("n_controls_t=%d\n",n_controls_t);

	amigo_model=allocate_AMIGO_model(n_states,n_observables,n_pars,
		n_opt_pars,n_times,n_opt_ics,n_controls, n_controls_t,exp_num);


	//privstruct.index_observables
	pa=mxGetField(privstruct_ptr,0, "index_observables");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_observables; i++){
		amigo_model->index_observables[i]=(int)mxGetPr(cell_ptr)[i]-1;
		//mexPrintf("%d\n",amigo_model->index_observables[i]);
	}

	//Simulation Pars
	pa=mxGetField(inputs_ptr,0, "model");
	pa=mxGetField(pa,0, "par");
	for (i = 0; i < n_pars; i++){
		amigo_model->pars[i]=(double)mxGetPr(pa)[i];
		//mexPrintf("%f\n",amigo_model->pars[i]);
	}

	//inputs.PEsol.index_global_theta
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "index_global_theta");
	for (i = 0; i < n_opt_pars; i++){
		amigo_model->index_opt_pars[i]=(int)mxGetPr(pa)[i] - 1;
		//mexPrintf("%f\n",amigo_model->index_opt_pars[i]);
	}


	//inputs.PEsol.global_theta_guess
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "global_theta_guess");
	for (i = 0; i < n_opt_pars; i++){
		amigo_model->pars_guess[i]=(double)mxGetPr(pa)[i];
		//mexPrintf("%f\n",amigo_model->pars_guess[i]);
	}

	//inputs.PEsol.global_theta_min
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "global_theta_min");
	for (i = 0; i < n_opt_pars; i++){
		amigo_model->pars_LB[i]=(double)mxGetPr(pa)[i];
		//mexPrintf("%f\n",amigo_model->pars_LB[i]);
	}

	//inputs.PEsol.global_theta_max
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0, "global_theta_max");
	for (i = 0; i < n_opt_pars; i++){
		amigo_model->pars_UB[i]=(double)mxGetPr(pa)[i];
		//mexPrintf("%f\n",amigo_model->pars_UB[i]);
	}

	//Simulation times
	//inputs.exps.t_int{iexp}
	pa=mxGetField(privstruct_ptr,0, "t_in");
	cell_ptr = mxGetCell(pa, exp_num);
	amigo_model->t0=(double)mxGetPr(cell_ptr)[0];
	//mexPrintf("%f\n",amigo_model->t0);

	//privstruct.t_f{iexp}
	pa=mxGetField(privstruct_ptr,0, "t_f");
	cell_ptr = mxGetCell(pa, exp_num);
	amigo_model->tf=(double)mxGetPr(cell_ptr)[0];
	//mexPrintf("%f\n",amigo_model->tf);


	//privstruct.t_int{iexp}
	pa=mxGetField(privstruct_ptr,0, "t_int");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_times; i++){
		amigo_model->t[i]=(double)mxGetPr(cell_ptr)[i];
		//   mexPrintf("%f\n",amigo_model->t[i]);
	}


	//Initial conditions
	//inputs.exps.exp_y0
	pa=mxGetField(inputs_ptr,0, "exps");
	pa=mxGetField(pa,0, "exp_y0");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_states; i++){
		amigo_model->y0[i]=(double)mxGetPr(cell_ptr)[i];
		//mexPrintf("%f\n",amigo_model->y0[i]);
	}
	//inputs.PEsol.index_local_theta_y0{1}
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0,"index_local_theta_y0");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_opt_ics; i++){
		amigo_model->index_opt_ics[i]=(int)mxGetPr(cell_ptr)[i]-1;
		//mexPrintf("%d\n",amigo_model->index_opt_ics[i]);
	}

	//inputs.PEsol.local_theta_y0_guess{1,2,...}
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0,"local_theta_y0_guess");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_opt_ics; i++){
		amigo_model->y0_guess[i]=(double)mxGetPr(cell_ptr)[i];
		//mexPrintf("%f\n",amigo_model->y0_guess[i]);
	}


	//inputs.PEsol.local_theta_y0_min{1,2,...}
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0,"local_theta_y0_min");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_opt_ics; i++){
		amigo_model->y0_LB[i]=(double)mxGetPr(cell_ptr)[i];
		//mexPrintf("%f\n",amigo_model->y0_LB[i]);
	}


	//inputs.PEsol.local_theta_y0_max{1,2,...}
	pa=mxGetField(inputs_ptr,0, "PEsol");
	pa=mxGetField(pa,0,"local_theta_y0_max");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_opt_ics; i++){
		amigo_model->y0_UB[i]=(double)mxGetPr(cell_ptr)[i];
		//mexPrintf("%f\n",amigo_model->y0_UB[i]);
	}

	//Controls
	//privstruct.t_con{iexp}
	pa=mxGetField(privstruct_ptr,0, "t_con");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i <n_controls_t; i++){
		amigo_model->controls_t[i]=(double)mxGetPr(cell_ptr)[i];
		//mexPrintf("%f\n",amigo_model->controls_t[i]);
	}

	//privstruct.u{iexp}
	counter=0;
	pa=mxGetField(privstruct_ptr,0, "u");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_controls; i++) {
		for (j= 0; j < n_controls_t-1; j++){
			// mexPrintf("%d %d\n",i,j);
			amigo_model->controls_v[i][j]=(double)mxGetPr(cell_ptr)[counter++];
			//mexPrintf("%f\t",amigo_model->controls_v[i][j]);
		}
		//mexPrintf("\n");
	}

	//inputs.exps.pend{iexp}
	counter=0;
	pa=mxGetField(privstruct_ptr,0, "pend");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i= 0; i < amigo_model->n_controls; i++){
		for (j= 0; j < amigo_model->n_controls_t-1; j++){
			amigo_model->slope[i][j]=(double)mxGetPr(cell_ptr)[counter++];
			//    mexPrintf("%f\t",amigo_model->slope[i][j]);
		}
		//mexPrintf("\n");
	}

	//Storing matrixes
	counter=0;
	pa=mxGetField(privstruct_ptr,0, "exp_data");
	cell_ptr = mxGetCell(pa, exp_num);
	if(mxGetNumberOfElements(cell_ptr)>0){
		//mexPrintf("\n%d\n",mxGetNumberOfElements(cell_ptr));
		for (i = 0; i < n_observables; i++){
			for (j = 0; j < n_times; j++){
				amigo_model->exp_data[i][j]=(double)mxGetPr(cell_ptr)[counter++];
				//    mexPrintf("%f\t",amigo_model->exp_data[i][j]);
			}
			//  mexPrintf("\n");
		}
		//mexPrintf("\n");

		counter=0;
		pa=mxGetField(inputs_ptr,0, "exps");
		pa=mxGetField(pa,0, "Q");
		cell_ptr = mxGetCell(pa, exp_num);
		for (i = 0; i < n_observables; i++){
			for (j = 0; j < n_times; j++){
				amigo_model->Q[i][j]=(double)mxGetPr(cell_ptr)[counter++];
				//    mexPrintf("%f\t",amigo_model->Q[i][j]);
			}
			//mexPrintf("\n");
		}
		//mexPrintf("\n");
	}

	counter=0;
	pa=mxGetField(privstruct_ptr,0, "w_obs");
	cell_ptr = mxGetCell(pa, exp_num);
	for (i = 0; i < n_observables; i++){}
	//amigo_model->w_obs[i]=(double)mxGetPr(cell_ptr)[counter++];

	//Simulation Related Parameter
	//inputs.ivpsol.rtol
	pa=mxGetField(inputs_ptr,0, "ivpsol");
	if(mxGetFieldNumber(pa, "rtol") != - 1){
		pa=mxGetField(pa,0, "rtol");
		amigo_model->reltol=(double)mxGetPr(pa)[0];
	}else amigo_model->reltol=1e-6;
	//mexPrintf("%e\n",amigo_model->reltol);

	//inputs.ivpsol.atol
	pa=mxGetField(inputs_ptr,0, "ivpsol");
	if(mxGetFieldNumber(pa, "atol") != - 1){
		pa=mxGetField(pa,0, "atol");
		amigo_model->atol=(double)mxGetPr(pa)[0];
	}else amigo_model->atol=1e-6;
	//mexPrintf("%e\n",amigo_model->atol);

	//privstruct.ivpsol.max_step_size
	pa=mxGetField(inputs_ptr,0, "ivpsol");
	if(mxGetFieldNumber(pa, "max_step_size") != - 1){
		pa=mxGetField(pa,0, "max_step_size");
		amigo_model->max_step_size=(double)mxGetPr(pa)[0];
	}else amigo_model->max_num_steps=DBL_MAX;
	//mexPrintf("%f\n",amigo_model->max_step_size);

	//privstruct.ivpsol.max_num_steps
	pa=mxGetField(inputs_ptr,0, "ivpsol");
	if(mxGetFieldNumber(pa, "max_num_steps") != - 1){
		pa=mxGetField(pa,0, "max_num_steps");
		amigo_model->max_num_steps=(int)mxGetPr(pa)[0];
	}else amigo_model->max_num_steps=10000;
	//mexPrintf("%d\n",amigo_model->max_num_steps);


	//privstruct.ivpsol.max_error_test_fails
	pa=mxGetField(inputs_ptr,0, "ivpsol");
	if(mxGetFieldNumber(pa, "max_error_test_fails") != - 1){
		pa=mxGetField(pa,0, "max_error_test_fails");
		amigo_model->max_error_test_fails=(int)mxGetPr(pa)[0];
	}else amigo_model->max_error_test_fails=50;
	//mexPrintf("%d\n",amigo_model->max_error_test_fails);

	pa=mxGetField(inputs_ptr,0, "PEsol");
	if(mxGetFieldNumber(pa, "mkl_tol") != - 1){
		pa=mxGetField(pa,0, "mkl_tol");
		amigo_model->mkl_tol=(double)mxGetPr(pa)[0];
	}


#ifdef AMIGO_JAC
	amigo_model->jac=amigoJAC;
	amigo_model->use_jacobian=1;
#else
	amigo_model->use_jacobian=0;
#endif

return(amigo_model);

}

void mxSolveAMIGOivp(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File){

	mxArray *outputs_ptr,*stats_ptr,*stats_value,*sim_data,*obs_data, *sim_cell,*stats_cell,*flag_cell,*observables_cell,*flag_value;
	const char *field_names[] = {"success", "simulation","observables","sim_stats"};
	const char *stats_names[] = {"flag","nst","nfe","nsetups","netf","nni","ncfn"};
	mwSize dims[2]={1,1};
	mwSize dims_cell[1];
	//MATFile *pmat;
	int i,j,k,flag,counter;

	outputs_ptr=mxCreateStructArray(1, dims,4,field_names);

	dims_cell[0]=amigo_problem->n_models;

	sim_cell=mxCreateCellArray(1,dims_cell);
	flag_cell=mxCreateCellArray(1,dims_cell);
	observables_cell=mxCreateCellArray(1,dims_cell);
	stats_cell=mxCreateCellArray(1,dims_cell);

	mxSetField(outputs_ptr, 0,"simulation", sim_cell);
	mxSetField(outputs_ptr, 0,"success", flag_cell);
	mxSetField(outputs_ptr, 0,"observables", observables_cell);
	mxSetField(outputs_ptr, 0,"sim_stats", stats_cell);


	for (i = 0; i < amigo_problem->n_models; i++){

		flag=simulate_AMIGO_model_observables(amigo_problem->amigo_models[i],0);

		flag_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(flag_value)[0]=(double)flag;

		sim_data=mxCreateDoubleMatrix(
			amigo_problem->amigo_models[i]->n_times ,
			amigo_problem->amigo_models[i]->n_states ,
			mxREAL
			);

		obs_data=mxCreateDoubleMatrix(
			amigo_problem->amigo_models[i]->n_times ,
			amigo_problem->amigo_models[i]->n_observables,
			mxREAL
			);

		stats_ptr=mxCreateStructArray(1, dims,7,stats_names);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->error_flag;
		mxSetField(stats_ptr, 0,"flag", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nst;
		mxSetField(stats_ptr, 0,"nst", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfe;
		mxSetField(stats_ptr, 0,"nfe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nsetups;
		mxSetField(stats_ptr, 0,"nsetups", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->netf;
		mxSetField(stats_ptr, 0,"netf", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nni;
		mxSetField(stats_ptr, 0,"nni", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->ncfn;
		mxSetField(stats_ptr, 0,"ncfn", stats_value);

		mxSetCell(stats_cell,i, stats_ptr);
		mxSetCell(flag_cell,i, flag_value);
		mxSetCell(sim_cell, i, sim_data);
		mxSetCell(observables_cell, i, obs_data);

		counter=0;
		for (j = 0; j < amigo_problem->amigo_models[i]->n_states; j++){
			for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
				mxGetPr(sim_data)[counter++]=amigo_problem->amigo_models[i]->sim_results[j][k];
			}
		}

		counter=0;
		for (j = 0; j < amigo_problem->amigo_models[i]->n_observables; j++){
			for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
				mxGetPr(obs_data)[counter++]=amigo_problem->amigo_models[i]->obs_results[j][k];
			}
		}
	}

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace)
		mexPutVariable("caller", "outputs", outputs_ptr);
	else mxDestroyArray(outputs_ptr);
}

void mxSolveAMIGO_FSA(AMIGO_problem* amigo_problem, int save2Workspace, char* save2File){

	mxArray *outputs_ptr,*sens_data,*stats_ptr,*stats_cell,*stats_value,*sim_data,*sens_cell, *sim_cell,*flag_cell,*flag_value;
	const char *field_names[] = {"flag", "simulation","sensitivities","cvodes_stats"};
	const char *stats_names[] = {"flag","nst","nfe","nsetups","netf","nni","ncfn"
		,"nfSe","nfeS","nsetupsS","netfS","nniS","ncfnS"};

	mwSize dims[2]={1,1};
	mwSize dims_sens[3];
	mwSize dims_cell[1];
	MATFile *pmat;
	int i,j,k,m,flag,counter;

	outputs_ptr=mxCreateStructArray(1, dims,4,field_names);

	dims_cell[0]=amigo_problem->n_models;

	flag_cell=mxCreateCellArray(1,dims_cell);
	sim_cell=mxCreateCellArray(1,dims_cell);
	sens_cell=mxCreateCellArray(1,dims_cell);
	stats_cell=mxCreateCellArray(1,dims_cell);

	mxSetField(outputs_ptr, 0,"flag", flag_cell);
	mxSetField(outputs_ptr, 0,"simulation", sim_cell);
	mxSetField(outputs_ptr, 0,"sensitivities", sens_cell);
	mxSetField(outputs_ptr, 0,"cvodes_stats", stats_cell);

	for (i = 0; i < amigo_problem->n_models; i++){

		flag=get_AMIGO_model_sens(amigo_problem->amigo_models[i],1,0);

		flag_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(flag_value)[0]=(double)flag;

		dims_sens[0]=amigo_problem->amigo_models[i]->n_times;
		dims_sens[1]=amigo_problem->amigo_models[i]->n_states;
		dims_sens[2]=amigo_problem->amigo_models[i]->n_total_x;

		sens_data=mxCreateNumericArray(3,dims_sens,mxDOUBLE_CLASS, mxREAL);

		sim_data=mxCreateDoubleMatrix(
			amigo_problem->amigo_models[i]->n_times ,
			amigo_problem->amigo_models[i]->n_states ,
			mxREAL
			);

		counter=0;
		for (j = 0; j < amigo_problem->amigo_models[i]->n_states; j++){
			for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
				mxGetPr(sim_data)[counter++]=amigo_problem->amigo_models[i]->sim_results[j][k];
			}
		}

		counter=0;
		for (m = 0; m < amigo_problem->amigo_models[i]->n_total_x; m++){
			for (j = 0; j < amigo_problem->amigo_models[i]->n_states; j++){
				for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
					mxGetPr(sens_data)[counter++]=amigo_problem->amigo_models[i]->sens_results[j][m][k];
				}
			}
		}

		stats_ptr=mxCreateStructArray(1, dims,13,stats_names);
		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->error_flag;
		mxSetField(stats_ptr, 0,"flag", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nst;
		mxSetField(stats_ptr, 0,"nst", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfe;
		mxSetField(stats_ptr, 0,"nfe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nsetups;
		mxSetField(stats_ptr, 0,"nsetups", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->netf;
		mxSetField(stats_ptr, 0,"netf", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nni;
		mxSetField(stats_ptr, 0,"nni", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->ncfn;
		mxSetField(stats_ptr, 0,"ncfn", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfSe;
		mxSetField(stats_ptr, 0,"nfSe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfeS;
		mxSetField(stats_ptr, 0,"nfeS", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->netfS;
		mxSetField(stats_ptr, 0,"netfS", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nniS;
		mxSetField(stats_ptr, 0,"nniS", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->ncfnS;
		mxSetField(stats_ptr, 0,"ncfnS", stats_value);

		mxSetCell(stats_cell,i, stats_ptr);
		mxSetCell(flag_cell,i, flag_value);
		mxSetCell(sim_cell, i, sim_data);
		mxSetCell(sens_cell, i, sens_data);
	}

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace)
		mexPutVariable("caller", "outputs", outputs_ptr);
	else mxDestroyArray(outputs_ptr);
}

void mxSolveAMIGO_MKLSENS(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File){

	mxArray *outputs_ptr,*sens_data,*sim_data,*sens_cell, *sim_cell,*flag_cell,*flag_value;
	const char *field_names[] = {"flag", "simulation","sensitivities"};


	mwSize dims[2]={1,1};
	mwSize dims_sens[3];
	mwSize dims_cell[1];
	// MATFile *pmat;
	int i,j,k,m,flag,counter;

	outputs_ptr=mxCreateStructArray(1, dims,3,field_names);

	dims_cell[0]=amigo_problem->n_models;

	flag_cell=mxCreateCellArray(1,dims_cell);
	sim_cell=mxCreateCellArray(1,dims_cell);
	sens_cell=mxCreateCellArray(1,dims_cell);

	mxSetField(outputs_ptr, 0,"flag", flag_cell);
	mxSetField(outputs_ptr, 0,"simulation", sim_cell);
	mxSetField(outputs_ptr, 0,"sensitivities", sens_cell);

	for (i = 0; i < amigo_problem->n_models; i++){

#ifdef MKL
		//last flag
		flag=get_AMIGO_model_sens(amigo_problem->amigo_models[i],0,1);
#else
		mexPrintf("Compile using MKL library and try again.");
		printf("Compile using MKL library and try again.");
		mexErrMsgTxt("Compile using MKL library and try again.");
#endif

		flag_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(flag_value)[0]=(double)flag;

		dims_sens[0]=amigo_problem->amigo_models[i]->n_times;
		dims_sens[1]=amigo_problem->amigo_models[i]->n_states;
		dims_sens[2]=amigo_problem->amigo_models[i]->n_total_x;

		sens_data=mxCreateNumericArray(3,dims_sens,mxDOUBLE_CLASS, mxREAL);

		sim_data=mxCreateDoubleMatrix(
			amigo_problem->amigo_models[i]->n_times ,
			amigo_problem->amigo_models[i]->n_states ,
			mxREAL
			);

		counter=0;
		for (j = 0; j < amigo_problem->amigo_models[i]->n_states; j++){
			for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
				mxGetPr(sim_data)[counter++]=amigo_problem->amigo_models[i]->sim_results[j][k];
			}
		}

		counter=0;
		for (m = 0; m < amigo_problem->amigo_models[i]->n_total_x; m++){
			for (j = 0; j < amigo_problem->amigo_models[i]->n_states; j++){
				for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
					mxGetPr(sens_data)[counter++]=amigo_problem->amigo_models[i]->sens_results[j][m][k];
				}
			}
		}

		mxSetCell(flag_cell,i, flag_value);
		mxSetCell(sim_cell, i, sim_data);
		mxSetCell(sens_cell, i, sens_data);
	}

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace)
		mexPutVariable("caller", "outputs", outputs_ptr);
	else mxDestroyArray(outputs_ptr);

}

void mxSolveAMIGOnl2sol(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File){


	mxArray *outputs_ptr,*xbest,*fbest,*nevals,*flag, *fails_cell, *stats_ptr, *stats_value,*observables_cell,*obs_data,*par;
	const char *field_names[] = {"par","xbest", "fbest", "nfevals","fail_stats","observables","flag"};
	const char *stats_names[] = {"flag","nst","nfe","nsetups","netf","nni","ncfn"
		,"nfSe","nfeS","nsetupsS","netfS","nniS","ncfnS"};

	mwSize dims[2]={1,1};
	int i,j,k,counter,total_par_size, n_model_pars, n_ics;
	mwSize dims_cell[1];
	mwSize dims_obs_cell[1];

	dims_obs_cell[0]=amigo_problem->n_models;

	if(amigo_problem->cvodes_gradient || amigo_problem->mkl_gradient)
		NL2SOL_AMIGO_pe(amigo_problem,1);
	else
		NL2SOL_AMIGO_pe(amigo_problem,0);

	outputs_ptr=mxCreateStructArray(1,dims,7,field_names);

	n_model_pars=amigo_problem->amigo_models[0]->n_pars;
	n_ics=amigo_problem->n_ics;
	total_par_size=n_model_pars+n_ics;

	xbest=mxCreateDoubleMatrix(1,amigo_problem->nx,mxREAL);
	par=mxCreateDoubleMatrix(1,total_par_size,mxREAL);
	fbest=mxCreateDoubleMatrix(1,1,mxREAL);
	nevals=mxCreateDoubleMatrix(1,1,mxREAL);
	observables_cell=mxCreateCellArray(1,dims_obs_cell);

	counter=0;

	for (i = 0;  i < amigo_problem->nx; i++){
		mxGetPr(xbest)[i]=amigo_problem->xbest[i];
	}

	for (i = 0;  i < amigo_problem->amigo_models[0]->n_pars; i++){
		mxGetPr(par)[i]=amigo_problem->amigo_models[0]->pars[i];
	}

	for (i = 0;  i < amigo_problem->n_pars; i++){
		mxGetPr(par)
			[amigo_problem->amigo_models[0]->index_opt_pars[i]]=amigo_problem->xbest[i];
		counter++;
	}

	for (i = n_model_pars;  i < total_par_size; i++){
		mxGetPr(par)[i]=amigo_problem->xbest[counter++];
	}

	mxGetPr(fbest)[0]=amigo_problem->local_fbest;
	mxGetPr(nevals)[0]=amigo_problem->local_nfeval;

	dims_cell[0]=amigo_problem->n_stored_fails;
	fails_cell=mxCreateCellArray(1,dims_cell);


	for (i = 0; i < amigo_problem->n_stored_fails; i++){


		stats_ptr=mxCreateStructArray(1, dims,13,stats_names);
		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->error_flag;
		mxSetField(stats_ptr, 0,"flag", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nst;
		mxSetField(stats_ptr, 0,"nst", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nfe;
		mxSetField(stats_ptr, 0,"nfe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nsetups;
		mxSetField(stats_ptr, 0,"nsetups", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->netf;
		mxSetField(stats_ptr, 0,"netf", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nni;
		mxSetField(stats_ptr, 0,"nni", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->ncfn;
		mxSetField(stats_ptr, 0,"ncfn", stats_value);

		if(amigo_problem->amigo_stats_containers[i]->sens){

			stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nfSe;
			mxSetField(stats_ptr, 0,"nfSe", stats_value);

			stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nfeS;
			mxSetField(stats_ptr, 0,"nfeS", stats_value);

			stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->netfS;
			mxSetField(stats_ptr, 0,"netfS", stats_value);

			stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->nniS;
			mxSetField(stats_ptr, 0,"nniS", stats_value);

			stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_stats_containers[i]->ncfnS;
			mxSetField(stats_ptr, 0,"ncfnS", stats_value);
		}

		mxSetCell(fails_cell,i, stats_ptr);

	}

	for (i = 0; i < amigo_problem->n_models; i++){

		obs_data=mxCreateDoubleMatrix(
			amigo_problem->amigo_models[i]->n_times ,
			amigo_problem->amigo_models[i]->n_observables,
			mxREAL
			);

		mxSetCell(observables_cell, i, obs_data);

		counter=0;
		for (j = 0; j < amigo_problem->amigo_models[i]->n_observables; j++){
			for (k = 0; k < amigo_problem->amigo_models[i]->n_times; k++){
				mxGetPr(obs_data)[counter++]=amigo_problem->amigo_models[i]->obs_results[j][k];
			}
		}

	}

	flag=mxCreateDoubleMatrix(1,1,mxREAL);

	mxGetPr(flag)[0]=(double)amigo_problem->local_flag;

	mxSetField(outputs_ptr, 0,"par", par);
	mxSetField(outputs_ptr, 0,"xbest", xbest);
	mxSetField(outputs_ptr, 0,"fbest", fbest);
	mxSetField(outputs_ptr, 0,"nfevals", nevals);
	mxSetField(outputs_ptr, 0,"fail_stats", fails_cell);
	mxSetField(outputs_ptr, 0,"observables", observables_cell);
	mxSetField(outputs_ptr, 0,"flag", flag);

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace){
		mexPutVariable("caller", "outputs", outputs_ptr);
	}
	else mxDestroyArray(outputs_ptr);

}


void mxEvalAMIGOlsq(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File){

	mxArray *outputs_ptr,*stats_ptr,*stats_value,*f_value,*stats_cell,*flag_cell,*flag_value;

	const char *field_names[] = {"success", "f","sim_stats"};

	const char *stats_names[] = {"flag","nst","nfe","nsetups","netf","nni","ncfn"};

	mwSize dims[2]={1,1};
	mwSize dims_cell[1];

	int i,j,k,flag,counter;

	outputs_ptr=mxCreateStructArray(1, dims,3,field_names);

	dims_cell[0]=amigo_problem->n_models;

	f_value=mxCreateDoubleMatrix(1,1,mxREAL);
	flag_cell=mxCreateCellArray(1,dims_cell);
	stats_cell=mxCreateCellArray(1,dims_cell);

	mxSetField(outputs_ptr, 0,"f", f_value);
	mxSetField(outputs_ptr, 0,"success", flag_cell);
	mxSetField(outputs_ptr, 0,"sim_stats", stats_cell);

	mxGetPr(f_value)[0]=(double)eval_AMIGO_problem_LSQ(amigo_problem);

	for (i = 0; i < amigo_problem->n_models; i++){

		flag_value=mxCreateDoubleMatrix(1,1,mxREAL);   

		stats_ptr=mxCreateStructArray(1, dims,7,stats_names);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->error_flag;
		mxSetField(stats_ptr, 0,"flag", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nst;
		mxSetField(stats_ptr, 0,"nst", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfe;
		mxSetField(stats_ptr, 0,"nfe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nsetups;
		mxSetField(stats_ptr, 0,"nsetups", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->netf;
		mxSetField(stats_ptr, 0,"netf", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nni;
		mxSetField(stats_ptr, 0,"nni", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->ncfn;
		mxSetField(stats_ptr, 0,"ncfn", stats_value);

		mxSetCell(stats_cell,i, stats_ptr);
		mxSetCell(flag_cell,i, flag_value);

	}

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace){
		mexPutVariable("caller", "outputs", outputs_ptr);
	}
	else mxDestroyArray(outputs_ptr);
}

void mxEvalAMIGOllk(AMIGO_problem* amigo_problem,int save2Workspace, char* save2File){

	mxArray *outputs_ptr,*stats_ptr,*stats_value,*f_value,*stats_cell,*flag_cell,*flag_value;

	const char *field_names[] = {"success", "f","sim_stats"};

	const char *stats_names[] = {"flag","nst","nfe","nsetups","netf","nni","ncfn"};

	mwSize dims[2]={1,1};
	mwSize dims_cell[1];

	int i,j,k,flag,counter;

	outputs_ptr=mxCreateStructArray(1, dims,3,field_names);

	dims_cell[0]=amigo_problem->n_models;

	f_value=mxCreateDoubleMatrix(1,1,mxREAL);
	flag_cell=mxCreateCellArray(1,dims_cell);
	stats_cell=mxCreateCellArray(1,dims_cell);

	mxSetField(outputs_ptr, 0,"f", f_value);
	mxSetField(outputs_ptr, 0,"success", flag_cell);
	mxSetField(outputs_ptr, 0,"sim_stats", stats_cell);

	mxGetPr(f_value)[0]=(double)eval_AMIGO_problem_LLK(amigo_problem);

	for (i = 0; i < amigo_problem->n_models; i++){

		flag_value=mxCreateDoubleMatrix(1,1,mxREAL);   

		stats_ptr=mxCreateStructArray(1, dims,7,stats_names);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->error_flag;
		mxSetField(stats_ptr, 0,"flag", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nst;
		mxSetField(stats_ptr, 0,"nst", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nfe;
		mxSetField(stats_ptr, 0,"nfe", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nsetups;
		mxSetField(stats_ptr, 0,"nsetups", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->netf;
		mxSetField(stats_ptr, 0,"netf", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->nni;
		mxSetField(stats_ptr, 0,"nni", stats_value);

		stats_value=mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr(stats_value)[0]=(double)amigo_problem->amigo_models[i]->amigo_model_stats->ncfn;
		mxSetField(stats_ptr, 0,"ncfn", stats_value);

		mxSetCell(stats_cell,i, stats_ptr);
		mxSetCell(flag_cell,i, flag_value);

	}

	mxAMIGOsave2File(outputs_ptr,save2File);

	if(save2Workspace){
		mexPutVariable("caller", "outputs", outputs_ptr);
	}
	else mxDestroyArray(outputs_ptr);
}

int mxAMIGOsave2File(mxArray *outputs_ptr,char* save2File){

	MATFile *pmat;

	if(strcmp(save2File,"")){
		pmat = matOpen(save2File, "w");
		if (pmat == NULL) {
			printf("Error creating file %s\n", save2File);
			printf("(Do you have write permission in this directory?)\n");
			return(EXIT_FAILURE);
		}

		matPutVariable(pmat, "outputs", outputs_ptr);

		if (matClose(pmat) != 0) {
			printf("Error closing file %s\n",save2File);
			return(EXIT_FAILURE);
		}
	}
	return(EXIT_SUCCESS);
}