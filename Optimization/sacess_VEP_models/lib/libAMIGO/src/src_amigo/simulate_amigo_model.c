#include "simulate_amigo_model.h"


#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define ZERO  RCONST(0.0)

static int check_flag(void *flagvalue,int opt, void* cvode_mem, AMIGO_model* amigo_model);


int set_simulation_stats(void* cvode_mem, AMIGO_model* amigo_model){
	  
	long int nst;
	long int nfe, nsetups, nni, ncfn, netf;
	long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
	long int nli, ncfl, npe, nps;
	int flag;
	
	flag=CVodeGetNumSteps(cvode_mem, &nst);
	if(flag>=0)amigo_model->amigo_model_stats->nst=nst;else amigo_model->amigo_model_stats->nst=-1;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	if(flag>=0)amigo_model->amigo_model_stats->nfe=nfe;else amigo_model->amigo_model_stats->nfe=-1;
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	if(flag>=0)amigo_model->amigo_model_stats->nsetups=nsetups;else amigo_model->amigo_model_stats->nsetups=-1;
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	if(flag>=0)amigo_model->amigo_model_stats->netf=netf;else amigo_model->amigo_model_stats->netf=-1;
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	if(flag>=0)amigo_model->amigo_model_stats->nni=nni;else amigo_model->amigo_model_stats->nni=-1;
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	if(flag>=0)amigo_model->amigo_model_stats->ncfn=ncfn;else amigo_model->amigo_model_stats->ncfn=-1;

	if (amigo_model->compute_sens) {
		amigo_model->amigo_model_stats->sens=1;
		flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
		if(flag>=0)amigo_model->amigo_model_stats->nfSe=nfSe;else amigo_model->amigo_model_stats->nfSe=-1;
		flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
		if(flag>=0)amigo_model->amigo_model_stats->nfeS=nfeS;else amigo_model->amigo_model_stats->nfeS=-1;
		flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
		if(flag>=0)amigo_model->amigo_model_stats->nsetupsS=nsetupsS;else amigo_model->amigo_model_stats->nsetupsS=-1;
		flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
		if(flag>=0)amigo_model->amigo_model_stats->netfS=netfS;else amigo_model->amigo_model_stats->netfS=-1;
		flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
		if(flag>=0)amigo_model->amigo_model_stats->nniS=nniS;else amigo_model->amigo_model_stats->nniS=-1;
		flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
		if(flag>=0)amigo_model->amigo_model_stats->ncfnS=ncfnS;else amigo_model->amigo_model_stats->ncfnS=-1;
	}else{
		amigo_model->amigo_model_stats->sens=0;
	}
	return(0);
}



int simulate_amigo_model(AMIGO_model* amigo_model,int verbose){
	int i,j,k,neq,NS,flag,stop_ti,temp,TS_index;
	booleantype err_con;
	realtype tout, tstop, ti, tf, t,previous_tstop;
	
	int* plist;

	N_Vector y;
	N_Vector *yS;

	void *cvode_mem;
	int counter;
	cvode_mem = NULL;
	y = NULL;
	neq=(*amigo_model).n_states;
	y = N_VNew_Serial(neq);
	NS=amigo_model->n_opt_pars+amigo_model->n_opt_ics;
	err_con=FALSE;
	ti=amigo_model->t0;
	tf=amigo_model->tf;



	if (check_flag((void *)y, 0,cvode_mem,amigo_model))return(0);

	counter=0;
	plist=(int*)malloc(sizeof(int)*NS);
	for (i = 0;  i < amigo_model->n_opt_pars; i++){
		plist[counter++]=amigo_model->index_opt_pars[i];
	}

	for (i = 0;  i < amigo_model->n_opt_ics; i++){
		plist[counter++]=i+amigo_model->n_opt_pars;
	}

	//Fill the first position of the outputs with y0
	if((*amigo_model).t[0]==(*amigo_model).t0){
		tout=ti;
		TS_index=1;
		for(i=0; i<(*amigo_model).n_states; i++){
			Ith(y,i) =(realtype)(*amigo_model).y0[i];
			(*amigo_model).sim_results[i][0]=(*amigo_model).y0[i];
		}
	}
	else{
		tout=ti;
		TS_index=0;
		for(i=0; i<(*amigo_model).n_states; i++){
			Ith(y,i) =(realtype)(*amigo_model).y0[i];
		}
	}

	//Initialize CVODE with BDF and NEWTON iterartion for stiff problem
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if(check_flag((void *)cvode_mem,  1,cvode_mem,amigo_model)) {
		N_VDestroy_Serial(y);
		return(0);
	}

	flag = CVodeInit(cvode_mem,(*amigo_model).rhs, ti, y);
	if (check_flag(&flag, 1,cvode_mem,amigo_model)){

		N_VDestroy_Serial(y);
		/* Free integrator memory */
		CVodeFree(&cvode_mem);
		return(0);
	}

	/* Set f_data */
	flag = CVodeSetUserData(cvode_mem, amigo_model);
	if(check_flag(&flag,  1,cvode_mem,amigo_model)){
		N_VDestroy_Serial(y);
		/* Free integrator memory */
		CVodeFree(&cvode_mem);
		return(0);
	}

	flag = CVodeSStolerances(cvode_mem,(realtype)(*amigo_model).reltol,(realtype)(*amigo_model).atol);

	if(check_flag(&flag,  1,cvode_mem,amigo_model)) return(0);

	flag = CVodeSetErrFile(cvode_mem, NULL);	

//	flag = CVLapackDense(cvode_mem, neq);
	flag = CVDense(cvode_mem, neq);
	if (check_flag(&flag, 1,cvode_mem,amigo_model)) return(0);

	/* Set maxnumsteps */
	flag = CVodeSetMaxNumSteps(cvode_mem, (*amigo_model).max_num_steps);
	if(check_flag(&flag,  1,cvode_mem,amigo_model)){
		N_VDestroy_Serial(y);
		/* Free integrator memory */
		CVodeFree(&cvode_mem);
		return(0);
	
	}

	CVodeSetMaxStep(cvode_mem,(realtype)(*amigo_model).max_step_size);

	CVodeSetMaxErrTestFails(cvode_mem, (*amigo_model).max_error_test_fails);

	//COnfigure sensibility Analysis*/
	if (amigo_model->compute_sens){

		yS = N_VCloneVectorArray_Serial(NS, y);
		if (check_flag((void *)yS, 0,cvode_mem,amigo_model)) return(0);

		//Initialize sensibilities for parameters with 0
		for (i=0;i<amigo_model->n_opt_pars;i++){
			N_VConst(ZERO, yS[i]);
			if((*amigo_model).t[0]==(*amigo_model).t0){
				for (j=0;j<(*amigo_model).n_states;j++){
					(*amigo_model).sens_results[j][i][0]=ZERO;
				}
			}
		}
		for (i=0;i<amigo_model->n_opt_ics;i++){
	
			for (j=0;j<(*amigo_model).n_states;j++){
				if(j==amigo_model->index_opt_ics[i]){
					Ith(yS[i+amigo_model->n_opt_pars],j)=1;						
					if((*amigo_model).t[0]==(*amigo_model).t0){
						(*amigo_model).sens_results[j][i+amigo_model->n_opt_pars][0]=1;
					}
				 }else{
					Ith(yS[i+amigo_model->n_opt_pars],j)=0;
					if((*amigo_model).t[0]==(*amigo_model).t0){
						(*amigo_model).sens_results[j][i+amigo_model->n_opt_pars][0]=0;
					}
				}
			}
		}
	
		flag = CVodeSensInit1(cvode_mem, NS, CV_STAGGERED, NULL, yS);
		if(check_flag(&flag,  1,cvode_mem,amigo_model)) return(0);

		flag = CVodeSensEEtolerances(cvode_mem);
		if(check_flag(&flag,  1,cvode_mem,amigo_model)) return(0);

		flag = CVodeSetSensErrCon(cvode_mem, err_con);
		if (check_flag(&flag,  1,cvode_mem,amigo_model)) return(0);

		flag = CVodeSetSensParams(cvode_mem,(realtype*)(*amigo_model).pars, NULL,plist);
		if (check_flag(&flag, 1,cvode_mem,amigo_model)) return(0);         

	} 

	//If the Jacobian matrix has been provided, then use it
	if((*amigo_model).use_jacobian){

		/* Attach linear solver */
		flag = CVDense(cvode_mem, neq);
       //         flag = CVLapackDense(cvode_mem, neq);
		if(check_flag(&flag, 1,cvode_mem,amigo_model)) return(0);

		flag = CVDlsSetDenseJacFn(cvode_mem, (*amigo_model).jac);

		if (check_flag(&flag, 1,cvode_mem,amigo_model)) return(0);
	}

	tstop=tout;
	stop_ti=0;

	while(stop_ti<(*amigo_model).n_controls_t && tstop<=tout){
		previous_tstop=tstop;
		tstop=(*amigo_model).controls_t[stop_ti++];
	}

	for (i = TS_index;  i < (*amigo_model).n_times; i++){

		//update tout 
		tout=(*amigo_model).t[i];
		(*amigo_model).tlast=previous_tstop;
		(*amigo_model).index_t_stim=stop_ti-2;

		//Try to update tstop
		while(tstop<=tout){
			//Find if there are still controls to be executed
			if(stop_ti<(*amigo_model).n_controls_t && (*amigo_model).n_controls_t>0){
				flag = CVode(cvode_mem, tstop, y, &t, CV_TSTOP_RETURN);

				if (check_flag(&flag,  1,cvode_mem,amigo_model)){

					N_VDestroy_Serial(y);
					if (amigo_model->compute_sens) {
						//Free yS vector 
						N_VDestroyVectorArray_Serial(yS, NS);   
					}
					CVodeFree(&cvode_mem);
					free(plist);
					return(0);
				}

				//We have reached TSTOP, reinitialize memory and update tstop
				amigo_model->changeYatTcon(amigo_model,tstop,y);
				CVodeReInit(cvode_mem, tstop, y); 
				previous_tstop=tstop;
				tstop=(*amigo_model).controls_t[stop_ti++];
				CVodeSetStopTime(cvode_mem, tstop);
				(*amigo_model).index_t_stim=stop_ti-2;
				(*amigo_model).tlast=previous_tstop;
				

			}
			//In case there are no more control times break the while
			else break;

		}	
		//Make sure tout is different from tstop. Dont integrate twice... causes an error
		if(tout!=previous_tstop){
			flag = CVode(cvode_mem, tout, y, &t, CV_TSTOP_RETURN);
			if (check_flag(&flag,1,cvode_mem,amigo_model)){	
				N_VDestroy_Serial(y);
				if (amigo_model->compute_sens) {
					//Free yS vector 
					N_VDestroyVectorArray_Serial(yS, NS);   
				}
				CVodeFree(&cvode_mem);	
				free(plist);
				return(0);
			}
		}
		for (j = 0; j < (*amigo_model).n_states; j++){
			(*amigo_model).sim_results[j][i]=(double)Ith(y,j);
		}
		//Get sensitivies
		if ((*amigo_model).compute_sens) {

			flag = CVodeGetSens(cvode_mem, &tf, yS);

			if (check_flag(&flag, 1,cvode_mem,amigo_model)) 
				break;

			for (j=0;j<(*amigo_model).n_states;j++){
				for (k=0;k<NS;k++){
					(*amigo_model).sens_results[j][k][i]=Ith(yS[k],j);
				}
			}
		}
	}
	N_VDestroy_Serial(y);

	if (amigo_model->compute_sens) {
		//Free yS vector 
		N_VDestroyVectorArray_Serial(yS, NS);   
	}

	/* Free integrator memory */
	set_simulation_stats(cvode_mem, amigo_model);
	amigo_model->amigo_model_stats->error_flag=flag;
	CVodeFree(&cvode_mem);
	free(plist);
	return(1);
}

static int check_flag(void *flagvalue,int opt, void* cvode_mem, AMIGO_model* amigo_model)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL){ 
		set_simulation_stats(cvode_mem, amigo_model);
		amigo_model->amigo_model_stats->error_flag=NULL;
		return(1);
	}

	/* Check if flag < 0 */
	else if (opt == 1){
		errflag = (int *) flagvalue;
		if (*errflag < 0){
			set_simulation_stats(cvode_mem,amigo_model);
			amigo_model->amigo_model_stats->error_flag=*errflag;
			return(1);
		}
	}

	
	return(0);
}
