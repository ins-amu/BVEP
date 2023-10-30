#include "AMIGO_model.h"

//! Allocate necessay memory for for an AMIGO_model
/*!
 Receives the variables necessary to know the dimensions of the the problem.
Includes model definitions, experimental design and data.
 */

AMIGO_model* allocate_AMIGO_model
(
	int n_states,		int n_observables,	int n_pars,	
	int n_opt_pars,		int n_times,		int n_opt_ics,		
	int n_controls,		int n_controls_t,	int exp_num
){

	int i,j;
	//Allocate the structure
	AMIGO_model new_amigo_model;
	AMIGO_model* amigo_model;
	amigo_model=(AMIGO_model*)malloc(sizeof(new_amigo_model));
	amigo_model->exp_num=exp_num;
	
	//States and observables
	amigo_model->n_states=n_states;
	amigo_model->n_observables=n_observables;
	amigo_model->index_observables=(int*)malloc(sizeof(int)*n_observables);
	for (i = 0; i < n_observables; i++)amigo_model->index_observables[i]=i;
	

	//Simulation Pars
	amigo_model->n_pars=n_pars;
	amigo_model->n_opt_pars=n_opt_pars;
        amigo_model->n_total_x=n_opt_pars+n_opt_ics;
	
	//->pars is larger because it will be used for sensitivity. It does not matter what values we put in here!!!
	amigo_model->pars=(double*)malloc(sizeof(double)*(n_pars+n_opt_ics));
	for (i = 0; i < n_pars; i++)amigo_model->pars[i]=0;
	
	amigo_model->index_opt_pars=(int*)malloc(sizeof(int)*n_opt_pars);
	
	for (i = 0; i < n_opt_pars; i++)amigo_model->index_opt_pars[i]=i;
	
	amigo_model->pars_guess=(double*)malloc(sizeof(double)*(n_opt_pars));
	amigo_model->pars_LB=(double*)malloc(sizeof(double)*(n_opt_pars));
	amigo_model->pars_UB=(double*)malloc(sizeof(double)*(n_opt_pars));
	

	//Simulation times
	amigo_model->t0=0;
	amigo_model->tf=(double)n_times-1;
	amigo_model->n_times=n_times;
	amigo_model->t=(double*)malloc(sizeof(double)*n_times);
	for (i = 0; i < n_times; i++)amigo_model->t[i]=(double)i;
	
	//Initial conditions
	amigo_model->y0=(double*)malloc(sizeof(double)*n_states);
	amigo_model->n_opt_ics=n_opt_ics;
	amigo_model->index_opt_ics=(int*)malloc(sizeof(int)*n_opt_ics);
	amigo_model->y0_guess=(double*)malloc(sizeof(double)*n_opt_ics);
	amigo_model->y0_LB=(double*)malloc(sizeof(double)*n_opt_ics);
	amigo_model->y0_UB=(double*)malloc(sizeof(double)*n_opt_ics);
	
	//Controls
	amigo_model->n_controls=n_controls;
	amigo_model->n_controls_t=n_controls_t;
	amigo_model->controls_t=(double*)malloc(n_controls_t*sizeof(double));

	for (i = 0; i < amigo_model->n_controls_t; i++) amigo_model->controls_t[i]=i;
	
	amigo_model->controls_v=(double**)malloc(amigo_model->n_controls*sizeof(double*));
	for (i = 0; i < amigo_model->n_controls; i++) {
		amigo_model->controls_v[i]=(double*)malloc((amigo_model->n_controls_t-1)*sizeof(double));
		
		for (j= 0; j < amigo_model->n_controls_t-1; j++)amigo_model->controls_v[i][j]=1;
	}
	
	amigo_model->slope=(double**)malloc((amigo_model->n_controls)*sizeof(double*));
	for (i= 0; i < amigo_model->n_controls; i++){	
		amigo_model->slope[i]=(double*)malloc((amigo_model->n_controls_t-1)*sizeof(double));
		for (j= 0; j < amigo_model->n_controls_t-1; j++){
			amigo_model->slope[i][j]=0;   			
		}
	}

	//Storing matrixes
	amigo_model->sim_results=(double**)malloc(sizeof(double*)*n_states);
	for (i = 0; i < n_states; i++){
		amigo_model->sim_results[i]=(double*)malloc(sizeof(double)*n_times);
	}	
	
	amigo_model->sens_results=(double***)malloc(sizeof(double**)*n_states);
	for (i = 0; i < (*amigo_model).n_states; i++){
		amigo_model->sens_results[i]=(double**)malloc(sizeof(double*)*(n_opt_pars+n_opt_ics));
		for (j = 0; j < (n_opt_pars+n_opt_ics); j++){
			amigo_model->sens_results[i][j]=(double*)malloc(sizeof(double)*n_times);
		}
	}

	amigo_model->sens_obs=(double***)malloc(sizeof(double**)*n_observables);
	for (i = 0; i < n_observables; i++){
		amigo_model->sens_obs[i]=(double**)malloc(sizeof(double*)*(n_opt_pars+n_opt_ics));
		for (j = 0; j < (n_opt_pars+n_opt_ics); j++){
			amigo_model->sens_obs[i][j]=(double*)malloc(sizeof(double)*n_times);
		}
	}
	
	amigo_model->obs_results=(double**)malloc(sizeof(double*)*n_observables);
	for (i = 0; i < n_observables; i++){
		amigo_model->obs_results[i]=(double*)malloc(sizeof(double)*n_times);
	}

	//Experimental Data
	amigo_model->Q=(double**)malloc(sizeof(double*)*n_observables);
	amigo_model->exp_data=(double**)malloc(sizeof(double*)*n_observables);
	for (i = 0; i < n_observables; i++){
		amigo_model->Q[i]=(double*)malloc(sizeof(double)*n_times);
		amigo_model->exp_data[i]=(double*)malloc(sizeof(double)*n_times);
		for (j = 0; j < n_times; j++){
			amigo_model->Q[i][j]=1;
			amigo_model->exp_data[i][j]=0;
		}
	}

	amigo_model->w_obs=(double*)malloc(sizeof(double)*n_observables);
	for (i = 0; i < n_observables; i++)amigo_model->w_obs[i]=1;
	
	
	//Simulation Related Parameter
	amigo_model->reltol=1e-6;
	amigo_model->atol=1e-6;
	amigo_model->max_step_size=DBL_MAX;
	amigo_model->max_num_steps=1000000;
	amigo_model->max_error_test_fails=50;
	amigo_model->use_jacobian=0;
	amigo_model->compute_sens=0;
	amigo_model->mkl_tol=1e-3;

	amigo_model->amigo_model_stats=(AMIGO_model_stats*)malloc(sizeof(AMIGO_model_stats));
	amigo_model->amigo_model_stats->npars=0;
        amigo_model->use_obs_func=0;
        amigo_model->use_sens_obs_func=0;

	return(amigo_model);
}

void free_AMIGO_model(AMIGO_model* amigo_model){

	int i,j;
	//States and observables
	
	free(amigo_model->index_observables);
	
	//Simulation Pars
	
	free(amigo_model->pars);
	free(amigo_model->index_opt_pars);
	free(amigo_model->pars_guess);
	free(amigo_model->pars_LB);
	free(amigo_model->pars_UB);

	//Simulation times
	free(amigo_model->t);
	
	//Initial conditions
	free(amigo_model->y0);
	free(amigo_model->index_opt_ics);
	free(amigo_model->y0_guess);
	free(amigo_model->y0_LB);
	free(amigo_model->y0_UB);

	//Controls
	free(amigo_model->controls_t);

	for (i = 0; i < amigo_model->n_controls; i++) free(amigo_model->controls_v[i]);
	free(amigo_model->controls_v);

	for (i= 0; i < amigo_model->n_controls; i++) free(amigo_model->slope[i]);
	free(amigo_model->slope);

	//Storing matrixes
	for (i = 0; i < amigo_model->n_states; i++) free(amigo_model->sim_results[i]);
	free(amigo_model->sim_results);
	
	
	for (i = 0; i < (*amigo_model).n_states; i++){
		for (j = 0; j < amigo_model->n_opt_pars+amigo_model->n_opt_ics; j++){
			free(amigo_model->sens_results[i][j]);
		}
		free(amigo_model->sens_results[i]);
	}
	free(amigo_model->sens_results);

	for (i = 0; i < (*amigo_model).n_observables; i++){
		for (j = 0; j < amigo_model->n_opt_pars+amigo_model->n_opt_ics; j++){
			free(amigo_model->sens_obs[i][j]);
		}
		free(amigo_model->sens_obs[i]);
	}
	free(amigo_model->sens_obs);
	
	for (i = 0; i < amigo_model->n_observables; i++)free(amigo_model->obs_results[i]);
	free(amigo_model->obs_results);
	
	//Error data

	for (i = 0; i < amigo_model->n_observables; i++){
		free(amigo_model->Q[i]);
		free(amigo_model->exp_data[i]);
	}
	free(amigo_model->Q);
	free(amigo_model->exp_data);
	free(amigo_model->w_obs);

	free_AMIGO_model_stats(amigo_model->amigo_model_stats);

	//structure
	free(amigo_model);
}

 int simulate_AMIGO_model(AMIGO_model* amigo_model,int sens){

	amigo_model->compute_sens=sens;
	simulate_amigo_model(amigo_model,0);
 	return 1;
}

int simulate_AMIGO_model_observables(AMIGO_model* amigo_model,int sens){
	
	int i,j,flag;
	amigo_model->compute_sens=sens;

	flag=simulate_amigo_model(amigo_model,0);
	if(amigo_model->use_obs_func){
		if(amigo_model->n_observables>0)
			amigo_model->obs_func(amigo_model);
	}else{
		for (i = 0; i < amigo_model->n_observables; ++i){
			for (j = 0; j < amigo_model->n_times; ++j){
                            //printf("amigo_model->index_observables[i] %d\n", amigo_model->index_observables[i]);
                            //printf("amigo_model->w_obs[i] %lf\n", amigo_model->w_obs[i]);
                            //printf("amigo_model->sim_results[amigo_model->index_observables[i]][j] %lf\n", amigo_model->sim_results[amigo_model->index_observables[i]][j]);
                            
				amigo_model->obs_results[i][j]=
					amigo_model->sim_results[amigo_model->index_observables[i]][j]*amigo_model->w_obs[i];
			}
		}
	}
	return(flag);
}

double lsq_AMIGO_model(AMIGO_model* amigo_model){
	int i,j;
	double value;
	value=0;
	for (i = 0; i < amigo_model->n_observables; ++i){
		for (j = 0; j < amigo_model->n_times; ++j){
			value+=(pow((amigo_model->exp_data[i][j]-
				amigo_model->obs_results[i][j]),2))/
				pow(amigo_model->Q[i][j],2);
		}
	}
	return(value);
}

#ifdef MKL
#include <mkl.h>

int get_mkl_sens(AMIGO_model* amigo_model) {

	int i, j, k,counter=0;
	int np=amigo_model->n_opt_pars+amigo_model->n_opt_ics;
	int ns=amigo_model->n_states*amigo_model->n_times; 
	double *x=(double*)malloc(sizeof(double)*np);
	
	extern void mkl_jacobian_function(MKL_INT*, MKL_INT*, double*, double*, void*);
	double* dxdt=(double*)malloc(np*ns*sizeof(double));

	for (i=0;i<amigo_model->n_opt_pars;i++){
		x[counter++]=(*amigo_model).pars[amigo_model->index_opt_pars[i]];			
	}

	//Initial conditions after
	for (i=0;i<amigo_model->n_opt_ics;i++){
		x[counter++]=amigo_model->y0[amigo_model->index_opt_ics[i]];
	}

	djacobix(mkl_jacobian_function, &np, &ns, dxdt, x, 
		&amigo_model->mkl_tol, amigo_model);

	counter=0;
	for (i = 0; i < np; ++i) {
		for (j = 0; j < amigo_model->n_states; ++j) {
			for (k = 0; k < amigo_model->n_times; ++k) {
				amigo_model->sens_results[j][i][k]=dxdt[counter++];
			}
		}
	}
	free(x);
	free(dxdt);
}

void mkl_jacobian_function(MKL_INT* m, MKL_INT* n, double* x, double* f, void* data){

	int i, j, counter=0;
	
	AMIGO_model* amigo_model=(AMIGO_model*) data;

	int np= amigo_model->n_opt_pars+amigo_model->n_opt_ics;
	
	//kinetic pars first
	for (i=0;i<amigo_model->n_opt_pars;i++){
		(*amigo_model).pars[amigo_model->index_opt_pars[i]]=x[counter++];			
	}

	//Initial conditions after
	for (i=0;i<amigo_model->n_opt_ics;i++){
		amigo_model->y0[amigo_model->index_opt_ics[i]]=x[counter++];
	}

	counter=0;
	//Simulate observables with cvodes sensitivtiy off sens=0
	amigo_model->compute_sens=0;

	simulate_AMIGO_model_observables(amigo_model,0);

	for (i = 0; i < amigo_model->n_states; ++i) {
		for (j = 0; j < amigo_model->n_times; ++j) {
				f[counter++]=
					amigo_model->sim_results[i][j];
		}
	}
}

#endif

int get_AMIGO_model_sens(AMIGO_model* amigo_model,int cvodes, int mkl){

	int flag=0,i,j,k;

	if(cvodes){
            flag=simulate_AMIGO_model_observables(amigo_model,1);
        }

	//If both are active and if cvodes fails the function also tries MKL
	#ifdef MKL
		if(mkl || !flag){
			flag=get_mkl_sens(amigo_model);
			flag=1;
		}
	#endif	

    if(amigo_model->use_sens_obs_func){
		if(amigo_model->n_observables>0){
			amigo_model->obs_sen_func(amigo_model);
        }
	}else{
		for (i = 0; i < amigo_model->n_observables; ++i) {
			for (j = 0; j < amigo_model->n_total_x; ++j) {
				for (k = 0; k < amigo_model->n_times; ++k) {
					amigo_model->sens_obs[i][j][k]=
						amigo_model->sens_results[amigo_model->index_observables[i]][j][k]*
																			amigo_model->w_obs[i];
				}
			}
		}
    }

	return flag;
}

