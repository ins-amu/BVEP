#include "AMIGO_problem.h"
#include "omp.h"

EXPORTIT AMIGO_problem* allocate_AMIGO_problem(int n_models, AMIGO_model** amigo_models){

	int i,j,nx,counter;
	AMIGO_problem new_amigo_problem;
	AMIGO_problem* amigo_problem=(AMIGO_problem*)malloc(sizeof(new_amigo_problem));
	amigo_problem->n_models=n_models;

	nx=amigo_models[0]->n_opt_pars;
	amigo_problem->n_pars=nx;
	amigo_problem->n_ics=0;

	amigo_problem->n_data=0;

	for (i = 0;  i < n_models; i++){
		nx+=amigo_models[i]->n_opt_ics;
		amigo_problem->n_ics+=amigo_models[i]->n_opt_ics;
		amigo_problem->n_data+=amigo_models[i]->n_observables*amigo_models[i]->n_times;
	}
	amigo_problem->amigo_models=amigo_models;

	amigo_problem->nx=nx;

	amigo_problem->x=(double*)malloc(sizeof(double)*nx);
	amigo_problem->xbest=(double*)malloc(sizeof(double)*nx);
	amigo_problem->x0=(double*)malloc(sizeof(double)*nx);
	amigo_problem->LB=(double*)malloc(sizeof(double)*nx);
	amigo_problem->UB=(double*)malloc(sizeof(double)*nx);

	counter=0;
	for (i = 0;  i < amigo_problem->n_pars; i++){
		amigo_problem->x0[i]=amigo_problem->amigo_models[0]->pars_guess[i];
		amigo_problem->xbest[i]=amigo_problem->amigo_models[0]->pars_guess[i];
		amigo_problem->x[i]=amigo_problem->amigo_models[0]->pars[i];
		amigo_problem->LB[i]=amigo_problem->amigo_models[0]->pars_LB[i];
		amigo_problem->UB[i]=amigo_problem->amigo_models[0]->pars_UB[i];
		counter++;
	}

	for (i = 0;  i < amigo_problem->n_models; i++){
		for (j = 0;  j < amigo_problem->amigo_models[i]->n_opt_ics; j++){
			amigo_problem->x0[counter]=amigo_problem->amigo_models[i]->y0_guess[j];
			amigo_problem->xbest[counter]=amigo_problem->amigo_models[i]->y0_guess[j];
			amigo_problem->x[counter]=amigo_problem->amigo_models[i]->y0[j];
			amigo_problem->LB[counter]=amigo_problem->amigo_models[i]->y0_LB[j];
			amigo_problem->UB[counter]=amigo_problem->amigo_models[i]->y0_UB[j];
			counter++;
		}

	}

	amigo_problem->n_fails=0;
	amigo_problem->n_stored_fails=0;
	amigo_problem->n_max_store_fails=1000;

	amigo_problem->local_fbest=DBL_MAX;
	amigo_problem->local_max_iter=500;
	amigo_problem->local_max_evals=1000;
	amigo_problem->nevals=0;
	amigo_problem->local_nfeval=0;

	amigo_problem->nevals=0;
	amigo_problem->verbose=0;

	amigo_problem->temp_min=DBL_MAX;
	amigo_problem->nthreads=1;


	return(amigo_problem);
}

EXPORTIT void free_AMIGO_problem(AMIGO_problem* amigo_problem){

	int i;

	for (i = 0;  i < amigo_problem->n_models; i++){
		free_AMIGO_model(amigo_problem->amigo_models[i]);
	}

	free(amigo_problem->x);
	free(amigo_problem->xbest);
	free(amigo_problem->x0);
	free(amigo_problem->LB);
	free(amigo_problem->UB);

	for (i = 0;  i < amigo_problem->n_stored_fails; i++){
		free_AMIGO_model_stats(amigo_problem->amigo_stats_containers[i]);
	}
	if(amigo_problem->n_stored_fails>0){
		free(amigo_problem->amigo_stats_containers);
	}

	free(amigo_problem);

}

EXPORTIT void set_AMIGO_problem_pars(double* x, AMIGO_problem* amigo_problem){

	int counter,i,j;
	counter=0;
	//update pars in all models

	for (i = 0;  i < amigo_problem->nx; i++){
		amigo_problem->x[i]=x[i];
	}

	for (i = 0;  i < amigo_problem->n_pars; i++){
		for (j = 0;  j < amigo_problem->n_models; j++){
			
			amigo_problem->amigo_models[j]->pars[amigo_problem->amigo_models[j]->index_opt_pars[i]]=x[i];
		}
		counter++;
	}

	for (i = 0;  i < amigo_problem->n_models; i++){
		for (j = 0;  j < amigo_problem->amigo_models[i]->n_opt_ics; j++){
			amigo_problem->amigo_models[i]->y0[amigo_problem->amigo_models[i]->index_opt_ics[j]]=x[counter++];
		}
	}
}




EXPORTIT  double eval_AMIGO_problem_LSQ(AMIGO_problem* amigo_problem){



	int i, j, k, n_exps, n_times, n_obs,flag;

	int counter=0;

	double res=0;

	




        for (i = 0; i < amigo_problem->n_models; ++i) {

		

			flag=simulate_AMIGO_model_observables(amigo_problem->amigo_models[i],0);

			if(!flag)res+=DBL_MAX;

	}	


	for (i = 0; i < amigo_problem->n_models; ++i) {

		n_obs=amigo_problem->amigo_models[i]->n_observables;

			for (j = 0; j < n_obs; ++j) {

			

				n_times=amigo_problem->amigo_models[i]->n_times;



				for (k = 0; k < n_times; ++k) {

					if( !isnan(amigo_problem->amigo_models[i]->exp_data[j][k]) &&

						!isnan(amigo_problem->amigo_models[i]->Q[j][k]) && 

						amigo_problem->amigo_models[i]->Q[j][k]!=0){

							res+=pow(

								(

									(amigo_problem->amigo_models[i]->obs_results[j][k]-

									amigo_problem->amigo_models[i]->exp_data[j][k])

									*amigo_problem->amigo_models[i]->Q[j][k]

								),

								2

							);

						}

				}

			}	

		}



	return(res);

}


EXPORTIT double eval_AMIGO_problem_LLK(AMIGO_problem* amigo_problem){

		int i, j, k, n_exps, n_times, n_obs,flag;
	int counter=0;
	double res=0.0;
	
//	#pragma omp parallel num_threads(amigo_problem->nthreads)
//    {
//        #pragma omp for schedule(dynamic,1) private(i,flag)
		for (i = 0; i < amigo_problem->n_models; ++i) {
		
			//0 For no sensitivity
			flag=simulate_AMIGO_model_observables(amigo_problem->amigo_models[i],0);
		}	
//	}


	for (i = 0; i < amigo_problem->n_models; ++i) {
		n_obs=amigo_problem->amigo_models[i]->n_observables;
			for (j = 0; j < n_obs; ++j) {
			
				n_times=amigo_problem->amigo_models[i]->n_times;

				for (k = 0; k < n_times; ++k) {
					if( !isnan(amigo_problem->amigo_models[i]->exp_data[j][k]) &&
						!isnan(amigo_problem->amigo_models[i]->Q[j][k]) && 
						amigo_problem->amigo_models[i]->Q[j][k]!=0){
							res+=pow(
								(
									(amigo_problem->amigo_models[i]->obs_results[j][k]-
									amigo_problem->amigo_models[i]->exp_data[j][k])
									/amigo_problem->amigo_models[i]->Q[j][k]
								),
								2
							);
						}
				}
			}	
		}

	return(res);
}



void handle_AMIGO_problem_stat_fails(int ith_model,AMIGO_problem* amigo_problem){

	amigo_problem->n_fails++;

	if(amigo_problem->n_stored_fails==0){
		//If it is the first then the container is not allocated at all. Allocate with size one
		amigo_problem->amigo_stats_containers=(AMIGO_model_stats**)malloc(sizeof(AMIGO_model_stats*));
		//Copy the value here
		amigo_problem->amigo_stats_containers[0]=amigo_problem->amigo_models[ith_model]->amigo_model_stats;
		////put new stats back in the modle
		amigo_problem->amigo_models[ith_model]->amigo_model_stats=(AMIGO_model_stats*)malloc(sizeof(AMIGO_model_stats));
		amigo_problem->amigo_models[ith_model]->amigo_model_stats->npars=0;
		amigo_problem->n_stored_fails=1;

	}else if(amigo_problem->n_stored_fails<=amigo_problem->n_max_store_fails){

		//Expand size
		amigo_problem->n_stored_fails++;
		amigo_problem->amigo_stats_containers=(AMIGO_model_stats**)
			realloc(amigo_problem->amigo_stats_containers,sizeof(AMIGO_model_stats*)*amigo_problem->n_stored_fails);
		//copy stats to the container
		amigo_problem->amigo_stats_containers[amigo_problem->n_stored_fails-1]=amigo_problem->amigo_models[ith_model]->amigo_model_stats;
		//put new stats back in the modle
		amigo_problem->amigo_models[ith_model]->amigo_model_stats=(AMIGO_model_stats*)malloc(sizeof(AMIGO_model_stats));
		amigo_problem->amigo_models[ith_model]->amigo_model_stats->npars=0;

	}
}

void set_AMIGO_problem_rhs(AMIGO_problem* amigo_problem, int(*rhs)(realtype,N_Vector, N_Vector, void*),void(change_y_func)(void*,realtype,N_Vector)){

	int i;

	for (i = 0; i < amigo_problem->n_models; ++i){

		amigo_problem->amigo_models[i]->rhs=rhs;
	}

	for (i = 0; i < amigo_problem->n_models; ++i){

		amigo_problem->amigo_models[i]->rhs=rhs;
		amigo_problem->amigo_models[i]->changeYatTcon=change_y_func;
	}
}

void set_AMIGO_problem_obs_function(AMIGO_problem* amigo_problem, void(*obs)(void*),void(*obs_sen_func)(void*)){

	int i;

	for (i = 0; i < amigo_problem->n_models; ++i){

		amigo_problem->amigo_models[i]->obs_func=obs;

		amigo_problem->amigo_models[i]->obs_sen_func=obs_sen_func;
		
	}
}
