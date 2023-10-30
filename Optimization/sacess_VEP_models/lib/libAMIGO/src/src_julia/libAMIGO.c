// julia.cpp : Defines the exported functions for the DLL application.
//






/*
EXPORTIT AMIGO_model* test_julia(
				int n_states,		int n_observables,	int n_pars,	
				int n_opt_pars,		int n_times,		int n_opt_ics,		
				int n_controls,		int n_controls_t,	int exp_num);


EXPORTIT int julia_get_num_AMIGO_opt_pars(void* data){
	
	int nx;
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	nx=amigo_problem->nx;
	return(nx);
}

EXPORTIT double* julia_get_AMIGO_problem_LB(void* data){
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	return(amigo_problem->LB);
}

EXPORTIT double* julia_get_AMIGO_problem_UB(void* data){
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	return(amigo_problem->UB);
}

EXPORTIT double* julia_get_AMIGO_problem_x(void* data){
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	return(amigo_problem->x);
}

EXPORTIT void julia_set_AMIGO_problem_x(void* data, double* x){
	
	int i;
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	set_AMIGO_problem_pars(x,amigo_problem);
}


EXPORTIT double julia_eval_lsq_AMIGO_model(void* data){
	
	double f;
	
	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	f=eval_AMIGO_problem_LSQ(amigo_problem);

	return(f);
}

*/

