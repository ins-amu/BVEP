#pragma once

#include <AMIGO_model.h>
#include "float.h"

#ifdef MKL
	#include <mkl.h>
#endif

#ifdef MATLAB
	#include "mex.h"
#endif

#define AMIGO_PI 3.141592653589793238462643

#ifdef _WIN32
#  define isnan _isnan 
#endif


typedef struct
{
	int n_models;
	AMIGO_model** amigo_models;
	
	int nx;
	int n_ics;
	int n_pars;
        int n_exp;
	int n_data;

	double* x;
	double* x0;		
	double* xbest;
	double* LB;
	double* UB;

	int use_gradient;
	int cvodes_gradient;
	int mkl_gradient;

	int n_fails;
	int n_stored_fails;
	int n_max_store_fails;
	AMIGO_model_stats** amigo_stats_containers;

	int nevals;
	double fbest;
	double temp_min;

	int local_flag;
	int local_nfeval;
	int local_niter;
	double local_fbest;
	int local_max_iter;
	int local_max_evals;

	int verbose;

	int nthreads;

	void* data;

}AMIGO_problem;

EXPORTIT AMIGO_problem* openMatFileAMIGO(const char* file);

EXPORTIT AMIGO_problem* allocate_AMIGO_problem(int n_models, AMIGO_model** amigo_models);

EXPORTIT void free_AMIGO_problem(AMIGO_problem* amigo_problem);

EXPORTIT int simulate_AMIGO_model_observables(AMIGO_model* amigo_model,int sens);

EXPORTIT void set_AMIGO_problem_pars(double* x, AMIGO_problem* amigo_problem);

EXPORTIT void handle_AMIGO_problem_stat_fails(int ith_model,AMIGO_problem* amigo_problem);

EXPORTIT double eval_AMIGO_problem_LSQ(AMIGO_problem* amigo_problem);

EXPORTIT double eval_AMIGO_problem_LLK(AMIGO_problem* amigo_problem);

EXPORTIT void set_AMIGO_problem_rhs(AMIGO_problem* amigo_problem, int(*rhs)(realtype,N_Vector, N_Vector, void*),void(change_y_func)(void*,realtype,N_Vector));

EXPORTIT double AMIGO_dummy(AMIGO_problem* amigo_problem);


EXPORTIT void set_AMIGO_problem_obs_function(AMIGO_problem* amigo_problem, void(*obs)(void*),void(*obs_sen_func)(void*));

EXPORTIT AMIGO_problem* openMatFileAMIGO(const char* file);