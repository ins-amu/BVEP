#pragma once 

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_dense.h>  
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <nvector/nvector_serial.h>/* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_math.h> 
#include <amigoRHS.h>
#include <simulate_amigo_model.h>
#include <math.h>
#include <stddef.h>
#include <AMIGO_model_stats.h>



#if defined WIN32 || defined __WIN32 || defined _WIN32 || defined _WIN64 || defined WIN64 || defined MSVC || defined win32 || defined _win32 || defined __win32 

	#if defined EXPORT
		#define EXPORTIT __declspec(dllexport) 
	#else
		#define EXPORTIT __declspec(dllexport) 
	#endif

	#if defined IMPORT
		#define EXPORTIT __declspec(dllimport) 
	#endif
#else
	#define EXPORTIT

#endif	



#ifdef _MSC_VER
	#define INFINITY (DBL_MAX+DBL_MAX)
	#define NAN (INFINITY-INFINITY)
#endif

#ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef  min
        #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif 

#ifndef fmax
	 #define fmax( a, b )  (max(a,b))
#endif 

#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

#ifdef AMIGO_OBS_STATES
	void amigo_OBS_states(void* data);
#endif

#ifdef AMIGO_OBS_SENS
	void amigo_OBS_sens(void* data);
#endif

typedef struct
{
	//States and observables
	int n_states;
	int n_observables;
	int* index_observables;
	int exp_num;

	//Simulation Pars
	int n_pars;
	int n_opt_pars;
    int n_total_x;
	double* pars;
	double* pars_guess;
	double* pars_LB;
	double* pars_UB;

	int* index_opt_pars;
	
	//Simulation times
	double t0;
	double tf;
	int n_times;
	double* t;
	
	//Initial conditions
	double* y0;
	int n_opt_ics;
	int* index_opt_ics;
	double* y0_guess;
	double* y0_LB;
	double* y0_UB;

	//Contols
	int n_controls;
	int n_controls_t;
	double* controls_t;
	double** controls_v;
	double** slope;
	double tlast;
	int index_t_stim;

	//Function
	int(*rhs)(realtype,N_Vector, N_Vector, void*);
	int(*jac)(int, realtype,N_Vector, N_Vector, DlsMat, void*, N_Vector, N_Vector, N_Vector);
	void (*jY)(int);
	void (*jS)(int);
	void (*obs_func)(void*);
	void (*obs_sen_func)(void*);
	void (*changeYatTcon)(void*, realtype, N_Vector);


	//Storing matrixes
	double**  sim_results;
	double*** sens_results;
	double**  obs_results;
	double*** sens_obs;

	//Experimental Data
	double** exp_data;
	double** Q;
	double* w_obs;

	//Simulation Related Parameter
	double reltol;
	double atol;
	double max_step_size;
	int max_num_steps;
	int max_error_test_fails;
	int use_jacobian;
	int compute_sens;
	double mkl_tol;
    
    int use_obs_func;
    int use_sens_obs_func;
	

	AMIGO_model_stats* amigo_model_stats;

	//Extra Data for RHS
	void* data;

}AMIGO_model;

EXPORTIT AMIGO_model* allocate_AMIGO_model
(
	int n_states,		int n_observarbles,	int n_pars,	
	int n_opt_pars,		int n_times,		int n_opt_ics,		
	int n_controls,		int n_controls_t,	int exp_num
);

EXPORTIT void free_AMIGO_model(AMIGO_model* amigo_model);

EXPORTIT int simulate_amigo_model(AMIGO_model* amigo_model,int verbose);

EXPORTIT int simulate_AMIGO_model_observables(AMIGO_model* amigo_model,int sens);

EXPORTIT double lsq_AMIGO_model(AMIGO_model* amigo_model);

EXPORTIT extern void dgemm(
    char   *transa,
    char   *transb,
    ptrdiff_t *m,
    ptrdiff_t *n,
    ptrdiff_t *k,
    double *alpha,
    double *a,
    ptrdiff_t *lda,
    double *b,
    ptrdiff_t *ldb,
    double *beta,
    double *c,
    ptrdiff_t *ldc
);

EXPORTIT int get_AMIGO_model_sens(AMIGO_model* amigo_model,int cvodes, int mkl);

