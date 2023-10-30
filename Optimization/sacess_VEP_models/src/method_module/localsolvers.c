/**
 * @file localsolvers.c
 * @author David R. Penas
 * @brief File containing function interfaces to use N2SOL local solver with 
 * eSS.
 */


#include <stdio.h>
#include <stdlib.h>
#include <structure_paralleltestbed.h>
#include <limits.h>
#include <AMIGO_problem.h>
#include <common_solver_operations.h>
#include <output.h>
#include <benchmark_functions_SystemBiology.h>
#include <string.h>
#include <AMIGO_pe.h>
#include <time.h>
#include "simulate_amigo_model.h"
#include <parallel_functions_cooperative_eSS.h>
#include <string.h>
#include <AMIGO_pe.h>
#ifdef MPI2
#include <mpi.h>
#endif
#ifdef OPENMP
#include <omp.h>
#endif

/**
  * @brief this method calls to fortran subroutine to run the local method nl2sol
  * @param double x0 is the initial point
  * @param gradient parameter specify what subroutine is called: dn2gb (with gradient) or dn2fb
  * @param NPROC is the number of parallel processors
  * @param nevals is the number of evaluations spent in the local solver
  * @param exp1 experiment_total struct, the main struct of program.
  * benchmark.
*/

void NL2SOL_pe(double *x0, int gradient, long *nevals, void *exp1, long maxevals, long iter, int *NPROC,  void *(*fitnessfunction)(double*, void*)) {
    int i, counter;
    int *ui;
    int liv, lty, lv, n, p, one, printLevel;
    int* IV;
    double *TY, *B, *X, *V;
    AMIGO_problem* amigo_problem;
    experiment_total *exp;
    int evalaux, nR;
    double *R;
    output_function *res;
    exp = (experiment_total *) exp1;
    amigo_problem = exp->amigo;

    one = 1;
    lty = 10;
    if ( (*exp).test.bench.use_amigo == 1 ) {
    	p = amigo_problem->nx;
    	n = amigo_problem->n_data;
    } else {
	p = (*exp).test.bench.dim; 
        res = (output_function *) fitnessfunction(x0, (void *) exp);
        if ( res->size_r > 0 ) {
           n = res->size_r;
        } else n=0;
        deallocateoutputfunction_(res,exp1);
    }
    printLevel = 0;
    //LIV GIVES THE LENGTH OF IV.  IT MUST BE AT LEAST 82 + 4*P.
    liv = 4 * p + 82;

    IV = (int*) malloc(liv * sizeof (int));
    for (i=0;i<liv;i++) IV[i]=0;
    
    TY = (double*) malloc(2 * sizeof (double));
    for (i=0;i<2;i++) TY[i]=0.0;
    
    //LV GIVES THE LENGTH OF V.  THE MINIMUM VALUE
    //FOR LV IS LV0 = 105 + P*(N + 2*P + 21) + 2*N
    lv = 105 + p * (n + 2 * p + 21) + 2 * n;

    V = (double*) malloc(lv * sizeof (double));
    for (i=0;i<lv;i++) V[i]=0.0;
    
    X = (double*) malloc(p * sizeof (double));
    for (i=0;i<p;i++) X[i]=0.0;
    
    ui = (int *) malloc(2 * sizeof (int));
    for (i=0;i<2;i++) ui[i]=0.0;
    
    //Set Options
    //DFAULT(iv,v);								//set defaults (OLD v2.2)

    //DIVSET(ALG, IV, LIV, LV, V) 
    divset_(&one, IV, &liv, &lv, V);

    IV[13] = IV[14] = 0; //no covar

    IV[16] = maxevals; //limit on fevals + gevals
    IV[17] = iter; //max iter

    IV[18] = 1; //no iteration printing
    IV[19] = 1; //no default printing
    IV[20] = 1; //no output unit printing
    IV[21] = 1; //no x printing
    IV[22] = 1; //no summary printing
    IV[23] = 1;
    //v[30] = fatol;   
    //v[31] = frtol;
    //V[31] = 1e-9f;
    //V[32] = 1e-9f;
    //MEX Options 

    //The circadian problem fails if you don't specify this
    V[33] = 1e-9f;

    V[41] = 0.0046; // V(DLTFDC) Default = MACHEP^1/3. error tolerance in CVODES 1e-7, 1e-7^-3=0.0046
    V[42] = 3.1623e-4; // V(DLTFDJ) Default = MACHEP^1/2  error tolerance in CVODES 1e-7  1e-7^-2=3.1623e-4

    ui[1] = printLevel;

    B = (double*) malloc(2 * p * sizeof (double));

    counter = 0;

    if ( (*exp).test.bench.use_amigo == 1 ) {    
    	for (i = 0; i < p; ++i) {
        	B[counter++] = amigo_problem->LB[i];
        	B[counter++] = amigo_problem->UB[i];
		X[i] = x0[i];
    	}
    } else {
        for (i = 0; i < p; ++i) {
                B[counter++] = exp->test.bench.min_dom[i];
                B[counter++] = exp->test.bench.max_dom[i];
                X[i] = x0[i];
        }
    }
  
    if ( (*exp).test.bench.use_amigo == 1 ) {
	evalaux = amigo_problem->nevals;
    } else {
    	(*exp).execution.inner_ls_evals=0;
    }
    if (gradient) {
        dn2gb_(&n, &p, X, B, IV, &liv, &lv, V, ui, TY, exp1, NPROC, fitnessfunction);
    } else {
        dn2fb_(&n, &p, X, B, IV, &liv, &lv, V, ui, TY, exp1, NPROC, fitnessfunction);
    }

    if ( (*exp).test.bench.use_amigo == 1 ) {
	    //Final Rnorm
	    amigo_problem->local_fbest = V[9]*2; //nl2sol uses 1/2 sum(resid^2)
	    //Save Status & Iterations
	    amigo_problem->local_flag = IV[0];
	    amigo_problem->local_niter = IV[30];
	    memmove(amigo_problem->xbest,X,p*sizeof(double));
	    *nevals = *nevals + (amigo_problem->nevals - evalaux) ;
	    (*exp).amigo = amigo_problem;

    } else {
	    *nevals = *nevals + exp->execution.inner_ls_evals;
    }

    memmove(x0,X,p*sizeof(double));
    
    free(B);
    free(X);
    free(IV);
    free(V);
    free(ui);
    free(TY);
    B = NULL;
    X = NULL;
    IV = NULL;
    V = NULL;
    ui = NULL;
    TY = NULL;
   
}


/**
 * @brief this method configures all the parameters to call nl2sol local solver
 * @param exp experiment_total struct, the main struct of program.
 * @param idp parallel identification number.
 * @return void function pointer of evaluation function of the selected
 * benchmark.
*/
void callnl2sol_(void * exp1_, void *(*fitnessfunction)(double*, void*),
         double *x0, long *neval, long *maxevals, long *iter, int *NPROC, int *opt) {
    experiment_total *exp1;
    int i, D;
    void *request_recep;
    AMIGO_problem* amigo_problem;

    exp1 = (experiment_total *) exp1_;
    amigo_problem = exp1->amigo;

    if ( exp1->test.bench.use_amigo == 1 ) {
        set_AMIGO_problem_pars(x0, amigo_problem);
        memmove(amigo_problem->x0,x0,amigo_problem->nx*sizeof(double));
    }
    exp1 = (experiment_total *) exp1_;
    NL2SOL_pe(x0, *opt, neval, exp1, *maxevals, *iter, NPROC,fitnessfunction);

}



int calcr_(int *n, int *p, double *x, int *nf, double *r__, int *lty, double *ty, void* exp1, int *stop, int *NPROC, void *(*fitnessfunction)(double*, void*)) {

    int nr, i;
    experiment_total *exp;
    exp = (experiment_total *) exp1;
    output_function *res;

    if  (exp->test.bench.translation == 1) {
        detranslation( x, exp->execution.transconst, exp->test.bench.dim );
    }

    if ( (*exp).test.bench.use_amigo == 1 ) { 
	calcramigo_(n,p,x,nf,r__,lty,ty,exp1,stop,NPROC);
    } else {
	(*exp).execution.inner_ls_evals++;
        res = (output_function *) fitnessfunction(x, (void *) exp);
        if ( res->size_r > 0 ) {
	   for (i=0;i<res->size_r;i++) {
		r__[i]=res->R[i];
	   }	
	}
        deallocateoutputfunction_(res,exp1);
    }

    if  (exp->test.bench.translation == 1) translation( x, exp->execution.transconst, exp->test.bench.dim );

    return 1;
}


int calcj_(int *n, int* p, double *x, int *nf, double *dr__, 
	int *lty, double *ty, void*exp1, int *stop, int *NPROC, void *(*fitnessfunction)(double*, void*)) {
    experiment_total *exp;
    exp = (experiment_total *) exp1;
    int ndr,i;
    output_function *res;
    if ( (*exp).test.bench.use_amigo == 1 ) {
        calcjamigo_(n,p,x,nf,dr__,lty,ty,exp1,stop,NPROC);
    } else {
	exp->execution.inner_ls_evals++;
        res = (output_function *) fitnessfunction(x, (void *) exp);
        if ( res->size_j > 0 ) {
           for (i=0;i<res->size_j;i++) {
                dr__[i]=res->J[i];
           }
        }
        deallocateoutputfunction_(res,exp1);
    }
    return 1;
}



void initlocalsolvervarsess_(void* exp1_){
    experiment_total *exp1;
    local_solver *local_s;

    exp1 = (experiment_total *) exp1_;    
    local_s = exp1->ls;
    local_s->counter = 0;
    local_s->total_local = 0.0;
    local_s->sucess_local = 0.0;
    local_s->state_ls = 0;
    local_s->sucess_interval = 0;
    local_s->num_ls = 0;
}


void addlocalscounteress_(void* exp1_){
    experiment_total *exp1;
    local_solver *local_s;

    exp1 = (experiment_total *) exp1_;    
    local_s = exp1->ls;
    local_s->counter++;
}


void setentervalueess_( void *exp1_, double *fval){
    experiment_total *exp1;
    local_solver *local_s;

    exp1 = (experiment_total *) exp1_;    
    local_s = exp1->ls;
    
    (*local_s).enter_value = *fval;
}
