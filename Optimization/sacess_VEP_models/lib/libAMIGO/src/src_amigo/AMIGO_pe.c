#include <AMIGO_pe.h>
#include <time.h>
#include <structure_paralleltestbed.h>
#include "simulate_amigo_model.h"
#include <parallel_functions_cooperative_eSS.h>
#include <string.h>
#ifdef MPI2
#include <mpi.h>
#endif
#ifdef OPENMP
#include <omp.h>
#endif

int calcramigo_(int *n, int *p, double *x, int *nf, double *r__,
        int *lty, double *ty, void* exp1, int *stop, int *NPROC) {

    int i, ii, j, k, n_exps, n_times, n_obs, flag;
    int counter = 0;
    double sum = 0;
    int dest, test;
    experiment_total *exp;
    exp = (experiment_total *) exp1;
    AMIGO_problem* amigo_problem = exp->amigo;
    int dist_criteria;

    getdistcriteria_(exp1,&dist_criteria);

    amigo_problem->nevals++;
    amigo_problem->local_nfeval++;
    set_AMIGO_problem_pars(x, amigo_problem);
    n_exps = amigo_problem->n_models;
    
    for (i = 0; i < n_exps; ++i) {
#ifdef MPI2
      if (dist_criteria == 1 ) {
        for (j=0;j<*NPROC;j++){
                dest=j;
                test = cooperativempitestess_(exp1, &dest );
                if ( test == 1 ) {
                     *stop = 1;
                     break;
                }
        }
      }
#endif
    if (*stop == 1) {
         break;
    }

    flag = simulate_AMIGO_model_observables(amigo_problem->amigo_models[i], 0);
    n_obs = amigo_problem->amigo_models[i]->n_observables;

    if (!flag) {
       handle_AMIGO_problem_stat_fails(i, amigo_problem);
       for (ii = 0; ii < n_exps; ++ii) {
         for (j = 0; j < n_obs; ++j) {
           n_times = amigo_problem->amigo_models[ii]->n_times;
           for (k = 0; k < n_times; ++k) {
              r__[counter++] = DBL_MAX;
           }
         }
       }
    }
    for (j = 0; j < n_obs; ++j) {
       n_times = amigo_problem->amigo_models[i]->n_times;
       for (k = 0; k < n_times; ++k) {
         if (!isnan(amigo_problem->amigo_models[i]->exp_data[j][k]) &&
              !isnan(amigo_problem->amigo_models[i]->Q[j][k]) &&
                amigo_problem->amigo_models[i]->Q[j][k] != 0) {
          	r__[counter++] =  (amigo_problem->amigo_models[i]->obs_results[j][k] -  amigo_problem->amigo_models[i]->exp_data[j][k]) *  amigo_problem->amigo_models[i]->Q[j][k];
                    sum += pow(r__[counter - 1], 2);
          } else {
                r__[counter++] = 0;
          }
       }
    }
    }
    return (0);
}

int calcjamigo_(int *n, int* p, double *x, int *nf, double *dr__, int *lty, double *ty, void*exp1, int *stop, int *NPROC) {

    int i, j, k, m, ii, n_exps, n_times, n_obs;
    int counter = 0;
    int n_ics = 0, count_ic = 0, flag;
    double t1,t2,t3,t4;
    int dest, test;

    experiment_total *exp;
    exp = (experiment_total *) exp1;
    AMIGO_problem* amigo_problem = exp->amigo;
    int dist_criteria;

    getdistcriteria_(exp1,&dist_criteria);

    amigo_problem->nevals++;
    amigo_problem->local_nfeval++;
    set_AMIGO_problem_pars(x, amigo_problem);
    n_exps = amigo_problem->n_models;
    t3 = clock();
    for (i = 0; i < n_exps; ++i) {
#ifdef MPI2
      if (dist_criteria == 1 ) {
        for (j=0;j<*NPROC;j++){
                dest=j;
                test = cooperativempitestess_(exp1, &dest );
                if ( test == 1 ) { 
		     *stop = 1;
		     break;
		}
        }
      }
#endif
        if (*stop == 1) {
            break;
        }
#ifdef MKL
        flag = get_AMIGO_model_sens(amigo_problem->amigo_models[i], amigo_problem->cvodes_gradient, amigo_problem->mkl_gradient);
#else
        flag = get_AMIGO_model_sens(amigo_problem->amigo_models[i], 1, 0);
#endif
        if (!flag) {
            handle_AMIGO_problem_stat_fails(i, amigo_problem);
            counter = 0;
            for (m = 0; m < *p; ++m) {
                for (ii = 0; ii < n_exps; ++ii) {
                    n_obs = amigo_problem->amigo_models[ii]->n_observables;
                    for (j = 0; j < n_obs; ++j) {
                        n_times = amigo_problem->amigo_models[ii]->n_times;
                        for (k = 0; k < n_times; ++k) {
                            dr__[counter++] = 0;
                        }
                    }
                }
            }
            if (amigo_problem->verbose) {
#ifdef MATLAB
                (int) mexPrintf("Gradient failed\n");
#else
                printf("Gradient failed\n");
#endif

            }

            return (0);
        }
    }
     t4 =clock();

    for (m = 0; m < *p; ++m) {
        n_ics = amigo_problem->n_pars;
        if (*stop >= 1) {
            break;
        }
        
        for (i = 0; i < n_exps; ++i) {

            n_obs = amigo_problem->amigo_models[i]->n_observables;

            if (m < amigo_problem->n_pars) {

                for (j = 0; j < n_obs; ++j) {


                    n_times = amigo_problem->amigo_models[i]->n_times;

                    for (k = 0; k < n_times; ++k) {

                        if (!isnan(amigo_problem->amigo_models[i]->exp_data[j][k]) &&
                                !isnan(amigo_problem->amigo_models[i]->Q[j][k]) &&
                                amigo_problem->amigo_models[i]->Q[j][k] != 0) {

                            dr__[counter++] =
                                    amigo_problem->amigo_models[i]->sens_obs[j][m][k]*
                                    (amigo_problem->amigo_models[i]->Q[j][k]);
                        } else {
                            dr__[counter++] = 0;

                        }
                    }

                }

            } else if (m >= amigo_problem->n_pars && m >= n_ics && m < n_ics + amigo_problem->amigo_models[i]->n_opt_ics) {

                for (j = 0; j < n_obs; ++j) {

                    n_times = amigo_problem->amigo_models[i]->n_times;

                    for (k = 0; k < n_times; ++k) {
                        if (!isnan(amigo_problem->amigo_models[i]->exp_data[j][k]) &&
                                !isnan(amigo_problem->amigo_models[i]->Q[j][k]) &&
                                amigo_problem->amigo_models[i]->Q[j][k] != 0) {

                            dr__[counter++] =
                                    (amigo_problem->amigo_models[i]->sens_obs[j][m - n_ics + amigo_problem->n_pars][k])*
                                    (amigo_problem->amigo_models[i]->Q[j][k]);
                        } else {
                            dr__[counter++] = 0;
                        }

                    }
                }
                count_ic++;
            } else {
                for (j = 0; j < n_obs; ++j) {

                    n_times = amigo_problem->amigo_models[i]->n_times;

                    for (k = 0; k < n_times; ++k) {
                        dr__[counter++] = 0;
                    }
                }
            }
            n_ics += amigo_problem->amigo_models[i]->n_opt_ics;
        }
    }
    return 0;
}

