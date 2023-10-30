#include <stdlib.h>
#include <string.h>
#include <structure_paralleltestbed.h>
#include <common_solver_operations.h>
#include <AMIGO_problem.h>
#include <configuration.h>
#include <ggn.h>
#include <amigoRHS_CIRCADIAN.h>
#include <amigoRHS_3step.h>
#include <amigoRHS_NFKB.h>
#include <amigoRHS_B1.h>
#include <amigoRHS_B2.h>
#include <amigoRHS_B3.h>
#include <amigoRHS_B4.h>
#include <amigoRHS_B5.h>
#include <math.h>
#include <input_AMIGO.h>
#include <hdf5.h>
#include <input_module.h>
#include <def_errors.h>

#define PI (3.141592653589793238462643383279)


int load_benchmark_test(experiment_total *exp, double *target, double *ftarget) {
    int D,i;
    double Up, Lo;
    
    
    D = (*exp).test.bench.dim;
    
    Up = exp->test.bench.max_dom[0]; 
    Lo = exp->test.bench.min_dom[0];
    
    exp->test.bench.max_dom = (double *) malloc(D* sizeof (double) );
    exp->test.bench.min_dom = (double *) malloc(D* sizeof (double) );
    
    for (i=0;i<D;i++) {
                exp->test.bench.max_dom[i] = Up;
                exp->test.bench.min_dom[i] = Lo;
    }
    
    *target = -100000;        
    *ftarget = -100000;
    
    return 1;
}





double fgeneric_SB(experiment_total *exp, double tolerance) {
    double ftarget;
    ftarget = eval_AMIGO_problem_LSQ(exp->amigo);

    return ftarget + tolerance;
}




const char* return_benchmark_SystemBiology(int i) {
    if (i == 0) {
        return "circadian";
    } else if (i == 1) {
        return "mendes";
    } else if (i == 2) {
        return "nfkb";
    } else if (i == 3) {
        return "B1";
    } else if (i == 4) {
        return "B2";
    } else if (i == 5) {
        return "B3";        
    } else if (i == 6) {
        return "B4";        
    } else if (i == 7) {
        return "B5";        
    } else if (i == 8) {
        return "B6";        
    } else {
        return "";
    }
}

char * getnameSB(int id){
    char * name;

    name = (char) calloc(500,sizeof(char) );
    if (id == 0) name = "Circadian problem";
    else if (id == 1) name = "3-step-pathway problem";
    else if (id == 2) name = "NfkB problem";
    else if (id == 3) name = "B1 problem";
    else if (id == 4) name = "B2 problem";
    else if (id == 5) name = "B3 problem";
    else if (id == 6) name = "B4 problem";
    else if (id == 7) name = "B5 problem";
    else if (id == 8) name = "B6 problem";
    else  name = "SB problem";


    return name;
}


int load_benchmark_SystemBiology(experiment_total *exp) {
    const char *path;
    int i,j;
    int type;
    int amigo_flag;
    int counter, exito, init_cond;
    double point;
    int *index_non_obs;
    type = -1;
    amigo_flag=1;
   
    (*exp).test.bench.estime_init_cond = 0; 
    if (exp->test.bench.current_bench == 0) {
        path = "benchmarks/systemsBiology/others/circadian/load_C.mat";
        type = 0;
        (*exp).test.VTR_default = 1e-5;
        (*exp).test.jf = 1e-5;
        
    } else if (exp->test.bench.current_bench == 1) {
        path = "benchmarks/systemsBiology/others/3-step_pathway/load_C.mat";
        type = 1;
        (*exp).test.VTR_default = 1e-5;
        (*exp).test.jf = 1e-5;
        
    } else if (exp->test.bench.current_bench == 2) {
        path = "benchmarks/systemsBiology/others/Nfkb/load_C.mat";
        type = 2;
        (*exp).test.VTR_default = 1e-2;
        (*exp).test.jf = 1e-2;
    } else if (exp->test.bench.current_bench == 3) {
        path = "benchmarks/systemsBiology/BioPredyn/B1/load_C.mat";
        type = 3;
        (*exp).test.VTR_default = 13753;
        (*exp).test.jf = 13753;
    }  else if (exp->test.bench.current_bench == 4) {
        path = "benchmarks/systemsBiology/BioPredyn/B2/load_C.mat";
        type = 4;
	(*exp).test.VTR_default = 250;
        (*exp).test.jf = 250;
    } 
    else if (exp->test.bench.current_bench == 5) {
        path = "benchmarks/systemsBiology/BioPredyn/B3/load_C.mat";
        type = 5;
        (*exp).test.VTR_default = 0.37029;
	(*exp).test.jf = 0.37029;
    } 
    else if (exp->test.bench.current_bench == 6) {
        path = "benchmarks/systemsBiology/BioPredyn/B4/load_C.mat";
        type = 6;
        (*exp).test.VTR_default = 55.0;
        (*exp).test.jf = 55.0;
    }  
    else if (exp->test.bench.current_bench == 7) {
        path = "benchmarks/systemsBiology/BioPredyn/B5/load_C.mat";
        type = 7;
        (*exp).test.jf = 4200;
        (*exp).test.VTR_default = 4200;
    } else if (exp->test.bench.current_bench == 8) {
        type = 8;
        amigo_flag=0;
        (*exp).test.VTR_default = 108330;
        (*exp).test.jf = 108330;
    } else {
        perror(error14);
        exit(14);
    }
    exp->test.bench.use_amigo=amigo_flag;
    exp->test.bench.idtype=type; 
    exp->param = NULL;
    exp->amigo = NULL;

// OPEN MAT CONDITIONAL 
    if (amigo_flag == 1){ 
    	exp->amigo =  openMatFileAMIGO(path); 
        
   	if (type == 0) {
        	    set_AMIGO_problem_rhs(exp->amigo, amigoRHS_CIRCADIAN, amigo_Y_at_tcon_CIRCADIAN);
       		    set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_CIRCADIAN, amigoRHS_get_sens_OBS_CIRCADIAN);
   	} else if (type == 1) {
        	    set_AMIGO_problem_rhs(exp->amigo, amigoRHS_MENDES, amigo_Y_at_tcon_MENDES);
        	    set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_MENDES, amigoRHS_get_sens_OBS_MENDES);
   	} else if (type == 2) {
        	    set_AMIGO_problem_rhs(exp->amigo, amigoRHS_NFKB, amigo_Y_at_tcon_NFKB);
        	    set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_NFKB, amigoRHS_get_sens_OBS_NFKB);    
    	} else if (type == 3) {
            	    set_AMIGO_problem_rhs(exp->amigo, amigoRHS_B1, amigo_Y_at_tcon_B1);
            	    set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_B1, amigoRHS_get_sens_OBS_B1);    
    	} else if (type == 4) {
	            set_AMIGO_problem_rhs(exp->amigo, amigoRHS_B2, amigo_Y_at_tcon_B2);
	            set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_B2, amigoRHS_get_sens_OBS_B2);   
    	} else if (type == 5) {
	            set_AMIGO_problem_rhs(exp->amigo, amigoRHS_B3, amigo_Y_at_tcon_B3);
	            set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_B3, amigoRHS_get_sens_OBS_B3);
    	} else if (type == 6) {
	            set_AMIGO_problem_rhs(exp->amigo, amigoRHS_B4, amigo_Y_at_tcon_B4);
	            set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_B4, amigoRHS_get_sens_OBS_B4);   
    	} else if (type == 7) {
	            set_AMIGO_problem_rhs(exp->amigo, amigoRHS_B5, amigo_Y_at_tcon_B5);
	            set_AMIGO_problem_obs_function(exp->amigo, amigoRHS_get_OBS_B5, amigoRHS_get_sens_OBS_B5);
  	}
        if ( (*exp).test.bench.estime_init_cond == 1) {
            init_cond = exp->amigo->amigo_models[0]->n_states - exp->amigo->amigo_models[0]->n_observables;
            exp->test.bench.dim = init_cond + exp->amigo->nx;
            index_non_obs = (int *) malloc(  init_cond * sizeof(int) );
            counter = 0;
            for (i=0;i<exp->amigo->amigo_models[0]->n_states;i++){
                exito=0;
                for (j=0;j<exp->amigo->amigo_models[0]->n_observables;j++) {
                    if (exp->amigo->amigo_models[0]->index_observables[j] == i){
                        exito = 1;
                        break;
                    }
                }
                if (exito == 0){
                    index_non_obs[counter]=i;
                    counter++;
                }
            }
            
            /*exp->test.bench.max_dom = (double *) malloc(exp->test.bench.dim * sizeof (double));
            exp->test.bench.min_dom = (double *) malloc(exp->test.bench.dim * sizeof (double));            
            for (i = 0; i < init_cond; i++) {
                exp->test.bench.max_dom[i] = 1.0;
                exp->test.bench.min_dom[i] = 0.0;
            }
            */
            counter = 0;
            for (i = init_cond; i < exp->test.bench.dim; i++) {

                exp->test.bench.max_dom[i] = exp->amigo->UB[counter];
                exp->test.bench.min_dom[i] = exp->amigo->LB[counter];
                counter++;
            }    
            
            free(index_non_obs);
            
        } 
	else {
            //exp->test.bench.dim = exp->amigo->nx;
            //exp->test.bench.max_dom = (double *) malloc(exp->test.bench.dim * sizeof (double));
            //exp->test.bench.min_dom = (double *) malloc(exp->test.bench.dim * sizeof (double));
           //for (i = 0; i < exp->test.bench.dim; i++) {
           //   exp->test.bench.max_dom[i] = exp->amigo->UB[i];
           //   exp->test.bench.min_dom[i] = exp->amigo->LB[i];
           //}
       }
    } else { 
    	if ( type == 8) {
           // exp->test.bench.BestSol = (double *) malloc( ( exp->test.bench.dim + 1) * sizeof (double));
           // exp->test.bench.dim = 37;
           // exp->test.bench.max_dom = (double *) malloc( exp->test.bench.dim * sizeof (double));
           // exp->test.bench.min_dom = (double *) malloc( exp->test.bench.dim * sizeof (double));
           //returnbounds( exp->test.bench.max_dom , exp->test.bench.min_dom );
    	}
    }

    return 1;
    
}

void manage_init_cond(experiment_total *exp,double *U, double *U_aux) {
    int n_IC, counter;
    int i, j, exito, *index_non_obs;
    
    n_IC = exp->amigo->amigo_models[0]->n_states - exp->amigo->amigo_models[0]->n_observables;
    index_non_obs = (int *) malloc(n_IC * sizeof (int));
    counter = 0;
    for (i = 0; i < exp->amigo->amigo_models[0]->n_states; i++) {
        exito = 0;
        for (j = 0; j < exp->amigo->amigo_models[0]->n_observables; j++) {
            if (exp->amigo->amigo_models[0]->index_observables[j] == i) {
                exito = 1;
                break;
            }
        }
        if (exito == 0) {
            index_non_obs[counter] = i;
            counter++;
        }
    }

    for (i=0;i<exp->amigo->n_exp;i++){
        for (j=0;j<n_IC;j++) {
            exp->amigo->amigo_models[i]->y0[index_non_obs[j]]=U[j];
        }
    }
    
    
    
    counter=0;
    for (i=n_IC;i<exp->test.bench.dim;i++){
        U_aux[counter]=U[i];
        counter++;
    }    
}


void repack_init_cond(experiment_total *exp,double *U, double *U_aux) {
    int n_IC,i,max,counter;
    
    n_IC = exp->amigo->amigo_models[0]->n_states - exp->amigo->amigo_models[0]->n_observables;
    max = exp->test.bench.dim;
    counter=0;
    for (i=n_IC; i<max ;i++){
        U[i]=U_aux[counter];
        counter++;
    }        
}


int destroySystemBiology(experiment_total *exp) {

    if (exp->amigo != NULL) {
        free_AMIGO_problem(exp->amigo);
        exp->amigo=NULL;
    }
    free(exp->execution.transconst);
    free(exp->test.bench.logindex);
    free(exp->test.bench.log_max_dom );
    free(exp->test.bench.log_min_dom );


    return 1;

}


void* evalSB_(double *U, void* data) {
    experiment_total *exp1;
    AMIGO_problem* amigo_problem;
    output_function *res;
    double *U_aux;

    exp1 = (experiment_total *) data;
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));

    if (  exp1->test.bench.use_amigo == 1 ) { 
        if (exp1->test.bench.estime_init_cond == 1) {
            U_aux = (double *) malloc(exp1->amigo->nx*sizeof(double));
            manage_init_cond(exp1, U, U_aux);
            amigo_problem = exp1[0].amigo;            
            set_AMIGO_problem_pars(U_aux, amigo_problem);
            res->value = eval_AMIGO_problem_LSQ(amigo_problem);
            amigo_problem->nevals++;    
            repack_init_cond(exp1, U, U_aux);

            free(U_aux);
        } else {
            amigo_problem = exp1[0].amigo;     
            set_AMIGO_problem_pars(U, amigo_problem);
            res->value = eval_AMIGO_problem_LSQ(amigo_problem);
            amigo_problem->nevals++;
        }
        
        
    } else {
	if (exp1->test.bench.idtype == 8) res->value = fitnessfunctionB6(U); 
    }

    return ((void *) res);

}




