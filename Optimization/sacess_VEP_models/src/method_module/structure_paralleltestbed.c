/**
 * @file structure_paralleltestbed.c
 * @author David R. Penas
 * @brief File containing setter and getter for all fields of the main struct
 * experiment_total.
 */


#include "structure_paralleltestbed.h"
#include <configuration.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <def_errors.h>
#include <configuration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>  
#include <string.h>

int returninitsol_( void *exp_ ){
    experiment_total *exp;

    exp = (experiment_total *) exp_;
    return  exp[0].test.bench.number_init_sol;
}

int getuseamigo_(void * exp_) {
    experiment_total *exp; 

    exp = (experiment_total *) exp_;
    return exp->test.bench.use_amigo;
}

void setdistcriteria_(void *exp_, int *dist) {
    experiment_total *exp;

    exp = (experiment_total *) exp_;
    exp->execution.dist_stopping_criteria= *dist;    	
}

void getdistcriteria_(void *exp_,  int *dist) {
    experiment_total *exp;

    exp = (experiment_total *) exp_;
    *dist = exp->execution.dist_stopping_criteria;
}


void loadinitpoints_( void *exp_, double *point, int *pos, double *fx ){
    experiment_total *exp;
    int i;
    exp = (experiment_total *) exp_;

    for (i=0;i<exp[0].test.bench.dim;i++){
        point[i] = exp[0].test.bench.X0[*pos][i];
    }

    fx[*pos]=exp[0].test.bench.F0[*pos];
}



void returnboundscfortran_(void *exp_, double *XU, double *XL, int *nvar) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp_;

    for (i=0;i<*nvar;i++) {
                XU[i] = exp1[0].test.bench.max_dom[i];
                XL[i] = exp1[0].test.bench.min_dom[i];
    }
}

void setvectorcfortran_(void *vector, double *in, int *nvar){
   double *aux;
   aux = (double *) vector;
   int i;

   for (i=0;i<*nvar;i++) {
          aux[i] = in[i];
   }


}


void returnvectorcfortran_(void *vector,double *out, int *nvar){
   double *aux;
   aux = (double *) vector;
   int i;

   for (i=0;i<*nvar;i++) {
          out[i] = aux[i];
   }
   

}


void sumfailevals_(void *exp_){
    experiment_total *exp;
    exp = (experiment_total *) exp_;    
    
    exp->execution.failevals = exp->execution.failevals + 1;

}

int returnfailevals_(void *exp_){
    experiment_total *exp;
    exp = (experiment_total *) exp_;    
    
    return exp->execution.failevals;
}

void familyslave_(void *exp1, int *option){
     experiment_total *exptotal;
     exptotal = (experiment_total *) exp1;
     exptotal->execution.family_slave = *option;
}

int init_argc_struct(argc_struct *arg){

arg->n_stuck= -1; 
arg->evals_threshold=-1;  
arg->mult_num_sendSol=-1;  
arg->minimum_num_sendSol=-1; 
arg->evalmax=-1; 

return 1;
}



int create_expetiment_struct(const char *file, experiment_total *exptotal, int NPROC, int id,  const char *path, int init, argc_struct arg){
    int error, counter;
	exptotal->test.bench.max_dom = NULL;
	exptotal->test.bench.min_dom = NULL;
	exptotal->test.bench.CL = NULL;
	exptotal->test.bench.CU = NULL;
	exptotal->methodDE = NULL;
	exptotal->methodScatterSearch = NULL;
	exptotal->par_st = NULL;
	exptotal->param = NULL;
	exptotal->amigo = NULL;
	exptotal->ls = NULL;
	exptotal->output = NULL;
	exptotal->test.output_graph=NULL;
	(*exptotal).test.output_path=NULL;
	(*exptotal).test.log_output = NULL;
	(*exptotal).test.log_percentage = NULL;
	(*exptotal).test.output_gant_log = NULL;
	(*exptotal).test.result_output = NULL;


    exptotal->execution.idp = id;
    exptotal->execution.NPROC = NPROC; 
    exptotal->execution.file = (char *) calloc(300,sizeof(char));
    strcpy(exptotal->execution.file,file);   
   
    exptotal->test.output_graph = (char *) calloc(1000, sizeof(char));
    strcpy(exptotal->test.output_graph,path);
    exptotal->test.output_path = (char *) calloc(1000, sizeof(char));
    
     
    if (id == 0) {
      error = -1;
      counter = 0;
      if (init==1) {
        while (error == -1) {
            strcpy(exptotal->test.output_path,exptotal->test.output_graph);
            error = mkdir( (const char *) exptotal->test.output_graph,0777);
            if ((error != -1)&&(error != 0)) {
                exit(EXIT_FAILURE);
            } else if (error == -1){
                counter++;
                sprintf(exptotal->test.output_graph, "%s_%d", path, counter);
	        strcpy(exptotal->test.output_path,exptotal->test.output_graph);
            }
        }
      }
    }
    
#ifdef MPI2
    MPI_Bcast(exptotal->test.output_graph, 1000, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
    if (load_configuration_XML(file, exptotal) == 0 ) {
        if (id == 0)
            printf(" Error in load configuration\n\n");
        exit(0);
    } 
        
    exptotal->execution.failevals = 0;
    exptotal->execution.initpath = 1;
    exptotal->execution.rep = NPROC + 10;
    exptotal->test.jfprint = 0;
    exptotal->ls = (local_solver *) malloc(sizeof (local_solver));
    exptotal->output = (output_struct *) malloc(sizeof (output_struct));
    exptotal->amigo=NULL;

    if (arg.n_stuck != -1 )  exptotal->methodScatterSearch->goptions->n_stuck=arg.n_stuck;


    return 1; 
}

void init_result_data(result_solver *result, int size) {
        // init result var
    int i;
    
        result->eval_value = 0;
        result->time_value = 0.0;
        result->time_vtr_value = 0.0;
        result->best_value = DBL_MAX;
        result->iterations_value = 0.0;
        result->totaltime = 0.0;
        result->paralleltime = 0.0;
        result->localsolvertime = 0.0;
        result->bestx_value = (double *) malloc(size * sizeof(double));
        
}


void destroy_result_data(result_solver *result) {
    free(result->bestx_value);
}
    



void updatenp_(void *exp_, int *dim_refset){
    experiment_total *exp;
    exp = (experiment_total *) exp_;
    
    exp->methodScatterSearch->_NP = *dim_refset;
    
}

const char* getname(experiment_total *exp) {
    const char *retchar;
    if (exp->methodDE != NULL) retchar = "Differential Evolution";
    else if (exp->methodScatterSearch != NULL) retchar = "ScatterSearch";
    else retchar= "algorithm";

    return retchar;

}

int getnumversion(experiment_total *exptotal){
    
    if (strcmp(exptotal->methodScatterSearch->eSSversion, "ScatterSearch") == 0) {
        return 0;
    } 
    else if (strcmp(exptotal->methodScatterSearch->eSSversion, "CeSS") == 0) {
        return 1;
    }
    else if (strcmp(exptotal->methodScatterSearch->eSSversion, "aCeSS_dist") == 0) {
        return 2;
    }
    else if (strcmp(exptotal->methodScatterSearch->eSSversion, "saCeSS") == 0) {
        return 3;
    }
    else if (strcmp(exptotal->methodScatterSearch->eSSversion, "eSSm") == 0) {
        return 4;
    }
}

const char* getversioness(int i){
    const char *retchar;
    
    if (i==0) {
        retchar="sequential eSS";
    } 
    else if (i==1) {
        retchar="CeSS - COOPERATIVE eSS";
    }
    else if (i==2) {
        retchar="aCeSS_dist - ASYNCHONOUS COOPERATIVE eSS";        
    }
    else if (i==3) {
        retchar="saCeSS - SELF-ADAPTED ASYNCHONOUS COOPERATIVE eSS";                
    }
    else if (i==4) {
        retchar="eSSm - eSS with multiple configuration and without communications";                
    }    
    
    return retchar;
}


const char * getlsess(experiment_total *exp){
    return exp->methodScatterSearch->loptions->solver;
}

const char* gettopologyess(int i){
    const char *retchar;
    
        
    
        if (i==0) {
            retchar="no communications";
        } 
        else if (i==1) {
            retchar="master-slave";
        }
        else if (i==2) {
            retchar="ring topology, without master";        
        }
        else if (i==3) {
            retchar="mixture between master-slave and star topology";                
        }
        else if (i==4) {
            retchar="no communications";                
        }   
    
    return retchar;
    
}

void destroyexp(experiment_total *exp) {
    const char *matlabproblem;

    matlabproblem="matlabproblem";
#ifdef MATLAB
     if (strcmp((*exp).test.bench.type, matlabproblem) == 0) {
#ifdef GNU
        __matlabproblem_MOD_closematlab(&exp->execution.ep_matlab);
#elif defined(INTEL)
        matlabproblem_mp_closematlab_(&exp->execution.ep_matlab);
#endif
    }     
#endif 

    if (exp->test.bench.max_dom != NULL) {
        free(exp->test.bench.max_dom);
        exp->test.bench.max_dom = NULL;
    }
    if (exp->test.bench.min_dom != NULL) {
        free(exp->test.bench.min_dom);
        exp->test.bench.min_dom = NULL;
    }
    if ( exp->test.ineq > 0 ) { 
    	if (exp->test.bench.CL != NULL) {
    	    free(exp->test.bench.CL);
    	    exp->test.bench.CL = NULL;
    	}
    	if (exp->test.bench.CU != NULL) {
    	    free(exp->test.bench.CU);
    	    exp->test.bench.CU = NULL;
    	}
    }
    if (exp->methodDE != NULL) {
        free(exp->methodDE);
        exp->methodDE = NULL;
    }

    if (exp->methodScatterSearch != NULL) {
        free(exp->methodScatterSearch);
        exp->methodScatterSearch = NULL;
    }
    if (exp->par_st != NULL) {
        free(exp->par_st);
        exp->par_st = NULL;
    }
    if (exp->param != NULL) {
        free(exp->param);
        exp->param = NULL;
    }
    if (exp->amigo != NULL) {
        free_AMIGO_problem(exp->amigo);
        exp->amigo = NULL;
    }
    if ((*exp).ls != NULL) {
        free((*exp).ls);
        (*exp).ls = NULL;
    }

    if ((*exp).output != NULL) {
        free((*exp).output);
        (*exp).output = NULL;
    }
    if ((*exp).test.output_graph!=NULL) {
        free((*exp).test.output_graph);
        (*exp).test.output_graph=NULL;
    }

    if ((*exp).test.output_path!=NULL) {
        free((*exp).test.output_path);
        (*exp).test.output_path=NULL;
    }
    if ((*exp).test.log_output != NULL) {
        free((*exp).test.log_output);
        (*exp).test.log_output=NULL;
    }
    
    if ((*exp).test.log_percentage != NULL) {
        free((*exp).test.log_percentage);
        (*exp).test.log_percentage=NULL;
    }
    if ((*exp).test.output_gant_log != NULL) {
        free((*exp).test.output_gant_log);
        (*exp).test.output_gant_log=NULL;
    }
    if ((*exp).test.result_output != NULL) {
        free((*exp).test.result_output);
        (*exp).test.result_output=NULL;
    }
    
}


int destroybenchmark(experiment_total *exp){
    int benchmark;
    
    benchmark = number_benchmark(exp);
    
    if (benchmark == 0)
    destroyBBOB(exp);
    if (benchmark == 1)
    destroySystemBiology(exp);

    return 1;
}

int is_asynchronous(experiment_total exp) {
    if (strcmp(exp.par_st->commu, "asynchronous") == 0) return 1;
    else return 0;
}

int is_parallel(experiment_total exp) {
    if (exp.par_st != NULL) return 1;
    else return 0;
}

int is_noise(experiment_total exp) {
    if (strcmp(exp.test.bench.type, "noiselessBBOB") == 0) return 0;
    else return 1;
}


int number_benchmark(experiment_total *exp) {
    char const *noiseBBOB;
    char const *noiselessBBOB;
    char const *system;
    char const *test;
    noiseBBOB = "noiseBBOB";
    noiselessBBOB = "noiselessBBOB";
    system = "systemBiology";
    test = "customized";

    if (strcmp((*exp).test.bench.type, test)) {
        return 0;
    } else if (strcmp((*exp).test.bench.type, system) == 0) {
        return 1;
    } else if (strcmp((*exp).test.bench.type, noiseBBOB) == 0 || strcmp((*exp).test.bench.type, noiselessBBOB) == 0) {
        return 2;
    }else {
        return -1;
    }

}

void check_fun(experiment_total exp, int current_bench, int noise_funcion) {
    char const *noiseBBOB;
    char const *noiselessBBOB;
    char const *system;
    char const *test;
    char const *LSGO;
    
    noiseBBOB = "noiseBBOB";
    noiselessBBOB = "noiselessBBOB";
    system = "systemBiology";
    test = "test";
    LSGO= "LSGO";


    if (strcmp(exp.test.bench.type, noiseBBOB) == 0 || strcmp(exp.test.bench.type, noiselessBBOB) == 0) {
        if (noise_funcion) {
            if (current_bench > 130) perror("ERROR FUNCTION");
            if (current_bench < 101) perror("ERROR FUNCTION");
        } else {
            if (current_bench > 1) perror("ERROR FUNCTION");
            if (current_bench < 24) perror("ERROR FUNCTION");
        }
    } else if (strcmp(exp.test.bench.type, system) == 0) {
        if (current_bench > 5) perror("ERROR FUNCTION");
        if (current_bench < 0) perror("ERROR FUNCTION");
    } else if (strcmp(exp.test.bench.type, test) == 0) {

    } else if (strcmp(exp.test.bench.type, LSGO) == 0) {
        if (current_bench > 20) perror("ERROR FUNCTION");
        if (current_bench < 1) perror("ERROR FUNCTION");
    }    
}

int chargedimension_(void *exp1_, int *D) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
            
    *D = exp1[0].test.bench.dim;
    
    return 1;
}

int chargedimensionopenmp_(void *exp1_, int *D) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    *D = exp1[0].test.bench.dim;
    
    return 1;
}



int chargeuseroptions_(void *exp_, double *maxtime, int *weight, double *tolc, double *prob_bound,
        int *nstuck_solution, int *strategy, int *inter_save, int *iterprint, int *plot, int *log_var, int *init_point, int *logvar) {
    experiment_total *exp;
    int i;

    exp = (experiment_total *) exp_;
    exp[0].test.output_stop = 0;  
    exp[0].test.jfprint = 0;
    if (exp[0].methodScatterSearch != NULL) {
          
        if (exp[0].methodScatterSearch->uoptions != NULL) {
                       
            *maxtime=exp[0].test.maxtime;
            *weight = exp[0].methodScatterSearch->uoptions->weight;
            *tolc = exp[0].methodScatterSearch->uoptions->tolc;
            *prob_bound = exp[0].methodScatterSearch->uoptions->prob_bound;
            *nstuck_solution = exp[0].methodScatterSearch->uoptions->nstuck_solution;
            *strategy = exp[0].methodScatterSearch->uoptions->strategy;
            *inter_save = exp[0].methodScatterSearch->uoptions->inter_save;
            *iterprint = exp[0].test.verbose;
            *plot = exp[0].test.output;
            if (exp[0].test._log == 1) {
		for (i=0;i<exp[0].test.bench.dim;i++) {
			if (exp[0].test.bench.logindex[i] == 1) {
				log_var[i] = 1;
			} else {
				log_var[i] = 0;
			}
		}
                *logvar=1;
            } else {
                for (i=0;i<exp[0].test.bench.dim;i++) {
                                log_var[i] = 0;
                }
		*logvar=0;
	    } 
            *init_point=0;
            return 1;
        } else return 0;
    } else return 0;
    return 1;
}

int chargeglobaloptions_( void *exp_, int *dim_ref, int *ndiverse, int *initiate, int *combination, int *regenerate, 
        char* delete, int *intens, double * tolf, int *diverse_criteria, double * tolx, int *n_stuck) {
    experiment_method_ScatterSearch *method;
    experiment_total *exp;
    const char *NULLchar;
    
    NULLchar = "\0";
    
    exp = (experiment_total *) exp_;
    method = exp[0].methodScatterSearch;
    if (method != NULL) {   
        if (method->goptions != NULL) {
            *dim_ref = method->goptions->dim_ref;
            *ndiverse = method->goptions->ndiverse;
            *initiate = method->goptions->initiate;
            *combination = method->goptions->combination;
            *regenerate = method->goptions->regenerate;
           
            if (method->goptions->delete1 == NULL ) {
                method->goptions->delete1 = strtok(method->goptions->delete1, NULLchar);
                sprintf(delete,"%s", (const char *) method->goptions->delete1);
            }
            *intens = method->goptions->intens; 
            *tolf = method->goptions->tolf;
            *diverse_criteria = method->goptions->diverse_criteria;
            *tolx = method->goptions->tolx;
            *n_stuck = method->goptions->n_stuck;
            return 1;
        } else return 0;
    } else return 0;
       
    return 1;     
}

int chargelsevalmax_( void *exp_, long *lsevals){
    experiment_total *exp;
    experiment_method_ScatterSearch *method;
    exp = (experiment_total *) exp_;
    
    if (exp[0].test.local_search == 1) {
        method = exp[0].methodScatterSearch;
        if (method != NULL) {
            *lsevals = method->loptions->evalmax;
        }
    }
    
}


int chargelocaloptions_( void *exp_, int *tol, int *iterprint, int *n1, int *n2, double * balance, 
        char * finish, int *bestx, int *merit_filter, int *distance_filter, double *thfactor, 
        double *maxdistfactor, int *wait_maxdist_limit, int *wait_th_limit, char *solver, 
        double *threshold_local) {
    experiment_method_ScatterSearch *method;
    experiment_total *exp;
    const char *NULLchar;
    NULLchar = "\0";
    
    exp = (experiment_total *) exp_;

    *balance=0.0;
    if (exp[0].test.local_search == 1) {
        method = exp[0].methodScatterSearch;
        if (method != NULL) {

            if (method->loptions != NULL) {

                * tol = method->loptions->tol;
                * iterprint = method->loptions->iterprint;
                * n1 = method->loptions->n1;
                * n2 = method->loptions->n2;
                * balance = method->loptions->balance;
                if (*n1 > -1) {
                    method->loptions->solver = strtok(method->loptions->solver, NULLchar);
                    if (method->loptions->solver != NULL) {
                        sprintf(solver, "%s", method->loptions->solver);
                    }
                    method->loptions->finish = strtok(method->loptions->finish, NULLchar);
                    if (method->loptions->finish != NULL) {
                        sprintf(finish, "%s", method->loptions->finish);
                    } else if (method->loptions->finish == NULL) {
                        finish = "";
                    }
                    *bestx = method->loptions->bestx;
                    *merit_filter = method->loptions->merit_filter;
                    *distance_filter = method->loptions->distance_filter;
                    *thfactor = method->loptions->thfactor;
                    *maxdistfactor = method->loptions->maxdistfactor;
                    *wait_maxdist_limit = method->loptions->wait_maxdist_limit;
                    *wait_th_limit = method->loptions->wait_th_limit;
                    return 1;
                } else {
                    return 0;
                }
            } else return 0;
        } else return 0;
    } else return 0;
}


int chargeproblemargs_(void *exp_, int *ineq, int *int_var, int *bin_var, int *neq) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    *ineq = exp1[0].test.ineq;
    *neq = exp1[0].test.neq;
    *int_var = exp1[0].test.int_var;
    *bin_var = exp1[0].test.bin_var;  
    
    return 1;
}


int getopenmpoption_(void *exp_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;

    return exp1->test.bench.openmp;
}

int chargeboundsnconst_(void *exp_, double *XU, double *XL, int *nvar , double *CU, double *CL, int *ineq) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp_;
    
    for (i=0;i<*nvar;i++) {
                XU[i] = exp1[0].test.bench.max_dom[i];
                XL[i] = exp1[0].test.bench.min_dom[i];
    }

    if (exp1[0].test.bench.CU != NULL && (*ineq > 0)) {
        for (i = 0; i<*ineq; i++) {
            CU[i] = exp1[0].test.bench.CU[i];
        }
    }
    
    if (exp1[0].test.bench.CL != NULL && (*ineq> 0)) {
        for (i = 0; i<*ineq; i++) {
            CL[i] = exp1[0].test.bench.CL[i];
        }
    }

    return 1;
}


int chargebounds_(void *exp_, double *XU, double *XL, int *nvar) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp_;
    
    
    for (i=0;i<*nvar;i++) {
                XU[i] = exp1[0].test.bench.max_dom[i];
                XL[i] = exp1[0].test.bench.min_dom[i];
    }
    
    return 1;
}



void getbench_(void *exp_, int *id) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    *id = exp1->test.bench.current_bench;
}

void setdblmax_(double *valor){
    *valor = DBL_MAX;
}

void setinfinity_(double *valor){
    *valor = INFINITY;
}

void setnan_(double *valor){
    *valor = NAN;
}


double gettotaltime_( result_solver *timestruct) {
    return timestruct->totaltime;
}


int ishete_(void *exp_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    char *homo;
    char *hete;
    
    homo = "homo";
    hete = "hete";
    
            
    if (strcmp(exp1->par_st->islandstrategy,homo)==0) {
        return 0;
    }
    else if (strcmp(exp1->par_st->islandstrategy,hete)==0 ){
        return 1;
    } else return 0;
    
}



int iscoop_(void *exp_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    char *coop;
    char *island;
    
    
    island="island";
    coop="cooperative";
            
    if (strcmp(exp1->par_st->type,island)==0) {
        return 0;
    }
    else if (strcmp(exp1->par_st->type,coop)==0 ){
        return 1;
    } else return 1;
    
}

double getparalleltime_( result_solver *timestruct) {
    return timestruct->paralleltime;
}

double getlocalsolvertime_( result_solver *timestruct) {
    return timestruct->localsolvertime;
    
}



void settotaltime_(result_solver *timestruct, double *time){
    timestruct->totaltime = timestruct->totaltime + *time;
}

void setparalleltime_(result_solver *timestruct, double *time){
    timestruct->paralleltime = timestruct->paralleltime + *time;
}

void setlocalsolvertime_(result_solver *timestruct, double *time){
    timestruct->localsolvertime = timestruct->localsolvertime + *time;
}


void setiteration_(void *exp_, int *ite) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    exp1->execution.iteration = *ite;
    
}


void setnumit_(void *exp_, int *ite) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    exp1->execution.num_it = *ite;
    
}



int getidp_(void *exp_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    
    return exp1->execution.idp;
}

void saveinittime_(void *exp_, double *starttime) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    exp1->execution.initTIME = *starttime;
}

double returninittime_(void *exp_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    return exp1->execution.initTIME;
}

void setamigolocalevals_(void *exp_, int *levals)  {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;
    
    exp1->amigo->local_max_evals = *levals;
}


void setparallelsacessfieldsmaster_(void *exp_, double *rth )  {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;

    if ( exp1->par_st->reception_threshold > 0 )
    	*rth = exp1->par_st->reception_threshold;
    else
	*rth = 0.0;
}

void setparallelsacessfieldsslaves_(void *exp_,  double *eth, int *numsend, int *minsend)  {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp_;

    *eth = (double) exp1->par_st->evals_threshold;
    *numsend = exp1->par_st->mult_num_sendSol;
    *minsend = exp1->par_st->minimum_num_sendSol;

}

void setinf_(double *value){
	*value = INFINITY;
}

void setminf_(double *value) {
	*value=-INFINITY;
}

