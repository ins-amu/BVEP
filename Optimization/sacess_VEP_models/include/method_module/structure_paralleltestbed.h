/* 
 * File:   structure_paralleltestbed.h
 * Author: david
 *
 * Created on 30 de mayo de 2013, 13:13
 */

#ifndef STRUCTURE_PARALLELTESTBED_H
#define	STRUCTURE_PARALLELTESTBED_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include <bbobStructures.h>  
#include <AMIGO_problem.h> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
    
#ifdef	MPI2
    #include <dynamic_list.h>
    #include <mpi.h> 
#endif    

#include <Python.h>
    typedef struct paramStruct ParamS;
    
    typedef struct {
        int *rigth;
        int num_r;
        int *left;
        int num_left;
        #ifdef MPI2
                MPI_Comm comunicator;
        #endif
    } topology_data; 

    typedef struct {
       int n_stuck;
       long evals_threshold;
       int mult_num_sendSol;
       int minimum_num_sendSol;
       int evalmax;
    } argc_struct;

    typedef struct {
        int current_bench;
        double *min_dom;
        double *max_dom;
        double *log_min_dom;
        double *log_max_dom;  
        int * logindex;
        int idtype;
        int use_amigo;      
        int estime_init_cond;
        double *CL;
        double *CU;           
        int dim;
        const char* type;   
        double *BestSol;
        int translation;
	int openmp;
	double **X0;
	double *F0;
	int number_init_sol;
    } experiment_benchmark;
    
    typedef struct {
        double _F;
        double _CR;
        int _NP;
        int cache;
	int ls_counter;
	double ls_threshold;
        const char *mutation_strategy;
        const char *name; 
        const char *solver;   
        int num_proc;
    } experiment_method_DE;
    
    
    typedef struct {
        int weight;
        double tolc;
        double prob_bound;
        int nstuck_solution;
        int strategy;
        int inter_save;
    } user_options;
    
    typedef struct {
        int dim_ref;
        int dim_ref_global;
        int ndiverse;
        int initiate;
        int combination;
        int regenerate;
        char* delete1;
        int intens;
        double tolf;
        int diverse_criteria;
        double tolx;
        int n_stuck;
    } global_options; 
    
    typedef struct {
        int *texp;
        double *yexp;
        double k1;
        double k2;
        double k3;
        double k4;
    } extraparametters;
    
    typedef struct {
        int tol;
        int iterprint;
        int n1;
        int n2;
        double balance;
        char * solver;
        char * finish;
        int bestx;
        int merit_filter;
        int distance_filter;
        double thfactor;
        double maxdistfactor;
        int wait_maxdist_limit;
        int wait_th_limit;
        long evalmax;        
        extraparametters extrap;
    } local_options;

    typedef struct {
        int _NP;
        int cache;
	int ls_counter;
	double ls_threshold;
        const char *name;  
        user_options *uoptions;
        global_options *goptions;
        local_options *loptions;
        int num_proc;
        int asynchronous;
        char *eSSversion;

    } experiment_method_ScatterSearch;  
    
    typedef struct {
        int instances;
        int init_point;
        experiment_benchmark bench;
        double   max_eval;
        double VTR;
        double VTR_default;
        double jf;
        int jfprint;
        double   repetitions;
        const char* type;   
        const char* namexml;   
        int _log;
	char* output_path; 
        char* output_graph;
        char* output_gant_log;
        int output;
        char *log_output;
        char *log_percentage;
        char *result_output;
        int verbose;
        int init_repetition;
        int local_search;
        int local_gradient;
        int output_stop;
        double maxtime;
        int int_var;
        int bin_var;
        int neq;
        int ineq;
    } experiment_testbed; 
    
    typedef struct {
        int NPROC;
        int NPROC_OPENMP;
        int migration_size;
        int migration_freq_ite;    
        double max_time_ite;
        const char* type;
        const char* commu;
        const char* islandstrategy;
        const char* topology;
        const char* SelectionPolicy;
        const char* ReplacePolicy;
        int evals_threshold;
        int mult_num_sendSol;
        int minimum_num_sendSol;
	double reception_threshold;
        double max_time;
    }parallelization_strategy;
    
    
    
    typedef struct {
        char *nameMatlab;
        topology_data topology;
        int size_send_buffer;
        int send_id;
        int *st;
        int idp;
        int num_it;
        double max_time_ite;
        double max_time_ite_last;
        int NPROC;
        int tam;
        int NP;
        int NM;
        int st_sent;
        int stuckcond_lowVar;
        int stuckcount;
        int migra_asin_wait;
        int migration;
        int contadorMigra;
        int enterMigrat;
        int update;
        int max_eval_ite_last;
        double minvarcondition;
        double ftarget;
        long maxfunevals;
        double initTIME;
        int rep;
        int initpath;
        int iteration;
        double *transconst;
#ifdef MPI2
        list l;
        MPI_Request *request_recept;
        MPI_Request *receptionrequestmaster; // recv in master
        MPI_Request *sendrequestmaster; // send to slave
        MPI_Request *receptionrequestslave;
        MPI_Request *sendrequestslave; 
        MPI_Request *adaptation_master_send;
        MPI_Request *adaptation_slave_send;
        MPI_Request *adaptation_master_recv;
        MPI_Request *adaptation_slave_recv;        
        int *adaptation_master_buffer_recv;
        double **adaptation_master_buffer_send;
        int *adaptation_slave_buffer_send;
        double *adaptation_slave_buffer_recv;
        double **receptionbuffermaster;
        double *receptionbufferslave; 
        double **sendbuffermaster;
        double *sendbufferslave;
        int *resetbuffer;
        long *masterlistrecp;
        double mastertime;
        double currentmastertime;
#endif  
        char *file; 
        long failevals;
	int dist_stopping_criteria;
	long inner_ls_evals;
        int ep_matlab;
        PyObject *pFunc;
        PyObject *pModule;
	PyGILState_STATE *gstate;
        const char *file_python_name;
        int python_active;
	int family_slave;
    } execution_vars;
    
    
    typedef struct {
        int counter;
        int state_ls;
        double total_local;
        double sucess_local;
        int num_ls;
        int sucess_interval;
        double enter_value;
    } local_solver;
    
    typedef struct {
        double st1;
        double st2;
        double st3;
        double point1;
        int point_counter;
        double oldbest;
    }output_struct;
       
    typedef struct {
        experiment_testbed test;
        experiment_method_DE  *methodDE;
        experiment_method_ScatterSearch *methodScatterSearch;
        parallelization_strategy *par_st;
        local_solver *ls;
        output_struct *output;
	void *fp;
        ParamS *param;
        AMIGO_problem *amigo;        
        execution_vars execution;
        const gsl_rng *random;
        double *seed;
        int contadorseed;
    }experiment_total;    
    
    typedef struct {
        double eval_avg; 
        double time_avg;
        double best_value_all;
        double restart;
        double value_avg;
        double porcenthits;
        int dim;
        int id;
    }result_values;
    
    
    typedef struct {
        long eval_value; 
        double time_value;
        double time_vtr_value;
        double iterations_value;
        double best_value; 
        double *bestx_value; 
        double totaltime;
        double paralleltime;
        double localsolvertime;        
    } result_solver;
    
    typedef struct {
        double value;
        double *g;
        double *R;
	int size_r;
	double *J;
	int size_j;
    } output_function;

void logandtranslation_(void *exp1);

int create_expetiment_struct(const char *file, experiment_total *exptotal, int NPROC, int id,  const char *path, int init,  argc_struct arg);

void init_result_data(result_solver *, int);
    
const char* getname(experiment_total *);
    
void updatenp_(void *, int *);    
    
void destroyexp(experiment_total *) ;

int is_asynchronous(experiment_total);

int is_parallel(experiment_total);    

int is_noise(experiment_total);  

void check_fun(experiment_total, int , int );

int chargedimension_(void *, int *);

int chargedimensionopenmp_(void *, int *);

void getbench_(void *, int *);

double gettotaltime_( result_solver *);

double getparalleltime_( result_solver *);

double getlocalsolvertime_(  result_solver *);

int getidp_(void *);

int ishete_(void *);

int iscoop_(void *);

void setiteration_(void *, int *);

void setnumit_(void *, int *);

int chargeproblemargs_(void *, int *, int *, int *, int *);

int chargeuseroptions_(void *, double *, int *, double *, double *,
        int *, int *, int *, int *, int *, int *, int *, int *);

int chargeglobaloptions_( void *, int *, int *, int *, int *, int *, 
        char* , int *, double * , int *, double * , int *);

int chargelocaloptions_( void *, int *, int *, int *, int *, double * , 
        char * , int *, int *, int *, double *, 
        double *, int *, int *, char *, double *) ;

int chargeboundsnconst_(void *, double *, double *, int * , double *, double *, int *);


int chargebounds_(void *, double *, double *, int *);    

int number_benchmark(experiment_total *);

void saveinittime_(void *, double *);

double returninittime_(void *);

int chargelsevalmax_( void *, long *);

int destroybenchmark(experiment_total *);

void setparallelsacessfieldsmaster_(void *, double *);

void setparallelsacessfieldsslaves_(void *,  double *, int *, int *);  

int getopenmpoption_(void *);

int returninitsol_(void *);

void loadinitpoints_(void *, double *, int *, double * );

void setdistcriteria_(void *, int *);

void getdistcriteria_(void *,  int *);


#ifdef	__cplusplus
}
#endif

#endif	/* STRUCTURE_PARALLELTESTBED_H */





