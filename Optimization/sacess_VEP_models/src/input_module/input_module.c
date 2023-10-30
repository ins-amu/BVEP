/**
 * @file input_module.c
 * @author David R. Penas
 * @brief File containing functions to parse the input XML, completing the main
 * struct of experiment_total data.
 */



#include <structure_paralleltestbed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <def_errors.h>
#include <hdf5.h>
#if  defined(OPENMP) 
    #include <omp.h>
#endif
#include <sys/stat.h> 


char* removeSpace(char *str) {
  char *p1 = str, *p2 = str;
  do 
    while (*p2 == ' ')
      p2++;
  while (*p1++ = *p2++);
  
  return str;
}

void removePoint(char *str) {
  char *p1 = str, *p2 = str;
  do 
    while (*p2 == '.')
      p2++;
  while (*p1++ = *p2++);
}


xmlNodePtr extract_init_node(xmlNodePtr cur, const char *name) {
    int exit;
    exit = 0;
    xmlNodePtr end_cur;
    
    end_cur = cur;
     while (end_cur != NULL && exit == 0)  {
        if (!xmlStrcmp(end_cur->name, (const xmlChar *) name )) {
            exit = 1;   
        } else {
            end_cur = end_cur->next;
        }
    }
    
    return end_cur;
} 

char*  extract_element_uniq(xmlDocPtr doc, xmlNodePtr cur, const char *name ) {
    xmlChar *key;
    xmlChar *output;
    
    
    cur = cur->xmlChildrenNode;

    output = NULL;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)name)))  {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            output = (xmlChar*) malloc(((xmlStrlen(key))+1) * sizeof (xmlChar));
            memmove(output,key,((xmlStrlen(key))+1)*sizeof(xmlChar));
            //strcpy((const char *) output, (const char *) key);
            xmlFree(key);
        }
        cur = cur->next;
    }
    return (char *) output;
}


int  count_element_multi(xmlDocPtr doc, xmlNodePtr cur, const char *name) {
    xmlChar *key;
    xmlNodePtr cur_pre;
    cur_pre = cur->xmlChildrenNode;
    cur = cur->xmlChildrenNode;
    int hit;
 
    hit=0;
    while (cur_pre != NULL) {
        if ((!xmlStrcmp(cur_pre->name, (const xmlChar *)name)))  {
            key = xmlNodeListGetString(doc,cur_pre->xmlChildrenNode, 1);
            hit++;
        }
        cur_pre = cur_pre->next;
    } 
    
    return hit;
}

void  extract_element_multi(xmlDocPtr doc, xmlNodePtr cur, const char *name,  char **vectorG) {
    xmlChar *key, *fxxml;
    xmlChar *output;
    int counter;

    cur = cur->xmlChildrenNode;

    counter=0;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)name)))  {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            memmove(vectorG[counter],key,((xmlStrlen(key))+1)*sizeof(xmlChar));
            xmlFree(key);
	    counter++;
        }

        cur = cur->next;
    }
}

void  extract_atribute_multi(xmlDocPtr doc, xmlNodePtr cur, const char *name, double *FX, const char *atrib) {
    xmlChar *key;
    xmlChar *output;
    int counter;

    cur = cur->xmlChildrenNode;

    counter=0;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)name)))  {
        key = xmlGetProp(cur, atrib);
        if ( key != NULL ) {
            FX[counter]= atof((char*) key);
        } else  {
            FX[counter]= DBL_MAX;
	}

        xmlFree(key);
        counter++;
	}
        cur = cur->next;

    }
}





int extract_element_test(xmlDocPtr doc, xmlNodePtr *root, experiment_testbed *test) {
    xmlNodePtr cur, init_cur, stc_cur; 
    const char *name_exp;
    const char *name_exptype;
    const char *name_stc;    
    const char *name_bench_type;
    const char *name_bench_dim;
    const char *name_bench_min;
    const char *name_bench_max;
    const char *name_eval;
    const char *name_tol;
    const char *name_rep;
    const char *name_minDom, *name_maxDom;
    const char *name_scale_log;
    char *vectorDim, *outputDim, *min, *max; 
    int *vectorsize;
    int i, size, hit;
    char const *noiseBBOB;
    char const *noiselessBBOB;
    char const *system;
    char const *name_output;    
    char const *name_verbose;
    char const *name_output_temp;
    char const *name_output_graph;  
    char const *name_init_point;
    char const *local_search;
    char const *local_gradient;
    char const * testname;
    char const * maxtimec;
    char const * LFSO;
    char const * idbench;
    char const * name_vtr;
    char const * python_problem;

    python_problem="python";
    testname = "test";
    noiseBBOB ="noiseBBOB";
    noiselessBBOB="noiselessBBOB";
    system="systemBiology";    
    LFSO = "LGSO";
    name_exp= "run";
    name_exptype= "exptype";
    name_stc= "stopping_criteria";
    name_bench_type= "typebench";    
    name_bench_dim= "dim";    
    name_bench_min= "min_bench";
    name_bench_max= "max_bench";
    name_eval= "maxevaluation";
    name_tol= "tolerance";
    name_rep= "repetitions";
    name_minDom = "minDomain";
    name_maxDom = "maxDomain";
    name_scale_log = "log_scale";
    name_verbose =    "verbose";
    name_output = "output";
    name_output_temp = "output_temp";
    name_output_graph = "output_graph";
    name_init_point = "init_point";
    name_vtr = "vtr";
    local_search= "local_search";
    local_gradient= "gradient";
    maxtimec= "maxtime";
    idbench= "id";
    
    init_cur = (*root)->xmlChildrenNode;
    cur = extract_init_node(init_cur, name_exp);
    
    if (cur == NULL) {
        perror(error28);
        error(28);
    }
    if (extract_element_uniq(doc,cur,name_bench_type) == NULL) {
        perror(error17);
        exit(17);
    }
    test->bench.type = removeSpace(extract_element_uniq(doc,cur,name_bench_type));
    
    if (extract_element_uniq(doc, cur, idbench) == NULL) {
        perror(error18);
        exit(18);
    }    
    test->bench.current_bench = atoi(extract_element_uniq(doc, cur, idbench)); 
    test->bench.use_amigo=0; 
    if (extract_element_uniq(doc, cur, name_scale_log) == NULL) {
        perror(error20);
        exit(20);
    }            
    test->_log = atoi(extract_element_uniq(doc, cur, name_scale_log));  
    
    if (extract_element_uniq(doc, cur, name_output) == NULL) {
        perror(error21);
        exit(21);
    }          
    test->output = atoi(extract_element_uniq(doc, cur, name_output));  
    
    if (extract_element_uniq(doc, cur, name_verbose) == NULL) {
        perror(error22);
        exit(22);
    }      
    test->verbose = atoi(extract_element_uniq(doc, cur, name_verbose)); 
    
    if (extract_element_uniq(doc,cur,maxtimec) == NULL) test->maxtime=-1; 
    else    test->maxtime = atof(extract_element_uniq(doc,cur,maxtimec));
    
    
    if (extract_element_uniq(doc, cur, local_search) == NULL) {
        perror(error23);
        exit(23);
    }          
    test->local_search = atoi(extract_element_uniq(doc, cur, local_search)); 
    
    stc_cur = cur->xmlChildrenNode;
    cur = extract_init_node(stc_cur, name_stc);
    if (cur == NULL) {
	perror(error42);
	exit(42);
    }

    if (extract_element_uniq(doc, cur, name_eval) == NULL) {
        perror(error19);
        exit(19);
    }
    test->max_eval = atof(extract_element_uniq(doc,cur,name_eval));

    if (extract_element_uniq(doc,cur,maxtimec) == NULL) {
	test->maxtime=-1;
    } else    {
	test->maxtime = atof(extract_element_uniq(doc,cur,maxtimec));
    } 

    if (extract_element_uniq(doc, cur, name_vtr) == NULL) {
        perror(error30);
        exit(30);
    }

    test->VTR = atof(extract_element_uniq(doc,cur,name_vtr));

    return 1;
}



int extract_element_method_ScatterSearch(xmlDocPtr doc, xmlNodePtr *root, experiment_method_ScatterSearch *method) {
    xmlNodePtr cur, init_cur,init_cur2, cur2; 

    const char *name_exp;
    const char *name_name;
    const char *int_var;
    const char *bin_var;
    const char *ineq;
    
    const char *user_optionsc;
    const char *weightc;
    const char *tolcc;
    const char *prob_boundc;
    const char *nstuck_solutionc;
    const char *strategyc;
    const char *inter_save;
    
    const char *global_optionsc;
    const char *dim_refc;
    const char *ndiversec;
    const char *initiatec;
    const char *combinationc;
    const char *regeneratec;
    const char *deletec;
    const char *intensc;
    const char *tolfc;
    const char *diverse_criteriac;
    const char *tolxc;
    const char *n_stuckc;

    const char * local_opstionsc;
    const char * tolc;
    const char * iterprintc;
    const char * n1c;
    const char * n2c;
    const char * balancec;
    const char * solverc;
    const char * finishc;
    const char * bestxc;
    const char * merit_filterc;
    const char * distance_filterc;
    const char * thfactorc;
    const char * maxdistfactorc;
    const char * wait_maxdist_limitc;
    const char * wait_th_limitc;
    const char * evalmaxc;
    const char * neq;
 
    evalmaxc = "evalmax";
    name_exp= "method";
    name_name = "name";
    
    user_optionsc = "user_options";
    weightc = "weight";
    tolcc = "tolc";
    prob_boundc = "prob_bound";
    nstuck_solutionc = "nstuck_solution";
    strategyc = "strategy";
    inter_save = "inter_save";
    
    global_optionsc="global_options";
    dim_refc="dim_ref";
    ndiversec="ndiverse"; 
    initiatec="initiate"; 
    combinationc="combination"; 
    regeneratec="regenerate"; 
    deletec="delete"; 
    intensc="intens"; 
    tolfc="tolf"; 
    diverse_criteriac="diverse_criteria"; 
    tolxc="tolx"; 
    n_stuckc="n_stuck"; 
    
    local_opstionsc = "local_options";
    tolc = "tol";
    iterprintc="iterprint";
    n1c="n1";
    n2c="n2";
    balancec="balance";
    finishc="finish";
    bestxc="bestx";
    merit_filterc="merit_filter";
    distance_filterc="distance_filter";
    thfactorc="thfactor";
    maxdistfactorc="maxdistfactor";
    wait_maxdist_limitc="wait_maxdist_limit";
    wait_th_limitc="wait_th_limit";
    solverc="solver";
    
    init_cur = (*root)->xmlChildrenNode;
    cur = extract_init_node(init_cur, name_exp);
    
    if (cur == NULL) {
        perror(error29);
        exit(29);
    }
      
    init_cur2 = cur->xmlChildrenNode;
    cur2 = extract_init_node(init_cur2, user_optionsc);
    if (cur2 != NULL ) {
        method->uoptions = (user_options *) malloc(sizeof(user_options));
        if ((extract_element_uniq(doc,cur2,weightc) != NULL ) &&
	   (strcmp(removeSpace(extract_element_uniq(doc,cur2,weightc)),"default")!=0))
	{                
                method->uoptions->weight = atoi(extract_element_uniq(doc,cur2,weightc));
        } 
        else  method->uoptions->weight = -1;

        if (( extract_element_uniq(doc,cur2,tolcc) != NULL ) &&
           (strcmp(removeSpace(extract_element_uniq(doc,cur2,tolcc)),"default")!=0))
        {
            method->uoptions->tolc = atof(extract_element_uniq(doc,cur2,tolcc));
        } else method->uoptions->tolc = -1.0;

        if (( extract_element_uniq(doc,cur2,prob_boundc) != NULL ) &&
           (strcmp(removeSpace(extract_element_uniq(doc,cur2,prob_boundc)),"default")!=0))
        {
            method->uoptions->prob_bound = atof(extract_element_uniq(doc,cur2,prob_boundc));
        } else  method->uoptions->prob_bound = -1.0;

        if ((  extract_element_uniq(doc,cur2,nstuck_solutionc) != NULL ) &&
	    (strcmp(removeSpace(extract_element_uniq(doc,cur2,nstuck_solutionc)),"default")!=0))
        {
            method->uoptions->nstuck_solution = atoi(extract_element_uniq(doc,cur2,nstuck_solutionc));
        } else method->uoptions->nstuck_solution = -1;

        if (( extract_element_uniq(doc,cur2,strategyc) != NULL  ) &&
	    (strcmp(removeSpace(extract_element_uniq(doc,cur2,strategyc)),"default")!=0))
        {
            method->uoptions->strategy = atoi(extract_element_uniq(doc,cur2,strategyc));
        } else method->uoptions->strategy = -1;

        if (( extract_element_uniq(doc,cur2,inter_save) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,inter_save)),"default")!=0))
        {
            method->uoptions->inter_save = atoi(extract_element_uniq(doc,cur2,inter_save));
        } else method->uoptions->inter_save = -1;
    } else {
        method->uoptions = (user_options *) malloc(sizeof(user_options));
        method->uoptions->weight = -1;
        method->uoptions->tolc = -1.0;
        method->uoptions->prob_bound = -1.0;
        method->uoptions->nstuck_solution = -1;
        method->uoptions->strategy = -1;
        method->uoptions->inter_save = -1;        
    }
    
    cur2 = extract_init_node(init_cur2, global_optionsc);
    if (cur2 != NULL ) {
        method->goptions = (global_options *) malloc(sizeof(global_options));
        if ((extract_element_uniq(doc,cur2,dim_refc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,dim_refc)),"default")!=0))
 	{
            method->goptions->dim_ref = atoi(extract_element_uniq(doc,cur2,dim_refc));
        }  else method->goptions->dim_ref = -1;

        if ((extract_element_uniq(doc,cur2,ndiversec) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,ndiversec)),"default")!=0))
	{
            method->goptions->ndiverse = atoi(extract_element_uniq(doc,cur2,ndiversec));
        } else method->goptions->ndiverse = -1;

        if (( extract_element_uniq(doc,cur2,initiatec) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,initiatec)),"default")!=0))
	{
            method->goptions->initiate = atoi(extract_element_uniq(doc,cur2,initiatec));
        }    else method->goptions->initiate = -1;

        if ( (extract_element_uniq(doc,cur2,combinationc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,combinationc)),"default")!=0))
	{
            method->goptions->combination = atoi(extract_element_uniq(doc,cur2,combinationc));
        } else method->goptions->combination = -1;

        if ( ( extract_element_uniq(doc,cur2,regeneratec) != NULL ) &&  
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,regeneratec)),"default")!=0))
	{
            method->goptions->regenerate  = atoi(extract_element_uniq(doc,cur2,regeneratec));  
        }

        if (( extract_element_uniq(doc,cur2,deletec) != NULL  )  &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,deletec)),"default")!=0))
	{
            method->goptions->delete1 = removeSpace(extract_element_uniq(doc,cur2,deletec));
        } else method->goptions->delete1 = "";

        if ( (extract_element_uniq(doc,cur2,intensc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,intensc)),"default")!=0))
	{
            method->goptions->intens  = atoi(extract_element_uniq(doc,cur2,intensc)); 
        } else method->goptions->intens = -1;

        if (( extract_element_uniq(doc,cur2,tolfc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,tolfc)),"default")!=0))
	{
            method->goptions->tolf = atof(extract_element_uniq(doc,cur2,tolfc));
        } else  method->goptions->tolf  = -1.0;

        if (( extract_element_uniq(doc,cur2,diverse_criteriac) != NULL   ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,diverse_criteriac)),"default")!=0))
	{
            method->goptions->diverse_criteria = atoi(extract_element_uniq(doc,cur2,diverse_criteriac));
        } else method->goptions->diverse_criteria = -1;

        if (( extract_element_uniq(doc,cur2,tolxc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,tolxc)),"default")!=0))
	{
            method->goptions->tolx = atof(extract_element_uniq(doc,cur2,tolxc));
        } else method->goptions->tolx = -1.0;

        if (( extract_element_uniq(doc,cur2,n_stuckc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,n_stuckc)),"default")!=0))
	{
            method->goptions->n_stuck = atoi(extract_element_uniq(doc,cur2,n_stuckc));
        } else method->goptions->n_stuck = -1;

    } else {
        method->goptions = (global_options *) malloc(sizeof(global_options));
        method->goptions->dim_ref = -1;
        method->goptions->ndiverse = -1;
        method->goptions->initiate = -1;
        method->goptions->combination = -1;
        method->goptions->delete1 = "";
        method->goptions->intens = -1;
        method->goptions->tolf  = -1.0;
        method->goptions->diverse_criteria = -1;
        method->goptions->tolx = -1.0;
        method->goptions->n_stuck = -1;        
    }

    cur2 = extract_init_node(init_cur2, local_opstionsc);

    if (cur2 != NULL ) {
        method->loptions = (local_options *) malloc(sizeof(local_options));
        if ((  extract_element_uniq(doc,cur2,tolc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,tolc)),"default")!=0))
        {
            method->loptions->tol = atoi(extract_element_uniq(doc,cur2,tolc)); 
        } else method->loptions->tol = -1;

        if (( extract_element_uniq(doc,cur2,iterprintc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,iterprintc)),"default")!=0)){
            method->loptions->iterprint = atoi(extract_element_uniq(doc,cur2,iterprintc)); 
        } else method->loptions->iterprint = -1;

        if ((  extract_element_uniq(doc,cur2,n1c) != NULL )  &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,n1c)),"default")!=0)){
            method->loptions->n1 = atoi(extract_element_uniq(doc,cur2,n1c)); 
        } else method->loptions->n1 = -1;

        if ((  extract_element_uniq(doc,cur2,n2c) != NULL )  &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,n2c)),"default")!=0)){
                method->loptions->n2 = atoi(extract_element_uniq(doc,cur2,n2c)); 
        } else  method->loptions->n2 = -1;

        if (( extract_element_uniq(doc,cur2,balancec) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,balancec)),"default")!=0)){
            method->loptions->balance = atof(extract_element_uniq(doc,cur2,balancec));
        } else method->loptions->balance = -1;

        if ((  extract_element_uniq(doc,cur2,evalmaxc) != NULL )  &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,evalmaxc)),"default")!=0)){
                method->loptions->evalmax = atoi(extract_element_uniq(doc,cur2,evalmaxc)); 
        } else  method->loptions->evalmax = 5000;
        
        if ( extract_element_uniq(doc,cur2,solverc) != NULL ) {
            method->loptions->solver = removeSpace(extract_element_uniq(doc,cur2,solverc));
            if (( strcmp(method->loptions->solver, "nl2sol.dn2gb")!=0) && 
		( strcmp(method->loptions->solver, "nl2sol.dn2fb")!=0) &&
                (strcmp(method->loptions->solver, "dhc")!=0) &&
                (strcmp(method->loptions->solver, "misqp")!=0) &&
                (strcmp(method->loptions->solver, "")!=0))
            {
                perror(error13);
                exit(13);
            }
        } else  method->loptions->solver = "";     
        
        if (( extract_element_uniq(doc,cur2,finishc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,finishc)),"default")!=0)){
            method->loptions->finish = removeSpace(extract_element_uniq(doc,cur2,finishc));
        } else  method->loptions->finish = NULL;

        if (( extract_element_uniq(doc,cur2,bestxc) != NULL )  &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,bestxc)),"default")!=0)){
            method->loptions->bestx = atoi(extract_element_uniq(doc,cur2,bestxc)); 
        } else method->loptions->bestx = -1;

        if ((  extract_element_uniq(doc,cur2,merit_filterc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,merit_filterc)),"default")!=0)){
            method->loptions->merit_filter = atoi(extract_element_uniq(doc,cur2,merit_filterc)); 
        } else method->loptions->merit_filter = -1;

        if (( extract_element_uniq(doc,cur2,distance_filterc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,distance_filterc)),"default")!=0)){
            method->loptions->distance_filter = atoi(extract_element_uniq(doc,cur2,distance_filterc)); 
        } else method->loptions->distance_filter = -1;

        if ((   extract_element_uniq(doc,cur2,thfactorc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,thfactorc)),"default")!=0)){
            method->loptions->thfactor = atof(extract_element_uniq(doc,cur2,thfactorc));
        } else method->loptions->thfactor = -1.0;

        if ((  extract_element_uniq(doc,cur2,maxdistfactorc) != NULL ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2, maxdistfactorc)),"default")!=0)){
            method->loptions->maxdistfactor = atof(extract_element_uniq(doc,cur2,maxdistfactorc));
        } else method->loptions->maxdistfactor = -1.0;

        if (( extract_element_uniq(doc,cur2,wait_maxdist_limitc) != NULL  ) &&
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,wait_maxdist_limitc)),"default")!=0)){
            method->loptions->wait_maxdist_limit = atoi(extract_element_uniq(doc,cur2,wait_maxdist_limitc)); 
        } else method->loptions->wait_maxdist_limit = -1;

        if (( extract_element_uniq(doc,cur2,wait_th_limitc) != NULL  ) && 
            (strcmp(removeSpace(extract_element_uniq(doc,cur2,wait_th_limitc)),"default")!=0)){
            method->loptions->wait_th_limit = atoi(extract_element_uniq(doc,cur2,wait_th_limitc)); 
        } else method->loptions->wait_th_limit = -1;

    } else {
        method->loptions = (local_options *) malloc(sizeof(local_options));
        method->loptions->finish = NULL;
        method->loptions->n1 = -1;
        method->loptions->n2 = -1;
        method->loptions->tol = -1;
        method->loptions->iterprint = -1;
        method->loptions->balance = -1;
        method->loptions->evalmax = 5000;
        method->loptions->solver = "";
        method->loptions->bestx = -1;
        method->loptions->merit_filter = -1;
        method->loptions->distance_filter = -1;
        method->loptions->tol = -1;
        method->loptions->thfactor = -1.0;
        method->loptions->maxdistfactor = -1.0;
        method->loptions->wait_maxdist_limit = -1;
        method->loptions->wait_th_limit = -1;
    }

    return 1;
}

int extract_element_parallelization(xmlDocPtr doc, xmlNodePtr *root, parallelization_strategy *parallel) {
    xmlNodePtr cur, init_cur; 
        
    const char *name_exp;    
    const char *name_evals_threshold;
    const char *name_mult_num_sendSol;
    const char *name_minimum_num_sendSol;
    const char *name_max_time_ite;
    const char *name_reception_threshold;
    name_exp= "parallelization";    
    name_evals_threshold= "evals_threshold";
    name_mult_num_sendSol= "mult_num_sendSol";
    name_minimum_num_sendSol= "minimum_num_sendSol";
    name_max_time_ite = "migration_max_time";
    name_reception_threshold= "reception_threshold";
    
    init_cur = (*root)->xmlChildrenNode;
    cur = extract_init_node(init_cur, name_exp);
        
    
    parallel->reception_threshold = atof(extract_element_uniq(doc,cur,name_reception_threshold));
    parallel->evals_threshold = atoi(extract_element_uniq(doc,cur,name_evals_threshold));
    parallel->mult_num_sendSol = atoi(extract_element_uniq(doc,cur,name_mult_num_sendSol));
    parallel->minimum_num_sendSol = atoi(extract_element_uniq(doc,cur,name_minimum_num_sendSol));    
    parallel->max_time_ite =  atof(extract_element_uniq(doc,cur,name_max_time_ite));
     
    return 1;
}



int load_configuration_XML(char *docname, experiment_total *exptotal){
    xmlDocPtr doc; 
    xmlNodePtr root,init_cur, cur;
    xmlChar *value;
    const char *name_experiment, *island, *cooperative;
    const char *method, *name, *parallelization, *problems;
    int exit1;
    int NPROC, id;
    
    id = exptotal->execution.idp;
    NPROC = exptotal->execution.NPROC;

    problems="problem";        
    method="method";
    name="name";
    parallelization="parallelization";
    island="island";
    cooperative="cooperative";

    doc = xmlParseFile((const char *) docname);
    if (doc == NULL ) printf("Document parsing failed. \n");
    root = xmlDocGetRootElement(doc);     
    if (root == NULL) {
        xmlFreeDoc(doc);
        printf("Document is Empty!!!\n");
        return 0;
    }    
    init_cur = root->xmlChildrenNode;
    // EXTRACT TEST ELEMENTS      
    extract_element_test(doc, &root, &exptotal->test);
    exptotal->execution.nameMatlab = (char *) malloc(100*sizeof(char));
    memmove( exptotal->execution.nameMatlab, strchr(exptotal->test.output_graph, '/') +1, strlen(strchr(exptotal->test.output_graph, '/'))  );
    exptotal->test.namexml = (const char*) docname;
    exptotal->test.init_repetition = 0;

       
    // EXTRACT METHOD ELEMENTS
    cur = extract_init_node(init_cur, method);
    if (cur != NULL) {
        exit1 = 0;
        value = NULL;
        while (cur != NULL || exit1 == 0) {
            if ((!xmlStrcmp(cur->name, (const xmlChar *) method))) {
                value = xmlGetProp(cur, (const xmlChar *) name);
                exit1 = 1;
            }
            cur = cur->next;
        }
        if ((strcmp((const char *) value, "ScatterSearch")==0) || (strcmp((const char *) value, "CeSS") == 0) 
                ||  (strcmp((const char *) value, "saCeSS") == 0) ||  (strcmp((const char *) value, "aCeSS_dist") == 0) ||
                (strcmp((const char *) value, "eSSm") == 0 ) ||  (strcmp((const char *) value, "coSHADE") == 0)) {

            exptotal->methodScatterSearch = (experiment_method_ScatterSearch *) malloc(sizeof (experiment_method_ScatterSearch));
            extract_element_method_ScatterSearch(doc, &root, exptotal->methodScatterSearch);                        
            exptotal->methodDE=NULL;
            exptotal->methodScatterSearch->eSSversion = (char *) calloc(500, sizeof(char));
            
            if (strcmp((const char *) value, "ScatterSearch") == 0) {
                exptotal->methodScatterSearch->eSSversion = "ScatterSearch";
                #ifdef MPI2
                perror(error31);
                exit(31);
                #endif
            } 
            else if (strcmp((const char *) value, "CeSS") == 0){
                exptotal->methodScatterSearch->eSSversion = "CeSS";
                #ifndef MPI2
                perror(error32);
                exit(32);
                #endif        
                #ifndef OPENMP
                perror(error33);
                exit(33);                
                #endif
            }  
            else if (strcmp((const char *) value, "saCeSS") == 0){
                exptotal->methodScatterSearch->eSSversion = "saCeSS";
                #ifndef MPI2
                perror(error32);
                exit(32);
                #endif    
                #ifndef OPENMP
                perror(error33);
                exit(33);                
                #endif
            }             
            else if (strcmp((const char *) value, "aCeSS_dist") == 0){
                exptotal->methodScatterSearch->eSSversion = "aCeSS_dist";
                #ifndef MPI2
                perror(error32);
                exit(32);
                #endif   
                #ifndef OPENMP
                perror(error33);
                exit(33);                
                #endif                
            }
            else if (strcmp((const char *) value, "eSSm") == 0){
                exptotal->methodScatterSearch->eSSversion = "eSSm";
                #ifndef MPI2
                perror(error32);
                exit(32);
                #endif       
                #ifndef OPENMP
                perror(error33);
                exit(33);                
                #endif                
            } 
            else if (strcmp((const char *) value, "coSHADE") == 0){
	        exptotal->methodScatterSearch->eSSversion = "eSSm";
                #ifndef MPI2
		perror(error32);
		exit(32);
		#endif
		#ifndef OPENMP
                perror(error33);
                exit(33);
                #endif
            }
            else {
                perror(error12);
                exit(12);        
            }

        }        
        
    } else return 0;
    
    
    // EXTRACT PARALLELIZATION ELEMENTS        
    cur = extract_init_node(init_cur, parallelization);

#if  defined(OPENMP) || defined(MPI2)
    if (cur != NULL) {
        exit1 = 0;
        
        while (cur != NULL || exit1 == 0) {
            if ((!xmlStrcmp(cur->name, (const xmlChar *) parallelization))) {
                value = xmlGetProp(cur, (const xmlChar *) name);
                exit1 = 1;
            }
            cur = cur->next;
        }
        if (strcmp((const char *) value, island) == 0 
                ||
                strcmp((const char *) value, cooperative) == 0) {
            exptotal->par_st = (parallelization_strategy *) malloc(sizeof(parallelization_strategy));
            
            if (strcmp((const char *) value, island) == 0) exptotal->par_st->type=(char *) island;
            else if (strcmp((const char *) value, cooperative) == 0) exptotal->par_st->type=(char *) cooperative;
            
            extract_element_parallelization(doc, &root, exptotal->par_st);
            
            exptotal->par_st->NPROC = NPROC;
                #if OPENMP
            #pragma omp parallel shared(NPROC) 
            {       
                exptotal->par_st->NPROC_OPENMP =  omp_get_max_threads();
            }
    #endif
        }
        else {
            perror(error9);
            exit(9);
        }
    } else {
        perror(error27);
        exit(27);
    }
    #else
        exptotal->par_st = NULL;
    #endif


    // EXTRACT ELEMENT PROBLEM  
    cur = extract_init_node(init_cur, problems);
    if (cur != NULL) {
            extract_element_problem(doc, &root, exptotal,  &exptotal->test);
    }
    xmlFree(value);       
    
    return 1;    
    
}


int extract_element_problem(xmlDocPtr doc, xmlNodePtr *root, experiment_total *exptotal,experiment_testbed *test ) {
    xmlNodePtr cur, init_cur, children_problem, children_points, curchildren;
    const char *name_exp,*atrib;
    const char *name_vtr;
    const char *name_dim, *int_var, *bin_var, *ineq, *neq;
    const char *name_sol, *name_bounds, *name_ub, *name_lb, *name_r;
    char *pt;
    int counter;
    char *max_dom, *min_dom;
    char **vector_char;
    int i;    
    const char *cl_s, *cu_s;
    char *max_const,*min_const;
    int SIZE_STRING;
    SIZE_STRING=1000000;
    cl_s="cl";
    cu_s="cu";
    name_dim="dim";
    name_exp= "problem";
    name_bounds= "bounds";
    name_sol="point";
    name_ub="ub";
    name_lb="lb";
    name_vtr="vtr";
    atrib="fx";
    int_var = "int_var";
    bin_var = "bin_var";
    ineq = "ineq";
    neq  = "neq";

    init_cur = (*root)->xmlChildrenNode;
    cur = extract_init_node(init_cur, name_exp);
    
    if (cur == NULL) {
        perror(error34);
        error(34);
    }

    if (extract_element_uniq(doc,cur,int_var) == NULL) {
        perror(error25);
        exit(25);
    }
    test->int_var =  atoi(extract_element_uniq(doc,cur,int_var));

    if (extract_element_uniq(doc,cur,bin_var) == NULL) {
        perror(error26);
        exit(26);
    }
    test->bin_var = atoi(extract_element_uniq(doc,cur,bin_var));

    if (extract_element_uniq(doc,cur,ineq) == NULL) {
        perror(error24);
        exit(24);
    }
    test->ineq = atoi(extract_element_uniq(doc,cur,ineq));

    if (extract_element_uniq(doc,cur,neq) == NULL) {
        perror(error24);
        exit(24);
    }
    test->neq = atoi(extract_element_uniq(doc,cur,neq));


        
//    if (extract_element_uniq(doc, cur, name_vtr) == NULL) {
//        perror(error30);
//       exit(30);
//    }
   
//    test->VTR = atof(extract_element_uniq(doc,cur,name_vtr));
	
    if (extract_element_uniq(doc, cur, name_dim) == NULL) {
        perror(error35);
        exit(35);
    }          
    test->bench.dim = atoi(extract_element_uniq(doc, cur, name_dim)); 
    if (extract_element_uniq(doc, cur, name_ub) == NULL) {
        perror(error36);
        exit(36);
    }

//MAX DOMAIN
    max_dom = (char *) malloc(SIZE_STRING*sizeof(char));
    max_dom=extract_element_uniq(doc, cur, name_ub);
    if (max_dom == NULL) {
        perror(error38);
        exit(38);
    }
    test->bench.max_dom = (double *) malloc(test->bench.dim* sizeof (double) );
    counter=0;
    pt = strtok (max_dom,",");
    while (pt != NULL) {
        test->bench.max_dom[counter] = atof(pt);
        pt = strtok (NULL, ",");
	counter++;
    }    		
    if (counter != test->bench.dim){
        perror(error37);
        exit(37);
    }
    free(max_dom);

// MIN DOMAIN
    min_dom = (char *) malloc(SIZE_STRING*sizeof(char));
    min_dom=extract_element_uniq(doc, cur, name_lb);
    if (min_dom == NULL) {
        perror(error38);
        exit(38);	
    }
    test->bench.min_dom = (double *) malloc(test->bench.dim* sizeof (double) );
    counter=0;
    pt = strtok (min_dom,",");
    while (pt != NULL) {
        test->bench.min_dom[counter] = atof(pt);
        pt = strtok (NULL, ",");
        counter++;
    }
    if (counter != test->bench.dim){
        perror(error37);
        exit(37);
    }
    free(min_dom);
// CONSTRAINTS MIN
    if (extract_element_uniq(doc, cur, cl_s) != NULL){
    min_const = (char *) malloc(SIZE_STRING*sizeof(char));
    min_const = extract_element_uniq(doc, cur, cl_s);
    if (min_const == NULL) {
        perror(error38);
        exit(38);
    }
    test->bench.CL = (double *) malloc((*exptotal).test.ineq * sizeof (double) );
    counter=0;
    pt = strtok (min_const,",");
    while (pt != NULL) {
        test->bench.CL[counter] = atof(pt);
        pt = strtok (NULL, ",");
        counter++;
    }
    if (counter != (*exptotal).test.ineq){
        perror(error37);
        exit(37);
    }
    free(min_const);
    }
// CONSTRAINTS MAX
    if (extract_element_uniq(doc, cur, cu_s) != NULL){
    max_const = (char *) malloc(SIZE_STRING*sizeof(char));
    max_const = extract_element_uniq(doc, cur, cu_s);
    if (max_const == NULL) {
        perror(error38);
        exit(38);
    }
    test->bench.CU = (double *) malloc((*exptotal).test.ineq * sizeof (double) );
    counter=0;
    pt = strtok (max_const,",");
    while (pt != NULL) {
        test->bench.CU[counter] = atof(pt);
        pt = strtok (NULL, ",");
        counter++;
    }
    if (counter != (*exptotal).test.ineq){
        perror(error37);
        exit(37);
    }
    free(max_const);
    }

// INIT POINT
    if (extract_element_uniq(doc, cur, name_sol) != NULL){
	test->bench.number_init_sol = count_element_multi(doc, cur, name_sol);
	
        test->bench.X0 = (double **) malloc( test->bench.number_init_sol * sizeof(double *) );
        test->bench.F0 = (double *) malloc( test->bench.number_init_sol * sizeof(double));	
        
        vector_char = (char **) malloc(test->bench.number_init_sol* sizeof(char *) );
	for (i=0;i<test->bench.number_init_sol;i++) {
	     vector_char[i] = (char *) malloc(100000* sizeof(char));
	     test->bench.X0[i] = (double *) malloc(test->bench.dim*sizeof(double));
	}
	extract_element_multi(doc,cur,name_sol,vector_char);
        for (i=0;i<test->bench.number_init_sol;i++) {
    		counter=0;
    		pt = strtok (vector_char[i],",");
    		while (pt != NULL) {
       	     		test->bench.X0[i][counter] = atof(pt);
	        	pt = strtok (NULL, ",");
	        	counter++;
	    	}
	}

        extract_atribute_multi(doc,cur,name_sol,test->bench.F0,atrib);

    } else  {
	test->bench.number_init_sol=0;
    }
  
    for (i=0;i<test->bench.number_init_sol;i++) {
	free(vector_char[i]);
    }
   // free(vector_char);
 
    return 1;
}











