#include <stdlib.h>
#include <string.h>
#include <structure_paralleltestbed.h>
#include <bbobStructures.h>
#include <configuration.h>


int initializeBBOB(experiment_total *exp, const char *namealg){
    int NN;
    
    exp->amigo = NULL;
    exp->param = (ParamStruct *) malloc(sizeof(ParamStruct));
    *(exp->param) = fgeneric_getDefaultPARAMS();
    strcpy(exp->param->dataPath, "TIME");
    strcpy(exp->param->algName, namealg);
    strcpy(exp->param->comments, "Differential_Evolution");
    
    exp->param->DIM = (unsigned int) (*exp).test.bench.dim;
    exp->param->precision = (*exp).test.VTR;
    NN = 1;
    
    (*exp).test.instances =9999; 
    
    return 1;
}

int updateFunctionsBBOB(experiment_total *exp, int ifun, int instance) {
    int i;
    exp->param->funcId = ifun;
    exp->param->instanceId = instance;
    exp->param->DIM = exp->test.bench.dim;
    exp->param->precision = (*exp).test.VTR;
    inicialize_functions((void *) exp);
    
//    exp->test.bench.max_dom = (double *) malloc(exp->test.bench.dim* sizeof (double) );
//    exp->test.bench.min_dom = (double *) malloc(exp->test.bench.dim* sizeof (double) );    
//    for (i=0;i<exp->test.bench.dim;i++) {
//                exp->test.bench.max_dom[i] = 5;
//                exp->test.bench.min_dom[i] = -5;
//    }
    

    return 1;
}

int destroyBBOB(experiment_total *exp){
    finalize_functions(exp);

    return 1;
}


void * fgeneric_noise(double * X, void *data) {
    output_function *res;
    
    res = (output_function *) malloc(1 * sizeof(output_function));
    res->value = fgeneric_evaluate_noise_without_writefile(X,data);
    
    return (void *)res;
}

void * fgeneric_noiseless(double * X, void *data) {
    output_function *res;
    
    res = (output_function *) malloc(1 * sizeof(output_function));
    res->value = fgeneric_evaluate_noiseless_without_writefile(X,data);
    
    return (void *)res;
}


