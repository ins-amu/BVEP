#ifdef MATLAB
#ifdef GNU


#include <structure_paralleltestbed.h>
#include <examplematlab.h>


void allocatevectorrg_(void *omt_p, int *size_R, double *R, int *size_G, double *g) {
    outputmatlab *omt;
    int i;

    omt = (outputmatlab *) omt_p;
    if (*size_R > 0) omt->residual = (double *) malloc(*size_R*sizeof(double));
    if (*size_G > 0) {
	omt->penalty  = (double *) malloc(*size_G*sizeof(double));
    }
  
    if (*size_R > 0) {
       for (i=0;i<*size_R;i++){
	   omt->residual[i]=R[i];
       }
    }
    if (*size_G > 0) {
       for (i=0;i<*size_G;i++){
           omt->penalty[i]=g[i];
       }
    }

    omt->size_residual = *size_R;
    omt->size_penalty  = *size_G;
}





// define the function
void* examplefunctionmatlab(double *x, void *data) {
    experiment_total *exp1;
    output_function *res;
    outputmatlab *omt;
    double y, resultado;
    double *residuo;
    double *g; 
    int i,D,error,size_r,size_g;

      
    exp1 = (experiment_total *) data;
    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));
    omt = (outputmatlab *) malloc(sizeof(outputmatlab));  
    
#ifdef GNU
    __matlabproblem_MOD_matlabobjfunc(x,&D,&resultado,omt,&exp1->execution.ep_matlab);
#elif defined(INTEL)
    matlabproblem_mp_matlabobjfunc_(x,&D,&resultado,residuo,omt,&exp1->execution.ep_matlab);
#endif

    res->value = resultado;
    res->g = (double *) malloc( exp1[0].test.ineq + exp1[0].test.neq * sizeof(double));
    for (i=0;i<(exp1[0].test.ineq + exp1[0].test.neq);i++) {
	res->g[i]=omt->penalty[i];
    }
    if (res->size_r > 0) {
	res->R = (double *) malloc( res->size_r * sizeof(double));
	for (i=0;i<res->size_r;i++) res->R[i] = omt->residual[i];
    }

    free(omt);

    return res;
}




#endif
#endif
