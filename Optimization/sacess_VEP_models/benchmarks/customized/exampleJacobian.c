#include <structure_paralleltestbed.h>


// define the function
void* examplefunction(double *x, void *data) {
    experiment_total *exp1;
    output_function *res;
    double y, sum; 
    int i,D;
    
    exp1 = (experiment_total *) data;
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));
    D=(*exp1).test.bench.dim;
    
    sum = 0.0;
    for (i=0;i<D;i++){
	sum = sum + x[i]*sin(sqrt(abs(x[i])));
    }

    y = 418.9829*D - sum;
    
    

    res->value = y;
    
    return res;
}



//void residuals(double *x, double *R, int *nr, void *exp) {
//    experiment_total *exp1;
//    int D;
//    exp1 = (experiment_total *) exp;
//    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM
	
    // CODE RESIDUALS
//}

//void jacobian(double *x, double *J, int *nJ, void *exp) {
//    experiment_total *exp1;
//    int D;
//    exp1 = (experiment_total *) exp;
//    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM

    // CODE JACOBIAN
//}

