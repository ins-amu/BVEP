#include <structure_paralleltestbed.h>


// define the function
// in this case: ROSENBROCK FUNCTION
// with DIMENSION=2 ---> FX = [X0, X1]

void* examplefunction(double *x, void *data) {
    experiment_total *exp1;
    output_function *res;
    double FX,FX2;
    int D;
    const char *error;

    exp1 = (experiment_total *) data;
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));
    D=(*exp1).test.bench.dim;

    if ( D != 2) {
	perror("Wrong dimension of the problem. This example was configured for D=2\n");	
	exit(0);
    }
    FX = 100.0 * pow( x[1] - pow(x[0],2), 2.0)+pow(x[0]-1,2.0); 
    res->value = FX;
    res->size_r=2;
    res->R = (double *) malloc(res->size_r*sizeof(double));
    res->R[0] = 10.0 * (x[1] - pow(x[0],2.0));
    res->R[1] = 1.0 - x[0];
    
    return res;
}

// FUNCTION TO COMPUIIN
void exampleJacobian(double *x, void *exp) {
    experiment_total *exp1;
    int D; 
    exp1 = (experiment_total *) exp;
    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM	

    // Modify next code
    res->size_j = 4;
    res->J = (double *) malloc(res->size_j * sizeof(double));
    J[0] = -20.0 * x[0]; // (0,0)
    J[2] = 10.0;         // (0,1)
    J[1] = -1.0;         // (1,0)
    J[3] = 0.0;          // (1,1)
    
}
