#include <structure_paralleltestbed.h>

void* examplefunction(double *x, void *data) {
    experiment_total *exp1;
    output_function *res;
    double x1, F, dim, *g;
    int x2,x3,x4;
 
    exp1 = (experiment_total *) data;

//  DIMENSION
    dim = 4;

//  SET THE VARS
    x1 = x[0];
    x2 = (int) x[1];
    x3 = (int) x[2];
    x4 = (int) x[3];
    
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));
    g = (double *) calloc(3,sizeof(output_function));
//  F = x(2)^2 + x(3)^2 + 2.0*x(1)^2 + x(4)^2 - 5.0*x(2) 
//  	- 5.0*x(3) - 21.0*x(1) + 7.0*x(4);
    F = (double) integer_power(x2,2); //x(2)^2 
    F = F + ((double)integer_power(x3,2)); //x(3)^2 
    F = F + 2.0*pow(x1,2.0); //2.0*x(1)^2
    F = F + ((double)integer_power(x4,2)); //x(4)^2
    F = F - ((double) 5*x2); //- 5.0*x(2)
    F = F - ((double) 5*x3); //- 5.0*x(3)	
    F = F - 21.0*x1; // - 21.0*x(1)
    F = F + ((double)7*x4); //+ 7.0*x(4)	
//  CONSTRAINTS :
//  g(1) = x(2)^2 + x(3)^2 + x(1)^2 + x(4)^2 + x(2) - x(3) + x(1) - x(4);
    g[0]= (double) integer_power(x2,2); // x(2)^2
    g[0]= g[0] + ((double) integer_power(x3,2)) + pow(x1,2.0); // x(3)^2+x(1)^2
    g[0]= g[0] + ((double) integer_power(x4,2)); //  x(4)^2
    g[0]= g[0] + ((double) x2) - ((double) x3);  //  + x(2) - x(3)
    g[0]= g[0] + x1 - ((double) x4); // + x(1) - x(4)

//  g(2) = x(2)^2 + 2.0*x(3)^2 + x(1)^2 + 2.0*x(4)^2 - x(2) - x(4);
    g[1]= (double) integer_power(x2,2); // x(2)^2
    g[1]=g[1]+((double) 2*integer_power(x3,2)); // + 2.0*x(3)^2
    g[1]=g[1]+  pow(x1,2.0); // + x(1)^2 
    g[1]=g[1]+((double) 2*integer_power(x4,2)); //+ 2.0*x(4)^2
    g[1]=g[1]+((double) - x2 - x4 ); // - x(2) - x(4);

//  g(3) = 2.0*x(2)^2 + x(3)^2 + x(1)^2 + 2.0*x(2) - x(3) - x(4);
    g[2]= (double) 2 * integer_power(x2,2); // 2.0*x(2)^2
    g[2]= g[2] + integer_power(x3,2); // + x(3)^2
    g[2]= g[2] + pow(x1,2.0); // + x(1)^2 
    g[2]= g[2] + ((double)  2*x2); //+ 2.0*x(2)
    g[2]= g[2] - x3 -x4; //- x(3) - x(4)

// PUT IN THE output_function     
    res->value = F;
    res->g = g;
    
    return res;
}

int integer_power(int base, unsigned int exp){

    if (exp == 0)
        return 1;
    int temp = integer_power(base, exp/2);
    if (exp%2 == 0)
        return temp*temp;
    else
        return base*temp*temp;
}


void residuals(double *x, double *R, int *nr, void *exp) {
    experiment_total *exp1;
    int D;
    exp1 = (experiment_total *) exp;
    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM
    // CODE RESIDUALS
}

void jacobian(double *x, double *J, int *nJ, void *exp) {
    experiment_total *exp1;
    int D;
    exp1 = (experiment_total *) exp;
    D=(*exp1).test.bench.dim; // DIMENSION OF THE PROBLEM
    // CODE JACOBIAN
}

