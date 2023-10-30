#include <structure_paralleltestbed.h>


void* examplefunction(double *x, void *data) {
    experiment_total *exp1;
    output_function *res;
    double F, dim, *g;
    double k1,k2,k3,k4;

    k1=0.09755988;
    k3=0.0391908;
    k2=0.99*k1;
    k4=0.9*k3;

    exp1 = (experiment_total *) data;
//  DIMENSION
    dim = 6;
//  SET THE VARS
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));
    g = (double *) calloc(4,sizeof(output_function));
//  f=-x(4);
    F = -x[3];
//  Equality constraints :
//  g(1)=x(4)-x(3)+x(2)-x(1)+k4*x(4).*x(6);
    g[0]=x[3]-x[2]+x[1]-x[0]+k4*x[3] *x[5];

//  g(2)=x(1)-1+k1*x(1).*x(5);
    g[1]=x[0]-1+k1*x[0]*x[4];
              
//  g(3)=x(2)-x(1)+k2*x(2).*x(6);
    g[2]=x[1]-x[0]+k2*x[1] *x[5];
                      
//  g(4)=x(3)+x(1)-1+k3*x(3).*x(5);
    g[3]=x[2]+x[0]-1+k3*x[2] *x[4];
                      		
//  Inequality constraint
//  g(5)=x(5).^0.5    +    x(6).^0.5;
    g[4]=pow(x[4],0.5)+pow(x[5],0.5);
                      				
//  PUT IN THE output_function
    res->value = F;
    res->g = g;
           				            
    return res;
}


void residuals(double *x, double *R, int *nr, void *exp) {
    experiment_total *exp1;
    int D;
    exp1 = (experiment_total *) exp;
    D=(*exp1).test.bench.dim; 
    // CODE RESIDUALS
}

void jacobian(double *x, double *J, int *nJ, void *exp) {
    experiment_total *exp1;
    int D;
    exp1 = (experiment_total *) exp;
    D=(*exp1).test.bench.dim; 
    // CODE JACOBIAN 
}
