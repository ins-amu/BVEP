/**
 * @file evaluationinterface.c
 * @author David R. Penas
 * @brief File containing functions about the evaluation of objective function.
 */


#include <structure_paralleltestbed.h>
#include <configuration.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <def_errors.h>
#include <configuration.h>
#include <time.h>  
#include <common_solver_operations.h>
#include <solversinterface.h>

void deallocateoutputfunction_(output_function *res, void * exp) {
    experiment_total *exp1;

    exp1 = (experiment_total *) exp;

    if (res->size_r > 0) free(res->R);
    if (res->size_j > 0) free(res->J);
    if ((exp1[0].test.ineq+exp1[0].test.neq)>0) free(res->g);
    free(res);

}


void DE_correction_bounds(double *X, int D, double* Xl, double* Xu) {
    register int i;

    for (i = 0; i < D; i++) {
    	if (Xl[i] == Xu[i]) {
            X[i] = Xl[i];
        } else {
            while (X[i] < Xl[i] || X[i] > Xu[i]) {
                if (X[i] < Xl[i]) X[i] = 2.0 * Xl[i] - X[i];
                if (X[i] > Xu[i]) X[i] = 2.0 * Xu[i] - X[i];
            }
        }
    }
 
}


void DE_correction_bounds2(double *X, int D, double* Xl, double* Xu) {
    register int i;

    for (i = 0; i < D; i++) {
    	if (Xl[i] == Xu[i]) {
            X[i] = Xl[i];
        } else {
            if (X[i] < Xl[i]) {
                X[i] = Xl[i];
            }
            else if (X[i] > Xu[i]) {
                X[i] = Xu[i];
            }
        }
    }
 
}


double callfitnessfunction_(void *(*fitnessfunction)(double*, void *), void *exp1_, double *U, int *D, double *Xl, double *Xm) {

    experiment_total *exp1;
    output_function *res;
    res = (output_function *) malloc(sizeof(output_function));
    exp1 = (experiment_total *) exp1_;

    
    
    DE_correction_bounds(U, *D, Xl, Xm);

    if (exp1[0].test._log == 1)
        converttonormal_(U, D);

    
    res = (output_function *) fitnessfunction(U,  &(exp1[0]) );
    U[*D] = res->value;
    
    if (exp1[0].test._log == 1)
        converttolog_(U, D);

    deallocateoutputfunction_(res,exp1);
 
    return U[*D]; 
}


void translation( double *X, double *transconst, int size ) {
	int i;

	for (i=0;i<size;i++) {
		if ( transconst[i] != 0.0 ) {
			X[i] = X[i] + transconst[i];
		}
	}
}

void detranslation( double *X, double *transconst, int size ) {
        int i;

        for (i=0;i<size;i++) {
                if ( transconst[i] != 0.0 ) {
                        X[i] = X[i] - transconst[i];
                }
        }
}

void detranslationinterface_( double *U, void *exp1_ ) {
    experiment_total *exp1;



    exp1 = (experiment_total *) exp1_;
    if  (exp1->test.bench.translation == 1) detranslation( U, exp1->execution.transconst, exp1->test.bench.dim );

}


double callfitnessfunctionfortran_(void *(*fitnessfunction)(double*, void *), void *exp1_, double *U, int *D, double *Xl, double *Xm,double*nlc) {
    experiment_total *exp1;
    int i;
    output_function *res;
    double value;
    
    res = NULL;    
    res = (output_function *) malloc(sizeof(output_function));
    
    exp1 = (experiment_total *) exp1_;

    DE_correction_bounds2(U, *D, Xl, Xm);

    if  (exp1->test.bench.translation == 1) detranslation( U, exp1->execution.transconst, exp1->test.bench.dim );

    res = (output_function *) fitnessfunction(U, & (exp1[0]) );
    
    if (res->g != NULL ) {
        for (i=0;i<(exp1[0].test.ineq+exp1[0].test.neq);i++) {
                nlc[i] = res->g[i];
        }
    }
            
    if  (exp1->test.bench.translation == 1) translation( U, exp1->execution.transconst, exp1->test.bench.dim );

    value = res->value;
    
    
    deallocateoutputfunction_(res,exp1);
 
    return value; 
}


double callfitnessfunctionfortranopenmp_(void *(*fitnessfunction)(double*, void *), void *exp1_, 
                double *U, int *D, double *Xl, double *Xm,double*nlc, int *idp) {
    
    experiment_total *exp1;
    int i;
    output_function *res;
    double value;
    
   
    exp1 = (experiment_total *) exp1_;

 
    DE_correction_bounds2(U, *D, Xl, Xm);


    
    if  (exp1->test.bench.translation == 1) {
	detranslation( U, exp1->execution.transconst, exp1->test.bench.dim );
    }

    res = (output_function *) fitnessfunction(U, &(exp1[*idp]));
    if (res->g != NULL ) {
        for (i=0;i<(exp1[0].test.ineq+exp1[0].test.neq);i++) {
                nlc[i] = res->g[i];
        }
    }

    if  (exp1->test.bench.translation == 1) translation( U, exp1->execution.transconst, exp1->test.bench.dim );

    value=res->value;

    deallocateoutputfunction_(res,exp1);
 
    return value; 
}






double callfitnessfunctionfortranopenmp2_(void *(*fitnessfunction)(double*, void *), void *exp1_, double *U,double*nlc, int *idp) {
    
    experiment_total *exp1;
    int i;
    output_function *res;
    double value;
    
    
    exp1 = (experiment_total *) exp1_;

    if  (exp1->test.bench.translation == 1)  detranslation( U, exp1->execution.transconst, exp1->test.bench.dim );


    res = (output_function *) fitnessfunction(U, &(exp1[*idp]));
   
    if (res->g != NULL ) {
        for (i=0;i<(exp1[0].test.ineq+exp1[0].test.neq);i++) {
                nlc[i] = res->g[i];
        }
    }

    if  (exp1->test.bench.translation == 1)  translation( U, exp1->execution.transconst, exp1->test.bench.dim );
    value=res->value;

    deallocateoutputfunction_(res,exp1);
   
    return value; 
}


double callfitnessfunctionopenmp_(void *(*fitnessfunction)(double*, void *), void *exp1_, double *U, int *D, double *Xl, double *Xm, int *idp) {
    output_function *res;
    
    res = (output_function *) malloc(sizeof(output_function) );
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    DE_correction_bounds(U, *D, Xl, Xm);
    if (exp1[*idp].test._log == 1)
        converttonormal_(U, D);
    res = (output_function *) fitnessfunction(U, &(exp1[*idp]));
    U[*D] = res->value;
    
    if (exp1[*idp].test._log == 1)
        converttolog_(U, D);
    
    deallocateoutputfunction_(res,exp1);
   
    return U[*D]; 
}

