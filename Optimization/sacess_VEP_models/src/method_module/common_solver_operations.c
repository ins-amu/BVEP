#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */


/**
 * @file common_solver_operations.c
 * @author David R. Penas
 * @brief File containing a set of functions with different purposes: arithmetic
 *  operations, management of random numbers, geometrical operations, search 
 * algorithms.
 */

#include <structure_paralleltestbed.h>
#include <evaluationinterface.h>
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
#include <limits.h>
#include "output.h"
#include <unistd.h>

void logandtranslation_(void *exp1) {
    int point, i;
    experiment_total *exp;

    exp = (experiment_total *) exp1;
    
    point = 1.0;
    exp->test.bench.logindex = (int *) malloc( exp->test.bench.dim * sizeof (int));
    for (i = 0; i < exp->test.bench.dim ; i++) exp->test.bench.logindex[i] = 0;

    exp->test.bench.translation = 0;
    if ((*exp).test._log == 1) {
        for (i = 0; i < exp->test.bench.dim; i++) {
            if ((exp->test.bench.min_dom[i] < 0.0 )) {
                exp->test.bench.translation = 1;
                exp->execution.transconst = (double *) malloc(exp->test.bench.dim * sizeof(double) );
                break;
            }
        }
        
        if ((exp->execution.idp == 0)&& (exp->test.bench.translation == 1)) 
            printtranslationmassage();
    }

    for (i = 0; i < exp->test.bench.dim; i++) {
       if (( exp->test.bench.translation == 1)&&(exp->test.bench.min_dom[i] < 0.0)) exp->execution.transconst[i] = point  - (exp->test.bench.min_dom[i]);
       else if (( exp->test.bench.translation == 1)&&(exp->test.bench.min_dom[i] >= 0.0))  exp->execution.transconst[i] = 0.0;

 
       if (exp->test.bench.translation == 1) {
           if ( exp->execution.transconst[i] != 0 ) {
               exp->test.bench.min_dom[i] = point; 
               exp->test.bench.max_dom[i] = exp->test.bench.max_dom[i] + exp->execution.transconst[i];
           }
       }
    }


    if ((*exp).test._log == 1) {
          exp->test.bench.log_max_dom = (double *) malloc( exp->test.bench.dim* sizeof (double));
          exp->test.bench.log_min_dom = (double *) malloc( exp->test.bench.dim* sizeof (double));

          for (i = 0; i < exp->test.bench.dim; i++) {
                if (exp->test.bench.min_dom[i] == 0.0) exp->test.bench.log_min_dom[i] = log(EPSILON_LOG);
                else exp->test.bench.log_min_dom[i] = log(exp->test.bench.min_dom[i]);

                if (exp->test.bench.max_dom[i] == 0.0) exp->test.bench.log_max_dom[i] = log(EPSILON_LOG);
                else exp->test.bench.log_max_dom[i] = log(exp->test.bench.max_dom[i]);

                (*exp).test.bench.logindex[i]=1;
          }
    }
}

double calcSD(double mean, double *population, int sizePOPUL) {
    int i;
    double total;

    total = 0.0;
    for (i = 0; i < sizePOPUL; i++) {
        total = total + powl(fabsl(population[i] - mean), 2.0);
    }
    total = total / (((double) sizePOPUL) - 1.0);
    return sqrt(total);
}

double calctime(void *exp, double starttime){
    double current;
    experiment_total *exp1;

    current = 0.0;
    exp1 = (experiment_total *) exp;
    
    if (exp1->test.output == 1) {    
        current = (double) clock();
    }
    
    return ((double) (current - starttime) / (double) CLOCKS_PER_SEC);
}

int returnseedcounter_(void *exp1_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;  
    
    return exp1[0].contadorseed;
}

void returnseed_(void *exp1_, double *seed) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;  
    memcpy(seed, exp1[0].seed, exp1[0].contadorseed * sizeof(double));
}


void initrngrandomserial_(void *exp1_) {
    experiment_total *exp1;
    long init_seed;
    exp1 = (experiment_total *) exp1_;

    exp1[0].seed = (double *) calloc(1,sizeof(double));
    exp1[0].contadorseed = 1;
    init_seed = time(NULL) * getpid();
    exp1[0].seed[0]  = (double) init_seed;
    printf("SERIAL SEED %lf\n", exp1[0].seed[0]  );
    exp1[0].random = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (exp1[0].random , (unsigned long int) init_seed );
}


void initrngrandomparallel_(void *exp1_, int *idp) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    struct timespec spec;
    long init_seed;
    exp1[0].contadorseed = 1;
    exp1[0].seed = (double *) malloc(1*sizeof(double));
    //clock_gettime(CLOCK_REALTIME, &spec);
    init_seed = time(NULL) * (getpid()*(*idp+1)) ;
    exp1[0].seed[0]  = (double) init_seed;
    exp1[0].random = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (exp1[0].random ,  (unsigned long int)  init_seed  );
}


void reinitrngrandomparallel_(void *exp1_, int *idp) {
    experiment_total *exp1;
    long init_seed;
    exp1 = (experiment_total *) exp1_;

    exp1[0].contadorseed = exp1[0].contadorseed  + 1;
    
    exp1[0].seed = (double *) realloc(exp1[0].seed, exp1[0].contadorseed*sizeof(double));
    init_seed = time(NULL) * ((*idp+1)) ;

    exp1[0].seed[exp1[0].contadorseed-1]  = init_seed;
    gsl_rng_free(exp1[0].random);
    exp1[0].random = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (exp1[0].random ,  (unsigned long int)  init_seed );
}

double getrngrandomserial_(void *exp1_) {
    experiment_total *exp1;
    double res;
    
    exp1 = (experiment_total *) exp1_;
    
    res = gsl_rng_uniform(exp1[0].random);
    return res;
}

int allocatematrixdouble_(double *vector, int *fil, int *col) {
    
    vector = (double *)malloc((*fil)*(*col)*sizeof(double));
    
    return 1;
}



void converttonormal_(double *vector, int *D) {
    int i;
    for (i = 0; i < *D; i++) {
        if (vector[i] == EPSILON_LOG) vector[i] = 0.0;
        else     vector[i] = exp(vector[i]);
    }
}

void converttonormal2_( double *vector, int *D, int *index) {
    int i;
    for (i = 0; i < *D; i++) {
          if ( index[i] == 1 ) {
                if (vector[i] == EPSILON_LOG) vector[i] = 0.0;
                else     vector[i] = exp(vector[i]);
        }
    }
    
}


void converttolog_(double *vector, int *D) {
    int i;

    for (i = 0; i < *D; i++) {
        if (vector[i] == 0.0) vector[i] = log(EPSILON_LOG);
        else vector[i] = log(vector[i]);
    }
}

void converttolog2_(double *vector, int *D, int *index) {
    int i;

    for (i = 0; i < *D; i++) {
	if ( index[i] == 1 ) {
                if (vector[i] == 0.0) vector[i] = log(EPSILON_LOG);
                else vector[i] = log(vector[i]);
        }
    }
}

void convtrans_(double *vector, int *D, void *exp1_) {

    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    int i;
 
    if ( exp1->test.bench.translation == 1) translation( vector, exp1->execution.transconst, exp1->test.bench.dim );
/*
    for (i = 0; i < *D; i++) {
                if (vector[i] == 0.0) vector[i] = log(EPSILON_LOG);
                else vector[i] = log(vector[i]);
    }
*/
}


void coord_transformation(double *OLD, double *New, int D, double *UB, double *LB) {
    int i;

    for (i = 0; i < D; i++) {
        if (OLD[i] == LB[i]) New[i] = 0.0;
        else if ((fabs(LB[i]) + fabs(UB[i])) == 0) New[i] = 0.0;
        else if (LB[i] == UB[i]) New[i] = 0.0;
        else New[i] = (OLD[i] - LB[i]) / (fabs(LB[i]) + fabs(UB[i]));
    }


}

int extract_worst_index(double *populLocal, int tam, int D) {
    double max_value;
    int i, index;

    // MIRAMOS O MELLOR
    max_value = populLocal[0 * (D + 1) + D];
    index = 0;
    for (i = 0; i < tam; i++) {
        if (populLocal[i * (D + 1) + D] > max_value) {
            max_value = populLocal[i * (D + 1) + D];
            index = i;
        }
    }


    return index;
}

void insert_matrix_tabu_list(double *U, double *matrix, int D, int NP, int contador) {
    int i, index;

    if (contador >= NP) {
        index = extract_worst_index(matrix, NP, D);
        for (i = 0; i < (D + 1); i++)
            matrix[index * (D + 1) + i] = U[i];
    } else {
        for (i = 0; i < (D + 1); i++)
            matrix[contador * (D + 1) + i] = U[i];
    }
}



double calc_euclidean_distance(double *U1, double *U2, int D, double *UB, double *LB) {
    int i;
    double sum;
    double *U1_NEW, *U2_NEW;
    double *tempU1, *tempU2;
    

    tempU1 = (double *) malloc(D* sizeof (double));
    tempU2 = (double *) malloc(D* sizeof (double));
    
    memmove(tempU1,U1,D*sizeof(double));
    memmove(tempU2,U2,D*sizeof(double));
    
    U1_NEW = (double *) malloc(D* sizeof (double));
    U2_NEW = (double *) malloc(D* sizeof (double));

    coord_transformation(tempU1, U1_NEW, D, UB, LB);
    coord_transformation(tempU2, U2_NEW, D, UB, LB);
    sum = 0.0;


    
    for (i = 0; i < D; i++) {
        sum = sum + pow(U1_NEW[i] - U2_NEW[i], 2.0);      
    }
    
    if (tempU1!=NULL) {
        free(tempU1);
        tempU1 = NULL;
    }
    
    if (tempU2!=NULL) {
        free(tempU2);
        tempU2 = NULL;
    }
    
    if (U1_NEW != NULL) {
        free(U1_NEW);
        U1_NEW = NULL;
    }
    
    if (U2_NEW != NULL) {
        free(U2_NEW);
        U2_NEW = NULL; 
    }
    
    return sqrt(fabs(sum));
}

int initializebenchmarks_(void *exp1_, double *Xl, double *Xm, int *D) {
    int j;
    char const *noiseBBOB;
    char const *noiselessBBOB;
    char const *system;
    char const *LSGO;

    noiseBBOB = "noiseBBOB";
    noiselessBBOB = "noiselessBBOB";
    system = "systemBiology";
    LSGO = "LSGO";
    
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    if ((strcmp(exp1[0].test.bench.type, noiseBBOB) == 0) || (strcmp(exp1[0].test.bench.type, noiselessBBOB) == 0)) {
        for (j = 0; j < *D; j++) {
            Xl[j] = (double) exp1[0].test.bench.min_dom[0];
            Xm[j] = (double) exp1[0].test.bench.max_dom[0];
        }
    } else if (strcmp(exp1[0].test.bench.type, system) == 0) {
        if ( exp1[0].test._log == 0 ) {
            for (j = 0; j < *D; j++) {
                Xl[j] = (double) exp1[0].test.bench.min_dom[j];
                Xm[j] = (double) exp1[0].test.bench.max_dom[j];
            }
        } else {
            for (j = 0; j < *D; j++) {
                Xl[j] = (double) exp1[0].test.bench.log_min_dom[j];
                Xm[j] = (double) exp1[0].test.bench.log_max_dom[j];
            }            
        }
    } else if ((strcmp(exp1[0].test.bench.type, LSGO) == 0) ) {
        for (j = 0; j < *D; j++) {
            Xl[j] = (double) exp1[0].test.bench.min_dom[0];
            Xm[j] = (double) exp1[0].test.bench.max_dom[0];
        }
    }
    
    return 1;
}


int initializebenchmarksopenmp_(void *exp1_, double *Xl, double *Xm, int *D, int *idp) {
    int j;
    char const *noiseBBOB;
    char const *noiselessBBOB;
    char const *system;

    noiseBBOB = "noiseBBOB";
    noiselessBBOB = "noiselessBBOB";
    system = "systemBiology";
    
    char const *LSGO;
    LSGO = "LSGO";
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    if ((strcmp(exp1[*idp].test.bench.type, noiseBBOB) == 0) || (strcmp(exp1[*idp].test.bench.type, noiselessBBOB) == 0)) {
        for (j = 0; j < *D; j++) {
            Xl[j] = (double) exp1[*idp].test.bench.min_dom[0];
            Xm[j] = (double) exp1[*idp].test.bench.max_dom[0];
        }
    } else if (strcmp(exp1[*idp].test.bench.type, system) == 0) {
        for (j = 0; j < *D; j++) {
            Xl[j] = (double) exp1[*idp].test.bench.min_dom[j];
            Xm[j] = (double) exp1[*idp].test.bench.max_dom[j];
        }
    } else if ((strcmp(exp1[0].test.bench.type, LSGO) == 0) ) {
        for (j = 0; j < *D; j++) {
            Xl[j] = (double) exp1[*idp].test.bench.min_dom[0];
            Xm[j] = (double) exp1[*idp].test.bench.max_dom[0];
        }
    }
    return 1;
}

/* SUM VAR */
double sumvar(double *Matriz, int sizeMatriz, int sizeVector) {
    int i, j;
    double *vector1, *vector_matrices;
    double avg, suma, varian;

    vector1 = (double *) malloc(sizeMatriz* sizeof (double));
    vector_matrices = (double *) malloc(sizeVector* sizeof (double));
    // VARIANZAS
    for (i = 0; i < sizeVector; i++) {
        suma = 0.0;
        for (j = 0; j < sizeMatriz; j++) {
            vector1[j] = Matriz[j * (sizeVector + 1) + i];
            suma = suma + vector1[j];
        }
        avg = suma / (double) sizeMatriz;
        suma = 0.0;
        for (j = 0; j < sizeMatriz; j++) {
            suma = suma + powl((vector1[j] - avg), 2);
        }
        varian = suma / (double) sizeMatriz;
        vector_matrices[i] = varian;
    }

    // SUMA
    suma = 0.0;
    for (j = 0; j < sizeVector; j++) {
        suma = suma + vector_matrices[j];
    }


    if (vector1 != NULL) {
        free(vector1);
        vector1 = NULL; 
    }
    
    if (vector_matrices != NULL) {
        free(vector_matrices);
        vector_matrices = NULL;
    }
    
    return fabsl(suma);
}

/* ALLOCATE QUICKSHORT */
int allocate_QuickShort(double *v, int b, int NP, int D) {
    int i; 
            //z;
    int pivote;
    double *valor_pivote;
    double *temp;

    valor_pivote = (double *) malloc((D + 1)* sizeof (double));
    temp = (double *) malloc((D + 1)* sizeof (double));

    pivote = b;
    
    memmove(valor_pivote, &v[pivote * (D + 1)], (D+1)*sizeof(double));

    for (i = b + 1; i <= NP; i++) {
        if (v[i * (D + 1) + D] < valor_pivote[D]) {
            pivote++;
            memmove(temp,&v[i * (D + 1)],(D+1)*sizeof(double));
            memmove(&v[i * (D + 1)],&v[pivote * (D + 1)],(D+1)*sizeof(double));
            memmove(&v[pivote * (D + 1)],temp,(D+1)*sizeof(double));
        }
    }

    memmove(temp,&v[b * (D + 1)], (D+1)*sizeof(double));
    memmove(&v[b * (D + 1)],&v[pivote * (D + 1)], (D+1)*sizeof(double));
    memmove(&v[pivote * (D + 1)],temp, (D+1)*sizeof(double));

    if (valor_pivote != NULL) {
        free(valor_pivote);
        valor_pivote=NULL;
    }
    
    if (temp != NULL) {
        free(temp);
        temp=NULL;
    }

    return pivote;
}

/* REORDER QUICK SHORT */
double* reorderVector_QuickShort(double* v, int b, int NP, int D) {
    int pivote;


    if (b < NP) {
        pivote = allocate_QuickShort(v, b, NP, D);
        reorderVector_QuickShort(v, b, pivote - 1, D);
        reorderVector_QuickShort(v, pivote + 1, NP, D);
    }

    return v;
}

/* REPLACE_WORST */
void replaceWorst(double *v, int NP, int D, double *best, int sizeBest) {
    int i;
    int z = 0;


    v = reorderVector_QuickShort(v, 0, NP - 1, D);

    
    
    for (i = (NP - sizeBest); i < NP; i++) {
        memmove(&v[i * (D + 1)], &best[z * (D + 1)], (D+1) * sizeof(double) );
        z++;
    }


}


void replaceWorstRecp(double *v, int NP, int D, double *best, int sizeBest, double *LB, double *UB, int rest) {
    int i,j,z,k;
    int *index_state;
    double *U1, *U2;

    index_state = (int *) malloc(sizeBest*sizeof(int));
    U1 = (double *) malloc((D+1)*sizeof(double));
    U2 = (double *) malloc((D+1)*sizeof(double));
            
    v = reorderVector_QuickShort(v, 0, NP - 1, D);
    
    for (i=0; i<sizeBest; i++) {
        memmove(U1,&best[i * (D + 1)],(D+1-rest)*sizeof(double));
        for (j=0; j<NP; j++) {
            memmove(U2,&v[j * (D + 1)],(D+1-rest)*sizeof(double));
            
            if ( calc_euclidean_distance(U1, U2, D, UB, LB) <= 1e-1 ) {
                index_state[i] = 0;
                break;
            } 
        }
        if (j==NP) {
            index_state[i] = 1;
        }
    }

    //printf("\nINDEX : ");
    //    for (i=0;i<sizeBest;i++)  printf(" %d ", index_state[i]);
    //        printf("\n");
                
    //               printf("\nBEST : ");
    //                   for (i=0;i<sizeBest;i++)  printf(" %lf ", best[i * (D + 1) + D]);
    //                       printf("\n");
            
    z=0;
    k=NP-1;
   
    for (i = (NP - sizeBest); i < NP; i++) {
        if (index_state[z] == 1 ) {
            memmove(&v[k * (D + 1)], &best[z * (D + 1)], (D+1) * sizeof(double) );
            k--;
        }
        z++;
    }
    //printf("\n");

       
   // printf("POPUL : ");
   //     for (i=0;i<NP;i++)  printf(" %lf ", v[i * (D + 1)+D]);
   //         printf("\n"); 
    free(index_state);
    index_state = NULL;
    free(U1);
    U1=NULL;
    free(U2);
    U2=NULL;
}

/* REPLACE RANDOM */
void replaceRandom(double *v, int NP, int D, double *best, int sizeBest) {
    int i, j, punto, valido;
    int *vector_puntos;


    vector_puntos = (int *) malloc(sizeBest* sizeof (int));
    for (i = 0; i < sizeBest; i++) {
        vector_puntos[i] = -1;
    }

    for (i = 0; i < sizeBest; i++) {
        valido = 0;
        while (valido == 0) {
            punto = (int) (NP * URAND_BBOB);
            for (j = 0; j < sizeBest; j++) {
                if (punto == vector_puntos[j]) break;
            }
            if (j == sizeBest) {
                valido = 1;
            }
        }

        vector_puntos[i] = punto;
        for (j = 0; j < (D + 1); j++) {
            v[punto * (D + 1) + j] = best[i * (D + 1) + j];

        }
    }

    if (vector_puntos != NULL) {
        free(vector_puntos);
        vector_puntos=NULL;
    }
}

/* REPLACE BEST */
void returnBest(double *v, int N, int D, double *best, int sizeBest) {
    int i,  cont = 0;
    double *vect;

    vect = (double *) malloc(N * (D + 1) * sizeof (double));

    vect = reorderVector_QuickShort(v, 0, N - 1, D);

   for (i = 0; i < N; i++) {       
       memmove(&v[i * (D + 1)],&vect[i * (D + 1)],(D+1)*sizeof(double));
   }
    
    cont = 0;
    for (i = 0; i < (N - 1); i++) {
        if (cont < sizeBest) {
            memmove(&best[cont * (D + 1)], &vect[i * (D + 1)], (D+1)*sizeof(double));
            cont++;
        } else break;
    }


   // free(vect);
   // vect = NULL;
}

/* REPLACE RANDOM */
void returnRandom(double *v, int NP, int D, double *best, int sizeBest) {
    int i, punto;

    for (i = 0; i < sizeBest; i++) {
        punto = (int) (NP * URAND_BBOB);
        memmove(&best[i * (D + 1)], &v[punto * (D + 1)], (D+1)*sizeof(double));
    }


}



/* R
void returnBestTabuList(experiment_total exp1, double *v, int N, int D, double *best, int sizeBest, double *UB, double *LB, int rest) {
    int i, j, cont = 0;
    double *vect, *U1, *U2;

    vect = (double *) malloc(N * (D + 1) * sizeof (double));
    U1 = (double *) malloc((D+1)*sizeof(double));
    U2 = (double *) malloc((D+1)*sizeof(double));
    
    
    vect = reorderVector_QuickShort(v, 0, N - 1, D);

    for (i = 0; i < N; i++) {       
       memmove(&v[i * (D + 1)],&vect[i * (D + 1)],(D+1)*sizeof(double));
    }
    
    cont = 0;
    for (i = 0; i < (N - 1); i++) {
        if (cont < sizeBest) {
            memmove(U1,&vect[i * (D + 1-rest)],(D+1)*sizeof(double));
            for (j=0; j<exp1.execution.tabusend.size; j++) {
                memmove(U2,&exp1.execution.tabusend.array[j * (D + 1-rest)],(D+1)*sizeof(double));
                if ( calc_euclidean_distance(U1, U2, D, UB, LB) <= 1e-1 ) {
                    break;
                }                 
            }
            if (j == cont) {
                memmove(&best[cont * (D + 1)], &vect[i * (D + 1)], (D+1)*sizeof(double));
                cont++;
            }
        } else break;
    }

    
    if (cont < sizeBest) {
        memmove(U1, best, (D+1)*sizeof(double)); 
        for (i=cont;i<sizeBest;i++) {
            memmove(&best[i * (D + 1)], U1, (D+1)*sizeof(double));
        }
    } else {
        for (i=0;i<sizeBest;i++) {
            insert_matrix_tabu_list(  
                    &best[i*(D+1)], 
                    exp1.execution.tabusend.array,
                    D,exp1.execution.tabusend.size,
                    exp1.execution.tabusend.counter);
        }
                
    }

    
    
        cont = 0;
    for (i = 0; i < (N - 1); i++) {
        if (cont < sizeBest) {
            memmove(&best[cont * (D + 1)], &vect[i * (D + 1)], (D+1)*sizeof(double));
            cont++;
        } else break;
    }
        
        
    free(U1);
    U1 = NULL;
    
    free(U2);
    U2 = NULL;
}
*/

void returnIndv(experiment_total exp1, double *v, int N, int D, double *best, int sizeBest, double *UB, double *LB, int rest) {
    const char *bestword;
    const char *random;

    bestword = "Best";
    random = "Random";
    
    
    if (strcmp(exp1.par_st->SelectionPolicy, bestword) == 0) {
	returnBest(v, N, D, best, sizeBest);
        //returnBestTabuList(exp1,v, N, D, best, sizeBest,UB,LB, rest);
    } else if (strcmp(exp1.par_st->SelectionPolicy, random) == 0) {
        returnRandom(v, N, D, best, sizeBest);
    } else {
        perror(error1);
        exit(2);
    }


}

void replaceIndv(experiment_total exp1, double *v, int N, int D, double *worst, int sizeBest, double *Xl, double *Xm, int rest) {
    const char *worstw;
    const char *random;

    worstw = "Worst";
    random = "Random";

    if (strcmp(exp1.par_st->ReplacePolicy, worstw) == 0) {
        replaceWorst(v, N, D, worst, sizeBest);
        //replaceWorstRecp(v, N, D, worst, sizeBest,Xl,Xm, rest);
    } else if (strcmp(exp1.par_st->ReplacePolicy, random) == 0) {
        replaceRandom(v, N, D, worst, sizeBest);
    } else {
        perror(error2);
        exit(2);
    }
}



int extract_best_index(double *populLocal, int tam, int D) {
    double min_value;
    int i, index;

    // MIRAMOS O MELLOR
    min_value = populLocal[0 * (D + 1) + D];
    index = 0;
    for (i = 0; i < tam; i++) {
        if (populLocal[i * (D + 1) + D] < min_value) {
            min_value = populLocal[i * (D + 1) + D];
            index = i;
        }
    }


    return index;
}



void reorder_best(double *v, int N, int D) {
    int i, j;
    double *vect;

    vect = (double *) malloc(N * (D + 1) * sizeof (double));

    vect = reorderVector_QuickShort(v, 0, N - 1, D);

    for (i = 0; i < N; i++) {
        for (j = 0; j < (D + 1); j++) {
            v[i * (D + 1) + j] = vect[i * (D + 1) + j];
        }
    }
}



long getmaxlong_(){
	return LONG_MAX;
}




void setpoint_( double  *x){
x[0]=0.033451;
x[1]=0.895183;
x[2]=0.755114;
x[3]=0.014721;
x[4]=0.999943;
x[5]=0.105286;
x[6]=0.081112;
x[7]=0.000082;
x[8]=0.982555;
x[9]=1.864532;
x[10]=1.497671;
x[11]=0.100000;
x[12]=2.999893;
x[13]=1.489494;
x[14]=3.000000;
x[15]=0.100000;
x[16]=0.100000;
x[17]=1.000000;
x[18]=0.201452;
x[19]=0.099913;
x[20]=3.000000;
x[21]=0.100000;
x[22]=1.206606;
x[23]=0.323483;
x[24]=1.692353;
x[25]=0.100000;
x[26]=1.004461;
x[27]=0.100000;
x[28]=0.100000;
x[29]=1.042373;
x[30]=0.364734;
x[31]=3.000000;
x[32]=0.405295;
x[33]=0.078115;
x[34]=2.999788;
x[35]=1.499896;
x[36]=0.006831;
x[37]=1.166474;
x[38]=0.104781;
x[39]=3.000000;
x[40]=0.100018;
x[41]=0.011443;
x[42]=2.190206;
x[43]=0.795857;
x[44]=0.098583;
x[45]=2.732832;
x[46]=1.302585;
x[47]=2.797659;
x[48]=1.338829;
x[49]=3.000000;
x[50]=0.101281;
x[51]=1.244177;
x[52]=1.352100;
x[53]=0.008033;
x[54]=2.987493;
x[55]=0.100000;
x[56]=1.000000;
x[57]=1.157116;
x[58]=2.586270;
x[59]=0.100000;
x[60]=0.099913;
x[61]=2.064953;
x[62]=0.504286;
x[63]=0.100000;
x[64]=2.999889;
x[65]=0.117774;
x[66]=2.958100;
x[67]=0.513129;
x[68]=2.575058;
x[69]=0.504680;
x[70]=2.646648;
x[71]=0.889934;
x[72]=0.090826;
x[73]=2.600779;
x[74]=0.172736;
x[75]=2.765911;
x[76]=1.185380;
x[77]=2.987153;
x[78]=1.491439;
x[79]=0.100000;
x[80]=1.002631;
x[81]=0.100000;
x[82]=1.320622;
x[83]=0.101281;
x[84]=2.987720;
x[85]=0.108629;
x[86]=2.996614;
x[87]=0.729074;
x[88]=2.995539;
x[89]=0.266913;
x[90]=0.084206;
x[91]=2.703964;
x[92]=0.420162;
x[93]=0.100000;
x[94]=1.004542;
x[95]=0.178126;
x[96]=2.257473;
x[97]=0.188348;
x[98]=2.494030;
x[99]=1.490460;
x[100]=0.097895;
x[101]=1.675567;
x[102]=0.100000;
x[103]=0.054985;
x[104]=2.999894;
x[105]=0.882348;
x[106]=1.000000;
x[107]=0.100000;
x[108]=0.015688;
x[109]=2.955229;
x[110]=0.317292;
x[111]=2.997790;
x[112]=0.140812;
x[113]=0.100000;
x[114]=2.936550;
x[115]=0.190130;
x[116]=1.503120;
x[117]=0.161634;
x[118]=2.495773;
x[119]=1.254092;
x[120]=0.100000;
x[121]=3.000000;
x[122]=0.305996;
x[123]=3.000000;
x[124]=0.262628;
x[125]=1.078525;
x[126]=0.120045;
x[127]=0.100000;
x[128]=1.228515;
x[129]=1.477455;
x[130]=0.008562;
x[131]=3.000000;
x[132]=0.260938;
x[133]=2.017894;
x[134]=0.126150;
x[135]=0.100000;
x[136]=2.867177;
x[137]=0.327551;
x[138]=1.196621;
x[139]=0.480832;
x[140]=0.005000;
x[141]=2.998170;
x[142]=0.184825;
x[143]=0.005000;
x[144]=1.000000;
x[145]=0.000000;
x[146]=1.000000;
x[147]=1.000000;
x[148]=1.000000;
x[149]=0.000000;
x[150]=0.000000;
x[151]=0.000000;
x[152]=0.000000;
x[153]=0.000000;
x[154]=0.000000;
x[155]=1.000000;
x[156]=0.000000;
x[157]=1.000000;
x[158]=1.000000;
x[159]=1.000000;
x[160]=0.000000;
x[161]=1.000000;
x[162]=0.000000;
x[163]=0.000000;
x[164]=1.000000;
x[165]=1.000000;
x[166]=1.000000;
x[167]=0.000000;
x[168]=0.000000;
x[169]=0.000000;
x[170]=1.000000;
x[171]=0.000000;
x[172]=0.000000;
x[173]=0.000000;
x[174]=0.000000;
x[175]=0.000000;
x[176]=0.000000;
x[177]=0.000000;
x[178]=0.000000;
x[179]=0.000000;
x[180]=1.000000;
x[181]=1.000000;
x[182]=0.000000;
x[183]=1.000000;
x[184]=1.000000;
x[185]=1.000000;
x[186]=1.000000;
x[187]=1.000000;
x[188]=0.000000;
x[189]=0.000000;
x[190]=0.000000;
x[191]=0.000000;
x[192]=1.000000;
x[193]=0.000000;
x[194]=1.000000;
x[195]=0.000000;
x[196]=0.000000;
x[197]=0.000000;
x[198]=0.000000;
x[199]=0.000000;
x[200]=0.000000;
x[201]=1.000000;
x[202]=1.000000;
x[203]=0.000000;
x[204]=0.000000;
x[205]=0.000000;
x[206]=0.000000;
x[207]=0.000000;
x[208]=0.000000;
x[209]=0.000000;
x[210]=1.000000;
x[211]=0.000000;
x[212]=1.000000;
x[213]=1.000000;
x[214]=1.000000;
x[215]=1.000000;
x[216]=1.000000;
x[217]=0.000000;
x[218]=1.000000;
x[219]=1.000000;
x[220]=0.000000;
x[221]=0.000000;
x[222]=0.000000;
x[223]=0.000000;
x[224]=1.000000;
x[225]=0.000000;
x[226]=1.000000;
x[227]=0.000000;
x[228]=0.000000;
x[229]=0.000000;
x[230]=0.000000;
x[231]=0.000000;
x[232]=1.000000;
x[233]=0.000000;
x[234]=0.000000;
x[235]=1.000000;
x[236]=1.000000;
x[237]=1.000000;
x[238]=1.000000;
x[239]=0.000000;
x[240]=0.000000;
x[241]=0.000000;
x[242]=0.000000;
x[243]=1.000000;
x[244]=0.000000;
x[245]=1.000000;
x[246]=0.000000;
x[247]=0.000000;
x[248]=1.000000;
x[249]=0.000000;
x[250]=1.000000;
x[251]=1.000000;
x[252]=1.000000;

}
