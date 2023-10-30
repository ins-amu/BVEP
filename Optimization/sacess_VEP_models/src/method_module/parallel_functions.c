/**
 * @file parallel_functions.c
 * @author David R. Penas
 * @brief File containing general functions about the parallelization of the 
 * algorithms with MPI.
 */

#ifdef MPI2

#include <structure_paralleltestbed.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <def_errors.h>
#include <AMIGO_problem.h>
#include <float.h>
#include <string.h>

void setnproc_(void *exp, int *NPROC){
    
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    *NPROC = exp1[0].par_st->NPROC;;
}


double initmpi_(){

    return (double) MPI_Wtime();
}

double calctimempi_(void *exp, double *starttime) {
    double current;
    experiment_total *exp1;

    exp1 = (experiment_total *) exp;

    current = (double) MPI_Wtime();


    return (current - *starttime);
}



int returnmigrationsize_(void *exp) {
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    
    return exp1[0].par_st->migration_size;
}




void charge_island_size(experiment_total *exp1, int *tam, int *NP, int *NM, int idp, int coop){
    int NPROC, size;

    NPROC = exp1[0].par_st->NPROC;

    if (coop == 1) {
        *tam = *NP;
        size = (int) floorl((double) (*tam / exp1[0].par_st->migration_size));
        if ((size > 1) && (*tam > size) ) {
            *NM = (int) floorl((double) (*tam / exp1[0].par_st->migration_size));
        } else {
            *NM = 1;
        }

        if (*tam < size) {
            if (idp == 0) {
                perror(error6);
            }
            exit(6);
        }
	*tam = 40;
        *NP = 40;
    }
    else if (coop == 0) {

        *tam = (int) ceil((*NP  / NPROC));
        size = (int) floorl((double) (*tam / exp1[0].par_st->migration_size));
        if ((size > 1) && (*tam > size) ) {
            *NM = (int) floorl((double) (*tam / exp1[0].par_st->migration_size));
        } else {
            *NM = 1;
        }
        if (*tam < size) {
            if (idp == 0) {
                perror(error6);
            }
            exit(6);
        }
    }
    

}



void chargeid_(void *exp1_, int *idp ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    *idp = exp1[0].execution.idp;
}



int test_stop (void *array, int *stop, int NPROC) {
    int flag, l;
    MPI_Request *request_recep = (MPI_Request *) array;
        if (*stop < 1) {
            for (l = 0; l < NPROC; l++) {
                MPI_Test(&request_recep[l], &flag, MPI_STATUS_IGNORE);
                if (flag == 1) {
                    *stop = 1;
                    break;
                }
            }
        }
    
    return *stop;
}


int createtopology_(void *exp_){
    int di[1], period[1],reorder,i;
    MPI_Comm ring;
    char const *ringw;
    char const *starw;
    experiment_total *exp;
    
    exp = (experiment_total *) exp_;
    ringw = "ring";
    starw = "star";
    exp->par_st->topology= "star";

    if (strcmp(exp->par_st->topology, ringw) == 0) {
        exp->execution.topology.left = (int *) malloc(sizeof (int));
        exp->execution.topology.rigth = (int *) malloc(sizeof (int));
        exp->execution.topology.num_left = 1;
        exp->execution.topology.num_r = 1;
        di[0] = exp->par_st->NPROC;
        period[0] = 1;
        reorder = 1;
        MPI_Cart_create(MPI_COMM_WORLD, 1, di, period, reorder, &ring);
        exp->execution.topology.comunicator = ring;
        MPI_Cart_shift(exp->execution.topology.comunicator, 0, 1, exp->execution.topology.left, exp->execution.topology.rigth);
    }
    else if (strcmp(exp->par_st->topology, starw) == 0) {
        //MPI_Cart_create(MPI_COMM_WORLD, 1, di, period, reorder, &exp->execution.topology.comunicator);
        if ((exp->execution.idp == 0) && (exp->execution.NPROC>1)) {
                exp->execution.topology.left =  (int *) malloc( (exp->par_st->NPROC - 1 )* sizeof (int));
                exp->execution.topology.rigth = (int *) malloc( (exp->par_st->NPROC - 1 )* sizeof (int));
                exp->execution.topology.num_left = exp->par_st->NPROC - 1;
                exp->execution.topology.num_r = exp->par_st->NPROC - 1;
                exp->execution.topology.comunicator = MPI_COMM_WORLD;    
                for (i=1;i<exp->par_st->NPROC;i++) {
                    exp->execution.topology.rigth[i-1] = i;
                    exp->execution.topology.left[i-1]  = i;
                }
        } 
        else {
                exp->execution.topology.left =  (int *) malloc(sizeof (int));
                exp->execution.topology.rigth = (int *) malloc(sizeof (int));
                exp->execution.topology.num_left = 1;
                exp->execution.topology.num_r = 1;
                exp->execution.topology.comunicator = MPI_COMM_WORLD;   
                exp->execution.topology.rigth[0] = 0; 
                exp->execution.topology.left[0]  = 0;                
        }
        
        
    }
    else {
        perror(error4);
        exit(4);
    }
    
    return 1;
}

int destroy_topology(topology_data *topology) {
    if (topology->left != NULL) {
        free(topology->left);
        topology->left = NULL;
    }
    
    if (topology->rigth != NULL) {
        free(topology->rigth);
        topology->rigth = NULL;
    }

    return 1;
}

void cooperativedist_(void *exp, double *matrix, int *D, double *vector) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Scatter(matrix, exp1[0].execution.tam * (*D + 1), MPI_DOUBLE, vector, exp1[0].execution.tam * (*D + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}

void cooperativedistelement_(void *exp, double *matrix, int *size, double *vector) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Scatter(matrix, *size, MPI_DOUBLE, vector, *size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}



void cooperativegathertelement_(void *exp, double *vect, int *size, double *mat) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    MPI_Gather(vect, *size, MPI_DOUBLE, mat, *size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}


void cooperativegathertelementint_(void *exp, int *vect, int *size, int *mat) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    MPI_Gather(vect, *size, MPI_INT, mat, *size, MPI_INT, 0, MPI_COMM_WORLD);
    
}

void cooperativebcastelement_(void *exp, double *matrix, int *size, int *idp) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Bcast(matrix, *size, MPI_DOUBLE, *idp, MPI_COMM_WORLD);
}


void cooperativebcastelementint_(void *exp, int *matrix, int *size) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Bcast(matrix, *size, MPI_INT,  0, MPI_COMM_WORLD);
}


void returnmaxelementint_(void *exp, int *value, int *valuemax) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Reduce(value, valuemax, 1,  MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
}



void returnminelementint_(void *exp, int *value, int *valuemax) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Reduce(value, valuemax, 1,  MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
}



void returnsumelementint_(void *exp, int *value, int *valuemax) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Reduce(value, valuemax, 1,  MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

void returnsumelementlong_(void *exp, long *value, long *valuemax) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Reduce(value, valuemax, 1,  MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
}

void returnminelement_(void *exp, double *value, double *valuemax) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
    MPI_Reduce(value, valuemax, 1,  MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
}



void returnavgelementintdouble_(void *exp, int *value, double *valueavg, int *NPROC) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    double castdouble, total;
    castdouble = (double) *value;
    MPI_Reduce( &castdouble, &total, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *valueavg  = total / ((double) *NPROC);
}



void returnminlocelement_(void *exp, double *value, double *valuemax, int *loc, int *idp) { 
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    int i;
    double *vals, aux;
    
    MPI_Barrier(MPI_COMM_WORLD);
    vals = (double *) malloc( exp1[0].execution.NPROC * sizeof(double));
    
    MPI_Gather(value, 1 , MPI_DOUBLE, vals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vals, exp1[0].execution.NPROC,MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    
    aux = DBL_MAX;
    for (i=0;i<exp1[0].execution.NPROC;i++) {
        if (vals[i]<aux) {
            aux = vals[i];
            *loc = i;
        }
    }
    *valuemax = aux;
    free(vals);
    vals = NULL;
}


void mpibarrieress_(){
    MPI_Barrier(MPI_COMM_WORLD);
}

void sleepmpi_(int *num){
    sleep(*num);
}

#endif
        
        
        

