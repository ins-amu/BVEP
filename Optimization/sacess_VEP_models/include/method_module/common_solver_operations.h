#include <structure_paralleltestbed.h>

void logandtranslation_(void *);

double calcSD(double , double *, int ) ;

double calctime(void*,double );

int returnseedcounter_(void *);

void returnseed_(void *, double *);

void initrngrandomserial_(void *);

void initrngrandomparallel_(void *, int *);

void reinitrngrandomparallel_(void *, int *);

double getrngrandomserial_(void *);

double sumvar(double *, int , int);

int allocate_QuickShort(double *, int , int, int);

double* reorderVector_QuickShort(double*, int, int, int);

void replaceWorst(double *, int, int,  double *, int);

void replaceWorstRecp(double *, int, int,  double *, int, double *, double *, int );

void replaceRandom(double *, int, int,  double *, int);

void returnBest(double *, int, int,  double *, int);

void returnBestTabuList(experiment_total , double *, int , int , double *, int , double *, double *, int );

void returnIndv(experiment_total, double *, int, int, double *, int, double *, double *, int );

void replaceIndv(experiment_total, double *, int, int, double *, int, double *, double *, int );

int extract_best_index(double *, int, int);

int initializebenchmarks_(void *, double *, double *, int *) ;

int initializebenchmarksopenmp_(void *, double *, double *, int *, int *) ;

double calc_euclidean_distance(double *, double *, int , double *, double *);

void insert_matrix_tabu_list(double *, double *, int, int, int);

void converttonormal_(double *, int *);

void converttonormal2_(double *, int *, int *);

void converttolog_(double *, int *);

void converttolog2_(double *, int *, int *);

int allocatematrixdouble_(double *, int *, int *);

void setinfinity_(double*);

void setnan_(double*);

void reorder_best(double *, int , int );

int extract_worst_index(double *, int , int );

void setpoint(double *);

