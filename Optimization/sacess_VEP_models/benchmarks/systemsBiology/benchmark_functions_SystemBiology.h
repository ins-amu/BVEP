#include <structure_paralleltestbed.h>
#include <AMIGO_problem.h>

void copylocalexperimentsb_(experiment_total *, experiment_total *, int);

void*  functiontest( double *,  void*);

void*  functiontest2( double *,  void*);

void* evalSB_( double *,  void*);


const char* return_benchmark_SystemBiology(int);

int load_benchmark_SystemBiology(experiment_total *);

int destroySystemBiology(experiment_total *);

void init_point_SystemBiology_Problems_matrix_(void * ,  double *, int, int, int);

int load_benchmark_test(experiment_total *, double *, double *);

int initializetest(experiment_total *);

