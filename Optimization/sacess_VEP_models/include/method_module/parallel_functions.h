#include <structure_paralleltestbed.h>

double calctimeMPI(void*,double);

void DE_correction_bounds(double *,int,double*,double*);

#ifdef MPI2
int create_topology(void *, topology_data *, int);
int destroy_topology( topology_data *);
#endif


void chargecooperativeparametersfortran_(void *, int *, int *, int *, double *, long * );

    
void charge_island_size(experiment_total *, int *, int *, int *, int );

int createtopology_(void *);


void chargeid_(void *, int *);


double initmpi_();



int returnmigrationsize_(void *);

void setnproc_(void *, int *);

void cooperativedist_(void *, double *, int *, double *) ;

void cooperativedistelement_(void *, double *, int *, double *);

void cooperativegathertelement_(void *, double *, int *, double *);

void cooperativegathertelementint_(void *, double *, int *, double *);

void cooperativebcastelement_(void *, double *, int *, int *);

void cooperativebcastelementint_(void *, int *, int *);

void returnmaxelementint_(void *, int *, int *);

void returnminelementint_(void *, int *, int *);

void returnminelement_(void *, double *, double *) ;

void returnminlocelement_(void *, double *, double *, int *, int *);

void returnsumelementint_(void *, int *, int *);

void returnsumelementlong_(void *, long *, long *);
