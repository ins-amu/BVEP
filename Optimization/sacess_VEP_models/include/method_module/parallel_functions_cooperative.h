#include <structure_paralleltestbed.h>

void setnproc_(void *, int *);

void asynchronousstoppingcriteria_(void *, int *, long *, double *, double *, double *, int *);

void gatherresultsserialize_(void *, double *, double *, int *, double *, double *, long *, long *);

void createcooperativetopology_(void *);

void createcooperativestructs_(void *, int *);

void cooperativedist_(void *, double *, int *, double *);

void cooperativedistelement_(void *, double *, int *, double *);

void cooperativebcastelement_(void *, double *, int *, int *);

void cooperativebcastelementint_(void *, int *, int *);

int destroycooperativestructures_(void *);

void initcooperativestoppingcriteria_(void * );

int checkcooperativemigrationcriteria_(void *);

void migrationcooperativesent_(void *, double *, int *, int *, double *, double *, int *);

void migrationcooperativereception_( void *, void*(*fitnessfunction)(double*, void *),double *, int *, int *, double *, double *, int *);

void asynchronousstoppingcriteriawithrestart_(void *,  int *, long *, double *,double *,int *);

void synchronousstoppingcriteriawithrestart_(void *, int *, long *, double *, double *, int *);

int cooperativempitest_(void *, int *);

void gatherresults_(void *, double *, double *, int *, double *, double *, long *, long *) ;

void migrationsynchcooperativesent_(void *, int *, double *, double *, double *, double *, double *);

void setcountstopsvar_(void *, int *, int *);

void destroysendbuffer_(void *);

double * returnssendbuffer_(void *,  int *, double *);

double * initsendbuffer_(void *, int *);

void setexpexecutionstuckcount(void *, int *);

int getexpexecutionstuckcount(void *);

void initcooperativestoppingcriteria( void *) ;

double * seriallizevector_ (double *, double *, double *, double *, double *, int *, int *, int *, double *);

int sizeseriallize_(int *,  int *, double *);

void returnmaxelementint_(void *, int *, int *);

void returnminelementint_(void *, int *, int *);

void returnminelement_(void *exp, double *, double *);

void returnminlocelement_(void *, double *, double *, int *, int *);