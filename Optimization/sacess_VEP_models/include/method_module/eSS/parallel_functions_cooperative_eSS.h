#ifdef MPI2

void asynchronousstoppingcriteriaessmaster_(void *, int *, long *, double *, double *, double *, int *);
void asynchronousstoppingcriteriaessslave_(void *, int *, long *, double *, double *, double *, int *);

void incrementmaxtimemig_(void *);

void chargecooperativeparametersfortran_(void *, int *, int *, int *, double *, long * ) ;

void createcooperativetopologyess_(void *);

void createcooperativestructsess_(void *, int *);

int destroycooperativestructuresess_(void *);

void initcooperativestoppingcriteriaess_( void *);

double * initsendbufferess_(void *, int *);

double * returnssendbufferess_(void *,  int *, double *);

void destroysendbufferess_(void *);

int checkcooperativemigrationcriteriaess_(void *);

int checkcooperativemigrationcriteriacessinner_(void *);

int checkcooperativemigrationcriteriacess_(void *, double *); 

int checkcooperativemigrationcriteriacesstao_(void *, long *, int *); 

void migrationcooperativesentess_(void *, double *, int *, int *, double *, double *, int *);

void migrationsynchcooperativesentess_(void *, int *, double *, double *, double *);

void migrationcooperativereceptioness_( void *, void*(*fitnessfunction)(double*, void *),double *, int *, int *, double *, double *, int *);

void asynchronousstoppingcriteriawithrestartess_(void *, int *, long *, double *,double *,int *);

void asynchronousstoppingcriteriaess_(void *, int *, long *, double *, double *, double *, int *) ;

void synchronousstoppingcriteriacess_(void *, int *, long *, double *, double *, double *, int *);

void synchronousstoppingcriteriacessmaster_(void *, int *, long *, double *, double *, double *, int *);

void setcountstopsvaressmaster_(void *, int *, int *);

void migrationsynchcooperativesentcess_(void *,int *, double *, double *, int *, int *);

void masterupdatesolutionslaves_( void *, double *, int *, int *, int *, double *,int *, double *, int *, int *  );

void returnwindowvalueslave_(void *, double *, int *, int *);

void setcountstopsvaress_(void *, int *, int *) ;

int cooperativempitestess_(void *, int *) ;

void gatherresultsess_(void *, double *, double *, int *, double *, double *, long *, long *) ;


void gatherresultsserializeess_(void *, void *, double *, double *, int *, double *, double *,  long *, long *);


void setexpexecutionstuckcountess(void *, int *) ;


int getexpexecutionstuckcountess(void *);


#endif
