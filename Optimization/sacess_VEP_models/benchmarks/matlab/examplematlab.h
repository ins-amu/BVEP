#ifdef MATLAB
#ifdef GNU
#include <structure_paralleltestbed.h>

typedef struct {
  double *residual;
  int size_residual;
  double *penalty;
  int size_penalty;
} outputmatlab;

void allocatevectorrg_(void *, int *, double *, int *, double *);
// define the function
void* examplefunctionmatlab(double *, void *);
#endif
#endif 
