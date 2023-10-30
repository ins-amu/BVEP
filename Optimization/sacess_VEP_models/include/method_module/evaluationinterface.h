void deallocateoutputfunction_(output_function *, void *);

void DE_correction_bounds(double *, int, double* , double* );

void DE_correction_bounds2(double *X, int D, double* Xl, double* Xu);

double callfitnessfunction_(void *(*fitnessfunction)(double*, void *), void *, double *, int *, double *, double *);

double callfitnessfunctionfortran_(void *(*fitnessfunction)(double*, void *,double*), void *, double *, int *, double *, double *, double *);

double callfitnessfunctionopenmp_(void *(*fitnessfunction)(double*, void *), void *, double *, int *, double *, double *, int *);

double callfitnessfunctionfortranopenmp2_(void *(*fitnessfunction)(double*, void *), void *, double *,double*, int *);
