typedef void*(*function)(double*,void*);

extern void* examplefunction(double*,void*);
extern void* fgeneric_noiseless(double*,void*);
extern void* fgeneric_noise(double*,void*);
extern void* evalSB_(double*,void*);


function setup_benchmark(experiment_total *, int, int *);
