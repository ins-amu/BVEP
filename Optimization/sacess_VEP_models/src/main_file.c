/**
 * @file main_file.c
 * @author David R. Penas
 * @brief This is the main file of the program.
 */

#include <stdio.h>
#include <configuration.h>
#include <string.h>
#include <stdlib.h>
#include <structure_paralleltestbed.h>
#include <input_module.h>
#include <output.h>
#include <example.h>
#include <setup_benchmarks.h>
#include <def_errors.h>

#ifdef MPI2
        #include <mpi.h> 
#endif

#ifdef OPENMP
        #include <omp.h>
#endif




int main(int argc, char** argv) {
    int id, NPROC, error, i, NPROC_OPENMP;
    const char *file, *path;
    experiment_total *exptotal;
    result_solver result;
    function func;
    int first, init;
    argc_struct arg; 
    
    #ifdef MPI2
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
        MPI_Comm_rank(MPI_COMM_WORLD, &id);          
    #else        
        id = 0;
        NPROC = 1;    
    #endif
        
    
// INPUT :
// path --> foreach RUN the path must be different
    if (argc < 3) {
        if (id == 0)
            printf(" ARG1 file_config PATH\n\n");
        return 0;
    }
        
    file = argv[1];
    if (argv[2] != NULL) path = argv[2];

    init_argc_struct(&arg);

    if (argv[3] != NULL) {
   	 for (int i = 3; i < argc; i++) {
         	if (strcmp(argv[i], "-nstuck") == 0) {
			if (i + 1 < argc) { 
				arg.n_stuck = atoi(argv[i + 1]); 
				i++;
			} else {
				fprintf(stderr, "Option -nstuck requires an argument.\n");
				return 1;
			}
        	 }
    	}
    }
   

// CREATE EXPTOTAL STRUCT
#ifdef OPENMP
#pragma omp parallel
    {
            NPROC_OPENMP = omp_get_num_threads();
    }

    exptotal = (experiment_total *) malloc( NPROC_OPENMP * sizeof(experiment_total));
    init=1;
    for (i = 0; i < NPROC_OPENMP; i++) {
        // PARSE THE OPTIONS OF THE SOLVER
    	create_expetiment_struct(file, &(exptotal[i]), NPROC, id, path, init, arg);
	init=0;
    }
    // INIT MESSAGE    
    if ( (id == 0) ) init_message(NPROC, &(exptotal[0]), NPROC_OPENMP );  
    // INIT BENCHMARK
    first=1;
    for (i = 0; i < NPROC_OPENMP; i++) {
        func = setup_benchmark(&(exptotal[i]), id, &first);
    }   

    init_result_data(&result, exptotal[0].test.bench.dim);

#else
    exptotal = (experiment_total *) malloc(sizeof(experiment_total));
// PARSE THE OPTIONS OF THE SOLVER    
    init=1;
    create_expetiment_struct(file, exptotal, NPROC, id, path, init, arg);
// INIT MESSAGE 
    NPROC_OPENMP=1;   
    if (id == 0) init_message(NPROC, &exptotal[0], NPROC_OPENMP ) ;
// INIT BENCHMARK
    first=1;
    func = setup_benchmark(exptotal,id,&first);  
    
        
// INIT RESULT STRUCT    
    init_result_data(&result, exptotal->test.bench.dim);
    
#endif    
    
#ifdef MPI2
    MPI_Barrier(MPI_COMM_WORLD);    
#endif
// RUNNING THE SOLVER
    if ( execute_Solver(exptotal, &result, func) == 0) {
	perror(error8);
    }
    
// PRINT RESULTS   
    if (id == 0) {
        printsolution_(exptotal, &result.bestx_value[0], &result.best_value);
        printresults_end(exptotal,result);
    }
// PLOT RESULTS
    plot(&exptotal[0] ); 
    if (id == 0) {
        graphs_message(exptotal);
    }
    
// DESTROY EXP DATA       
    destroyexp(exptotal);

    #ifdef MPI2
        MPI_Finalize();
    #endif
 
    return (EXIT_SUCCESS);
}




