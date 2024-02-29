/**
 * @file setup_benchmarks.c
 * @author David R. Penas
 * @brief File containing functions about the load for different benchmarks.
 * When users typical want to put their benchmarks, they should modify this
 * file.
 */

#include <stdio.h>
#include <configuration.h>
#include <string.h>
#include <stdlib.h>
#include <structure_paralleltestbed.h>
#include <input_module.h>
#include <output.h>
#include <benchmark_functions_SystemBiology.h>
#include <benchmark_functions_BBOB.h>
#include <python_interface.h>
#include <example.h>
#ifdef MATLAB
#include <examplematlab.h>
#endif	
#include <def_errors.h>

#ifdef MPI2
        #include <mpi.h> 
#endif

#ifdef OPENMP
        #include <omp.h>
#endif

#define MINNUMSB 0
#define MAXNUMSB 11

typedef void*(*function)(double*,void*);


/**
 * @brief this the setup function where the usar.
 * @param exp experiment_total struct, the main struct of program.
 * @param idp parallel identification number.
 * @param first this variable is 1, when it is the first time that the program
 * enters in this function.
 * @return void function pointer of evaluation function of the selected.
 * benchmark.
 */
function setup_benchmark(experiment_total *exp, int idp, int *first) {
    const char *namealg, *custom, *noiseBBOB, *noiselessBBOB, *systemsBiology, *matlabproblem, *pythonNLP;
    char *name;
    function func;
     
    systemsBiology="systemBiology";
    custom="customized";
    noiselessBBOB="noiselessBBOB";
    noiseBBOB="noiseBBOB";
    matlabproblem="matlabproblem";
    pythonNLP = "pythonNLP";

    name = (char *) calloc(500,sizeof(char));
    namealg = getname(exp);
    
    /// CUSTOMIZED BENCHMARK
    if (strcmp((*exp).test.bench.type, custom) == 0) {
        logandtranslation_(exp);
        func = examplefunction;
        name = "user benchmarks";
        (*exp).test.VTR_default=(*exp).test.VTR;
        (*exp).test.bench.openmp=1;
    }    
    /// SYSTEM BIOLOGY PROBLEMS
    else if (strcmp((*exp).test.bench.type, systemsBiology) == 0) {
        /// select a specific system biology problem:
        /// current_bench = 0 --> CIRCADIAN PROBLEM
        /// current_bench = 1 --> step PATHWAY PROBLEM
        /// current_bench = 2 --> NFKB PROBLEM
        /// current_bench = 3 --> B1 BIOPREDYN PROBLEM
        /// current_bench = 4 --> B2 BIOPREDYN PROBLEM
        /// current_bench = 5 --> B3 BIOPREDYN PROBLEM
        /// current_bench = 6 --> B4 BIOPREDYN PROBLEM
        /// current_bench = 7 --> B5 BIOPREDYN PROBLEM
        /// current_bench = 8 --> B6 BIOPREDYN PROBLEM
        manage_options(1,exp->test.bench.current_bench, (*exp).methodScatterSearch->loptions->solver);
        load_benchmark_SystemBiology(exp);
	if (exp->test.bench.current_bench == 8) (*exp).test.bench.openmp=0; // B6 is not compatible with openmp
	else  (*exp).test.bench.openmp=1;
        logandtranslation_(exp);
        func = evalSB_;
        name = (char *) getnameSB(exp->test.bench.current_bench );
    }
    /// BBOB NOISELESS BENCHMARKS
    else if (strcmp((*exp).test.bench.type, noiselessBBOB) == 0) {
        manage_options(2,exp->test.bench.current_bench, (*exp).methodScatterSearch->loptions->solver);
        initializeBBOB(exp, namealg);
        updateFunctionsBBOB(exp, exp->test.bench.current_bench, 1);
        logandtranslation_(exp);
        func = fgeneric_noiseless; 
        (*exp).test.bench.openmp=0;
 
        sprintf(name,"bbob noiseless benchmark ID=%d\n", exp->test.bench.current_bench);
        (*exp).test.VTR_default = fgeneric_ftarget_tol(1e-8);
    }    
    /// BBOB NOISE BENCHMARKS    
    else if (strcmp((*exp).test.bench.type, noiseBBOB) == 0) {
        manage_options(3,exp->test.bench.current_bench, (*exp).methodScatterSearch->loptions->solver);
        initializeBBOB(exp, namealg);
        updateFunctionsBBOB(exp, exp->test.bench.current_bench, 1);
        logandtranslation_(exp);
        func = fgeneric_noise; 
        (*exp).test.bench.openmp=0;
        
        sprintf(name,"bbob noise benchmark ID=%d\n", exp->test.bench.current_bench);
        (*exp).test.VTR_default = fgeneric_ftarget_tol(1e-8);
    }
#ifdef MATLAB
    // MATLAB PROBLEM
    else if (strcmp((*exp).test.bench.type, matlabproblem) == 0) {
#ifdef GNU
        logandtranslation_(exp);
        func = examplefunctionmatlab;
        name = "matlab problem";
        (*exp).test.VTR_default=(*exp).test.VTR;
        (*exp).test.bench.openmp=1;
#ifdef GNU
        __matlabproblem_MOD_openmatlab(&exp->execution.ep_matlab);
#elif defined(INTEL)
        matlabproblem_mp_openmatlab_(&exp->execution.ep_matlab);
#endif
    	
#else
	perror(error43);	
#endif
    }
#endif
    // Python
    else if (strcmp((*exp).test.bench.type, pythonNLP) == 0) {

       logandtranslation_(exp);
       name="python cost function NLP problem";
       (*exp).test.VTR_default=(*exp).test.VTR;
       (*exp).test.bench.openmp=1;
       manage_options(4, exp->test.bench.current_bench, (*exp).methodScatterSearch->loptions->solver);

       const char *name_file;
       //sprintf(name_file, "problem%d", exp->test.bench.current_bench);
       if ( exp->test.bench.current_bench == 1)       name_file = "P1_ForwardSimulation_VEPmodel_Sourcelevel";
       if ( exp->test.bench.current_bench == 2)       name_file = "P2_ForwardSimulation_VEPmodel_Sensorlevel";
       if ( exp->test.bench.current_bench == 3)       name_file = "P3_ForwardSimulation_VEPmodel_Sourcelevel_nopropagation";
       if ( exp->test.bench.current_bench == 4)       name_file = "P4_ForwardSimulation_VEPmodel_Sensorlevel_nopropagation";
       if ( exp->test.bench.current_bench == 5)       name_file = "P4_ForwardSimulation_VEPmodel_stochastic_estimateTau_Sensorlevel";
       if ( exp->test.bench.current_bench == 6)       name_file = "P5_ForwardSimulation_VEPmodel_Sourcelevel_SDE_stiff_withpropagation";
       if ( exp->test.bench.current_bench == 7)       name_file = "P6_ForwardSimulation_VEPmodel_Sensorlevel_SDE_stiff_withpropagation";
       if ( exp->test.bench.current_bench == 8)       name_file = "P6_ForwardSimulation_VEPmodel_Sensorlevel_SDE_stiff_withpropagation_100sim";
       if ( exp->test.bench.current_bench == 9)       name_file = "Problem1_ForwardSimulation_VEPmodel_Sourcelevel_nopropagation_NN42";
       if ( exp->test.bench.current_bench == 10)      name_file = "Problem1_ForwardSimulation_VEPmodel_Sourcelevel_nopropagation_NN84";
       if ( exp->test.bench.current_bench == 11)      name_file = "Problem1_ForwardSimulation_VEPmodel_Sourcelevel_nopropagation_NN162";
       if ( exp->test.bench.current_bench == 12)      name_file = "Problem1_ForwardSimulation_VEPmodel_Sourcelevel_nopropagation_NN400";
       const char *string2 = "cost_function";

       int NPROC_OPENMP=1;
       int i;

//#ifdef OPENMP
//#pragma omp parallel
//{
//       NPROC_OPENMP = omp_get_num_threads();
//}
//        (*exp).test.bench.openmp=1;
//	for (i = 0; i < NPROC_OPENMP; i++) {
//		init_python(name_file, string2, &exp[i]);
//	}
//	func = call_to_obj_function_python;
//# else

       (*exp).test.bench.openmp=0;
       init_python(name_file, string2, exp);
       func = call_to_obj_function_python;
//#endif
       exp->execution.python_active = 1;

    }
    else {
        perror(error10);
        exit(10);
    }

    
    if ((idp == 0)&&(*first==1)) {
        bechmark_message(exp,name );
        *first=0;
    }
    /// to order the printf messages
    sleep(1);
    
    return func;
}



/**
 * @brief this is a control function: check if exist the id range for type of 
 * benchmark selected.
 * @param id type of the benchmark.
 * @param i id of the benchmark selected.
 * @param solver name of the selected solver.
 * @return void function pointer of evaluation function of the selected 
 * benchmark.
 */
void manage_options(int id, int i, const char* solver){
    
    if (id==1){
        if ((i<MINNUMSB)||(i>MAXNUMSB)) {
            perror(error14);
            exit(14);
        }
        
        
    }
    else if (id==2) {
        if ((i<1)||(i>24)) {
            perror(error15);
            exit(15);            
        }
    }
    else if (id==3) {
        if ((i<101)||(i>130)) {
            perror(error16);
            exit(16);            
        }        
    }
    
    if ((id==2)||(id==3)){
        if (strcmp(solver,"n2sol")==0) {
            perror(error11);
            exit(11);
        }
    }
	
//    if (id==4){
//        if ((i<1)||(i>7)) {
//            perror(error14);
//            exit(14);
//        }
//    }
}



