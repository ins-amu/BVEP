/**
 * @file output.c
 * @author David R. Penas
 * @brief File containing all functions about the different outputs of the 
 * program: print screens, files logs, matlab files generation etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <structure_paralleltestbed.h>
#include <configuration.h>
#include <read_a_file_line_by_line.h>
#include <common_solver_operations.h>
#include <hdf5.h>
#include <benchmark_functions_SystemBiology.h>
#ifdef MPI2
        #include <mpi.h>
#endif

void printtranslationmassage() {
    printf("* TRANSLATION ENABLED: logarithmic space is enabled and there are negative values in the bounds\n");
}

void printganthyper_(void *exp1_, double *currenttime, int *code) {


    experiment_total *exp1;
    FILE *p3;


    exp1 = (experiment_total *) exp1_;

//    if (exp1[0].test.output == 1) {
//        p3 = fopen((char *) exp1[0].test.hyper_log, "a");
//        if (p3 == NULL) exit(0);

 //       fprintf(p3, "%.20lf,%ld\n", *currenttime, *code);

  //      fclose(p3);
  //      p3 = NULL;
 //   }

}


void printsolution_(void *exp1_, double *xbest, double *fbest) {


    experiment_total *exp1;
    output_struct *output;
    FILE *p3;
    int i;

    exp1 = (experiment_total *) exp1_;
    output = exp1->output;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
            if (p3 == NULL) exit(0);

            fprintf(p3, "SOLUTION[%d] - fx = %0.15f\n", exp1[0].execution.idp, *fbest );
            fprintf(p3, "XSOLUTION ");
            for (i=0; i< exp1[0].test.bench.dim; ++i) {
                   if ( i != ( exp1[0].test.bench.dim - 1) ) {
                           fprintf(p3, "%0.15f, ", exp1[0].execution.idp, xbest[i]);
                   }
                   else {
                           fprintf(p3, "%0.15f", exp1[0].execution.idp, xbest[i]);
                   }
            }
            fprintf(p3, "\n");

            fclose(p3);
            p3 = NULL;
   }

}



void updateresultsandprint2_(void *exp1_, void *result_, void *output_, void *local_s_, double *totaltime, long *evaluation,
        int *D, double *xbest, double *best, int *idp, int *NPROC ) {
    int i, indexvector;
    double auxbest;
    result_solver *result;
    output_struct *output;
    local_solver *local_s;
    experiment_total *exp1;
    double st3;
    double *vectorbestx,*vectorbest,*bestglobalx,bestglobal;

    st3 = 0.0;
    exp1 = (experiment_total *) exp1_;
    result = (result_solver *) result_;
    output = (output_struct *) output_;
    local_s = (local_solver *) local_s_;
    double totalL, sucessL;



#ifdef MPI2
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&(local_s->total_local), &totalL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(local_s->sucess_local), &sucessL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(output->st3), &st3, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);

    if (*idp == 0) {
            local_s->total_local = totalL / (double) *NPROC;
           local_s->sucess_local = sucessL/ (double) *NPROC;
            output->st3 = st3 / (double) *NPROC;
    }

    MPI_Bcast(&(local_s->total_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(local_s->sucess_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(output->st3), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);

    vectorbestx = (double *) malloc(*NPROC*(*D)*sizeof(double));
    vectorbest = (double *) malloc(*NPROC*sizeof(double));
    bestglobalx = (double *) malloc((*D)*sizeof(double));
    MPI_Gather(xbest,*D,MPI_DOUBLE,vectorbestx,*D,MPI_DOUBLE,0,exp1->execution.topology.comunicator);
    MPI_Gather(best,1,MPI_DOUBLE,vectorbest,1,MPI_DOUBLE,0,exp1->execution.topology.comunicator);

    if (*idp == 0) {
            auxbest = vectorbest[0];
            indexvector = 0;
            for (i=0;i<*NPROC;i++) {
                if ( vectorbest[i] < auxbest ) {
                    auxbest = vectorbest[i];
                    indexvector = i;
                }
            }

            bestglobal=auxbest;
            memmove(bestglobalx,&vectorbestx[indexvector*(*D)],*D*sizeof(double));
    }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(bestglobalx, *D, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
        MPI_Bcast(&bestglobal, 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
        if ((*result).best_value >= bestglobal) {
            (*result).best_value = bestglobal;
        }
        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;

        print_end_file_(exp1,bestglobalx,&bestglobal,result,D,output,local_s,idp);

        free(vectorbestx);
        vectorbestx=NULL;
        free(vectorbest  );
        vectorbest=NULL;
        free(bestglobalx  );
        bestglobalx=NULL;

        MPI_Barrier(MPI_COMM_WORLD);
#else

        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;

        print_end_file_(exp1,xbest,best,result,D,output,local_s,idp);


        if ((*result).best_value >= *best) {
            (*result).best_value = *best;
        }
#endif
}


/**
 * @brief this procedure initialize the differents output files of the program.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param best fx value of the best solution of the processor up to now.
 * @param par binary parametter: 1 if the program is parallel, 0 in the 
 * opposite case.
 * @param idp identification of the parallel processor.
 * @param currenttime time spent from the beginning.
 * @param iterations current iterations of the solver.
 * @param master binary parametter: 1 if the processors is the master, 0 in the 
 * opposite case.
 */
void initprintfile_(void *exp1_, double *best, int *par, int *idp, double *currenttime, long *evals, int *master) {
    experiment_total *exp1;
    output_struct *output;
    char *cadea2;
    int initp;
    FILE *p3, *p4, *p5, *p6;
    
    
    exp1 = (experiment_total *) exp1_;
    output = exp1->output;
    
    cadea2 = (char *) malloc(400* sizeof (char));

    (*output).st3 = 0.0;

    output->point1 = *currenttime;
    initp=0;
    p5 = NULL;
    if (exp1[0].test.output == 1) {
        output->oldbest = *best;
        if (*par == 0) {
            exp1[0].test.log_output = (char *) calloc(500, sizeof(char));
            exp1[0].test.result_output = (char *) calloc(500, sizeof(char));
	    char *temp = (char *) calloc(1000, sizeof(char));
	    strcpy(temp, exp1[0].test.output_graph);
            sprintf(exp1[0].test.output_graph, "%s/convergence_id%d.csv", temp,*idp);
            sprintf(exp1[0].test.log_output,     "%s/logfile_id%d", temp,*idp);
            p3 = fopen( (const char *) exp1[0].test.output_graph, "a");
	    free(temp);
        } else {

            if (exp1[0].execution.initpath == 1) {
                initp=1;
                exp1[0].execution.initpath = 0;
                exp1[0].test.log_output = (char *) calloc(500, sizeof(char));
                exp1[0].test.log_percentage = (char *) calloc(100, sizeof(char));
                exp1[0].test.output_gant_log = (char *) calloc(100, sizeof(char));
                exp1[0].test.result_output = (char *) calloc(500, sizeof(char));     
		char *temp = (char *) calloc(1000, sizeof(char));
                strcpy(temp, exp1[0].test.output_graph);
                sprintf(exp1[0].test.result_output,     "%s/result", temp);
                sprintf(exp1[0].test.log_output,     "%s/logfile_id%d", temp,*idp);
                if (*idp == 0)  sprintf(exp1[0].test.log_percentage, "%s/percentage_id%d.csv", temp,*idp);
                sprintf(exp1[0].test.output_gant_log,"%s/gantt_id%d.csv", temp,*idp);
                if (*master == 1) {
                    if (*idp != 0) 
                        sprintf(exp1[0].test.output_graph,   "%s/convergence_id%d.csv", temp, *idp);  
                } else {
                    sprintf(exp1[0].test.output_graph,   "%s/convergence_id%d.csv", temp, *idp); 
                }
            }

            
            p3 = fopen( (const char *) exp1[0].test.output_graph, "a");
            p4 = fopen( (const char *) exp1[0].test.log_output, "a");
            p5 = fopen( (const char *) exp1[0].test.log_percentage, "a");
            p6 = fopen( (const char *) exp1[0].test.output_gant_log, "a");

        }

        if (*master == 1) {
            if (*idp != 0) {
                if ((p3 != NULL)) {
                    fprintf(p3, "Time,BestFx,N_evals,N_Call_Local_solver,N_Migration\n");
                    fprintf(p3, "%.20lf,%.20lf,%ld,0,0\n", *currenttime, *best, *evals);
                    fclose(p3);
                    p3 = NULL;
                }
            } 
        } else {
                if ((p3 != NULL)) {
                    fprintf(p3, "Time,BestFx,N_evals,N_Call_Local_solver,N_Migration\n");
                    fprintf(p3, "%.20lf,%.20lf,%ld,0,0\n", *currenttime, *best, *evals);
                    fclose(p3);
                    p3 = NULL;
                }            
        }
       


#ifdef MPI2        
        if (*idp == 0) {
            if (p5 != NULL) {
                fprintf(p5, "Time,BestFx,Percentage\n");
                fclose(p5);
                p5 = NULL;
            }       
        }

        if (p6 != NULL) {
            if (p6 != NULL) {
                fprintf(p6, "Time,State\n");		
        	fclose(p6);
        	p6 = NULL;
	    }
        }  
        if (p4 != NULL) {
                fclose(p4);
        }

#endif        
    }

    
}


/**
 * @brief print current results in the logfile
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param best fx value of the best solution of the processor up to now.
 * @param par binary parametter: 1 if the program is parallel, 0 in the 
 * opposite case.
 * @param evaluation_local  function evaluations spent by the processor up to 
 * now.
 * @param currenttime time spent from the beginning.
 */
void printiterationcesslog_(void *exp1_,  double *best, long *evaluation_local, double *currenttime,
        int *ite) {


    experiment_total *exp1;
    output_struct *output;
    FILE *p3;


    exp1 = (experiment_total *) exp1_;
    output = exp1->output;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);

        fprintf(p3, "ID[%d] - Iteration: %d NFunEvals: %ld Bestf: %lf CURRENT TIME: %lf s\n",exp1[0].execution.idp, *ite, *evaluation_local,*best, *currenttime );

        fclose(p3);
        p3 = NULL;
    }

}


void printdescartsolution_(void *exp1_,  double *fpen, int  *idp) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION - discart solution  %lf from %d processor in reception subroutine\n", exp1[0].execution.idp, *fpen, *idp );
        fclose(p3);
        p3 = NULL;
    }
}



void printinititeration_(void *exp1_,  int *iter, double *cputime3) {
    experiment_total *exp1;
    FILE *p3;


    exp1 = (experiment_total *) exp1_;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\n+++ID[%d] --- INIT ITERATION %d - CURRENT TIME: %lf s\n",exp1[0].execution.idp, *iter+1,  *cputime3 );

        fclose(p3);
        p3 = NULL;
    }

}

void printlsinitlog_(void *exp1_, double *best ) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);

        fprintf(p3, "\tID[%d] - LOCAL SOLVER - Call local solver: Initial point function value: [%lf] \n", exp1[0].execution.idp, *best );

        fclose(p3);
        p3 = NULL;
    }

}

void printlocalsolverinsert_(void *exp1_, double *valor, double *localmatrix, int *sizelocal, int *exito, double *comp1, double *comp2 ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        if (*exito == 0) {
	  fprintf(p3, "\tID[%d] - LOCAL SOLVER - Point %lf is NOT put in refset. [NEW %lf compare with OLD %lf]. It is near to the set LOCAL_SOLVER_TABU_LIST : [ ", exp1[0].execution.idp, *valor, *comp1, *comp2 );
//          for (i=0;i<*sizelocal;i++){
//		fprintf(p3," %lf ", localmatrix[i] );
//	  }      
	  fprintf(p3,"***]\n");
	} else {
          fprintf(p3, "\tID[%d] - LOCAL SOLVER - Point %lf is put in refset. [NEW %lf compare with OLD %lf]. LOCAL_SOLVER_TABU_LIST :  [ ", exp1[0].execution.idp, *valor, *comp1, *comp2);
//          for (i=0;i<*sizelocal;i++){
//                fprintf(p3," %lf ", localmatrix[i] );
//          }
          fprintf(p3,"***]\n");	
        }

        fclose(p3);
        p3 = NULL;
    }

}

void printinitlocalsolverinsert_(void *exp1_, double *valor, double *localmatrix, int *sizelocal ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - LOCAL SOLVER - Init point X0 %lf is put in INITIAL_POINTS_TABU_LIST : [ ", exp1[0].execution.idp, *valor );
//        for (i=0;i<*sizelocal;i++){
//                fprintf(p3," %lf ", localmatrix[i] );
//        }
        fprintf(p3,"***]\n");
         

        fclose(p3);
        p3 = NULL;
    }

}

void printlsendlog_(void *exp1_, double *best, double *currenttime, long *evals ) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);

        fprintf(p3, "\tID[%d] - LOCAL SOLVER - OUTPUT: Local solution function value: [%lf] - %lf seconds spent in local solver - %ld evals \n",exp1[0].execution.idp,*best, *currenttime, *evals );

        fclose(p3);
        p3 = NULL;
    }

}

void printadaptation_(void *exp1_, double *balance, int *dim, int *counter, double *time) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");

        fprintf(p3, "\tID[%d] - ADAPTING THREAD -- new balance %lf, new size of refset %d, new ncounter %d - TIME[%lf]\n",
                exp1[0].execution.idp,  *balance, *dim, *counter, *time );

        fclose(p3);
        p3 = NULL;
    }

}

void printadaptationslave_(void *exp1_, int *recp, int *send, long *diff_evals, double *time) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");

        fprintf(p3, "\tID[%d] - SEND ADAPTING SIGNAL TO MASTER -- recv sol %d send sol %d diff-evals %ld - TIME[%lf]\n",
                exp1[0].execution.idp,  *recp, *send, *diff_evals, *time );

        fclose(p3);
        p3 = NULL;
    }

}

void printadaptationmaster_(void *exp1_, int *origin, int *dest, double *balance, int *dim, int *counter, int *accept, int *adap, double *time) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (*accept == 1) {
	    if (*adap == 1) {
            fprintf(p3, "\tID[%d] - SEND ADAPTING THREAD %d [RECV>SEND] -- new balance %lf, new size of refset %d, new ncounter %d [best id %d] [ACCEPT] - TIME[%lf]\n",
                exp1[0].execution.idp,  *origin, *balance, *dim, *counter, *dest, *time );
	    }
            if (*adap == 2) {
            fprintf(p3, "\tID[%d] - SEND ADAPTING THREAD %d [MAXEVALSx2] -- new balance %lf, new size of refset %d, new ncounter %d [best id %d] [ACCEPT] - TIME[%lf]\n",
                exp1[0].execution.idp,  *origin, *balance, *dim, *counter, *dest, *time );
            }
            if (*adap == 3) {
            fprintf(p3, "\tID[%d] - SEND ADAPTING THREAD %d [MAXEVALS&SEND=0] -- new balance %lf, new size of refset %d, new ncounter %d [best id %d] [ACCEPT] - TIME[%lf]\n",
                exp1[0].execution.idp,  *origin, *balance, *dim, *counter, *dest, *time );
            }

        } else {
            fprintf(p3, "\tID[%d] - SEND ADAPTING THREAD %d  -- new balance %lf, new size of refset %d, new ncounter %d [best id %d] [REJECT] - TIME[%lf]\n",
                exp1[0].execution.idp,  *origin, *balance, *dim, *counter, *dest, *time );            
        }
        fclose(p3);
        p3 = NULL;
    }

}

void printadaptationscores_(void *exp1_, double *scores, int *NPROC) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");
            
        fprintf(p3, "\tID[%d] - ADAPTING SCORES:  ", exp1[0].execution.idp );
	for (i=0;i<*NPROC;i++) {
		fprintf(p3, " proc(%d)= %lf ",i+1,scores[i]);	
	}
	fprintf(p3, "\n");

        fclose(p3);
        p3 = NULL;
    }

}

void printdesadaptationmaster_(void *exp1_, int *origin, int *dest, double *balance, int *dim, int *counter, double *time) {


    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_output, "a");

        fprintf(p3, "\tID[%d] - SEND (DES)ADAPTING THREAD %d  -- new balance %lf, new size of refset %d, new ncounter %d [best id %d] - TIME[%lf]\n",
                exp1[0].execution.idp,  *origin, *balance, *dim, *counter, *dest, *time );

        fclose(p3);
        p3 = NULL;
    }

}

#ifdef MPI2
/**
 * @brief this function prints in the logfile a summary about the configuration
 *  and the contribution of the different slaves.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param vector this vector contains the contributions of each processor.
 * @param size this variable represents the number of the slaves of the system.
 * @param masterslave binary parametter: 1 if the program has a master-slave
 * scheme, 0 in the opposite case.
 * @param dim_refset dimension of the Reference Set
 * @param n2 iteration to enter in the local solver.
 * @param balance balance variable used to select the solution to enter in the 
 * local solver.
 */
void printmasterocurrencesend_(void *exp1_, int *vector, int  *size, int *masterslave, int *dim_refset, int *n2, double *balance ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i, cont;
    int *n2_vector;
    int *dim_refset_vector;
    double *balance_vector;
    int sizelocal;

	
    n2_vector = (int *) malloc(*size*sizeof(int));
    dim_refset_vector = (int *) malloc(*size*sizeof(int));
    balance_vector = (double *) malloc(*size*sizeof(double));

    sizelocal=1;
    cooperativegathertelementint_(exp1, dim_refset, &sizelocal, dim_refset_vector);
    cooperativegathertelementint_(exp1, n2, &sizelocal, n2_vector);
    cooperativegathertelement_(exp1, balance, &sizelocal, balance_vector);

    if ( exp1[0].execution.idp == 0 ) { 

    	if (exp1[0].test.output == 1) {
        	p3 = fopen((char *) exp1[0].test.log_output, "a");
        	if (p3 == NULL) exit(0);
        	fprintf(p3, "\nID[%d] - MASTER :  contributions of solutions per procesor \n", exp1[0].execution.idp  );
		for ( i=0; i<*size-1; i++ ) {
			if (*masterslave == 1 ) cont = i + 1;
			else cont = i;
			fprintf(p3, "\tPROCESOR %d ---> %d submitted solutions\n",cont, vector[cont] );
		}

        	fprintf(p3, "\n\n" );

		for ( i=0; i<*size-1; i++ ) {
        		if (*masterslave == 1 ) cont=i+1;
			else cont=0;
			fprintf(p3, "\tCONFIGURATION PROCESOR %d", cont);
			fprintf(p3, " --- DIM_REFSET %d   ITERATIONS BEFORE TO ENTER IN LOCAL SOLVER %d   BALANCE %lf\n", dim_refset_vector[cont], n2_vector[cont], balance_vector[cont]  );
		}

	        fclose(p3);
	        p3 = NULL;
    	}
    }
	
    free(n2_vector);
    free(dim_refset_vector);
    free(balance_vector);

}
#endif

void printputmasterlog_(void *exp1_, double *fpen, double *time, int *id ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - MASTER - the new data %lf (from processor %d) is put in the slaves  --- CURRENT TIME %lf \n", exp1[0].execution.idp, *fpen, *id, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printputcheckmasterlog_(void *exp1_, double *new, double *old ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "ID[%d] - MASTER -check if new value %lf is better than %lf \n", exp1[0].execution.idp, *new, *old  );
        fclose(p3);
        p3 = NULL;
    }

}

void printreturnselectsolution_(void *exp1_,  double *fpen, double *oldvalue, double *currenttime, int *idp ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "ID[%d] - RECEPTION - New solution chooses in reception  %lf from %d processor compares with current best %lf --- CURRENT TIME %lf \n", exp1[0].execution.idp, *fpen, *idp, *oldvalue, *currenttime  );
        fclose(p3);
        p3 = NULL;
    }

}

void printcheckvtr_(void *exp1_,  double *best, long *evals,  double *currenttime, double *vtr ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if ((exp1[0].test.VTR >= *best)  &&  ( exp1[0].test.jfprint == 0)) {
        *vtr = *currenttime;
        exp1[0].test.jfprint = 1;
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        printf("[SLAVE ID=%d] fx = %.10lf < VTR =%.10lf -- CURRENT TIME %lf s\n", 
                exp1[0].execution.idp,
                *best,
                exp1[0].test.VTR,
                *currenttime);
        
        fprintf(p3, "\nID[%d] - FJ OVERCOME - value %lf - evals %ld  --- CURRENT TIME %lf s\n", exp1[0].execution.idp, *best, *evals, *currenttime  );
        fclose(p3);
        p3 = NULL;
    }

}

void printrecvmasterlog_(void *exp1_,  double *fpen, double *oldvalue, double *currenttime, int *idp ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - MASTER - New solution received %lf from %d processor compares with current best %lf --- CURRENT TIME %lf s\n", exp1[0].execution.idp, *fpen, *idp, *oldvalue, *currenttime  );
        fclose(p3);
        p3 = NULL;
    }

}

void printputmasterendlog_(void *exp1_,  double *time, double *currenttime, double *valuecandidate, int *idp ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "ID[%d] - MASTER ACCEPT CANDIDATE FROM %d - fx candidate %lf - %lf seconds spent in the last put --- CURRENT TIME %lf s\n", exp1[0].execution.idp, *idp, *valuecandidate, *time, *currenttime  );
        fclose(p3);
        p3 = NULL;
    }

}

void printcooperative_(void *exp1_, int *cooperative, int *size, double *time, int *iter ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - COOPERATIVE FLAG - Content of the cooperative flag iter %d - CURRENT TIME %lf s: \n", exp1[0].execution.idp, *iter, *time  );
        for (i=0;i<*size;i++) {
                fprintf(p3,"\t\t %d ", cooperative[i] );
                if (i%10 == 9) fprintf(p3,"\n");
        }
        if (i%10 != 0)  fprintf(p3,"\n");


        fclose(p3);
        p3 = NULL;
    }
}

void printrefset_(void *exp1_, double *refset, int *size, double *time, int *iter, int *index ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;
	
    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - REFSET FVAL - Content of the Refset in the end of iteration %d - CURRENT TIME %lf s: \n", exp1[0].execution.idp, *iter, *time  );
	for (i=0;i<*size;i++) {
		fprintf(p3,"\t\t %lf ", refset[i] );
		if (i%10 == 9) fprintf(p3,"\n");
	}
	if (i%10 != 0)  fprintf(p3,"\n");
        fprintf(p3, "\tID[%d] - INDEX CHANGE COUNTER - end iteration %d - CURRENT TIME %lf s: \n", exp1[0].execution.idp, *iter,
                *time  );
	for (i=0;i<*size;i++) {
		fprintf(p3,"\t\t %d ", index[i] );
		if (i%10 == 9) fprintf(p3,"\n");
	}
	if (i%10 != 0)  fprintf(p3,"\n");        

	
        fclose(p3);
        p3 = NULL;
    }
}

void printrefsetinner_(void *exp1_, double *refset, int *size, double *balance ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - REFSET FVAL (balance %lf)- Content of the Refset: \n", exp1[0].execution.idp, *balance );
        for (i=0;i<*size;i++) {
                fprintf(p3,"\t\t %lf ", refset[i] );
                if (i%10 == 9) fprintf(p3,"\n");
        }
        if (i%10 != 0)  fprintf(p3,"\n");


        fclose(p3);
        p3 = NULL;
    }
}

void printchilds_(void *exp1_, double *refset, double *childset, int *sizechild, int *childsetindex, int *childsetindex2 ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - CHILDSET FVAL - Content of the childset: \n", exp1[0].execution.idp );
        for (i=0;i<*sizechild;i++) {
                fprintf(p3,"\t\t %lf (%d - %lf <> %d - %lf)\n", childset[i], childsetindex[i]-1, refset[childsetindex[i]-1], childsetindex2[i]-1, refset[childsetindex2[i]-1] );
        }
        if (i%10 != 0)  fprintf(p3,"\n");


        fclose(p3);
        p3 = NULL;
    }
}

void printenditeration_(void *exp1_ ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;
	
    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        
        fprintf(p3, "\n+++ID[%d] - END ITERATION\n\n\n",  exp1[0].execution.idp );
	
        fclose(p3);
        p3 = NULL;
    }
}

void printrestart_(void *exp1_,  double *time ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RESTART SENT SIGNAL - CURRENT TIME %lf seconds\n", exp1[0].execution.idp, *time  );
        fclose(p3);
        p3 = NULL;
    }
}

void printrestartslave_(void *exp1_,  double *time ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RESTART BEGINS - CURRENT TIME %lf seconds\n", exp1[0].execution.idp, *time  );
        fclose(p3);
        p3 = NULL;
    }
}

void printreplaceslavelog_(void *exp1_, double *fpen, double *time, int *position ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION - ACCEPTED -  received data is introduced in refset - %lf - in position %d --- CURRENT TIME %lf seconds\n", exp1[0].execution.idp, *fpen, *position, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printreceivedslave_(void *exp1_, double *fpen, double *time ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION -   Candidate_from_Master %lf is received --- CURRENT TIME %lf sec\n", exp1[0].execution.idp, *fpen,  *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printdiscardreceivedslave_(void *exp1_, double *fpen, double *time, double *fbest ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION NOT ACCEPTED -   Candidate_from_Master %lf is discard : there is a better solution in REFSET %lf --- CURRENT TIME %lf sec\n", exp1[0].execution.idp, *fpen,  *fbest, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printdiscardreceivedslave3_(void *exp1_, double *fpen, double *time, double *fbest ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION NOT ACCEPTED -   Candidate_from_Master %lf is discard : this solution is the same than the last sent %lf --- CURRENT TIME %lf sec\n", exp1[0].execution.idp, *fpen,  *fbest, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printstuck_(void *exp1_, double *solution, int *ocurrence ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;
    int i;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RESTART SOLUTION --> %lf - iterations without changes --> %d \n", 
                exp1[0].execution.idp, *solution, *ocurrence  );

        fclose(p3);
        p3 = NULL;
    }
}

void printdiscardreceivedslave2_(void *exp1_, double *fpen, double *time, double *porc ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - RECEPTION NOT ACCEPTED -   Candidate_from_Master %lf is discard - it does not improve enough (%lf) --- CURRENT TIME %lf sec\n", exp1[0].execution.idp, *fpen, *porc, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printcomparenewsolutionslavelog_(void *exp1_, double *new, double *old, double *dist, double *time, double *th ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - CHECK TO SEND  - compared BKS_S %lf with BKS_M %lf ---> improving %lf %% (threshold %lf) -- CURRENT TIME %lf s\n", 
		exp1[0].execution.idp, *new, *old, *dist, *th, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printreceptioncompare_(void *exp1_, double *new, double *old, double *dist, double *time, double *maxthreshold ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - COMPARATION SOLUTIONS IN RECEPTION - compared candidate value %lf with fbest value %lf ---> improving %lf %% (min %lf%%) -- CURRENT TIME %lf s\n", 
		exp1[0].execution.idp, *new, *old, *dist, *maxthreshold, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printcomparenewsolutionmasterlog_(void *exp1_, double *new, double *old, double *dist, double *time, double *maxthreshold ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\tID[%d] - COMPARATION SOLUTIONS - compared candidate value %lf with best know value %lf ---> improving %lf %% (min %lf%%) -- CURRENT TIME %lf s\n", 
		exp1[0].execution.idp, *new, *old, *dist, *maxthreshold, *time  );
        fclose(p3);
        p3 = NULL;
    }

}

void printfinalsendslavelog_(void *exp1_, double *new, double *time ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    FILE *p3;

    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.log_output, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "\n\tID[%d] - ASYNCH SEND - slave final send ---> %lf -- CURRENT TIME %lf seconds\n", 
		exp1[0].execution.idp, *new, *time  );
        fclose(p3);
        p3 = NULL;
    }

}



/**
 * @brief print iteration in output files. This files are used to generate the
 * final Matlab graphs.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param k iteration spent from the beginning.
 * @param best best solution up to now.
 * @param evaluation_local function evaluations spent by the processor up to 
 * now.
 * @param currenttime  time spent from the beginning.
 */
void printiteration_(void *exp1_,
        int *k, double *best, long *evaluation_local, double *currenttime) {

    experiment_total *exp1;
    output_struct *output;
    FILE *p3;
        
    exp1 = (experiment_total *) exp1_;
    output = exp1->output;
    
    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.output_graph, "a");
        if (p3 == NULL) exit(0);
  
        if (*best < output->oldbest ) {	
            if ( exp1[0].test.output_stop != 1 )
            	fprintf(p3, "%.20lf,%.20lf,%ld,0,0\n", *currenttime, *best, *evaluation_local);
            if ( exp1[0].test.VTR >= *best ) exp1[0].test.output_stop = 1;
            output->oldbest = *best;
        }
        
        fclose(p3);
        p3 = NULL;
    }
    
    

}

/**
 * @brief print results in gant files. This files are used to generate the
 * Gant char.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param currenttime  time spent from the beginning.
 * @param code inner code to represent the kind of operation: (1) global search,
 * (2) local search, (3) communication.
 */
void printgant_(void *exp1_, double *currenttime, int *code) {

   
    experiment_total *exp1;
    FILE *p3;
        
        
    exp1 = (experiment_total *) exp1_;
    
    if (exp1[0].test.output == 1) {
        p3 = fopen((char *) exp1[0].test.output_gant_log, "a");
        if (p3 == NULL) exit(0);
   
        fprintf(p3, "%.20lf,%ld\n", *currenttime, *code);
        
        fclose(p3);
        p3 = NULL;
    }
    
    


}

/**
 * @brief print results in percentage files. This files are used to generate 
 * the percentages char.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param currenttime  time spent from the beginning.
 * @param value value of the solution sent.
 * @param percentage percentage of improvement.
 */
void printpercentage_(void *exp1_, double *currenttime, double *value, double *percentage) {

   
    experiment_total *exp1;
    FILE *p3;
        
        
    exp1 = (experiment_total *) exp1_;
    
    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.log_percentage, "a");
        if (p3 == NULL) exit(0);
        fprintf(p3, "%.20lf,%lf,%lf\n", *currenttime, *value, *percentage);
        
        fclose(p3);
        p3 = NULL;
    }
    
    


}


/**
 * @brief print results in percentage files. This files are used to generate 
 * the percentages char.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param best  best solution up to now.
 * @param evaluation_local function evaluations spent by the processor up to 
 * now.
 * @param currenttime  time spent from the beginning.
 * @param locals number of local solver calls spent from the beginning.
 * @param migration number of migration performed up to now.
 */
void printiterationcess_(void *exp1_, double *best, long *evaluation_local, double *currenttime,
        int *locals, int *migration) {

   
    experiment_total *exp1;
    output_struct *output;
    FILE *p3;
        
        
    exp1 = (experiment_total *) exp1_;
    output = exp1->output;
    
    if (exp1[0].test.output == 1) {

        p3 = fopen((char *) exp1[0].test.output_graph, "a");
        if (p3 == NULL) exit(0);

        if ((*best < output->oldbest ) && ((*locals == 0) && (*migration == 0))) {
       
            if ( exp1[0].test.output_stop != 1 )
            fprintf(p3, "%.20lf,%.20lf,%ld,0,0\n", *currenttime, *best, *evaluation_local);
            if ( exp1[0].test.VTR >= *best ) exp1[0].test.output_stop = 1;
            output->oldbest = *best;
        } 
        else if ((*locals > 0) || (*migration > 0)) {
            if ( exp1[0].test.output_stop != 1 )
            fprintf(p3, "%.20lf,%.20lf,%ld,%d,%d\n", *currenttime, *best, *evaluation_local, *locals, *migration);
            if ( exp1[0].test.VTR >= *best ) exp1[0].test.output_stop = 1;
        }
        
        fclose(p3);
        p3 = NULL;
    }


}

void printverboselocaloutput_(void *exp1_, int *D, double *U, double *fval,int *idp) {
    experiment_total *exp1;
    output_struct *output;
    
    exp1 = (experiment_total *) exp1_;
    output = exp1->output;
    
    if (exp1->test.verbose != 0) {
#ifdef MPI2
            (*output).st1 = MPI_Wtime();
#else
            (*output).st1 = clock();            
#endif

    }
}

void printverboselocaloutput2_(void *exp1_, double *fval, int *idp) {
    double time_local_solver;
    experiment_total *exp1;
    output_struct *output;
    
    exp1 = (experiment_total *) exp1_;
    output = exp1->output;
    
    time_local_solver = 0.0;
    
    if (exp1->test.verbose != 0) {
#ifdef MPI2
            (*output).st2 = MPI_Wtime();
            time_local_solver = (double) ((*output).st2 - (*output).st1);
#else
            (*output).st2 = clock();
            time_local_solver = ((double) ((*output).st2 - (*output).st1) / (double) CLOCKS_PER_SEC);        
#endif

        (*output).st3 = (*output).st3 + time_local_solver;


    }
}

void improvelocalsolver_(void *exp1_, double *fval, double *fvalold ) {
    local_solver *local_s;
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp1_;   
    local_s = exp1->ls;
    
        if ( *fval < *fvalold) {
            (*local_s).state_ls = 1;
            (*local_s).total_local =(*local_s).total_local  +1;
            (*local_s).sucess_local=(*local_s).sucess_local +1;
            (*local_s).num_ls++;
            (*local_s).sucess_interval++;

        } else {
            (*local_s).state_ls = 2;
            (*local_s).num_ls++;
            (*local_s).total_local++;
        }
}

void print_verbose_local_success(experiment_total exp1, int success) {
    if (exp1.test.verbose != 0) {
        if (success == 1) {
            printf(" Exito.\n");
        } else {
            printf(" falla\n");
        }
    }
}

void print_end_file_(experiment_total *exp1, double *U, double *best, result_solver *result) {
    int i, D, idp;
    output_struct *output;
    local_solver *local_s;
    
    D = exp1->test.bench.dim;
    idp = exp1->execution.idp;
    output = exp1->output;
    local_s = exp1->ls;

    if (exp1[0].test.verbose == 1) {
        //if (idp == 0)
        //    printf("Numero de chamadas solver local %.5lf (Exitos: %.5lf) tempo gastado --> %.20lf\n",
        //        local_s->total_local, local_s->sucess_local, (*output).st3);
        if (idp == 0) printf("\nBEST SOLUTION:\n");

        // FOR DIFERENTIAL EVOLUTION
        if (exp1->methodDE != NULL) {
            if (exp1[0].test._log  == 1) converttonormal_(U,&D); 
        }
        for (i = 0; i < D; i++) {
            if (idp == 0) printf("\t%.30lf\n", U[i]);
        }
        if (idp == 0) printf("f(x)=%.30lf\n", *best);


    }
}


void verboseiteration_(void *exp1_, int k,  double starttime, double ftarget, double best, long evaluation_local, int idp) {
    double mediumTime;
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp1_;
    
    if (exp1->test.verbose != 0) {
#ifdef MPI2
        mediumTime = (double) MPI_Wtime();
        printf("%d-iteracion %d best-coste %.10lf evals %.1ld -- TIME %.10lf s -- ftarget %.20lf\n\n", idp, k, best, evaluation_local,
                (double) (mediumTime - starttime), ftarget);  
#else 
        mediumTime = clock();
        
        printf("%d-iteracion %d best-coste %.10lf evals %.1ld -- TIME %.10lf s -- ftarget %.20lf\n\n", idp, k, best, evaluation_local,
                ((double) (mediumTime - starttime) / (double) CLOCKS_PER_SEC), ftarget);        
#endif


    }
}

void verboseiterationfortran_(void *exp1_, int *k, double *time,  double *ftarget, double *best, long *evaluation_local, int *par, int *idp) {
    experiment_total *exp1;

    exp1 = (experiment_total *) exp1_;
    
    if (exp1->test.verbose != 0) {
#ifdef MPI2
        printf("%d-iteracion %d best-coste %.20lf evals %ld -- TIME %lf s -- ftarget %.20lf\n\n", *idp, *k, *best, *evaluation_local,
                *time, *ftarget);  
#else 
        printf("%d-iteracion %d best-coste %.20lf evals %ld -- TIME %lf s -- ftarget %.20lf\n\n", *idp, *k, *best, *evaluation_local,
                *time, *ftarget);        
#endif


    }
}

char * concat_char(int *BUFFER, char *string_eval_total, char * string_eval) {
    char *string_eval_total_aux;
    int len1, len2;
    
    
    len1 = (int) strlen(string_eval_total);
    len2 = (int) strlen(string_eval);

    if (*BUFFER <  (len1 + len2 + 1)) {
        *BUFFER = len1*2 + len2 + 1;
        string_eval_total_aux = (char *) calloc((len1*2 + len2 + 1), sizeof(char));
        memmove( string_eval_total_aux , string_eval_total , len1);
        memmove( string_eval_total_aux+len1 , string_eval , len2);                        
                                
        string_eval_total = string_eval_total_aux;
    } else {
        string_eval_total_aux = strcat(string_eval_total,(const char *) string_eval);
        string_eval_total = string_eval_total_aux;
    }
    
    return string_eval_total;
    
}


char* delete_substring(char* str, const char* substr) {
    size_t len_str = strlen(str);
    size_t len_substr = strlen(substr);
    char* p = str;
    while ((p = strstr(p, substr)) != NULL) {
       memmove(p, p + len_substr, len_str - (p - str) - len_substr + 1);
       len_str -= len_substr;
   }
   return str;
}


/**
 * @brief this function generates a convergence MATLAB graph with the files of 
 * all slaves.
 * @param exp1_ void pointer for the main struct experiment_total.
 * @param string path of the convergence output graph of the master.
 * @param string2 path of the convergence output graph of the spacific slave.
 * @param color char variable for color lines.
 * @param ancho_linea char variable for line width.
 * @param marca char variable for the mark graph.
 * @param par binary parametter: 1 if the program is parallel, 0 in the 
 * opposite case.
 * @param idp  identification number of the processors.
 * @param NPROC total number of processors.
 * @param func binary parametter: 1 if it is the fisrt time to enter in this 
 * function, 0 in the opposite case.
 * @param end binary parametter: 1 if the specific slave is the last, 0 in the 
 * opposite case.
 * @param initp binary parametter: 1 indicates the id, where the program are 
 * going to begin to print the stairs matlab function; 0 in the opposite case.
 */
void matlab_plot_file(experiment_total exp1,  char *string,  char *string2, const char *color, int ancho_linea,  
        const char *marca, int par, int idp, int NPROC, int func, int end, int initp ) {

        char *path;
        int len,salir,init, contador_puntos;
        const char *mfile;
        const char *auxchar;
        double eval, fx;
        int iter, ls, mi;
        char *string_eval_total, *string_fx_total, *string_evalr_total, *string_ls_total, *string_mi_total;
        char *string_eval, *string_fx, *string_evalr, *string_mi, *string_ls;
        struct line_reader lr;
        /*int BIG_BUFFER = 10000;*/
        FILE *p3, *p4;
        int BUFFER1;
        int BUFFER2;
        int BUFFER3,BUFFER4,BUFFER5;
        int SIZE_PATH;
        int i,j;
        char *line;
        int contadorXXX;
        char *numberid;

        line = NULL;
        BUFFER1=1000;
        BUFFER2=1000;
        BUFFER3=1000;
        BUFFER4=1000;
        BUFFER5=1000;        
        SIZE_PATH = 5000;
        
        p3 = fopen(string2, "r");
        if (p3 == NULL) exit(0);
       
        
        mfile = ".m";
        path = (char *) malloc(SIZE_PATH*sizeof(char));
	numberid = (char *) malloc(100*sizeof(char));
	sprintf(numberid,"%d",idp);
        //strcpy(path,string,strlen(string)-4-3-strlen(numberid));
        strcpy(path, string);
        const char* substr = ".csv";
	delete_substring(path, substr);
        strcat(path, mfile);
        p4 = fopen(path, "a");  
        if (p4 == NULL) exit(0);
        
        string_fx_total = ( char *) malloc((BUFFER2+1)*sizeof(char));
        string_fx   = (char *) malloc((BUFFER2+1)*sizeof(char));     
        
        string_eval = (char *) malloc((BUFFER1+1)*sizeof(char));
        string_eval_total = ( char *) malloc((BUFFER1+1)*sizeof(char));

        string_evalr = (char *) malloc((BUFFER3+1)*sizeof(char));
        string_evalr_total = ( char *) malloc((BUFFER3+1)*sizeof(char));
        
        string_mi = (char *) malloc((BUFFER3+1)*sizeof(char));
        string_mi_total = ( char *) malloc((BUFFER3+1)*sizeof(char));
        
        string_ls = (char *) malloc((BUFFER3+1)*sizeof(char));
        string_ls_total = ( char *) malloc((BUFFER3+1)*sizeof(char));        
        
        sprintf(string_eval_total, " ");
        sprintf(string_fx_total, " ");
        sprintf(string_evalr_total, " ");        
        sprintf(string_eval, " ");
        sprintf(string_fx, " ");  
        sprintf(string_evalr, " "); 
        sprintf(string_ls_total, " ");        
        sprintf(string_ls, " ");
        sprintf(string_mi_total, " ");        
        sprintf(string_mi, " ");
        
        salir = 0;
        lr_init(&lr, p3);

    if (func == 1) {
        fprintf(p4, "function [ timeline linebest ] = convergence()\n");
        if (par == 1) {
            fprintf(p4, "function new_p = prev_interpol(x1,p1,xx);\n");
            fprintf(p4, "    new_p = ones(length(xx), 1 );\n");
            fprintf(p4, "    for i=1:length(xx)\n");
            fprintf(p4, "        ind = find(x1==xx(i));\n");
            fprintf(p4, "        if (isempty(ind)) \n");
            fprintf(p4, "            ind = max(find(x1<xx(i)));\n");
            fprintf(p4, "            if (isempty(ind) )\n");
            fprintf(p4, "                new_p(i) = unique(p1(1));\n");
            fprintf(p4, "            else\n");
            fprintf(p4, "                new_p(i) = unique(p1(ind));\n");
            fprintf(p4, "            end\n");
            fprintf(p4, "        else\n");
            fprintf(p4, "            new_p(i) = unique(p1(ind));\n");
            fprintf(p4, "        end\n");
            fprintf(p4, "    end\n");
            fprintf(p4, "    return\n");
            fprintf(p4, "end\n");
        }
    }

     contador_puntos = 0;
    eval = -1;
    fx = -1;
    init = 0;
    while (salir == 0) {
            
            line = next_line(&lr, &len);
            
            if (line != NULL ) {
                    auxchar = "%lf,%lf,%d,%d,%d";
                    sscanf(line, auxchar, &eval, &fx, &iter, &ls, &mi);
            } else {
                salir = 1;
            }
                if ((eval == -1 && fx == -1) || (salir == 1)) {
                    if (init == 1) {
                        sprintf(string_eval, " ];");
                        sprintf(string_fx, " ];");
                        sprintf(string_evalr, " ];");
                        sprintf(string_ls, " ];");
                        sprintf(string_mi, " ];");
                        
                        string_eval_total = concat_char(&BUFFER1, string_eval_total, string_eval);
                        string_ls_total = concat_char(&BUFFER1, string_ls_total, string_ls);
                        string_mi_total = concat_char(&BUFFER1, string_mi_total, string_mi);
                        string_fx_total = concat_char(&BUFFER2, string_fx_total, string_fx);
                        string_evalr_total = concat_char(&BUFFER3, string_evalr_total, string_evalr);
                        
                        fprintf(p4, "%s\n", string_eval_total);
                        fprintf(p4, "%s\n", string_mi_total);
                        fprintf(p4, "%s\n", string_ls_total);
                        fprintf(p4, "%s\n", string_fx_total);
                        fprintf(p4, "%s\n", string_evalr_total);
                        sprintf(string_eval, " ");
                        sprintf(string_fx, " ");
                        sprintf(string_evalr, " ");
                        contador_puntos++;
                    }
                    init = 1;
                    sprintf(string_eval_total, " ");
                    sprintf(string_mi_total, " ");
                    sprintf(string_ls_total, " ");
                    sprintf(string_fx_total, " ");
                    sprintf(string_eval, " ");
                    sprintf(string_mi, " ");
                    sprintf(string_ls, " ");                    
                    sprintf(string_fx, " ");   
                    sprintf(string_evalr_total, " ");
                    sprintf(string_evalr, " ");                      
                    if (par == 0) {
                        sprintf(string_fx_total, "p%d = [ ", contador_puntos+1);
                        sprintf(string_eval_total, "x%d = [ ", contador_puntos+1);
                        sprintf(string_mi_total, "mi%d = [ ", contador_puntos+1);
                        sprintf(string_ls_total, "ls%d = [ ", contador_puntos+1);
                        sprintf(string_evalr_total, "i%d = [ ", contador_puntos+1);
                    } else {
                        sprintf(string_fx_total, "p%d_%d = [ ", contador_puntos+1, idp+1);
                        sprintf(string_eval_total, "x%d_%d = [ ", contador_puntos+1, idp+1);
                        sprintf(string_mi_total, "mi%d_%d = [ ", contador_puntos+1, idp+1); 
                        sprintf(string_ls_total, "ls%d_%d = [ ", contador_puntos+1, idp+1); 
                        sprintf(string_evalr_total, "i%d_%d = [ ", contador_puntos+1,idp+1);                        
                    }
                } else {
                    sprintf(string_fx, " %.20lf ", fx);
                    sprintf(string_evalr, " %d ", iter);
                    if (ls == -1) sprintf(string_ls, " 0 "); else sprintf(string_ls, " %d ", ls);
                    if (ls == -1) sprintf(string_mi, " 0 "); else sprintf(string_mi, " %d ", mi);
                    sprintf(string_eval, " %.20lf ", eval);
                    string_eval_total = concat_char(&BUFFER1, string_eval_total, string_eval);
                    string_ls_total = concat_char(&BUFFER5, string_ls_total, string_ls);
                    string_mi_total = concat_char(&BUFFER4, string_mi_total, string_mi);
                    string_fx_total = concat_char(&BUFFER2, string_fx_total, string_fx);
                    string_evalr_total = concat_char(&BUFFER3, string_evalr_total, string_evalr);
                    eval = -1; fx = -1; iter=-1; ls = -1; mi=-1;
                }
        }
        
        
        if (!feof(p3)) {
		perror("next_line");
		exit(1);
	}
        
	lr_free(&lr);
        
        if (end == 1) {
#if MPI2        
            j=0;
            fprintf(p4, "vector_TIME_RUN%d = [", j+1);
            for (i = initp; i < NPROC+1; i++) {
                fprintf(p4, " x%d_%d(length(x%d_%d)) ", j + 1,i, j + 1,i);
            }
            fprintf(p4, "];\n");
            fprintf(p4, "vector_EVAL_RUN1 = [", j+1);
            for (i = initp; i < NPROC+1; i++) {
                fprintf(p4, " p%d_%d(length(p%d_%d)) ", j + 1,i, j + 1,i);
            }
            fprintf(p4, "];\n");
            fprintf(p4, "VTR=%.20lf;\n", exp1.test.VTR);
            fprintf(p4, "filter_FX = find(vector_EVAL_RUN%d<=VTR);\n", j+1);
            fprintf(p4, "if (isempty(filter_FX))\n");
            fprintf(p4, "[MIN_R, ii_R] = min (vector_TIME_RUN%d);\n", j+1);
            fprintf(p4, "else\n");
            fprintf(p4, "[MIN_R, ii_R] = min (vector_TIME_RUN%d(filter_FX));\n", j+1);
            fprintf(p4, "end\n");
            fprintf(p4, "index = find(vector_TIME_RUN%d == MIN_R);\n", j+1);
            fprintf(p4, "MIN_RUN1 = vector_TIME_RUN%d(index(1));\n", j+1);
            fprintf(p4, "xx1_RUN%d = [",j+1);
            for (i = initp; i < NPROC+1; i++) {
                fprintf(p4, " x%d_%d ", j+1,i);
            }
            fprintf(p4, "];\n");
            fprintf(p4, "xx2_RUN%d = unique(sort(xx1_RUN%d));\n",j+1,j+1);
            fprintf(p4, "indice_RUN%d = find(xx2_RUN%d==MIN_RUN%d);\n",j+1,j+1,j+1);
            fprintf(p4, "xx_RUN%d= xx2_RUN%d(1:indice_RUN%d);\n",j+1,j+1,j+1);
            fprintf(p4, "pp_RUN%d = ones(%d,length(xx_RUN%d));\n", j+1,NPROC,j+1);
            contadorXXX=1;
            for (i = initp; i < NPROC+1; i++) {
                fprintf(p4, " pp_RUN%d(%d,:) = prev_interpol(x%d_%d,p%d_%d,xx_RUN%d);\n", j+1,contadorXXX, j+1,i,j+1,i,j+1);
                contadorXXX++;
            }
            fprintf(p4, "if (size(pp_RUN1) > 1 )\n");
            fprintf(p4, "\t[ linpp_RUN%d, ind_RUN%d ] = min(pp_RUN%d);\n",j+1,j+1,j+1);
            fprintf(p4, "else\n");
            fprintf(p4, "\tlinpp_RUN1 = pp_RUN1;\n");
            fprintf(p4, "end;\n");
            
        fprintf(p4, "timeline = xx_RUN1;\n");
        fprintf(p4, "linebest = linpp_RUN1;\n");
        fprintf(p4, "stairs(timeline,linebest,'-r','LineWidth',2,'Markersize',4); hold on;\n");
        for (i = initp; i < NPROC+1; i++) {
            fprintf(p4, "stairs(x%d_%d,p%d_%d,'-k','LineWidth',2,'Markersize',4); hold on;\n",j+1,i,j+1,i);
        }
        fprintf(p4, "stairs(timeline,linebest,'-r','LineWidth',2,'Markersize',4); hold on;\n");
                   
        if ( exp1.test.VTR > 0  ) { 
        	fprintf(p4,"set(gca, 'YScale', 'log','XScale', 'log','FontSize',12,'TickDir','out','Box','On','LineWidth',2.0);\n");
	} else {
	        fprintf(p4,"set(gca,'FontSize',12,'TickDir','out','Box','On','LineWidth',2.0);\n");
	}
        fprintf(p4,"title('Convergence curves');\n");
        fprintf(p4,"xlabel('Wall-time (s)','FontSize',12);  ylabel('f(x)','FontSize',12);\n");
        fprintf(p4,"hleg1 = legend('minimum value in all runs', 'slave runs');\n");
        fprintf(p4,"\nend");  
#else    
        fprintf(p4, "timeline1 = x1;\n",NPROC);
        fprintf(p4, "linebest1 = p1;\n",NPROC);
        fprintf(p4,"\n\nend");        
#endif
        } 

        
        free(string_eval);
        free(string_fx);
        free(string_fx_total);
        free(string_eval_total);
        free(string_evalr);
        free(string_evalr_total);     
        free(path);
	free(numberid);
        string_eval = NULL;
        string_fx = NULL;
        string_fx_total = NULL;
        string_evalr = NULL;
        string_evalr_total = NULL;        
        string_eval_total = NULL;
        path = NULL;
        fclose(p3);
        fclose(p4);
        p3 = NULL;
        p4 = NULL;
	numberid = NULL;
}

/**
 * @brief this function generates a gantt char with the gantt files of all slaves.
 * @param exp1 main struct experiment_total.
 * @param init1 binary parametter: 1 if it is the fisrt time to enter in this 
 * function, 0 in the opposite case.
 * @param end binary parametter: 1 if the specific slave is the last, 0 in the 
 * opposite case.
 * @param stringpath path to print the matlab code.
 * @param maxsize  number of rows of the gantt file.
 */
void matlab_plot_file_gant(experiment_total exp1,  int init1, int end, char* stringpath, int maxsize) {

        char *path;
        const char *mfile;
        const char *auxchar;
        double time,OLDtime;
        int len, state, currentstate, contador;
        struct line_reader lr;
        /*int BIG_BUFFER = 10000;*/
        FILE *p3, *p4;
        int salir,i;
        char *line;
        line = NULL;
        int init;
     
        
        p3 = fopen(stringpath, "r");
        if (p3 == NULL) exit(0);
        mfile = ".m";
        path = (char *) calloc(1000,sizeof(char));
        memcpy(path,exp1.test.output_gant_log,strlen(exp1.test.output_gant_log)-4-4);
        strcat( path, mfile);
        p4 = fopen(path, "a");
        if (p4 == NULL) exit(0);
        
        lr_init(&lr, p3);


        salir = 0;
        init = 0;
        contador = 0;
        time = -1;
        state = -1;

 
        while (salir == 0) {
            
            line = next_line(&lr, &len);
            
            if (line != NULL ) {
                    auxchar = "%lf,%d";
                    sscanf(line, auxchar, &time, &state);
            } else {
                salir = 1;
            }
                
            if ((time != -1 && state != -1)) {
                if (init == 0) {
                    if (init1 == 1) {
                        fprintf(p4,"clear all;\n");
                        fprintf(p4,"close all;\n\n");
                        fprintf(p4,"data = [ %lf ", time);
                        contador++;                        
                    }
                    else {
                        fprintf(p4,"\n   %lf ", time);
                        contador++;
                    }
                    OLDtime= time;
                    currentstate=1;
                    init=1;
                } else { 
                
                if (currentstate==1) {
                    if (state == 1) {
                        fprintf(p4," %lf ", time);
                        OLDtime= time;
                        currentstate=1;
                        contador++;                        
                    }
                    if (state == 2) {
                        fprintf(p4," %lf ", time);
                        OLDtime= time;
                        currentstate=2;
                        contador++;                        
                    }                    
                    if (state == 3) {
                        fprintf(p4," %lf ", OLDtime);
                        contador++;  
                        fprintf(p4," %lf ", time);                        
                        OLDtime= time;
                        currentstate=3;
                        contador++;                        
                    }                     
                }
                
                if (currentstate==2) {
                    if (state == 1) {
                        fprintf(p4," %lf ", OLDtime);
                        contador++;  
                        fprintf(p4," %lf ", time);
                        OLDtime= time;
                        currentstate=1;
                        contador++;                        
                    }
                }  
                
                if (currentstate==3) {
                    if (state == 1) {
                        fprintf(p4," %lf ", time);
                        OLDtime= time;
                        currentstate=1;
                        contador++;                        
                    }
                }                 
                }
            }
        }
        
        if ( (contador -  (maxsize)) < 0 ) {
            contador = (maxsize)-contador;
            for (i=0;i<contador;i++) {
                fprintf(p4, " 0 ");
            }
        }
        
        if (end == 1) {
            fprintf(p4, "];\n");
            fprintf(p4, "slots = data(:,1:end);\n");
            fprintf(p4, "linex= [];\n");
            fprintf(p4, "liney= [];\n");
            
            fprintf(p4, "for i = length(slots(1,:)):-1:2\n");
            fprintf(p4, "    line=[ ];\n");
            fprintf(p4, "    for j = 1:1:length(slots(:,1))\n");
            fprintf(p4, "        if slots(j,i) > 0 \n");
            fprintf(p4, "           valueprev = slots(j,i);\n");
            fprintf(p4, "           slots(j,i) = slots(j,i)-slots(j,i-1); \n");
            fprintf(p4, "           if slots(j,i) > 0\n");            
            fprintf(p4, "                if mod(i,3) == 1 \n");
            fprintf(p4, "                    if i ~= 2\n");
            fprintf(p4, "                        linex = [linex; valueprev];\n");            
            fprintf(p4, "                        liney = [liney; j];\n");   
            fprintf(p4, "                    end\n");   
            fprintf(p4, "                end\n");   
            fprintf(p4, "           end\n");   
            fprintf(p4, "        end\n");   
            fprintf(p4, "    end\n");               
            fprintf(p4, "end\n\n");

            fprintf(p4, "linex=linex'\n");
            fprintf(p4, "liney=liney'\n");
            
            fprintf(p4, "figure\n");            
            fprintf(p4, "h = barh(slots,'stacked');\n");
            fprintf(p4, "title('Visual Logs');\n");
            fprintf(p4, "xlabel('Time','FontSize',25);\n");
            fprintf(p4, "ylabel('Processes','FontSize',25);\n");
            fprintf(p4, "set(gca, 'YDir', 'reverse','FontSize',25);\n");
            fprintf(p4, "set(h(1), 'facecolor', 'blue', 'EdgeColor', 'blue');\n");
            fprintf(p4, "for i = 2:3:(length(slots(1,:))-2)\n");
            fprintf(p4, "   set(h(i), 'facecolor', 'cyan', 'EdgeColor', 'cyan');\n");
            fprintf(p4, "   set(h(i+1), 'facecolor', 'green', 'EdgeColor', 'green');\n");
            fprintf(p4, "   set(h(i+2), 'facecolor', 'red', 'EdgeColor', 'red');\n");
            fprintf(p4, "end\n");
            fprintf(p4, "legend(h,'init stage','global search','local search', 'communications');\n");
            fprintf(p4, "hold on;\n");
            fprintf(p4, "plot(linex,liney,'rs','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor', [0.5 0 0]);\n");
            fprintf(p4, "set(gcf,'PaperPosition',[0 0 80 60]);\n");

                        
        } else {
            fprintf(p4, ";\n");
        }
     
        free(path);
        path = NULL;
        fclose(p3);
        fclose(p4);
        p3 = NULL;
        p4 = NULL;


}

/**
 * @brief this function generates a percentage char with the percentege output
 *  files of all slaves.
 * @param exp1 main struct experiment_total.
 * @param end binary parametter: 1 if the specific slave is the last, 0 in the 
 * opposite case.
 * @param stringpath path to print the matlab code.
 * @param idppar  identification number of the processor.
 */
void matlab_plot_file_porcentage(experiment_total exp1, int end, char* stringpath, int idppar) {
        int i;
        char *path;
        int init1;
        const char *mfile;
        const char *auxchar;
        double time,fx,porcentage;
        int len, contador;
        struct line_reader lr;
        /*int BIG_BUFFER = 10000;*/
        FILE *p3, *p4;
        int salir;
        char *line;
        line = NULL;
     
        init1 = 1;
        p3 = fopen(stringpath, "r");
        if (p3 == NULL) exit(0);
        mfile = ".m";
        path = (char *) calloc(500,sizeof(char));
        memcpy(path,exp1.test.log_percentage,strlen(exp1.test.log_percentage)-4);
        strcat( path, mfile);
        p4 = fopen(path, "a");
        if (p4 == NULL) exit(0);
        
        lr_init(&lr, p3);


        salir = 0;
        contador = 0;
        time = -1.0;
        porcentage = -1.0;
        fx = -1.0;
        
        while (salir == 0) {
            
            line = next_line(&lr, &len);
            
            if (line != NULL ) {
                    auxchar = "%lf,%lf,%lf";
                    sscanf(line, auxchar, &time, &fx, &porcentage);
            } else {
                salir = 1;
            }
                
            if ((time != -1.0 && porcentage != -1.0 && fx != -1.0)) {
                    if (init1 == 1) {
                        fprintf(p4,"porcent%d = [ %lf %lf %lf ;", idppar, time, fx, porcentage);
                        contador++;      
                        init1=0;
                    }
                    else {
                        fprintf(p4,"\n %lf  %lf  %lf;", time, fx, porcentage);
                        contador++;
                    }
            }
        }
        
                fprintf(p4, " ];");
        
                fprintf(p4, "\n\ntime%d = porcent%d(:,1)';\n", idppar, idppar);
                fprintf(p4, "fx%d = porcent%d(:,2)';\n", idppar, idppar);
                fprintf(p4, "sizeporcent = size(porcent%d(:,2)')';\n", idppar);
                fprintf(p4, "size%d=sizeporcent(2);\n", idppar);
                fprintf(p4, "vector%d=1:size%d;\n", idppar, idppar);
                fprintf(p4, "porc%d = porcent%d(:,3)';\n", idppar, idppar);
                fprintf(p4, "figure;\n");                
                fprintf(p4, "bar(vector%d,porc%d); hold on;\n", 
                        idppar, idppar);
                fprintf(p4, "title('Evolution of the improvement in the procesor %d','FontSize',22);\n", 
                        idppar);
                fprintf(p4, "xlabel('Number of migration','FontSize',22); \n" );
                fprintf(p4, "ylabel('Porcentage of improvement','FontSize',22);\n");
                fprintf(p4, "hold off; \n\n");
        
        
     
        free(path);
        path = NULL;
        fclose(p3);
        fclose(p4);
        p3 = NULL;
        p4 = NULL;

}

int row_count(char* stringpath) {

        const char *auxchar;
        double time;
        int len, state, currentstate, contador;
        struct line_reader lr;
        /*int BIG_BUFFER = 10000;*/
        FILE *p3;
        int salir;
        char *line;
        line = NULL;
        int init;
     
        
        p3 = fopen(stringpath, "r");
        if (p3 == NULL) exit(0);

        
        lr_init(&lr, p3);


        salir = 0;
        init = 0;
        contador = 0;
        time = -1.0;
        state = -1;
        
        while (salir == 0) {
            
            line = next_line(&lr, &len);
            
            if (line != NULL ) {
                    auxchar = "%lf,%d";
                    sscanf(line, auxchar, &time, &state);
            } else {
                salir = 1;
            }
                
            if ((time != -1 && state != -1)) {
                if (init == 0) {
                    contador++;
                    currentstate=1;
                    init=1;
                } else {
                
                  if (currentstate==1) {
                    if (state == 1) {
                        currentstate=1;
                        contador++;                        
                    }
                    if (state == 2) {
                        currentstate=2;
                        contador++;                        
                    }                    
                    if (state == 3) {
                        currentstate=3;
                        contador++;contador++;                          
                    }                     
                  }
                
                  if (currentstate==2) {
                    if (state == 1) {
                        currentstate=1;
                        contador++;contador++;                          
                    }
                  }  
                
                  if (currentstate==3) {
                    if (state == 1) {
                        currentstate=1;
                        contador++;                        
                    }
                  }  
                }
            }
        }
        
     
        fclose(p3);
        p3 = NULL;

        if ((contador-1) % 3 != 0) {
            contador = contador + (3-(contador-1)%3);
        }
        return contador;
}

double varianze(double *array, double avg, int size) {
    int i;
    double suma;
    suma = 0;
    for (i=0;i<size;i++){
        suma = pow(array[i]-avg,2);
    }
    
    suma = suma / (double) (size-1);
    
    return suma;
}

char select_color(int i) {
    if (i == 0) {
        return 'c';
    } else if (i == 1) {
        return 'm';
    } else if (i == 2) {
        return 'y';
    } else if (i == 3) {
        return 'k';
    } else if (i == 4) {
        return 'r';
    } else if (i == 5) {
        return 'g';
    } else if (i == 6) {
        return 'b';
    } else if (i == 7) {
        return 'w';
    } else {
        return 'c';
    }
    
}

char select_marca(int i) {
    int result;

    result = (int) (i / 8);

    if (result == 0) {
        return 'e';
    } else if (result == 1) {
        return 'x';
    } else if (result == 2) {
        return '*';
    } else if (result == 3) {
        return 'o';
    } else if (result == 4) {
        return '.';
    } else if (result == 5) {
        return '+';
    } else if (result == 6) {
        return 's';
    } else if (result == 7) {
        return 'd';
    } else {
        return 'e';
    }
}

/**
 * @brief this function generate a convergence MATLAB graph with the files of 
 * all slaves.
 * 
 */
void plot_file(experiment_total exp1, int par, int idp, int NPROC) {

    char *string,*string2, *string3, *string4, *string5, *string6;
    char *vectstring, *vectstring2, *vectstring3;
    int f,i,end,init,maxsize_aux,maxsize;
    
    if (par == 1) {
#ifdef MPI2
        char color;
        char marca;

        color = select_color(idp % 8);
        marca = select_marca(idp);

        string = (char *) malloc(100*sizeof(char));
        vectstring = (char *) malloc(100*NPROC*sizeof(char));
        string2 = (char *) malloc(100*sizeof(char));
        
        sprintf(string,"%s",exp1.test.output_graph);
        MPI_Gather(string, 100, MPI_CHAR, vectstring, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        if (idp == 0 ) {
            for (i=0;i<NPROC;i++) {
                memmove(string2, &vectstring[i*100], 100* sizeof(char));
                if (i == 0 ) f=1; else f=0;
                if (i == NPROC-1) end=1; else end=0;
                matlab_plot_file(exp1,string, string2, &color, 2, &marca, 1, i-1, NPROC-1,f,end, 0);
            }
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        free(string);
        free(string2);
        free(vectstring);
    
        //remove(exp1.test.log_percentage);

        
        string3 = (char *) malloc(100*sizeof(char));
        vectstring2 = (char *) malloc(100*NPROC*sizeof(char));
        sprintf(string3,"%s",exp1.test.output_gant_log);
        string4 = (char *) malloc(100*sizeof(char));
        

        MPI_Gather(string3, 100, MPI_CHAR, vectstring2, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

        maxsize = 0;
        if (idp == 0 ) {
            
            for (i=0;i<NPROC;i++) {
                memmove(string4, &vectstring2[i*100], 100* sizeof(char));
                maxsize_aux = row_count(string4);
                if ( maxsize_aux > maxsize) maxsize = maxsize_aux;
            }

            for (i=0;i<NPROC;i++) {
                memmove(string4, &vectstring2[i*100], 100* sizeof(char));
                if (i == 0) init=1; else init=0;
                if (i == (NPROC-1)) end=1; else end=0;
                matlab_plot_file_gant(exp1, init, end, string4,maxsize);
            }        
        }

        
        MPI_Barrier(MPI_COMM_WORLD);
        free(string3);
        free(string4);
        free(vectstring2);        
        
#endif

    } else {
        char color;
        char marca;
        string = (char *) malloc(100*sizeof(char));

        sprintf(string,"%s",exp1.test.output_graph);

        matlab_plot_file(exp1,exp1.test.output_graph, string,&color,2, &marca, 0, 0, NPROC,1,1,0);
	free(string);
	string=NULL;
    }
    
    
    
}

/**
 * @brief this function manages the generation of the differents graphs.
 * @param exp1 main struct experiment_total.
 * @param par binary parametter: 1 if the program is parallel, 0 in the 
 * opposite case.
 * @param idp identification number of the processors.
 * @param NPROC total number of processors.
 */
void plot_file_cess(experiment_total exp1, int par, int idp, int NPROC) {
    char *string,*string2, *string3, *string4, *string5, *string6;
    char *vectstring, *vectstring2, *vectstring3;
    int maxsize_aux, maxsize; 
    int f,i,end,init;
    
    if (par == 1) {
#ifdef MPI2

        char color;
        char marca;

        color = select_color(idp % 8);
        marca = select_marca(idp);

        string = (char *) malloc(100*sizeof(char));
        vectstring = (char *) malloc(100*NPROC*sizeof(char));
        sprintf(string,"%s",exp1.test.output_graph);
        string2 = (char *) malloc(100*sizeof(char));
        
        MPI_Gather(string, 100, MPI_CHAR, vectstring, 100, MPI_CHAR, 1, MPI_COMM_WORLD);
        
            
        if (idp == 1 ) {
            for (i=1;i<NPROC;i++) {
                memmove(string2, &vectstring[i*100], 100* sizeof(char));
                if (i == 1 ) f=1; else f=0;
                if (i == NPROC-1) end=1; else end=0;
                matlab_plot_file(exp1,string, string2, &color, 2, &marca, 1, i-1, NPROC-1,f,end, idp);
//                remove(string2);                
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        free(string);
        free(string2);
        free(vectstring);
        string=NULL; string2=NULL; vectstring=NULL;
        
        string5 = (char *) malloc(100*sizeof(char));
        vectstring3 = (char *) malloc(100*NPROC*sizeof(char));
        sprintf(string5,"%s",exp1.test.log_percentage);
        string6 = (char *) malloc(100*sizeof(char));
        MPI_Gather(string5, 100, MPI_CHAR, vectstring3, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        
        if (idp == 0 ) {
		i=0;
                memmove(string6, &vectstring3[i*100], 100* sizeof(char));
                end=1;
                matlab_plot_file_porcentage(exp1, end, string6, 0);
//                remove(string6);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        

        free(string5);
        free(string6);
        free(vectstring3);          
        
        string3 = (char *) malloc(100*sizeof(char));
        vectstring2 = (char *) malloc(100*NPROC*sizeof(char));
        sprintf(string3,"%s",exp1.test.output_gant_log);
        string4 = (char *) malloc(100*sizeof(char));
        
        MPI_Gather(string3, 100, MPI_CHAR, vectstring2, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
        maxsize = 0;
        
        if (idp == 0 ) {
            for (i=0;i<NPROC;i++) {
                memmove(string4, &vectstring2[i*100], 100* sizeof(char));
                maxsize_aux = row_count(string4);
                if ( maxsize_aux > maxsize) maxsize = maxsize_aux;
            }
            
            
            for (i=0;i<NPROC;i++) {
                memmove(string4, &vectstring2[i*100], 100* sizeof(char));
                if (i==0) init=1; else init=0;
                if (i == (NPROC-1)) end=1; else end=0;
                matlab_plot_file_gant(exp1, init, end, string4, maxsize); 
 //               remove(string4);
            }        
        }

        
        MPI_Barrier(MPI_COMM_WORLD);
        free(string3);
        free(string4);
        free(vectstring2);
        
#endif

    } else {
        char color;
        char marca;
        
        color = select_color(idp % 8);
        marca = select_marca(idp);
        
        matlab_plot_file(exp1,exp1.test.output_graph, exp1.test.output_graph,&color,2, &marca, 0, 0, NPROC,1,1,0);


    }
    
    
    
}


void initoutputvars_(void *exp1_){
    experiment_total *exp1;        
    exp1 = (experiment_total *) exp1_;

    output_struct *output;
    
    output = exp1->output;
    output->point_counter = 0;
    output->point1=0;
    output->oldbest=DBL_MAX;
    output->st1=0.0;
    output->st2=0.0;
    output->st3=0.0;
    
    
    
}



void updateresultsandprint_(void *exp1_, void *result_, double *totaltime, long *evaluation, double *fx, double *xbest) {
    int idp, NPROC, D; 
    int i, indexvector;
    double auxbest;
    result_solver *result;
    output_struct *output;
    local_solver *local_s;
    experiment_total *exp1;
    double st3;
    double *vectorbestx,*vectorbest,*bestglobalx,bestglobal;
    exp1 = (experiment_total *) exp1_;


    st3 = 0.0;
    idp = exp1->execution.idp;
    NPROC = exp1->execution.NPROC;
    D = exp1->test.bench.dim;
    result = (result_solver *) result_;
    output =  exp1->output;
    local_s = exp1->ls;
    double totalL, sucessL;



#ifdef MPI2
    MPI_Reduce(&(local_s->total_local), &totalL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(local_s->sucess_local), &sucessL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(output->st3), &st3, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
       
    if (idp == 0) {
            local_s->total_local = totalL / (double) NPROC;
           local_s->sucess_local = sucessL/ (double) NPROC;
            output->st3 = st3 / (double) NPROC;
    }
        
    MPI_Bcast(&(local_s->total_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(local_s->sucess_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(output->st3), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
            
    vectorbestx = (double *) malloc(NPROC*(D)*sizeof(double));
    vectorbest = (double *) malloc(NPROC*sizeof(double));
    bestglobalx = (double *) malloc((D)*sizeof(double));
    MPI_Gather(xbest,D,MPI_DOUBLE,vectorbestx,D,MPI_DOUBLE,0,exp1->execution.topology.comunicator);
    MPI_Gather(fx,1,MPI_DOUBLE,vectorbest,1,MPI_DOUBLE,0,exp1->execution.topology.comunicator);
        
    if (idp == 0) {
            auxbest = DBL_MAX;
            indexvector = -1;
            for (i=0;i<NPROC;i++) {
                printf("vectorbest[%d]- %.20lf\n",i,vectorbest[i]);
                if ( vectorbest[i] < auxbest ) {
                    auxbest = vectorbest[i];
                    indexvector = i;
                }
            }
          
            bestglobal=auxbest;
	    for (i=0;i++;i<D)
                bestglobalx[i] = vectorbestx[indexvector*(D)+i];     
    }
    MPI_Bcast(bestglobalx, D, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&bestglobal, 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    
    
    if ((*result).best_value >= bestglobal) {
            (*result).best_value = bestglobal;
    }
        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;  
        
        print_end_file_(exp1,bestglobalx,&bestglobal,result);
        
        free(vectorbestx);
        vectorbestx=NULL;
        free(vectorbest  );
        vectorbest=NULL;
        free(bestglobalx  );
        bestglobalx=NULL;
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        
#else 
    
        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;

        print_end_file_(exp1,xbest,fx,result);
        
        
        if ((*result).best_value >= *fx) {
            (*result).best_value = *fx;
        }
#endif
        
        
    
}



/**
 * @brief this function gathers the different results of all processors in 
 * the master, and then the master diffuses this information to the slaves.
 * @param exp1_ main struct experiment_total.
 * @param result_ result struct.
 * @param totaltime total time spents by the program.
 * @param evaluation total number of evaluations spents by the program.
 * @param fx value of the objective function of the best solution 
 * for each slave.
 * @param xbest vector of the best solution for each slave.
 * @param totaliter otal iteration spents by the solver.
 * @param timevtr time to reach to VTR.
 */
void updateresultsess_(void *exp1_, void *result_, double *totaltime, long *evaluation, double *fx, double *xbest,
        double *totaliter, double *timevtr) {
    int idp, NPROC, D; 
    int i, indexvector;
    double auxbest;
    result_solver *result;
    output_struct *output;
    local_solver *local_s;
    experiment_total *exp1;
    double st3;
    double *vectorbestx,*vectorbest,*bestglobalx,bestglobal;
    long fails;
    exp1 = (experiment_total *) exp1_;

    sleep(2);
    st3 = 0.0;
    idp = exp1->execution.idp;
    NPROC = exp1->execution.NPROC;
    D = exp1->test.bench.dim;
    result = (result_solver *) result_;
    output =  exp1->output;
    local_s = exp1->ls;
    double totalL, sucessL;



#ifdef MPI2
    MPI_Reduce(&(local_s->total_local), &totalL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(local_s->sucess_local), &sucessL, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(output->st3), &st3, 1, MPI_DOUBLE, MPI_SUM, 0, exp1->execution.topology.comunicator);
    MPI_Reduce(&(exp1->execution.failevals), &fails, 1, MPI_LONG, MPI_SUM, 0, exp1->execution.topology.comunicator);   
    exp1->execution.failevals = fails;
    if (idp == 0) {
            local_s->total_local = totalL / (double) NPROC;
           local_s->sucess_local = sucessL/ (double) NPROC;
            output->st3 = st3 / (double) NPROC;
    }
        
    MPI_Bcast(&(local_s->total_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(local_s->sucess_local), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(output->st3), 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&(exp1->execution.failevals), 1, MPI_LONG, 0, exp1->execution.topology.comunicator);
    
    vectorbestx = (double *) malloc(NPROC*(D)*sizeof(double));
    vectorbest = (double *) malloc(NPROC*sizeof(double));
    bestglobalx = (double *) malloc((D)*sizeof(double));
    MPI_Gather(xbest,D,MPI_DOUBLE,vectorbestx,D,MPI_DOUBLE,0,exp1->execution.topology.comunicator);
    MPI_Gather(fx,1,MPI_DOUBLE,vectorbest,1,MPI_DOUBLE,0,exp1->execution.topology.comunicator);
       
    (*result).iterations_value =  *totaliter;
    (*result).time_vtr_value = *timevtr;
    if (idp == 0) {
            auxbest = vectorbest[0];
            indexvector = 0;
            for (i=0;i<NPROC;i++) {
                if ( vectorbest[i] < auxbest ) {
                    auxbest = vectorbest[i];
                    indexvector = i;
                }
            }
          
            bestglobal=auxbest;
            memmove(bestglobalx,&vectorbestx[indexvector*(D)],D*sizeof(double));     
    }
    MPI_Bcast(bestglobalx, D, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    MPI_Bcast(&bestglobal, 1, MPI_DOUBLE, 0, exp1->execution.topology.comunicator);
    
    
    if ((*result).best_value >= bestglobal) {
            (*result).best_value = bestglobal;
    }
    
    
        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;  
        
        MPI_Barrier(MPI_COMM_WORLD);        
        print_end_file_(exp1,bestglobalx,&bestglobal,result);
        sleep(2);
        MPI_Barrier(MPI_COMM_WORLD);

        memmove(result->bestx_value, bestglobalx, D*sizeof(double));
        
        free(vectorbestx);
        vectorbestx=NULL;
        free(vectorbest);
        vectorbest=NULL;
        free(bestglobalx);
        bestglobalx=NULL;
        
        
        
#else 
        if  (exp1->test.bench.translation == 1)  detranslation( xbest, exp1->execution.transconst, exp1->test.bench.dim );

        (*result).time_value = (*result).time_value + *totaltime;
        (*result).eval_value = (*result).eval_value + *evaluation;
        (*result).iterations_value = (double) *totaliter;
        (*result).time_vtr_value = (int) *timevtr;
    
        print_end_file_(exp1,xbest,fx,result);
	
        for (i=0;i<D;i++) (*result).bestx_value[i] = xbest[i];

        if ((*result).best_value >= *fx) {
            (*result).best_value = *fx;
        }
#endif
        
        
    
}




int updateresultsrandomsearch_(void *exp1_, void *result_, double *totaltime, double *evaluation, double *best ) {
    result_solver *result;

    result = (result_solver *) result_;
    (*result).time_value = (*result).time_value + *totaltime;
    (*result).eval_value = (*result).eval_value + *evaluation;
    (*result).best_value = *best;

   
    return 1;
}




void init_message(int NPROC, experiment_total *exp, int NTHREAD) {
    const char *namealg,*nameversionalg,*nametopologyalg,*namels;
    int idsolver;
    
    idsolver = getnumversion(exp);
    namealg  = getname(exp);
    nameversionalg = getversioness(idsolver);
    nametopologyalg = gettopologyess(idsolver);
    namels = getlsess(exp);
    
    printf("\n\n=====================================================\n");
    printf("== PARALLEL SACESS OPTIMIZATION LIBRARY =============\n");
    printf("=====================================================\n\n");  
    
    printf("************ GENERAL INFORMATION ********************\n\n");

    printf("* SOLVER :: %s\n", namealg);
    printf("    Version solver : %s\n", nameversionalg);
    printf("    Topology       : %s\n", nametopologyalg);    
    printf("    Max.Evals      : %.1lf\n", (*exp).test.max_eval );
    printf("    Max.Time(s)    : %.1lf\n", (*exp).test.maxtime );
    
    if ( exp[0].test.local_search == 0) {
        printf("    Local Solver   : none\n\n", namels);
    } else {
        printf("    Local Solver   : %s\n\n", namels);        
    }
    printf("* TOTAL PROCESSORS USED :: %d\n", NPROC*NTHREAD);
    printf("    Number of MPI proc.                : %d\n", NPROC);
    printf("    Number of openMP threads per proc. : %d\n\n", NTHREAD);
}



void bechmark_message(experiment_total *exp, const char *name) {
    printf("* TYPE OF BENCKMARK :: %s\n", (*exp).test.bench.type );
    printf("    Name                     : %s\n", name);
    printf("    Number of parameters     : %d\n", (*exp).test.bench.dim);
    printf("    Value To Reach           : %.10lf\n", (*exp).test.VTR);
    printf("    [DEFAULT] Value To Reach : %.10lf\n\n", (*exp).test.VTR_default);
    
    printf("************ SLAVE INFORMATION ********************\n");
}



void graphs_message(experiment_total *exp) {
    int idsolver;
    char *m1,*m2,*m3;
 
    m1=(char *) calloc(2000,sizeof(m1));
    m2=(char *) calloc(2000,sizeof(m2));
    m3=(char *) calloc(2000,sizeof(m3));

    printf("\n************ GRAPHS GENERATED ***********************\n\n");
    strcpy(m1,exp->test.output_path);
    printf(" GRAPH PATH               : %s\n", m1);
    printf(" CONVERGENCE GRAPH MATLAB : %s/convergence_id*.csv.m\n", m1);
    printf(" CONVERGENCE GRAPHs CSV   : %s/convergence_id*.csv\n", m1);
#ifdef MPI2    
    idsolver = getnumversion(exp);
    memcpy(m2,exp->test.output_gant_log,strlen(exp->test.output_gant_log)-5);
    printf(" GANTT CHART MATLAB       : %s/gantt.m\n", m1);
    printf(" GANTT CHARTs CSV         : %s*.csv\n", m2);
    
    if (idsolver == 3) {
        memcpy(m3,exp->test.log_percentage,strlen(exp->test.log_percentage)-4);
    printf(" PERCENTAGE GRAPH MATLAB  : %s.m\n", m3);
    printf(" PERCENTAGE GRAPH CSV     : %s.csv\n\n\n", m3);
    }
#endif    
    free(m1);
    free(m2);
    free(m3);
}


void printresults_(void *exp1_, double *fx, long *evals, double *totaltime, double *timeVTR, double *iter) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
        
   // printf("- fx evals totaltime timeVTR iter evalfails\n");
   // printf("tab %.20f %ld %lf %lf %lf %ld\n", *fx, *evals, *totaltime, *timeVTR, *iter, exp1->execution.failevals);	

}



void printstatemaster_(double *BKS, double *time){
    printf("[MASTER] Best Known Solution fx = %.10lf -- CURRENT TIME %lf\n", *BKS, *time);
}


void printresults_end(experiment_total *exptotal, result_solver result) {
    double time;
    if (result.time_vtr_value == DBL_MAX )
        printf("\nTOTAL TIME TO REACH VTR : not reach\n");
    else 
        printf("\nTOTAL TIME TO REACH VTR : %lf\n", result.time_vtr_value);
        printf("TOTAL EVALUATIONS       : %ld\n", result.eval_value);
        printf("AVERAGE ITERATIONS      : %lf\n", result.iterations_value);
        printf("EVALUATION FAILED       : %ld\n", exptotal->execution.failevals);
    if (result.time_vtr_value == DBL_MAX ){
        printf("\nf-%d DIM %d fbest: %.20lf < ftarget %.20lf numevals (%ld) -- time VTR -> not reach -- iter %lf\n",
            exptotal[0].test.bench.current_bench,
            exptotal[0].test.bench.dim,
                        result.best_value,
                        exptotal[0].test.VTR,
                        result.eval_value,
                        result.iterations_value);
    } else {
        printf("\nf-%d DIM %d fbest: %.20lf < ftarget %.20lf numevals (%ld) -- time VTR -> %lf -- iter %lf\n",
            exptotal[0].test.bench.current_bench,
            exptotal[0].test.bench.dim,
                        result.best_value,
                        exptotal[0].test.VTR,
                        result.eval_value,
                        result.time_vtr_value,
                        result.iterations_value);        
    }
}



/**
 * @brief this function generate a convergence MATLAB graph with the files of 
 * all slaves.
 */
void plot(experiment_total *exptotal) {
    int NPROC,id;
    int idsolver;
    
    if (exptotal->test.output == 1) {
        idsolver = getnumversion(exptotal);
        NPROC=1;
        id=0;
#ifdef MPI2
        MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
        MPI_Comm_rank(MPI_COMM_WORLD, &id);   

        if ((*exptotal).methodScatterSearch != NULL) {
            if ((idsolver == 1 ) ||
                (idsolver == 3)     
                    ) {
                plot_file_cess(*exptotal, 1, id, NPROC);	
            } else  {
                plot_file(*exptotal, 1, id, NPROC);
            }

        } else {
	    plot_file(*exptotal, 1, id, NPROC);
        }


#else 
        plot_file(*exptotal, 0, id,1);
#endif
    }
}
    




