/**
 * @file parallel_functions_cooperative_eSS.c
 * @author David R. Penas
 * @brief File containing specific function for the parallelization of the eSS.
 */


#ifdef MPI2

#include <math.h>
#include <structure_paralleltestbed.h>
#include <parallel_functions.h>
#include <mpi.h>
#include <string.h>



// Function no que cada nodo manda AO ANILLO a solucion seleccionada 
void sendbestsolutiondist_( void *exp1_, double *serializedata, int *sizeser, int *init ) {
   experiment_total *exp1;
   exp1 = (experiment_total *) exp1_;
   int i, flag;
   MPI_Status status;

// Comunicar al de atras en el anillo
   flag=0;
   if (*init >= 2) {
         MPI_Test( &(exp1[0].execution.sendrequestmaster[0]) , &flag, MPI_STATUS_IGNORE);
   } else {
         flag = 1;
         *init = *init + 1;
   }
 
   if (flag == 1) {
      	memmove(exp1[0].execution.sendbuffermaster[0], serializedata, (*sizeser)*sizeof(double));
	if(exp1[0].execution.idp>0)
         	MPI_Isend( exp1[0].execution.sendbuffermaster[0], *sizeser, MPI_DOUBLE,exp1[0].execution.idp-1, 100, 
			MPI_COMM_WORLD, &(exp1[0].execution.sendrequestmaster[0]));
	else
 		MPI_Isend( exp1[0].execution.sendbuffermaster[0], *sizeser, MPI_DOUBLE,exp1[0].execution.NPROC-1, 100, 
 				MPI_COMM_WORLD, &(exp1[0].execution.sendrequestmaster[0]));


   }   
   
// Comunicar al de delante en el anillo
   flag=0;
   if (*init >= 2) {
               MPI_Test( &(exp1[0].execution.sendrequestmaster[1]) , &flag, MPI_STATUS_IGNORE);
   } else {
        flag = 1;
        *init = *init + 1;
   }

   if (flag == 1) {
        memmove( exp1[0].execution.sendbuffermaster[1], serializedata, (*sizeser)*sizeof(double));
	if(exp1[0].execution.idp<(exp1[0].execution.NPROC-1))
          	MPI_Isend(exp1[0].execution.sendbuffermaster[1], *sizeser, MPI_DOUBLE,exp1[0].execution.idp+1, 100, 
 				MPI_COMM_WORLD, &(exp1[0].execution.sendrequestmaster[1]));
	else
 		MPI_Isend(exp1[0].execution.sendbuffermaster[1], *sizeser, MPI_DOUBLE,0, 100, 
 				MPI_COMM_WORLD, &(exp1[0].execution.sendrequestmaster[1]));

   }


}

 // Funcion no cada nodo consulta si hai soluci�ns difundidas neste.
void returnwindowvaluesdist_(void *exp1_, double *serializedata, int *sizeser, int *returnflag, int *idsent ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    int flag;
    MPI_Status status;
    int i, id, size, idproc;
    double auxSOLUTION;
    

    auxSOLUTION = DBL_MAX;
    flag = 1;
    *returnflag = 0;
    size = 2;
 
    while (flag != 0) {
       flag = 0;
       id = -1;
       // comprobar si alguno de los dos (el de atras o el de delante en el anillo) han comunicado
       MPI_Testany( size, &(exp1[0].execution.receptionrequestmaster[0]) ,&id,&flag,MPI_STATUS_IGNORE);
 
       if (flag == 1) {
          
          if ( id == 0 ) {
		if ( exp1[0].execution.idp == 0 ) idproc = exp1[0].execution.NPROC-1;
                else idproc = exp1[0].execution.idp-1;
	  } else {
                if ( exp1[0].execution.idp == (exp1[0].execution.NPROC-1) ) idproc = 0;
                else idproc = exp1[0].execution.idp+1;
 	  }

          if ( exp1[0].execution.receptionbuffermaster[id][*sizeser-1] < auxSOLUTION  ) {
              auxSOLUTION = exp1[0].execution.receptionbuffermaster[id][*sizeser-1];
 	      memmove( serializedata, exp1[0].execution.receptionbuffermaster[id] , (*sizeser)*sizeof(double));
              *returnflag = 1;
              *idsent = idproc;
          } 
          else {
              if (id == 1)   
 	      printdescartsolution_(exp1, &(exp1[0].execution.receptionbuffermaster[id][*sizeser-1]), &idproc);
          }
          MPI_Irecv(exp1[0].execution.receptionbuffermaster[id],*sizeser,MPI_DOUBLE,idproc, 100,MPI_COMM_WORLD,
                                                                     &(exp1[0].execution.receptionrequestmaster[id]));
       }
   }

}


// Iniciamos as variables Request e Buffers necesarios na version distribuida sen maestro do acess
void asynchinitasynchronousvars_(void *exp1_, int *sizeser ) {
	 experiment_total *exp1;
	 exp1 = (experiment_total *) exp1_;
	 int i;

 // Solo se comunica con el anterior y el posterior en el anillo
	 exp1[0].execution.receptionbuffermaster = (double **) malloc(2* sizeof (double));
	 exp1[0].execution.sendbuffermaster = (double **) malloc(2* sizeof (double));
	 exp1[0].execution.receptionrequestmaster = (MPI_Request *) malloc(2 * sizeof (MPI_Request));
	 exp1[0].execution.sendrequestmaster = (MPI_Request *) malloc(2 * sizeof (MPI_Request));
 
 // solo hay dos con los que comunicar
         exp1[0].execution.receptionbuffermaster[0] = (double *) malloc(*sizeser *sizeof(double));
         exp1[0].execution.receptionbuffermaster[1] = (double *) malloc(*sizeser * sizeof(double));
         exp1[0].execution.sendbuffermaster[0] = (double *) malloc(*sizeser *sizeof(double));
         exp1[0].execution.sendbuffermaster[1] = (double *) malloc(*sizeser *sizeof(double));     

 	 if(exp1[0].execution.idp>0)
      	     MPI_Irecv(exp1[0].execution.receptionbuffermaster[0], *sizeser, MPI_DOUBLE,exp1[0].execution.idp-1,
                        100,MPI_COMM_WORLD, &(exp1[0].execution.receptionrequestmaster[0]));
         else
 	     MPI_Irecv(exp1[0].execution.receptionbuffermaster[0], *sizeser, MPI_DOUBLE,exp1[0].execution.NPROC-1,
            		100,MPI_COMM_WORLD, &(exp1[0].execution.receptionrequestmaster[0]));
	

	 if(exp1[0].execution.idp<(exp1[0].execution.NPROC-1))
      	     MPI_Irecv(exp1[0].execution.receptionbuffermaster[1], *sizeser, MPI_DOUBLE,exp1[0].execution.idp+1,
             		100,MPI_COMM_WORLD, &(exp1[0].execution.receptionrequestmaster[1]));
	 else
 	     MPI_Irecv(exp1[0].execution.receptionbuffermaster[1], *sizeser, MPI_DOUBLE,0,
                        100,MPI_COMM_WORLD, &(exp1[0].execution.receptionrequestmaster[1]));



}
	
void destroyasynchvars_(void *exp1_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    int i;

    for (i = 0; i < 2; i++) {
            free(exp1[0].execution.receptionbuffermaster[i] );
            free(exp1[0].execution.sendbuffermaster[i] );
    }

    free(exp1[0].execution.receptionrequestmaster );
    free(exp1[0].execution.sendrequestmaster );

}



void mpitest_(void *exp1_, int *init){
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    int flag,i;
    for (i = 0; i < (exp1[0].execution.NPROC); i++) {
        if (i != exp1[0].execution.idp){
            if (init[i] == 1)
                MPI_Test( &exp1[0].execution.sendrequestmaster[i] , &flag, MPI_STATUS_IGNORE);
        }
    }
}


void checktimemaster_(void *exp1_, double *currenttime, int *flagoutput, int *init) {
	experiment_total *exp1;
    	exp1 = (experiment_total *) exp1_;
        double tao;
        	
        if ( *init == 0 ) {
		exp1[0].execution.currentmastertime =  *currenttime;
                *init=1;
        } else {
        	tao = *currenttime - exp1[0].execution.currentmastertime;
		if ( tao >= exp1[0].execution.mastertime) {
                        exp1[0].execution.currentmastertime =  *currenttime;
			*flagoutput=1;
		} else {
			*flagoutput=0;
		}
        }

}

// Function de recepcion dos mellores valores por parte do master
void  receivesolutionsmaster_ (void *exp1_, int *dest, double *serializedata, int *sizeser, int *result)  {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;   
    int flag;
    MPI_Status status;
    flag = 0;

    MPI_Test( &exp1[0].execution.receptionrequestmaster[*dest-1], &flag, &status);
    
    if (flag == 1) {
        memmove(serializedata, exp1[0].execution.receptionbuffermaster[(*dest-1)], (*sizeser)*sizeof(double));
//      free(exp1[0].execution.receptionbuffermaster[(*dest-1)]);    
//      exp1[0].execution.receptionbuffermaster[(*dest-1)] = (double *) malloc( *sizeser * sizeof(double));
        MPI_Irecv(exp1[0].execution.receptionbuffermaster[*dest-1], *sizeser, MPI_DOUBLE,*dest, 150, MPI_COMM_WORLD, &exp1[0].execution.receptionrequestmaster[*dest-1]);
            
    } else flag=0;
    *result=flag;
}



// Function no que o master actualiza con un valor os slaves
void masterupdatesolutionslaves_( void *exp1_, double *serializedata, int *sizeser, int *sender, int *init) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;   
    int i, flag;
    MPI_Status status; 
    
    for (i = 0; i < (exp1[0].execution.NPROC-1); i++) {
	       memmove(exp1[0].execution.sendbuffermaster[i], serializedata, *sizeser*sizeof(double));
	       MPI_Isend(exp1[0].execution.sendbuffermaster[i], *sizeser, MPI_DOUBLE,i+1, 150, MPI_COMM_WORLD,&exp1[0].execution.sendrequestmaster[i]);
    }
}


// Funcion no que o esclavo obten o valor que hai no buffer de recepcion se existe
void returnwindowvalueslave_(void *exp1_, double *serializedata, int *sizeser, int *returnflag) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;   
    int flag;
    MPI_Status status;
 

    flag = 1;
    *returnflag = 0;
    while (flag != 0) {
        flag = 0;
	MPI_Test( exp1[0].execution.receptionrequestslave, &flag, MPI_STATUS_IGNORE);
    	if (flag == 1) {
    		memmove( serializedata,exp1[0].execution.receptionbufferslave, *sizeser*sizeof(double));
        	MPI_Irecv(exp1[0].execution.receptionbufferslave, *sizeser, MPI_DOUBLE, 0, 150, MPI_COMM_WORLD, exp1[0].execution.receptionrequestslave);
        	*returnflag = 1;
    	}  
    }
}


// funcion no que o esclavo envia a sua mellor solucion o master
void sendbestsolutionslave_ (void *exp1_, double *serializedata, int *sizeser, int *init) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;   
    MPI_Status status;
    
    if (*init == 0) {
        *init = 1;
    } 
    else {
        MPI_Wait(exp1[0].execution.sendrequestslave, &status);
    }
    
    memmove ( exp1[0].execution.sendbufferslave , serializedata, *sizeser * sizeof(double) );    
    MPI_Isend( exp1[0].execution.sendbufferslave , *sizeser, MPI_DOUBLE, 0, 150,
                MPI_COMM_WORLD, exp1[0].execution.sendrequestslave);  

}



// Funcion no que os esclavos e o master inicializan os seus buffers e requests de recepcion e envio
void asynchinitmasterandwindows_(void *exp1_, int *sizeser ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;       
    int i;
    
     exp1[0].execution.sendrequestslave = (MPI_Request *) malloc(sizeof(MPI_Request));
     exp1[0].execution.sendbufferslave = (double  *) malloc((*sizeser)*  sizeof(double) );

     exp1[0].execution.receptionrequestslave= (MPI_Request *) malloc(sizeof(MPI_Request));
     exp1[0].execution.receptionbufferslave = (double  *) malloc((*sizeser)*  sizeof(double) );
     
     exp1[0].execution.adaptation_slave_recv = (MPI_Request  *) malloc( sizeof(MPI_Request) );
     exp1[0].execution.adaptation_slave_buffer_recv = (double  *) malloc( 4 * sizeof(double) );
     
     exp1[0].execution.adaptation_slave_send = (MPI_Request  *) malloc( sizeof(MPI_Request) );
     exp1[0].execution.adaptation_slave_buffer_send = (int  *) malloc(  sizeof(int) );

     for (i=0;i<(*sizeser);i++) {
        exp1[0].execution.receptionbufferslave[i] = DBL_MAX;
     }

// MASTER
     if (exp1[0].execution.idp == 0) {
     	 exp1[0].execution.receptionbuffermaster = (double **) malloc((exp1[0].execution.NPROC-1)* sizeof (double));
     	 exp1[0].execution.sendbuffermaster = (double **) malloc((exp1[0].execution.NPROC-1)* sizeof (double));
         
         exp1[0].execution.adaptation_master_buffer_recv = (int *) malloc( (exp1[0].execution.NPROC-1) * sizeof(int) );
         exp1[0].execution.adaptation_master_buffer_send = (double **) malloc( (exp1[0].execution.NPROC-1) * sizeof(double) );
                  
     	 exp1[0].execution.receptionrequestmaster = (MPI_Request *) malloc((exp1[0].execution.NPROC-1) * sizeof (MPI_Request));
         exp1[0].execution.sendrequestmaster      = (MPI_Request *) malloc((exp1[0].execution.NPROC-1) * sizeof (MPI_Request)); 

         exp1[0].execution.adaptation_master_recv = (MPI_Request *) malloc((exp1[0].execution.NPROC-1) * sizeof (MPI_Request)); 
         exp1[0].execution.adaptation_master_send = (MPI_Request *) malloc((exp1[0].execution.NPROC-1) * sizeof (MPI_Request)); 
         
         for (i = 0; i < (exp1[0].execution.NPROC-1); i++) {
	      exp1[0].execution.receptionbuffermaster[i] = (double *) malloc((*sizeser)* sizeof (double));
              exp1[0].execution.adaptation_master_buffer_send[i] = (double *) malloc( 4 * sizeof (double));
              
	      exp1[0].execution.sendbuffermaster[i] = (double *) malloc((*sizeser)* sizeof (double));

              MPI_Irecv(exp1[0].execution.receptionbuffermaster[i], *sizeser, MPI_DOUBLE,i+1, 150,
                 MPI_COMM_WORLD, &exp1[0].execution.receptionrequestmaster[i]);
              
              MPI_Irecv(&exp1[0].execution.adaptation_master_buffer_recv[i], 1, MPI_INT,i+1, 155,
                 MPI_COMM_WORLD, &exp1[0].execution.adaptation_master_recv[i]);              
         }

    } else {
// SLAVES
       MPI_Irecv(exp1[0].execution.receptionbufferslave, *sizeser, MPI_DOUBLE,0, 150,
                 MPI_COMM_WORLD, exp1[0].execution.receptionrequestslave );
       
       MPI_Irecv(exp1[0].execution.adaptation_slave_buffer_recv, 4, MPI_DOUBLE,0, 155,
                 MPI_COMM_WORLD, exp1[0].execution.adaptation_slave_recv );       
       
    } 
     
}


void adaptationcheck_(void *exp1_, int *resultflag_adapt, int *idslave_to_adapt) {
    experiment_total *exp1;    
    int i, flag;

    exp1 = (experiment_total *) exp1_;  
    *resultflag_adapt = 0;   
    *idslave_to_adapt = -1;
    
    for (i=0;i<(exp1[0].execution.NPROC-1);i++) {
        MPI_Test(&(exp1[0].execution.adaptation_master_recv[i]),&flag,MPI_STATUS_IGNORE);
        
        if ((flag == 1)){
            *resultflag_adapt = exp1[0].execution.adaptation_master_buffer_recv[i];
            *idslave_to_adapt = i;
            
            MPI_Irecv(&exp1[0].execution.adaptation_master_buffer_recv[i], 1, MPI_INT,i+1, 155,
                 MPI_COMM_WORLD, &exp1[0].execution.adaptation_master_recv[i]);              
            break;
        }
    }
    
    
}

void adaptationsendmaster_(void *exp1_, int *idsend, int *size_dim, int *ncounter, double *balance, int *accept) {
    experiment_total *exp1;    
    double *buffersend;
    exp1 = (experiment_total *) exp1_;  
    
    buffersend = exp1[0].execution.adaptation_master_buffer_send[*idsend-1];

    buffersend[0] = (double) *size_dim;
    buffersend[1] = (double) *ncounter;
    buffersend[2] = *balance;
    buffersend[3] = *accept;
    
    MPI_Isend(buffersend, 4, MPI_DOUBLE, *idsend, 155, 
            MPI_COMM_WORLD, &exp1[0].execution.adaptation_master_send[*idsend-1] );
}


void sendslaveadaptsignal_(void *exp1_, int *code) {
    experiment_total *exp1;    
    exp1 = (experiment_total *) exp1_;  
    int adaptation_slave_buffer_send;
    
    adaptation_slave_buffer_send = *code; 
    MPI_Isend(&adaptation_slave_buffer_send, 1, MPI_INT, 0, 155, 
            MPI_COMM_WORLD, exp1[0].execution.adaptation_slave_send );    
}


void checkadaptsettings_(void *exp1_, int *size_dim, int *ncounter, double *balance, int *adaptflag, int *accept) {
    experiment_total *exp1;    
    int flag;
    exp1 = (experiment_total *) exp1_;  
    
    *adaptflag = 0;
    MPI_Test(exp1[0].execution.adaptation_slave_recv, &flag, MPI_STATUS_IGNORE);
    
    if (flag == 1) {
        *size_dim = (int)       exp1[0].execution.adaptation_slave_buffer_recv[0];  
        *ncounter = (int)       exp1[0].execution.adaptation_slave_buffer_recv[1];  
        *balance  = (double)    exp1[0].execution.adaptation_slave_buffer_recv[2]; 
        *accept = (int) exp1[0].execution.adaptation_slave_buffer_recv[3]; 
        *adaptflag = 1;
        
                
        MPI_Irecv(exp1[0].execution.adaptation_slave_buffer_recv, 4, MPI_DOUBLE,0, 155,
                 MPI_COMM_WORLD, exp1[0].execution.adaptation_slave_recv );  
    }
    
}









void destroyasynchinitmasterandwindows_(void *exp1_) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;        
    int i;

    // MASTER
    free(exp1[0].execution.sendbufferslave);
    exp1[0].execution.sendbufferslave=NULL;
    free(exp1[0].execution.receptionbufferslave);
    exp1[0].execution.receptionbufferslave=NULL;   
    
    if (exp1[0].execution.idp == 0) {
        for (i = 0; i < (exp1[0].execution.NPROC-1); i++) {
            free(exp1[0].execution.receptionbuffermaster[i] );
            free(exp1[0].execution.sendbuffermaster[i] );
            free(exp1[0].execution.adaptation_master_buffer_send[i] );
        }
        free(exp1[0].execution.receptionrequestmaster );
        free(exp1[0].execution.sendrequestmaster );
        free(exp1[0].execution.adaptation_master_send );
        free(exp1[0].execution.adaptation_master_recv );      
    }
    free( exp1[0].execution.adaptation_slave_send );    
    free( exp1[0].execution.adaptation_slave_recv );          
    free( exp1[0].execution.receptionrequestslave );
    free( exp1[0].execution.sendrequestslave );

}

void chargecooperativeparametersfortran_(void *exp1_, int *NP, int *tam, int *idp, double *ftarget, long *maxfunevals ) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp1_;
    
    exp1[0].execution.num_it = exp1[0].par_st->migration_freq_ite;
    exp1[0].execution.max_time_ite = exp1[0].par_st->max_time_ite;
    exp1[0].execution.NPROC = exp1[0].par_st->NPROC;
    exp1[0].execution.st_sent = 0;
    exp1[0].execution.stuckcond_lowVar = (int) 30 * exp1[0].execution.NPROC ;
    exp1[0].execution.minvarcondition = 1e-10 / exp1[0].execution.NPROC ;
    exp1[0].execution.stuckcount = 0;
    exp1[0].execution.migra_asin_wait = 0;
    exp1[0].execution.contadorMigra = 0;
    exp1[0].execution.enterMigrat = 0;  
    exp1[0].execution.migration = 0;
    exp1[0].execution.ftarget = *ftarget;
    exp1[0].execution.maxfunevals = *maxfunevals;
    exp1[0].execution.max_time_ite_last = 0.0;
    exp1[0].execution.max_eval_ite_last = 0.;
    *idp = exp1[0].execution.idp;
    exp1[0].execution.tam = *NP;
    *tam = exp1[0].execution.tam;
    exp1[0].execution.NM = exp1[0].par_st->migration_size;
    exp1[0].execution.NP = *NP;

}

void createcooperativetopologyess_(void *exp){
    
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    createtopology_(exp1);
    MPI_Comm_rank(exp1[0].execution.topology.comunicator, &(exp1[0].execution.idp));
    
    
}



void initcooperativestoppingcriteriaessmaster_( void *exp) {
    int i;
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;

    exp1[0].execution.st = (int *) malloc((exp1[0].execution.NPROC-1) * sizeof (int));
    exp1[0].execution.request_recept = (MPI_Request *) malloc((exp1[0].execution.NPROC-1) * sizeof (MPI_Request));

    for (i = 0; i < exp1[0].execution.NPROC-1; i++) {
        exp1[0].execution.st[i] = 0;
        MPI_Irecv(&(exp1[0].execution.st[i]), 1, MPI_INT,i+1, 200, MPI_COMM_WORLD, &(exp1[0].execution.request_recept[i]));
    }
}


void initcooperativestoppingcriteriaessslave_( void *exp) {
    int i;
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;

    exp1[0].execution.st = (int *) malloc(sizeof (int));
    exp1[0].execution.request_recept = (MPI_Request *) malloc( sizeof (MPI_Request));

    exp1[0].execution.st[0] = 0;
    MPI_Irecv(exp1[0].execution.st, 1, MPI_INT, 0, 200, MPI_COMM_WORLD, &exp1[0].execution.request_recept[0] );
}


void initcooperativestoppingcriteriaess_( void *exp) {
    int i;
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;

    exp1[0].execution.st = (int *) malloc(exp1[0].execution.NPROC * sizeof (int));
    exp1[0].execution.request_recept = (MPI_Request *) malloc(exp1[0].execution.NPROC * sizeof (MPI_Request));

    for (i = 0; i < exp1[0].execution.NPROC; i++) {
        exp1[0].execution.st[i] = 0;
        MPI_Irecv(&(exp1[0].execution.st[i]), 1, MPI_INT,
                i, 200, MPI_COMM_WORLD, &(exp1[0].execution.request_recept[i]));
    }
}

double * initsendbufferess_(void *exp, int *D) {
    experiment_total *exp1;    
    exp1 = (experiment_total *) exp;
    double *sendB;
    
    
    exp1[0].execution.size_send_buffer = 100;
    exp1[0].execution.send_id = 0;
    sendB = (double *) malloc(exp1[0].execution.size_send_buffer * exp1[0].execution.NM * (*D+1) * sizeof(double) );
    
    return sendB;
}

double * returnssendbufferess_(void *exp,  int *D, double *sendB) {
    experiment_total *exp1;    
    int id;
    exp1 = (experiment_total *) exp;

    id = exp1[0].execution.send_id;
    exp1[0].execution.send_id = exp1[0].execution.send_id + 1;
    if (exp1[0].execution.send_id > exp1[0].execution.size_send_buffer) {
        
        exp1[0].execution.send_id = 0;
    }
    
    return &sendB[id*exp1[0].execution.NM*(*D+1)];
}

void destroysendbufferess_(void *vect) {
    if (vect != NULL) {
        free(vect);
        vect = NULL;
    }
}

int checkcooperativemigrationcriteriaess_(void *exp) {
    experiment_total *exp1;
    exp1 = (experiment_total *) exp;
    
        if ( exp1[0].execution.NPROC != 1 ) {
            
            if (exp1[0].execution.num_it != -1) {
                exp1[0].execution.contadorMigra= exp1[0].execution.contadorMigra+1;
                if ( exp1[0].execution.contadorMigra >=  exp1[0].execution.num_it) {
                    exp1[0].execution.enterMigrat = 1;
                    exp1[0].execution.contadorMigra = 0;
                    return 1;
                } else return 0;
            } else return 0;
                    
        } else return 0;
    
}

int checkcooperativemigrationcriteriacessinner_(void *exp) {
    experiment_total *exp1;
    double time, currenttime, inittime;
    int return_mig;
    
    exp1 = (experiment_total *) exp;
    exp1[0].execution.enterMigrat = 0;
    return_mig = 0;
        if ( exp1[0].execution.NPROC != 1 ) {
            if (exp1[0].execution.max_time_ite != -1.0 ) {
                    inittime= returninittime_(exp1);
                    currenttime = calctimempi_(exp1, &inittime);
                    time = currenttime - exp1[0].execution.max_time_ite_last;
                    if ( time >=  exp1[0].execution.max_time_ite) {
                        return_mig = 1;
                    } else {
                        return_mig = 0;
                    }      
            }  else return_mig = 0;   
        } else return_mig = 0;
    
    return return_mig;
}

int checkcooperativemigrationcriteriacess_(void *exp, double *currenttime) {
    experiment_total *exp1;
    double time;
    
    exp1 = (experiment_total *) exp;
    exp1[0].execution.enterMigrat = 0;
        if ( exp1[0].execution.NPROC != 1 ) {

            if (exp1[0].execution.max_time_ite > 0.0 ) {
                if (exp1[0].execution.idp > 0) {
                    time = *currenttime - exp1[0].execution.max_time_ite_last;
                    if ( time >=  exp1[0].execution.max_time_ite) {
                        exp1[0].execution.enterMigrat = 1;
                    }

                } else {
                    exp1[0].execution.enterMigrat = 1;
                }               
            }     
        }
    
    return exp1[0].execution.enterMigrat;
}

void incrementmaxtimemig_(void *exp) {
    experiment_total *exp1;
    double inittime, currenttime;
    
    exp1 = (experiment_total *) exp;
    
    inittime= returninittime_(exp1);
       
    exp1[0].execution.max_time_ite_last = MPI_Wtime() - inittime;
    
}





void asynchronousstoppingcriteriawithrestartess_(void *exp, int *stop, long *evaluation_local, double *best,double *matrix,int *D) {
    experiment_total *exp1;
    int i;
    MPI_Request mpir;
    exp1 = (experiment_total *) exp;
    
               
        if (*stop < 1) {
            if ( *evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
                *stop = 2;
            } else if ((*best < exp1[0].execution.ftarget)
                    || (exp1[0].execution.stuckcount > exp1[0].execution.stuckcond_lowVar && 
                    (sumvar(matrix, exp1[0].execution.tam, *D) / *D) < exp1[0].execution.minvarcondition)
                    ) {
                if (*stop < 1) {
                    *stop = 1;
                    for (i = 0; i < exp1[0].execution.NPROC; i++) {
                        if (exp1[0].execution.idp != i) {
                            exp1[0].execution.st_sent = 1;
                            MPI_Isend(&(exp1[0].execution.st_sent), 1, MPI_INT, i, 200, MPI_COMM_WORLD, &mpir);
                        }
                    }
                }

            }
        }
}

void asynchronousstoppingcriteriaess_(void *exp, int *stop, long *evaluation_local, double *best, double *cputime, double *maxtime, int *stopOptimization) {
    experiment_total *exp1;
    int i;
    MPI_Request mpir;
    exp1 = (experiment_total *) exp;
    
    if (*stop < 1) {

        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 1; 
        }
        else if (*cputime >= *maxtime) {
            *stop = 2; 
        } else if (*stopOptimization == 1) {
            *stop = 4; 
        } else if (*best < exp1[0].execution.ftarget) {
            *stop = 3;
        }
    } 
    if (*stop > 0) {
                for (i = 0; i < exp1[0].execution.NPROC; i++) {
                    if (exp1[0].execution.idp != i) {
                        exp1[0].execution.st_sent = 1;
                        MPI_Isend(&(exp1[0].execution.st_sent), 1, MPI_INT, i, 200, MPI_COMM_WORLD, &mpir);
                    }
                }
    }
        
}


void asynchronousstoppingcriteriaessmaster_(void *exp, int *stop, long *evaluation_local, double *best, double *cputime, double *maxtime, int *stopOptimization) {
    experiment_total *exp1;
    int i;
    MPI_Request mpir;
    exp1 = (experiment_total *) exp;


    if (*stop < 1) {

        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 1;
        }
        else if (*cputime >= *maxtime) {
            *stop = 2;
        } else if (*stopOptimization == 1) {
            *stop = 4;
        } else if (*best < exp1[0].execution.ftarget) {
            *stop = 3;
        }
     }
        if ((*stop == 4) || (*stop == 3)) {
		exp1[0].execution.st_sent = 1;
                for (i = 0; i < (exp1[0].execution.NPROC-1); i++) {
                        MPI_Isend(&(exp1[0].execution.st_sent), 1, MPI_INT, i+1, 200, MPI_COMM_WORLD, &mpir);
                }
        }

}


void asynchronousstoppingcriteriaessslave_(void *exp, int *stop, long *evaluation_local, double *best, double *cputime, double *maxtime, int *stopOptimization) {
    experiment_total *exp1;
    int i;
    MPI_Request mpir;
    exp1 = (experiment_total *) exp;


    if (*stop < 1) {

        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 1;
        }
        else if (*cputime >= *maxtime) {
            *stop = 2;
        } else if (*stopOptimization == 1) {
            *stop = 4;
        } else if (*best < exp1[0].execution.ftarget) {
            *stop = 3;
        }
    }

    if ((*stop == 4) || (*stop == 3) || (*stop == 2) || (*stop == 1) ) {
            exp1[0].execution.st_sent = 1;
            MPI_Isend(&(exp1[0].execution.st_sent), 1, MPI_INT, 0, 200, MPI_COMM_WORLD, &mpir);
    }

}


void synchronousstoppingcriteriawithrestartess_(void *exp, int *stop, long *evaluation_local, double *best, double *matrix, int *D) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp;


    if (*stop < 1) {
        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 2;
        } else if ((*best < exp1[0].execution.ftarget)
                || (exp1[0].execution.stuckcount > exp1[0].execution.stuckcond_lowVar &&
                (sumvar(matrix, exp1[0].execution.tam, *D) / *D) < exp1[0].execution.minvarcondition)
                ) {
            if (*stop < 1) {
                *stop = 1;
                for (i = 0; i < exp1[0].execution.NPROC; i++) {
                    if (exp1[0].execution.idp != i) {
                        exp1[0].execution.st_sent = 1;
                        MPI_Send(&(exp1[0].execution.st_sent), 1, MPI_INT, i, 200, MPI_COMM_WORLD);
                    }
                }
            }

        }
    }

}

void synchronousstoppingcriteriacess_(void *exp, int *stop, long *evaluation_local, double *best, double *cputime, double *maxtime, int *stopOptimization) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp;
    

    if (*stop < 1) {
        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 1; 
        }
        else if (*cputime >= *maxtime) {
            *stop = 2; 
        } else if (*stopOptimization == 1) {
            *stop = 4; 
        } else if (*best < exp1[0].execution.ftarget) {
            *stop = 3;
        }
     
        if ((*stop == 4) || (*stop == 3)) {
                for (i = 0; i < exp1[0].execution.NPROC; i++) {
                        exp1[0].execution.st_sent = 1;
                        MPI_Send(&(exp1[0].execution.st_sent), 1, MPI_INT, i, 200, MPI_COMM_WORLD);
                }
        }
        
    }
}

void synchronousstoppingcriteriacessmaster_(void *exp, int *stop, long *evaluation_local, double *best, double *cputime, double *maxtime, int *stopOptimization) {
    experiment_total *exp1;
    int i;
    exp1 = (experiment_total *) exp;
    

    if (*stop < 1) {
        if (*evaluation_local * exp1[0].execution.NPROC >= exp1[0].execution.maxfunevals) {
            *stop = 1; 
        }
        else if (*cputime >= *maxtime) {
            *stop = 2; 
        } else if (*stopOptimization == 1) {
            *stop = 4; 
        } else if (*best < exp1[0].execution.ftarget) {
            *stop = 3;
        }
        
        if (exp1[0].execution.idp != 0) {
            if (*stop > 0) {
                exp1[0].execution.st_sent = 1;
                MPI_Send(&(exp1[0].execution.st_sent), 1, MPI_INT, 0, 200, MPI_COMM_WORLD);
            }
        }
        
    }
}

void setcountstopsvaressmaster_(void *exp, int *fin, int *exitcounter) {
    experiment_total *exp1;
    int i, counter, flag;
    exp1 = (experiment_total *) exp;
    
    
    if (exp1[0].execution.NPROC != 1) {
        MPI_Reduce(fin,exitcounter,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
        MPI_Bcast(exitcounter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        if (*fin > 0)
            *exitcounter = 1;
    }
}

int cooperativempitestess_(void *exp, int *l) {
    int flag;
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    flag = 0;
      
    MPI_Test(&(exp1[0].execution.request_recept[*l]), &flag, MPI_STATUS_IGNORE);
    return flag;
}





void gatherresultsess_(void *exp, double *matrixlocal, double *matrix, int *D, double *starttime, double *time_total, long *evaluation_local, long *eval_total) {
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;  
    double endtime;
    double localtotaltime;
    int l,flag;

    for (l = 0; l < exp1[0].execution.NPROC; l++) {
        MPI_Test( &(exp1[0].execution.request_recept[l]), &flag, MPI_STATUS_IGNORE);
        if (flag == 0) {
            MPI_Cancel( &(exp1[0].execution.request_recept[l]));
        }
    }    
    
    if (exp1[0].execution.request_recept != NULL) {
        free(exp1[0].execution.request_recept);
        exp1[0].execution.request_recept = NULL;
    }
    
    
    MPI_Gather(matrixlocal, exp1[0].execution.tam * (*D + 1), MPI_DOUBLE, matrix, exp1[0].execution.tam * (*D + 1), MPI_DOUBLE, 0, exp1[0].execution.topology.comunicator);
    endtime = MPI_Wtime();
    localtotaltime = endtime - *starttime;
    MPI_Reduce(&localtotaltime, time_total, 1, MPI_DOUBLE, MPI_MAX, 0, exp1[0].execution.topology.comunicator);
    MPI_Reduce(evaluation_local, eval_total, 1, MPI_LONG, MPI_SUM, 0, exp1[0].execution.topology.comunicator);
    MPI_Bcast(time_total,1,MPI_DOUBLE,0,exp1[0].execution.topology.comunicator);
    MPI_Bcast(eval_total,1,MPI_LONG, 0,exp1[0].execution.topology.comunicator);
}


void gatherresultsserializeess_(void *exp, void *local_s, double *smatrixlocal, double *smatrix, int *sizetotal, double *time1, double *time_total, 
        long *evaluation_local, long *eval_total) {
    experiment_total *exp1;
    local_solver *local1;
    int totallocal;
    int l,flag;
    exp1 = (experiment_total *) exp;  
    local1 = (local_solver *) local_s;
    totallocal = 0;
    
    for (l = 0; l < exp1[0].execution.NPROC; l++) {
        MPI_Test( &(exp1[0].execution.request_recept[l]), &flag, MPI_STATUS_IGNORE);
        if (flag == 0) {
            MPI_Cancel( &(exp1[0].execution.request_recept[l]));            
        }
    }    
    
    if (exp1[0].execution.request_recept != NULL) {
        free(exp1[0].execution.request_recept);
        exp1[0].execution.request_recept = NULL;
    }
    
    
    MPI_Gather(smatrixlocal, *sizetotal, MPI_DOUBLE, smatrix, *sizetotal, MPI_DOUBLE, 0, exp1[0].execution.topology.comunicator);
    MPI_Bcast(smatrix, *sizetotal*exp1->execution.NPROC, MPI_DOUBLE, 0,exp1[0].execution.topology.comunicator);
    MPI_Reduce(time1, time_total, 1, MPI_DOUBLE, MPI_MAX, 0, exp1[0].execution.topology.comunicator);
    MPI_Reduce(evaluation_local, eval_total, 1, MPI_LONG, MPI_SUM, 0, exp1[0].execution.topology.comunicator);
    MPI_Reduce(&(local1->counter), &totallocal, 1, MPI_INT, MPI_SUM, 0, exp1[0].execution.topology.comunicator);
    local1->counter = totallocal;
    MPI_Bcast(time_total,1,MPI_DOUBLE,0,exp1[0].execution.topology.comunicator);
    MPI_Bcast(eval_total,1,MPI_LONG, 0,exp1[0].execution.topology.comunicator);
    MPI_Bcast(&(local1->counter),1,MPI_INT, 0,exp1[0].execution.topology.comunicator);

}




void setexpexecutionstuckcountess(void *exp, int *valor) {
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    
    exp1[0].execution.stuckcount = *valor;
}


int getexpexecutionstuckcountess(void *exp) {
    experiment_total *exp1;
    
    exp1 = (experiment_total *) exp;
    
    return exp1[0].execution.stuckcount;
}


void migrationsynchcooperativesentcess_(void *exp,int *sizeD, double *matrix, double *concatmatrix, int *dim_size, int *value) {
    experiment_total *exp1;
    int t,j,contador,ssi;
    double *matrix_local;
    MPI_Status status2;

    exp1 = (experiment_total *) exp;

    if (exp1[0].execution.NPROC != 1) {
        ssi = *sizeD;
        if (exp1[0].execution.idp != 0) {
            MPI_Send(matrix, ssi*dim_size[exp1[0].execution.idp], MPI_DOUBLE, 0, 250, exp1[0].execution.topology.comunicator);
        }

        if (exp1[0].execution.idp == 0) {
            contador = 0;
            for (t = 1; t < exp1[0].execution.NPROC; t++) {
                matrix_local = (double *) malloc( ssi*dim_size[t] * sizeof(double));

                MPI_Recv(matrix_local, ssi*dim_size[t], MPI_DOUBLE, t, 250 , exp1[0].execution.topology.comunicator, &status2);
                for (j=0;j<dim_size[t];j++) {
                        memmove (&concatmatrix[ssi*j+contador], &matrix_local[j*ssi], ssi* sizeof(double));
                }


                contador=contador + ssi*(dim_size[t]);
                free(matrix_local);
                matrix_local = NULL;

            }
        }

        MPI_Bcast(concatmatrix, *value*ssi, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }
}


#endif
