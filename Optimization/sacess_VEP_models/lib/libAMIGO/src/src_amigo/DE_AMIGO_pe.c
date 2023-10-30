
#include <AMIGO_pe.h>
#include <de.h>

double DE_AMIGO_pe_objective(double *X,int D, double* Xl, double* Xu, void* data)
{
   register int i;
   double sum;
   AMIGO_problem* amigo_problem=(AMIGO_problem*) data;

/* Correction of boundary constraint violations, violating variable
   values are reflected back from the violated boundary        */

   for (i=0; i<D; i++)
   {
      while (X[i] < Xl[i] || X[i] > Xu[i])
      {
         if (X[i] < Xl[i]) X[i] = 2.0*Xl[i] - X[i];
         if (X[i] > Xu[i]) X[i] = 2.0*Xu[i] - X[i];
      }
   }

	set_AMIGO_problem_pars(X,amigo_problem);
	 
	sum=eval_AMIGO_problem_LSQ(amigo_problem);

	amigo_problem->nevals++;

	if(amigo_problem->verbose){

		if(sum<amigo_problem->temp_min){

			amigo_problem->temp_min=sum;
		
			#ifdef MATLAB
				mexPrintf("Eval %d f=%f\n",amigo_problem->nevals,sum);
			#else
				printf("Eval %d f=%f\n",amigo_problem->nevals,sum);
			#endif
		}
	}
	return(sum);
}


DE_opts* get_default_DE_opts(){
	
	DE_opts* de_opts=(DE_opts*)malloc(sizeof(DE_opts));
	
	de_opts->NP=100;
	de_opts->Gmax=300;

	return(de_opts);
}

int DE_AMIGO_pe(AMIGO_problem* amigo_problem){
	
	int i;
	int D=amigo_problem->nx;
	
	DE_opts* de_opts=get_default_DE_opts();
	
	de(D,de_opts,amigo_problem->LB,amigo_problem->UB,DE_AMIGO_pe_objective,amigo_problem);
	
	free(de_opts);
	
	return(0);
	
}

