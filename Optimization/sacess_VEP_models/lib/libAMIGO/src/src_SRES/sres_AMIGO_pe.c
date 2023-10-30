#include <stdio.h>
#include <stdlib.h>

#include "sharefunc.h"
#include "ESSRSort.h"
#include "ESES.h"
#include <AMIGO_pe.h>
#include <AMIGO_SRES_opts.h>

/*********************************************************************
** fitness(double *x, double *f, double *g)                        **
** change the function name as you want                            **
** x[dim], g[constraint]                                           **
*********************************************************************/
void SRES_fitness4amigo(double *, double *, double *,void* data);

/*********************************************************************
** double transform(double x)                                      **
** transform x: for coding                                         **
*********************************************************************/
double transform(double);

AMIGO_SRES_opts* get_default_sres_opts(){
	
	AMIGO_SRES_opts* amigo_sres_opts=(AMIGO_SRES_opts*) malloc(sizeof(AMIGO_SRES_opts));
	/*
	NP: 100
    itermax: 3000
         mu: 15
         pf: 0.4500
     varphi: 1
        var: 1
     vareta: 1
    cvarmax: 1.0000e-010*/
	return(amigo_sres_opts);
};

int sres_AMIGO_pe(AMIGO_problem* amigo_problem){

	int i;
	ESParameter *param;
	ESPopulation *population;
	ESStatistics *stats;
	ESfcnTrsfm *trsfm;
	unsigned int seed;
	int es;
	int constraint, dim;
	double *ub, *lb;
	int miu, lambda, gen;
	double gamma, alpha, varphi;
	int retry;
	double pf;


	/*********************************************************************
	** change the parameters here as you want                          **
	** random seed, gamma, alpha, varphi, retry number, pf,            **
	** if you dont know how to set, keep default settings              **
	**                                                                 **
	** constraint number, x dimension, miu:parent number,              **
	** lambda:offspring number, generation number,                     **
	** up and low bounds on x,                                         **
	*********************************************************************/

	seed = shareDefSeed;
	gamma = esDefGamma;
	alpha = esDefAlpha;
	varphi = esDefVarphi;
	retry = esDefRetry;
	pf = essrDefPf;
	es = esDefESSlash;
	dim = amigo_problem->nx;
	miu = 15;
	lambda = 100;
	gen = 500;

	constraint = 0;
	
	ub = NULL;
	lb = NULL;
	ub = ShareMallocM1d(dim);
	lb = ShareMallocM1d(dim);

	trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));

	/*********************************************************************
	** set the up and low bounds on x here                             **
	** lb[dim] and ub[dim]                                             **
	*********************************************************************/

	for(i=0;i<dim;i++){
		lb[i]=amigo_problem->LB[i];
		ub[i]=amigo_problem->UB[i];
	}

	/********************************************************************
	** end of user parameter setting                                  **
	** get started                                                    **
	********************************************************************/
	ESInitial(seed, &param, trsfm, SRES_fitness4amigo, es,
		constraint, dim, ub, lb,  miu, lambda, gen, gamma,
		alpha, varphi,  retry, &population, &stats,amigo_problem);

	while(stats->curgen < param->gen)
		ESStep(population, param, stats, pf,amigo_problem);

	//ESDeInitial(param, population, stats,amigo_problem);

	//Copy the results
	//amigo_problem->fbest=stats->bestindvdl->f;

	for(i=0;i<dim;i++){
		amigo_problem->xbest[i]=stats->bestindvdl->op[i];
	}

	ShareFreeM1d(ub);
	ub = NULL;
	ShareFreeM1d(lb);
	lb = NULL;
	ShareFreeM1c((char*)trsfm);
	trsfm = NULL;

	return 0;
}


/*********************************************************************
** set the fitness function here                                   **
** x[dim], g[constraint]                                           **
*********************************************************************/
void SRES_fitness4amigo(double *x, double *f, double *g,void* data){

	AMIGO_problem* amigo_problem=(AMIGO_problem*) data;
	set_AMIGO_problem_pars(x,amigo_problem);
	*f=eval_AMIGO_problem_LSQ(amigo_problem);
	amigo_problem->nevals++;
	if(amigo_problem->verbose){
		if(*f<amigo_problem->temp_min){
			amigo_problem->temp_min=*f;
		
			#ifdef MATLAB
				mexPrintf("Eval %d f=%f\n",amigo_problem->nevals,*f);
			#else
				printf("Eval %d f=%f\n",amigo_problem->nevals,*f);
			#endif
		}
	}
}


/*********************************************************************
** double transform(double x)                                      **
** transform x: for coding                                         **
*********************************************************************/
double transform(double x){

	double y;

	y = x;

	return y;
}
