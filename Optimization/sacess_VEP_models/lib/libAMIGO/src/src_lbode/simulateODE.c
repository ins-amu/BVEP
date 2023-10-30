/*$Id: simulateODE.c 1372 2012-01-24 11:16:50Z cokelaer $*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Sundials Header Files */

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include "mex.h"

#include "CNOStructure.h"

#define Ith(v,i) ( NV_DATA_S(v)[i] )

#define ZERO  RCONST(0.0)

typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data);

static int check_flag(void *flagvalue, char *funcname, int opt,int verbose);

int simulateODE
(
		CNOStructure* data,		int exp_num, 			int verbose,
		double reltol,			double atol,			double maxStepSize,
		int maxNumSteps,		int maxErrTestFails,    int sensi
)
{
	int i,j,k,neq,counter,flag;
    int is;
    booleantype err_con;
	realtype tout, ti, tf;
    
    //int plist[(*data).nPars];
	int* plist=(int*)malloc(sizeof(int)*data->nPars);
    
	N_Vector y;
    N_Vector *yS;
    
    int NS;
    
	void *cvode_mem;
	
    
	cvode_mem = NULL;
	y = NULL;
    

	neq=(*data).nStates;

    y = N_VNew_Serial(neq);
    
    NS=(*data).nPars;
    err_con=FALSE;
     
	if (check_flag((void *)y, "N_VNew_Serial", 0,verbose))
    {
		if(verbose)printf("\nSolver failed in N_VNew_Serial(neq). . .\n");
		return(0);
	}

    /* Initialize y */
	for(i=0; i<(*data).nRows; i++)
	{
		(*data).state_array[i] =(*data).unknown_ICs;
		(*data).inhibitor_array[i]=0;
	}

	for(i=0; i<(*data).nRows; ++i)
	{
		if((*data).isState[i])
		{
			for (j = 0; j < (*data).nSignals; j++)
			{
			/*	The passed indexes are from 1 to N intead from 0 to N-1*/
				if((*data).indexSignals[j]==i+1)
				{
					(*data).state_array[i] = (*data).valueSignals[exp_num][j];
				}
			}

			for (j = 0; j < (*data).nInhibitors; j++)
			{
				/*	The passed indexes are from 1 to N intead from 0 to N-1*/
				if((*data).indexInhibitors[j]==i+1)
				{
					(*data).inhibitor_array[i] = (*data).valueInhibitors[exp_num][j];
				}
			}
		}
		else
		{
			for (j = 0; j < (*data).nStimuli; j++)
			{
				/*	The passed indexes are from 1 to N intead from 0 to N-1*/
				if((*data).indexStimuli[j]==i+1)
				{
					(*data).state_array[i] = (*data).valueStimuli[exp_num][j];
				}
			}
		}
	}

	counter=0;
	for(i=0; i<(*data).nRows; i++)
	{
		if((*data).isState[i])
		{
			Ith(y,(*data).state_index[i]) =
					(*data).state_array[i];
			(*data).sim_results[exp_num][0][(*data).state_index[i]]=
					(*data).state_array[i];
			/*
			if(verbose)
			{
				printf("species number initial value %d",i);
				printf("\t%f\n",Ith(y,(*data).state_index[i]));
			}
			*/
		}
	}

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if(check_flag((void *)cvode_mem, "CVodeCreate", 0,verbose)) {
		if(verbose)printf("\nSolver failed in CVodeCreate(CV_BDF, CV_NEWTON) . . .\n");
		N_VDestroy_Serial(y);
		free(plist);
		return(0);
	}

	 ti=(*data).timeSignals[0];
	 tf=(*data).timeSignals[(*data).nTimes-1];
	  flag = CVodeInit(cvode_mem,*rhsODE, ti, y);
	 if (check_flag(&flag, "CVodeMalloc", 1,verbose))
	 {
		if(verbose)printf("\nSolver failed in CVodeInit(cvode_mem,*rhsODE, ti, y). . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		  free(plist);
	 	return(0);
	 }

	if(verbose)printf("Solver Memory Allocated\n");

	/* Set f_data */
	 flag = CVodeSetUserData(cvode_mem, data);
	 if(check_flag(&flag, "CVodeSetFdata", 1,verbose))
	 {
		 if(verbose)printf("\nSolver failed in flag = CVodeSetUserData(cvode_mem, data). . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		  free(plist);
		 return(0);
	 }

	 flag = CVodeSStolerances(cvode_mem,(realtype)reltol,(realtype)atol);
	  if(check_flag(&flag, "CVodeSStolerances", 1,verbose)) return(1);

      
     flag = CVodeSetErrFile(cvode_mem, stdout);
	  if(!verbose)
	  {
		  flag = CVodeSetErrFile(cvode_mem, NULL);
	  }

	  flag = CVDense(cvode_mem, neq);
	  if (check_flag(&flag, "CVDense", 1,verbose)) return(1);
	 if(verbose)printf("CVDENSE Solver Initiated\n");

	 /* Set maxnumsteps */
	 flag = CVodeSetMaxNumSteps(cvode_mem, maxNumSteps);
	 if(check_flag(&flag, "CVodeSetMaxNumSteps", 1,verbose))
	 {
		if(verbose)printf("\nSolver failed in CVodeSetMaxNumSteps(cvode_mem, maxNumSteps). . .\n");
		return(0);
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		  free(plist);
	 }
	 if(verbose)printf("Max number of steps: %i\n", maxNumSteps);

    CVodeSetMaxStep(cvode_mem,(realtype)maxStepSize);
    CVodeSetMaxErrTestFails(cvode_mem, maxErrTestFails);

    /*Experimental sensibility Analysis*/
    if (sensi) {
        
        for (i=0;i<(*data).nPars;i++) plist[i]=i; 
        
        yS = N_VCloneVectorArray_Serial(NS, y);
        if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0,0)) return(1);
        for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);
        
        flag = CVodeSensInit(cvode_mem, NS, CV_STAGGERED, NULL, yS);
        if(check_flag(&flag, "CVodeSensInit", 1,verbose)) return(1);
        flag = CVodeSensEEtolerances(cvode_mem);
        
        if(check_flag(&flag, "CVodeSensEEtolerances", 1,verbose)) return(1);
        flag = CVodeSetSensErrCon(cvode_mem, err_con);
        
        if (check_flag(&flag, "CVodeSetSensErrCon", 1,verbose)) return(1);
        flag = CVodeSetSensParams(cvode_mem,(*data).odeParameters, NULL, plist);
        
        if (check_flag(&flag, "CVodeSetSensParams", 1,verbose)) return(1);         
    } 
    
    for (i = 1; i < (*data).nTimes; ++i)
    {
    	tout=(*data).timeSignals[i];
    	flag = CVode(cvode_mem, tout, y, &tf, CV_NORMAL);
        if(verbose)printf("sim\n");
    	for (j = 0; j < (*data).nStates; j++)
    	{
    		(*data).sim_results[exp_num][i][j]= (double)Ith(y,j);
    		if(verbose)printf("%f\t",Ith(y,j));
    	}

    	if (check_flag(&flag, "CVode", 1,verbose))
    	{
    		if(verbose)fprintf(stderr,"\nSolver failed at flag = CVode(cvode_mem, tout, y, &tf, CV_NORMAL);. . .\n");
    		N_VDestroy_Serial(y);
    		CVodeFree(&cvode_mem);
			free(plist);
    		return(0);
    	}
    	if(verbose)printf("\n");
        
        if (sensi) {
            flag = CVodeGetSens(cvode_mem, &tf, yS);
            if (check_flag(&flag, "CVodeGetSens", 1,verbose)) break;
            for (j=0;j<(*data).nStates;j++){
                for (k=0;k<(*data).nPars;k++){
                    (*data).sensResults[exp_num][i][j][k]=NV_DATA_S(yS[k])[j];
                }
            }
        }
    }
    
    if(verbose)printf("\n");
    
    N_VDestroy_Serial(y);
    
    if (sensi) {
        N_VDestroyVectorArray_Serial(yS, NS);  /* Free yS vector */
    }
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    free(plist);
    return(1);
}

static int check_flag(void *flagvalue, char *funcname, int opt,int verbose)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		if(verbose)fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1); }

	/* Check if flag < 0 */
	else if (opt == 1)
	{
		errflag = (int *) flagvalue;
		if (*errflag < 0)
		{	if(verbose)fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",funcname, *errflag);
			return(1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL)
	{
		if(verbose)fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1);
	}

	return(0);
}




