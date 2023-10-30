



/*
Loosely inspired by fgeneric.m, the matlab version

%    Example: Optimize function f11 with MATLABS FMINUNC:
%
%       DIM = 5;
%       ftarget = fgeneric('initialize', 11, 1, 'testfminunc');
%       opt = optimset('TolFun', 1e-11);
%       X = fminunc('fgeneric', 8*rand(DIM,1) - 4, opt);
%       disp(fgeneric(X) - ftarget);
%       fgeneric('finalize');
%
%    This will create folder testfminunc. In this folder, the info file
%    bbobexp_f11.info will provide meta-information on the different
%    optimization runs obtained. Data of these runs are located in folder
%    testfminunc/data_f11 which will contain the file bbobexp_f11_DIM5.dat.
%
%fge
    fgeneric.m
    Author: Raymond Ros, Nikolaus Hansen, Steffen Finck, Marc Schoenauer
        (firstName.lastName@lri.fr)
    Version = 'Revision: $Revision: 646 $'
    Last Modified: $Date: 2009-02-18 13:54:40 +0100 (Wed, 18 Feb 2009) $
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>

/*int strcasecmp (const char *s1, const char *s2);*/
/*int strncasecmp (const char *s1, const char *s2, size_t n);*/

/* the specific includes for BBOB */
#include <bbobStructures.h>
#include "benchmarksdeclare.h"
#include "benchmarkshelper.h"
#include "benchmarks.h"
#include "benchmarksnoisy.h"


/* and the forward declarations of helper functions at the end of this file */
void writeNewIndexEntry(ParamStruct PARAMS);
void addIndexEntry(ParamStruct PARAMS);
void addDatIndexEntry(ParamStruct PARAMS);
void writeFinalData(ParamStruct PARAMS, LastEvalStruct BestFEval, double lastWriteEval, LastEvalStruct LastEval, double Fopt);
void writeBestF(char * dataFile, double BestFEval, double Fopt);
void sprintData(FILE* fout, double evals, double F, double bestF, double Fnoisy, double bestFnoisy, double * x, double Fopt);
void writeDataHeader(char * dataFile, double Fopt);
ParamStruct setNextDataFile(ParamStruct PARAMS, unsigned int isAllParamsMatching);

/* The following are some of the global static variables of the Matlab program */
  /*
  int DIM;
  unsigned int isInitDone;
  unsigned int trialid;
  double Fopt; */
  double fTrigger;
  double evalsTrigger;
  unsigned int idxEvalsTrigger, idxDIMEvalsTrigger;
  int idxFTrigger;

  double maxFunEvalsFactor = 1e6;
  unsigned int nbPtsEvals = 20;
  unsigned int nbPtsF = 5;
  int initDone = 0;

  bbobFunction actFunc = NULL;
  LastEvalStruct LastEval, BestFEval;
  double lastWriteEval;

/* The default values for the LastEvalStruct structures */
LastEvalStruct lastEvalInit = {0, DBL_MAX, DBL_MAX, DBL_MAX, {0}, 0};

/* If your compiler cries, you might need to comment the above, 
and uncomment the following
But if later you need to change DIM_MAX 
you will also need to modify this */

/*
LastEvalStruct lastEvalInit = {0, DBL_MAX, DBL_MAX, DBL_MAX, 
{
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
},
0};
*/

/* The default values for the ParamStruct structures */
/* UGLY cut-and-paste, don't forget to modify all if you modify one */
ParamStruct DefaultParam = {
    "not-specified", /* char algName[DefaultStringLength];  */
    "", /* char comments[LongDefaultStringLength]; */
    ".", /* char dataPath[DefaultStringLength]; */
    "bbobexp", /* char filePrefix[DefaultStringLength]; */
    0, /* unsigned int DIM; */
    1e-3, /* double precision; AKA deltaFTarget */
    9999, /* unsigned int funcId; set to 9999 to detect forgotten first init */
    9999, /* unsigned int itrial; set to 9999 to detect forgotten first init */
    "", /* char indexFilePrefix[DefaultStringLength]; */
    "", /* char indexFile[DefaultStringLength]; */
    "", /* char indexFileNameOnly[DefaultStringLength]; */
    "", /* char dataFile[DefaultStringLength]; */
    "", /* char dataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFile[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char rdataFile[DefaultStringLength]; */
    "", /* char rdataFileNameOnly[DefaultStringLength]; */
    0, /* brute force suffix (no need to extract it from file name!) */
    0 /* unsigned int runCounter; */
};

ParamStruct PreviousPARAMS = {
    "not-specified", /* char algName[DefaultStringLength];  */
    "", /* char comments[LongDefaultStringLength]; */
    ".", /* char dataPath[DefaultStringLength]; */
    "bbobexp", /* char filePrefix[DefaultStringLength]; */
    0, /* unsigned int DIM; */
    1e-3, /* double precision; AKA deltaFTarget */
    9999, /* unsigned int funcId; set to 9999 to detect forgotten first init */
    9999, /* unsigned int itrial; set to 9999 to detect forgotten first init */
    "", /* char indexFilePrefix[DefaultStringLength]; */
    "", /* char indexFile[DefaultStringLength]; */
    "", /* char indexFileNameOnly[DefaultStringLength]; */
    "", /* char dataFile[DefaultStringLength]; */
    "", /* char dataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFile[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char rdataFile[DefaultStringLength]; */
    "", /* char rdataFileNameOnly[DefaultStringLength]; */
    0, /* brute force suffix (no need to extract it from file name!) */
    0 /* unsigned int runCounter; */
};

ParamStruct CurrentPARAMS = {
    "not-specified", /* char algName[DefaultStringLength];  */
    "", /* char comments[LongDefaultStringLength]; */
    ".", /* char dataPath[DefaultStringLength]; */
    "bbobexp", /* char filePrefix[DefaultStringLength]; */
    0, /* unsigned int DIM; */
    1e-3, /* double precision; AKA deltaFTarget */
    9999, /* unsigned int funcId; set to 9999 to detect forgotten first init */
    9999, /* unsigned int itrial; set to 9999 to detect forgotten first init */
    "", /* char indexFilePrefix[DefaultStringLength]; */
    "", /* char indexFile[DefaultStringLength]; */
    "", /* char indexFileNameOnly[DefaultStringLength]; */
    "", /* char dataFile[DefaultStringLength]; */
    "", /* char dataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFile[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char hdataFileNameOnly[DefaultStringLength]; */
    "", /* char rdataFile[DefaultStringLength]; */
    "", /* char rdataFileNameOnly[DefaultStringLength]; */
    0, /* brute force suffix (no need to extract it from file name!) */
    0 /* unsigned int runCounter; */
};



/* so the user does not have to fill all params one by one */
ParamStruct fgeneric_getDefaultPARAMS(void)
{
   return DefaultParam; 
}

/*
   returns true or false depending whether the given FUNC_ID
    number is part of the testbed.
*/
int fgeneric_exist(unsigned int FUNC_ID)
{
    if ( (FUNC_ID <= handlesLength ) ||
         ( (100 < FUNC_ID) && (FUNC_ID <= 100+handlesNoisyLength) )
       )
        return 1;
    return 0;
}


/*
Memory aLlocations to be done only if DIM changes

initbenchmarkshelper : malloc for the variables used in both benchmark and benchamrknoisy
initbenchmarks and initbenchmarksnoisy : specific variables for benchmarks and benchmarksnoisy

the corresponding "free" are in
    finibenchmarks();
    finibenchmarksnoisy();
    finibenchmarkshelper();

For each benchmark function, the specific init are in the function itself (test on isInitDone)
*/


/* initialisation: 
   DIM, funcId, instanceId, and dataPath *must* be set in PARAMS 

@return FTARGET is the target function value of the specified fitness

For the format and folder/file structure of the output files see the
    technical documentation.
*/
double fgeneric_initialize(void *data)
{
    experiment_total *exp = (experiment_total *) data;
    ParamStruct PARAMS  = *(exp->param);
    char sLoc[1024];   /* general purpose buffer, mainly for sprintf */
    double * X;
    //TwoDoubles res;

    /* if nothing important has changed */
    if (!( (PARAMS.DIM == CurrentPARAMS.DIM) &&
         (PARAMS.funcId == CurrentPARAMS.funcId) &&
         (PARAMS.instanceId == CurrentPARAMS.instanceId))
       )
    {
        isInitDone = 0;
    }

    /* checks if no current run is still up */
    if (actFunc != NULL) {
        WARNING("Calling fgeneric_initialize while an experiment is still running (DIM %d, function %d, instance %d)\nCalling fgeneric_finalize to close it properly", CurrentPARAMS.DIM, CurrentPARAMS.funcId, CurrentPARAMS.instanceId);
        fgeneric_finalize(exp);
    }
    /* OK, here actFunc is NULL, and we can start a new experiment */
    /* First check the dimenstion, and (re)allocate memory if necessary */
    /************************************************************/
    if (PARAMS.DIM == 0)
        ERROR("You need to set the dimension of the problem greater than 0");

    if (PARAMS.DIM  > DIM_MAX)
        ERROR("You need to recompile the program to increase DIM_MAX if you want to run with dimension %d\nPlease also pay attention at the definition of lastEvalInit", PARAMS.DIM);

    /* Once the dimension is known, do the memory allocation if needed */
    if (CurrentPARAMS.DIM != PARAMS.DIM)
    {
        if (CurrentPARAMS.DIM != 0) /* i.e. the memory is already allocated, but not the right size */
        {
            finibenchmarks();
            finibenchmarksnoisy();
            finibenchmarkshelper();
        }
       /* sets the global variable (historical - and practical in loops?) */
        DIM = PARAMS.DIM;
        /* now allocate */
        initbenchmarkshelper(data);
        initbenchmarks();
        initbenchmarksnoisy();
    }

    /* tests the func ID and sets the global variable */
    /**************************************************/
    if (PARAMS.funcId-1 < handlesLength )
        actFunc = handles[PARAMS.funcId-1];
    else
    {
        if ( (100 < PARAMS.funcId) && (PARAMS.funcId-101 < handlesNoisyLength) )
            actFunc = handlesNoisy[PARAMS.funcId - 101];
        else
            ERROR("funcId in PARAMS is %d which is not a valid function identifier", PARAMS.funcId);
    }

    trialid = PARAMS.instanceId;

    /* the target value and book-keeping stuffs related to stopping criteria */
    /************************************************************************/
    Fopt = computeFopt(PARAMS.funcId, PARAMS.instanceId);

    LastEval = lastEvalInit; /* default values */
    BestFEval = LastEval;
    BestFEval.isWritten = 1;  /* ??? */

    lastWriteEval = 0;

    idxFTrigger = INT_MAX;
    fTrigger = DBL_MAX;   /* because 10^DBL_MAX will generate an error */
    idxEvalsTrigger = 0;
    evalsTrigger = floorl(powl(10.,(double)idxEvalsTrigger/(double)nbPtsEvals)); /* = 1 ! */
    idxDIMEvalsTrigger = 0;

    PARAMS.runCounter = 1;

    /* now the paths, file names, and different headers to write */
    /************************************************************/
    if ( strlen(PARAMS.dataPath) == 0 ) 
        ERROR("PARAMS.DATAPATH is expected to be a non-empty string. To set DATAPATH to the current working directory, input '.'");

    dirOK(PARAMS.dataPath); /* never returns if not OK :-) */

    /* Keeping both the complete file name with path, and the file name alone
       to avoid problems in function 
    */
    /*PARAMS.filePrefix is used as the prefix for both the data files and the index files.
    The problem is the data file prefix can change (we musn't have two different entries
    in the index that refer to the same data file).
    */
    sprintf(sLoc, "%s_f%d.info", PARAMS.filePrefix, PARAMS.funcId);
    createFullFileName(PARAMS.indexFile, PARAMS.dataPath, sLoc);

    sprintf(sLoc, "data_f%d", PARAMS.funcId);
    createFullFileName(PARAMS.dataFilePrefixNameOnly, sLoc, PARAMS.filePrefix);
    /* We prepend a folder name to the prefix: data_f%d is the folder name. So what we call
       xxxNameOnly actually contain a path.*/

    createFullFileName(sLoc, PARAMS.dataPath, sLoc);
    dirOK(sLoc); /* Create the directory data_fX if necessary. */

    createFullFileName(PARAMS.dataFilePrefix, PARAMS.dataPath, PARAMS.dataFilePrefixNameOnly);

    /* HERE: check that this PARAMS.dataFilePrefix is valid, ie. that we do not need to append
       to PARAMS.dataFilePrefix a -X where X is a number. We can write in the same file if 
       these attributes are the same as previously: algName, comments, indexFile, DIM, precision, funcId
     */
    if ((PARAMS.DIM == CurrentPARAMS.DIM) && (PARAMS.funcId == CurrentPARAMS.funcId) &&
        (PARAMS.precision == CurrentPARAMS.precision) && !strcmp(PARAMS.algName, CurrentPARAMS.algName) &&
        !strcmp(PARAMS.comments, CurrentPARAMS.comments)
       )
    {
        PARAMS.dataFileSuffix = CurrentPARAMS.dataFileSuffix;
        PARAMS = setNextDataFile(PARAMS, 1);
        if (!strcmp(PARAMS.dataFile, CurrentPARAMS.dataFile))
            addIndexEntry(PARAMS);
        else
            addDatIndexEntry(PARAMS);
    }
    else
    {
        PARAMS = setNextDataFile(PARAMS, 0);
        writeNewIndexEntry(PARAMS);
    }
    writeDataHeader(PARAMS.dataFile, Fopt);
    writeDataHeader(PARAMS.hdataFile, Fopt);
    writeDataHeader(PARAMS.rdataFile, Fopt);
    CurrentPARAMS = PARAMS;

    /* These lines are used to align the call to myrand with the ones in Matlab.
     * Don't forget to declare X (double * X;) and res (TwoDoubles res).*/
    X = (double*)malloc(DIM * sizeof(double));
    //res = (*actFunc)(X);
    free(X);
    return Fopt + PARAMS.precision; /* minimization */

}

void inicialize_functions2(void *data) {
    experiment_total *exp = (experiment_total *) data;  
    ParamStruct PARAMS = *(exp->param);
    isInitDone = 0;
    DIM = exp->test.bench.dim;

    initbenchmarkshelper(data);

    
    if (strcmp(exp->test.bench.type, "noiselessBBOB") == 0) {
        initbenchmarks();
        actFunc = handles[PARAMS.funcId-1];
    }  
    else {
        initbenchmarksnoisy();
        actFunc = handlesNoisy[PARAMS.funcId - 101];
    }

    trialid = PARAMS.instanceId;
    CurrentPARAMS = PARAMS;
    Fopt = computeFopt(PARAMS.funcId, PARAMS.instanceId);
  
}




void inicialize_functions(void *data) {
    experiment_total *exp = (experiment_total *) data  ;  
    ParamStruct PARAMS = *(exp->param);
    /* if nothing important has changed */

    if (!( (PARAMS.DIM == CurrentPARAMS.DIM) &&
         (PARAMS.funcId == CurrentPARAMS.funcId) &&
         (PARAMS.instanceId == CurrentPARAMS.instanceId))
       )
    {
        isInitDone = 0;
    }

    /* checks if no current run is still up */
    if (actFunc != NULL) {
        WARNING("Calling fgeneric_initialize while an experiment is still running (DIM %d, function %d, instance %d)\nCalling fgeneric_finalize to close it properly", CurrentPARAMS.DIM, CurrentPARAMS.funcId, CurrentPARAMS.instanceId);
        finalize_functions(data);
    }
    /* OK, here actFunc is NULL, and we can start a new experiment */
    /* First check the dimenstion, and (re)allocate memory if necessary */
    /************************************************************/
    if (PARAMS.DIM == 0)
        ERROR("You need to set the dimension of the problem greater than 0");

    if (PARAMS.DIM  > DIM_MAX)
        ERROR("You need to recompile the program to increase DIM_MAX if you want to run with dimension %d\nPlease also pay attention at the definition of lastEvalInit", PARAMS.DIM);    
    /* Once the dimension is known, do the memory allocation if needed */
    if (CurrentPARAMS.DIM != PARAMS.DIM)
    {
        if (CurrentPARAMS.DIM != 0) /* i.e. the memory is already allocated, but not the right size */
        {
            finibenchmarks();
            finibenchmarksnoisy();
            finibenchmarkshelper();
        }
        
       /* sets the global variable (historical - and practical in loops?) */
        DIM = PARAMS.DIM;
        /* now allocate */
        initbenchmarkshelper(data);
        initbenchmarks();
        initbenchmarksnoisy();
    }
        
    if (PARAMS.funcId-1 < handlesLength )
        actFunc = handles[PARAMS.funcId-1];
    else
    {
        if ( (100 < PARAMS.funcId) && (PARAMS.funcId-101 < handlesNoisyLength) )
            actFunc = handlesNoisy[PARAMS.funcId - 101];
        else
            ERROR("funcId in PARAMS is %d which is not a valid function identifier", PARAMS.funcId);
    }
    trialid = PARAMS.instanceId;
    CurrentPARAMS = PARAMS;
    Fopt = computeFopt(PARAMS.funcId, PARAMS.instanceId);
  
}

/* --------------------------------------------------------------------*/

/*  should be called to indicate the end
%    of a single run.  It writes data of the best-ever fitness value
%    and of the final function evaluation. It closes the data files.
%    RETURN the best true fitness value ever obtained.
*/
double fgeneric_finalize(void *data)
{
    //experiment_total *exp = (experiment_total *) data;
    //ParamStruct PARAMS;
    //PARAMS = *(exp->param);    
    if (actFunc == NULL)
        ERROR("Finalization process of fgeneric is called before the initialization process has occurred.");

    writeFinalData(CurrentPARAMS,BestFEval,lastWriteEval, LastEval, Fopt);

    CurrentPARAMS.runCounter = CurrentPARAMS.runCounter + 1;
    PreviousPARAMS = CurrentPARAMS;

    /* actFunc is the marker for finalized runs - see fgeneric_initialize */
    actFunc = NULL;
    return BestFEval.F;
}

void finalize_functions2(void *data) {
    experiment_total *exp = (experiment_total *) data;
    /* actFunc is the marker for finalized runs - see fgeneric_initialize */
    CurrentPARAMS.runCounter = CurrentPARAMS.runCounter + 1;
    PreviousPARAMS = CurrentPARAMS;

    /* actFunc is the marker for finalized runs - see fgeneric_initialize */
    actFunc = NULL;    
    if (strcmp(exp->test.bench.type, "noiselessBBOB") == 0) {
            finibenchmarks();
            finibenchmarkshelper();
    } 
    else {
            finibenchmarksnoisy();
            finibenchmarkshelper();            
    }
}

void finalize_functions(void *data) {
    /* actFunc is the marker for finalized runs - see fgeneric_initialize */
    CurrentPARAMS.runCounter = CurrentPARAMS.runCounter + 1;
    PreviousPARAMS = CurrentPARAMS;

    /* actFunc is the marker for finalized runs - see fgeneric_initialize */
    actFunc = NULL;    

    
}

/* returns the number of function
    evaluations since FGENERIC('initialize', ...) was called.
*/
double fgeneric_evaluations(void)
{
    if (actFunc == NULL) 
        WARNING("fgeneric_evaluations: fgeneric_initialized has not been called.");
    return  LastEval.num;
}


/*  returns the best function value obtained. 
*/
double fgeneric_best(void)
{
    if (actFunc == NULL) 
        WARNING("fgeneric_best: fgeneric_initialized has not been called.");
    return BestFEval.F;
}


/* returns the target function value.
*/
double fgeneric_ftarget(void)
{
    if (actFunc == NULL) {
        WARNING("fgeneric_ftarget: fgeneric_initialized has not been called.");
        return 0.0;
    }
    else
        return Fopt + CurrentPARAMS.precision;
}


/* returns the target function value.
*/
double fgeneric_ftarget_tol(double tolerance)
{
    if (actFunc == NULL) {
        WARNING("fgeneric_ftarget: fgeneric_initialized has not been called.");
        return 0.0;
    }
    else
        return Fopt + tolerance;
}


/*  returns the maximum number
     of function evaluations.
*/
double fgeneric_maxevals(unsigned int DIM)
{
    return (double)maxFunEvalsFactor * DIM;
}

/* Sets the seed of the noise (in benchmarkshelper) with a modulo 1e9 (constraint
 * on the seed provided to benchmarkshelper, this can happen when using time(NULL)).
 */
void fgeneric_noiseseed(unsigned long seed)
{
    seed = seed % 1000000000;
    setNoiseSeed(seed, seed);
}


/* Adds an output line to the restart-log .rdat. Call this if restarts occur within run_(your)_optimizer.
 * 'restart_reason' can be any string characterizing why the restart occured.
 */
void fgeneric_restart(char * restart_reason){
  FILE * rdataFileId;
  rdataFileId = bbobOpenFile(CurrentPARAMS.rdataFile);
  fprintf(rdataFileId, "%% restart: '%s'\n",  restart_reason);
  sprintData(rdataFileId, LastEval.num, LastEval.F, BestFEval.F, LastEval.Fnoisy, LastEval.bestFnoisy, LastEval.x, Fopt);
  fclose(rdataFileId);
}


/* -----------------------------------------------------------------
Now the computation functions
*/

/*   XX is an array of decision variable vectors of size howMany,
     STORED AS AN ARRAY OF DOUBLES
      and the corresponding fitness values are put in result
     This is because C cannot return a vector.
WARNING: the memory has been allowed before calling this function, no check is possible.
*/
void fgeneric_evaluate_vector(double * XX, unsigned int howMany, double * result)
{
    unsigned int i;
    double * XXtmp=XX;    /* moving vector of design variable */
    double * resultTmp = result ; /* moving result */
    /* a simple loop over the arrays */
    for (i=0; i<howMany; i++)
    {
        *resultTmp = fgeneric_evaluate(XXtmp);
        resultTmp++;
        XXtmp += DIM;
    }
}

/*   returns the fitness value of X,
%    X being the decision variable vector 
*/
double fgeneric_evaluate(double * X)
{
    int i;
    unsigned int boolImprovement = 0;
    double Fvalue, Ftrue;
    double evalsj; /*copia das evaluacióque se levan*/
    FILE * dataFileId;
    FILE * hdataFileId;
    TwoDoubles res;

    if (actFunc == NULL)
        ERROR("fgeneric has not been initialized. Please call 'fgeneric_initialize' first.");

    /* Substitue os valores de X na funcióbxectivo */
    res = (*actFunc)(X);
    Fvalue = res.Fval; /* Valor con ruido */
    Ftrue = res.Ftrue; /* valor sin ruido */

    /* should we print ? 2 possible conditions, # evals or fitness value */
    
    /* 
     ou ben: Si o nú de avaliacióque se levaba na úa evalació+1 e superor o umbral definido por evalTrigger
     ou ben: A diferencia entre fx sen ruido menos o resultado ómo da funcióopt éenor que o umbral Ftrigger
     *
     */
    if ( (LastEval.num+1 >= evalsTrigger) || (Ftrue-Fopt < fTrigger) ) 
    {
        evalsj = LastEval.num + 1;

        if (Fvalue < LastEval.bestFnoisy) /* minimization*/
            LastEval.bestFnoisy = Fvalue;

        /* Míse si o valor de Ftrue éellor que o valor o mellor valor gardado de fx en BestEval*/
        
        
        if (Ftrue < BestFEval.F) { /* minimization*/
            boolImprovement = 1; /* */
            BestFEval.F = Ftrue;
            BestFEval.bestFnoisy = LastEval.bestFnoisy;
            BestFEval.isWritten = 0;
        }

        /* CODIGO PARA IMPRIMIR NO .dat*/
        /* should we print something? First based on # evals */
        /* Se imprime solo si se chega o umbral de impresió*/
        if (evalsj >= evalsTrigger)
        {
            lastWriteEval = evalsj; /* se almacena o úo valor que foi impreso */
            BestFEval.isWritten = 1;
            dataFileId = bbobOpenFile(CurrentPARAMS.dataFile);
            /* IMPRESSION */
            sprintData(dataFileId, evalsj,Ftrue, BestFEval.F, Fvalue, LastEval.bestFnoisy, X,Fopt);
            fclose(dataFileId);
            /* update of next print triggers based on # evals */
            /* SE RECALCULA EVALSTRIGGERS */
            while (evalsj >= floorl(powl(10., (double)idxEvalsTrigger/(double)nbPtsEvals)) )
                idxEvalsTrigger = idxEvalsTrigger + 1;

            while ( evalsj >= DIM * powl(10., (double)idxDIMEvalsTrigger))
                idxDIMEvalsTrigger = idxDIMEvalsTrigger + 1;

            evalsTrigger = fmin(floorl(powl(10., (double)idxEvalsTrigger/(double)nbPtsEvals)), DIM * powl(10., (double) idxDIMEvalsTrigger));
        }

        /* CODIGO PARA IMPRIMIR NO .tdat*/
        /* now based on fitness values */
        if (Ftrue - Fopt < fTrigger)
        {
            hdataFileId = bbobOpenFile(CurrentPARAMS.hdataFile);
            sprintData(hdataFileId, evalsj,Ftrue,BestFEval.F, Fvalue, LastEval.bestFnoisy, X, Fopt);
            fclose(hdataFileId);

            if (Ftrue-Fopt <= 0)
                fTrigger = -DBL_MAX;
            else
            {
                if (idxFTrigger == INT_MAX)
                    idxFTrigger = ceill(log10(Ftrue-Fopt))*nbPtsF;

                while ( (Ftrue-Fopt) <= powl(10., (double)idxFTrigger/(double)nbPtsF) )
                    idxFTrigger = idxFTrigger - 1;

                fTrigger = fmin(fTrigger, powl(10., (double)idxFTrigger/(double)nbPtsF));
            }
        }

        
        
        /* Si o mellor valor non estáscrito e boolImprovement = 1 (a Ftrue éellor que a gardada en BEstEval)*/
        if ( ! BestFEval.isWritten && boolImprovement )
        {        
            BestFEval.num = LastEval.num+1; 
            BestFEval.Fnoisy = Fvalue;
            for (i=0; i<DIM; i++)
                BestFEval.x[i] = X[i];
        }
    }
    else
    {
        if (Ftrue < BestFEval.F)
        {
            BestFEval.num = LastEval.num+1;
            BestFEval.Fnoisy = Fvalue;
            for (i=0; i<DIM; i++)
                BestFEval.x[i] = X[i];
            BestFEval.F = Ftrue;
            BestFEval.bestFnoisy=fmin(LastEval.bestFnoisy, Fvalue);
            BestFEval.isWritten = 0;
        }
        LastEval.bestFnoisy = fmin(LastEval.bestFnoisy,Fvalue);
    } /* if (LastEval.num+POPSI >= evalsTrigger || bestFtrue-Fopt <= fTrigger) */

    
    LastEval.num = LastEval.num + 1;
    LastEval.F = Ftrue;
    LastEval.Fnoisy = Fvalue;
    for (i=0; i<DIM; i++)
        LastEval.x[i] = X[i];

    
    return Fvalue;
}




int write_in_file(TwoDoubles res, double * X){
    int i;
    int error=0;
    unsigned int boolImprovement = 0;
    double Fvalue, Ftrue;
    double evalsj;
    FILE * dataFileId;
    FILE * hdataFileId;
    

    if (actFunc == NULL)
        ERROR("fgeneric has not been initialized. Please call 'fgeneric_initialize' first.");

    Fvalue = res.Fval;
    Ftrue = res.Ftrue;
    //printf("%lf\n",Ftrue);
    /* should we print ? 2 possible conditions, # evals or fitness value */
    if ( (LastEval.num+1 >= evalsTrigger) || (Ftrue-Fopt < fTrigger) ) 
    {
        evalsj = LastEval.num + 1;
        if (Fvalue < LastEval.bestFnoisy) /* minimization*/
            LastEval.bestFnoisy = Fvalue;

        if (Ftrue < BestFEval.F) { /* minimization*/
            boolImprovement = 1;
            BestFEval.F = Ftrue;
            BestFEval.bestFnoisy = LastEval.bestFnoisy;
            BestFEval.isWritten = 0;
        }
    
    
    
    
        
        /* should we print something? First based on # evals */
        if (evalsj >= evalsTrigger)
        {
            lastWriteEval = evalsj;
            BestFEval.isWritten = 1;
            
            dataFileId = bbobOpenFile(CurrentPARAMS.dataFile);
            /* IMPRESSION */
            //if (dataFileId==NULL) printf("NULO\n");else printf("NO NULO\n");sleep(1);
            sprintData(dataFileId, evalsj,Ftrue, BestFEval.F, Fvalue, LastEval.bestFnoisy, X,Fopt);
            
            //sleep(1);
            
            fclose(dataFileId);
            /* update of next print triggers based on # evals */
            while (evalsj >= floorl(powl(10., (double)idxEvalsTrigger/(double)nbPtsEvals)) )
                idxEvalsTrigger = idxEvalsTrigger + 1;
            while ( evalsj >= DIM * powl(10., (double)idxDIMEvalsTrigger))
                idxDIMEvalsTrigger = idxDIMEvalsTrigger + 1;
            evalsTrigger = fmin(floorl(powl(10., (double)idxEvalsTrigger/(double)nbPtsEvals)), DIM * powl(10., (double) idxDIMEvalsTrigger));
        }
    
    
    
    
        /* now based on fitness values */
        if (Ftrue - Fopt < fTrigger)
        {
            // printf(" %lf - %lf = lf\n", Ftrue, Fopt, (Ftrue - Fopt));
           // printf("INIT TRIGEER %lf \n", fTrigger);
            hdataFileId = bbobOpenFile(CurrentPARAMS.hdataFile);
            sprintData(hdataFileId, evalsj,Ftrue,BestFEval.F, Fvalue, LastEval.bestFnoisy, X, Fopt);
            fclose(hdataFileId);

            if (Ftrue-Fopt <= 0)
                fTrigger = -DBL_MAX;
            else
            {
                if (idxFTrigger == INT_MAX)
                    idxFTrigger = ceill(log10(Ftrue-Fopt))*nbPtsF;

                while ( (Ftrue-Fopt) <= powl(10., (double)idxFTrigger/(double)nbPtsF) )
                    idxFTrigger = idxFTrigger - 1;

                fTrigger = fmin(fTrigger, powl(10., (double)idxFTrigger/(double)nbPtsF));
            }
            
            
            //printf("FIN TRIGEER %lf \n", fTrigger);sleep(1);
        }
        if ( ! BestFEval.isWritten && boolImprovement )
        {
            BestFEval.num = LastEval.num+1;
            BestFEval.Fnoisy = Fvalue;
            for (i=0; i<DIM; i++)
                BestFEval.x[i] = X[i];
        }
    }
    else
    {
        if (Ftrue < BestFEval.F)
        {
            BestFEval.num = LastEval.num+1;
            BestFEval.Fnoisy = Fvalue;
            for (i=0; i<DIM; i++)
                BestFEval.x[i] = X[i];
            BestFEval.F = Ftrue;
            BestFEval.bestFnoisy=fmin(LastEval.bestFnoisy, Fvalue);
            BestFEval.isWritten = 0;
        }
        LastEval.bestFnoisy = fmin(LastEval.bestFnoisy,Fvalue);
        
    } /* if (LastEval.num+POPSI >= evalsTrigger || bestFtrue-Fopt <= fTrigger) */
    
    LastEval.num = LastEval.num + 1;
    LastEval.F = Ftrue;
    LastEval.Fnoisy = Fvalue;
    for (i=0; i<DIM; i++)
        LastEval.x[i] = X[i];

    return error;    
}


/*   returns the fitness value of X,
%    X being the decision variable vector 
*/
 double fgeneric_evaluate_noise_without_writefile(double * X, void *data)
{
    TwoDoubles res;
    if (actFunc == NULL)
        ERROR("fgeneric has not been initialized. Please call 'fgeneric_initialize' first.");

    res = (*actFunc)(X);
    return res.Ftrue;
}


/*   returns the fitness value of X,
%    X being the decision variable vector 
*/
 double fgeneric_evaluate_noiseless_without_writefile( double * X, void *data)
{
    TwoDoubles res;
    
    if (actFunc == NULL)
        ERROR("fgeneric has not been initialized. Please call 'fgeneric_initialize' first.");

    res = (*actFunc)(X);
    
    return res.Fval ;
}



/* %--------------------------------------------------------------------------
%
%  Subfunctions
%
%--------------------------------------------------------------------------*/

/* write the comment line header in the data files */
void writeDataHeader(char * dataFile, double Fopt)
{
    FILE * dataFileId = bbobOpenFile(dataFile);

    fprintf(dataFileId, "%% function evaluation | noise-free fitness - Fopt (%13.12e) | best noise-free fitness - Fopt | measured fitness | best measured fitness | x1 | x2...\n", Fopt);
    fclose(dataFileId);
}

/*----------------------------------------------------------------- */


/* open the index file and write a new index entry */
void writeNewIndexEntry(ParamStruct PARAMS)
{
    FILE * indexFileId; /* historical name */
    int newline = 1;
    if ( !existFile(PARAMS.indexFile) ) /* will be created when opened */
        newline = 0;

    indexFileId = fopen(PARAMS.indexFile,"a");
    if (indexFileId == NULL)
        ERROR("Could not open %s", PARAMS.indexFile);

    if (newline == 1)
        fprintf(indexFileId,"\n");

    fprintf(indexFileId, "funcId = %d, DIM = %d, Precision = %.3e, algId = '%s'\n", PARAMS.funcId, PARAMS.DIM, PARAMS.precision, PARAMS.algName);
    fprintf(indexFileId,"%% %s\n%s, %d", PARAMS.comments, PARAMS.hdataFileNameOnly, PARAMS.instanceId);
    fclose(indexFileId);
}

/* --------------------------------------------------------------- */

/* open the index file and write a new index entry */
void addIndexEntry(ParamStruct PARAMS)
{
    FILE * indexFileId; /* historical name */

    if ( !existFile(PARAMS.indexFile) ) /* will be created when opened */
        ERROR("Could not find %s", PARAMS.indexFile);

    indexFileId = fopen(PARAMS.indexFile,"a");
    if (indexFileId == NULL)
        ERROR("Could not open %s", PARAMS.indexFile);

    fprintf(indexFileId,", %d", PARAMS.instanceId);
    fclose(indexFileId);
}

/* --------------------------------------------------------------- */

/* open the index file and write a new index entry */
void addDatIndexEntry(ParamStruct PARAMS)
{
    FILE * indexFileId; /* historical name */

    if ( !existFile(PARAMS.indexFile) ) /* will be created when opened */
        ERROR("Could not find %s", PARAMS.indexFile);

    indexFileId = fopen(PARAMS.indexFile,"a");
    if (indexFileId == NULL)
        ERROR("Could not open %s", PARAMS.indexFile);

    fprintf(indexFileId,", %s, %d", PARAMS.hdataFileNameOnly, PARAMS.instanceId);
    fclose(indexFileId);
}

/* --------------------------------------------------------------- */

/* complete the data file with unwritten information */
void writeFinalData(ParamStruct PARAMS, LastEvalStruct BestFEval, double lastWriteEval, LastEvalStruct LastEval, double Fopt)
{
    FILE * dataFileId; /* historical name */
    FILE * indexFileId; /* historical name */

    dataFileId = bbobOpenFile(PARAMS.dataFile);

    if ( ! BestFEval.isWritten )
    {
        if (BestFEval.num > lastWriteEval)
        {
            lastWriteEval = BestFEval.num;
            sprintData(dataFileId, BestFEval.num,BestFEval.F, BestFEval.F, BestFEval.Fnoisy, BestFEval.bestFnoisy, LastEval.x, Fopt);
        }
        else
        {   /* here, need to rewind dataFileId to write best at correct position */
            fclose(dataFileId);
            writeBestF(PARAMS.dataFile, BestFEval.num, Fopt);
            dataFileId = bbobOpenFile(PARAMS.dataFile);
        }
    }
    if (LastEval.num > lastWriteEval)
        sprintData(dataFileId, LastEval.num, LastEval.F, BestFEval.F, LastEval.Fnoisy, LastEval.bestFnoisy, LastEval.x, Fopt);

    fclose(dataFileId);

    /* now the index file */
    indexFileId = bbobOpenFile(PARAMS.indexFile);
    fprintf(indexFileId, ":%.0f|%.1e", LastEval.num, BestFEval.F - Fopt - PARAMS.precision);
    fclose(indexFileId);
}


/* ------------------------------------------------------------------ */
/* rewrite the data file with the information about the best F. */
void writeBestF(char * dataFile, double BestFEval, double Fopt)
{
    printf("Calling writeBestF with %s, %.0f, %g\n", dataFile, BestFEval, Fopt);
}

    /* %TODO: this function lacks efficiency.
%   disp(sprintf(['Rewriting %s to include evaluation of the best ' ...
%                'fitness value obtained at the function evaluation ' ...
%                'number %d.'],dataFile,BestFEval.num));
*/
  
/*
buffertmp = '';
  buffer = '';

  dataFileId = fopen(dataFile,'r');
  if dataFileId < 0
    warning('MATLAB:CouldNotOpen','Could not open %s: %s', dataFile,msg);
  else
    % Store the data of the last run in the variable buffer
   while feof(dataFileId) == 0
      line = fgets(dataFileId);

      buffertmp = [buffertmp line];

      if strcmp(line(1),'%')
        buffer = [buffer buffertmp];
        buffertmp = '';
      end
    end
    % At this point the last run (delimited by the header line
    % starting with a %) is stored in buffertmp.
    fclose(dataFileId);

    dataFileId = fopen(dataFile,'w');
    fprintf(dataFileId,'%s',buffer);

    % Split buffer into lines
    % lines = regexpl(buffertmp,'[^\n]*\n','match');
    lines = regexpl(buffertmp, ['[^' sprintf('\n') ']*' ...
                               sprintf('\n')], 'match');
    %This is proposed to make octave 2.9.14 work :S

    for j = 1:length(lines)
      % Browse through the lines to find the insertion point.
      if sscanf(lines{j},'%d',1) > BestFEval.num
        fprintf(dataFileId,'%s', ...
               sprintData(BestFEval.num,BestFEval.F,BestFEval.F, ...
                     BestFEval.Fnoisy, BestFEval.bestFnoisy,BestFEval.x,Fopt));
        break;
      end
      fprintf(dataFileId,'%s',lines{j});
    end
    fprintf(dataFileId,'%s',lines{j:end});

    fclose(dataFileId);
  end
}
*/

/* ------------------------------------------------------------------- */
/* write a formatted line into a data file */
void sprintData(FILE* fout, double evals, double F, double bestF, double Fnoisy, double bestFnoisy, double * x, double Fopt)
{
    int i;
    fprintf(fout, "%.0f", evals);
    fprintf(fout, " %+10.9e %+10.9e %+10.9e %+10.9e", F-Fopt, bestF-Fopt, Fnoisy, bestFnoisy);
    if (DIM < 22)
        for (i=0; i<DIM; i++)
            fprintf(fout, " %+5.4e",x[i]);
    fprintf(fout, "\n");
}

/* -------------------------------------------------------------------*/
/* return a data file prefix that does not exist yet  */
ParamStruct setNextDataFile(ParamStruct PARAMS, unsigned int isAllParamsMatching)
{

    if (PARAMS.dataFileSuffix == 0)
    {
        sprintf(PARAMS.dataFile, "%s_f%d_DIM%d.tdat", PARAMS.dataFilePrefix,
                PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.dataFileNameOnly, "%s_f%d_DIM%d.tdat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.hdataFile, "%s_f%d_DIM%d.dat", PARAMS.dataFilePrefix,
                PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.hdataFileNameOnly, "%s_f%d_DIM%d.dat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.rdataFile, "%s_f%d_DIM%d.rdat", PARAMS.dataFilePrefix,
                PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.rdataFileNameOnly, "%s_f%d_DIM%d.rdat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.funcId, PARAMS.DIM);
    }
    else
    {
        sprintf(PARAMS.dataFile, "%s-%02d_f%d_DIM%d.tdat", PARAMS.dataFilePrefix,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.dataFileNameOnly, "%s-%02d_f%d_DIM%d.tdat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.hdataFile, "%s-%02d_f%d_DIM%d.dat", PARAMS.dataFilePrefix,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.hdataFileNameOnly, "%s-%02d_f%d_DIM%d.dat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.rdataFile, "%s-%02d_f%d_DIM%d.rdat", PARAMS.dataFilePrefix,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        sprintf(PARAMS.rdataFileNameOnly, "%s-%02d_f%d_DIM%d.rdat", PARAMS.dataFilePrefixNameOnly,
                PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
    }

    if (!isAllParamsMatching)
    {
        PARAMS.dataFileSuffix = 0;
        while (existFile(PARAMS.dataFile) || existFile(PARAMS.hdataFile))
        {
            PARAMS.dataFileSuffix ++;
            sprintf(PARAMS.dataFile, "%s-%02d_f%d_DIM%d.tdat", PARAMS.dataFilePrefix,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
            sprintf(PARAMS.dataFileNameOnly, "%s-%02d_f%d_DIM%d.tdat", PARAMS.dataFilePrefixNameOnly,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
            sprintf(PARAMS.hdataFile, "%s-%02d_f%d_DIM%d.dat", PARAMS.dataFilePrefix,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
            sprintf(PARAMS.hdataFileNameOnly, "%s-%02d_f%d_DIM%d.dat", PARAMS.dataFilePrefixNameOnly,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
            sprintf(PARAMS.rdataFile, "%s-%02d_f%d_DIM%d.rdat", PARAMS.dataFilePrefix,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
            sprintf(PARAMS.rdataFileNameOnly, "%s-%02d_f%d_DIM%d.rdat", PARAMS.dataFilePrefixNameOnly,
                    PARAMS.dataFileSuffix, PARAMS.funcId, PARAMS.DIM);
        }
    }
    return PARAMS;
}


