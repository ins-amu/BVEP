#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include "benchmarksdeclare.h"
#include "bbobStructures.h"


double TOL;

static int seed = -1;
static int seedn = -1;

static double * gval;
static double * gval2;
static double * gvect;
static double * uniftmp;
static double * tmpvect;


 double round1( double a){
        return (double) floorl(a + 0.5);
}
 double fmin1( double a,  double b){
 return b < a ? b : a;
}
 double fmax1( double a,  double b){
 return b > a ? b : a;
}

/* set the seed for the noise.
 * If the seeds are larger than 1e9 they are set back to 1 in randn and myrand.
 */
void setNoiseSeed(unsigned int _seed, unsigned int _seedn)
{
    seed = _seed;
    seedn = _seedn;
}

void unif(double* r, int N, int inseed)
{
    /* generates N uniform numbers with starting seed*/
    int aktseed;
    int tmp;
    int *rgrand;
    rgrand = (int *) malloc(32*sizeof(int));
    //int rgrand[32];
    int aktrand;
    int i;

    if (inseed < 0)
        inseed = -inseed;
    if (inseed < 1)
        inseed = 1;
    aktseed = inseed;

    for (i = 39; i >= 0; i--)
    {
        tmp = (int)floorl((double)aktseed/(double)127773);
        aktseed = 16807  * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        if (i < 32)
           rgrand[i] = aktseed;
    }
    aktrand = rgrand[0];

    for (i = 0; i < N; i++)
    {
        tmp = (int)floorl((double)aktseed/(double)127773);
        aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        tmp = (int)floorl((double)aktrand / (double)67108865);
        aktrand = rgrand[tmp];
        rgrand[tmp] = aktseed;
        r[i] = (double)aktrand/2.147483647e9;
        if (r[i] == 0.)
        {
            //printf("Warning: zero sampled(?), set to 1e-99.\n");
            r[i] = 1e-99;
        }
    }
    free(rgrand);
    return;
}

void gauss(double * g, int N, int seed)
{
    /* samples N standard normally distributed numbers
       being the same for a given seed.*/
    int i;

    unif(uniftmp, 2*N, seed);

    for (i = 0; i < N; i++)
    {
        g[i] = sqrtl(-2*logl(uniftmp[i])) * cosl(2*M_PI*uniftmp[N+i]);
        if (g[i] == 0.)
            g[i] = 1e-99;
    }
    return;
}

void computeXopt(int seed, int _DIM) {
    int i;

    unif(tmpvect, _DIM, seed);
    for (i = 0; i < _DIM; i++)
    {
        Xopt[i] = 8 * floorl(1e4 * tmpvect[i])/1e4 - 4;
        if (Xopt[i] == 0.0)
            Xopt[i] = -1e-5;
    }
}

void monotoneTFosc(double* f) {
    double a = 0.1;
    int i;
    for (i = 0; i < DIM; i++)
    {
        if (f[i] > 0)
        {
            f[i] = logl(f[i])/a;
            f[i] = powl(expl(f[i] + 0.49*(sinl(f[i]) + sinl(0.79*f[i]))), a);
        }
        else if (f[i] < 0)
        {
            f[i] = logl(-f[i])/a;
            f[i] = -powl(expl(f[i] + 0.49*(sinl(0.55 * f[i]) + sinl(0.31*f[i]))), a);
        }
    }
}

void freeStarStar(double** M, int m)
{
    int i;
    for (i = 0; i < m; i++)
    {
        free(M[i]);
    }
    free(M);
    return;
}

double** reshape(double** B, double* vector, int m, int n)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = vector[j * m + i];
        }
    }
    return B;
}


void computeRotation(double ** B, int seed, int _DIM)
{
    double prod;
    int i, j, k; /*Loop over pairs of column vectors*/

    gauss(gvect, _DIM * _DIM, seed);
    reshape(B, gvect, _DIM, _DIM);
    /*1st coordinate is row, 2nd is column.*/

    for (i = 0; i < _DIM; i++)
    {
        for (j = 0; j < i; j++)
        {
            prod = 0;
            for (k = 0; k < _DIM; k++)
            {
                prod += B[k][i] * B[k][j];
            }
            for (k = 0; k < _DIM; k++)
            {
                B[k][i] -= prod * B[k][j];
            }
        }
        prod = 0;
        for (k = 0; k < _DIM; k++)
        {
            prod += B[k][i] * B[k][i];
        }
        for (k = 0; k < _DIM; k++)
        {
            B[k][i] /= sqrtl(prod);
        }
    }
}

double myrand(void) {
    /*Adaptation of myrand*/
    if (seed == -1)
        seed = time(NULL) % 1000000000; /* cannot be larger than 1e9 */

    seed ++;
    if (seed > 1e9)
        seed = 1;
    unif(uniftmp, 1, seed);
    return uniftmp[0];
}

double randn(void) {
    /*Adaptation of myrandn*/
    if (seedn == -1)
        seedn = time(NULL) % 1000000000; /* cannot be larger than 1e9 */

    seedn ++;
    if (seedn > 1e9)
        seedn = 1;
    gauss(uniftmp, 1, seedn);
    return uniftmp[0];
}

double FGauss(double Ftrue, double beta)
{
    double Fval = Ftrue * expl(beta * randn());
    Fval += 1.01 * TOL;
    
    
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

double FUniform(double Ftrue, double alpha, double beta)
{
    double Fval = powl(myrand(), beta) * Ftrue * fmax1(1., powl(1e9/(Ftrue+1e-99), alpha * myrand()));
    Fval += 1.01 * TOL;
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

double FCauchy(double Ftrue, double alpha, double p)
{
    double Fval;
    double tmp = randn()/fabsl(randn()+1e-199);
    /* tmp is so as to actually do the calls to randn in order for the number
     * of calls to be the same as in the Matlab code.
     */
    if (myrand() < p)
        Fval = Ftrue + alpha * fmax1(0., 1e3 + tmp);
    else
        Fval = Ftrue + alpha * 1e3;

    Fval += 1.01 * TOL;
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

int compare_doubles (const void *a, const void *b)
{
    double temp = peaks[*(const int*)a] - peaks[*(const int*)b];
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}

double computeFopt(int _funcId, int _trialId) {
    int rseed, rrseed;
    if (_funcId == 4)
        rseed = 3;
    else if (_funcId == 18)
        rseed = 17;
    else if (_funcId == 101 || _funcId == 102 || _funcId == 103 || _funcId == 107 || _funcId == 108 || _funcId == 109)
        rseed = 1;
    else if (_funcId == 104 || _funcId == 105 || _funcId == 106 || _funcId == 110 || _funcId == 111 || _funcId == 112)
        rseed = 8;
    else if (_funcId == 113 || _funcId == 114 || _funcId == 115)
        rseed = 7;
    else if (_funcId == 116 || _funcId == 117 || _funcId == 118)
        rseed = 10;
    else if (_funcId == 119 || _funcId == 120 || _funcId == 121)
        rseed = 14;
    else if (_funcId == 122 || _funcId == 123 || _funcId == 124)
        rseed = 17;
    else if (_funcId == 125 || _funcId == 126 || _funcId == 127)
        rseed = 19;
    else if (_funcId == 128 || _funcId == 129 || _funcId == 130)
        rseed = 21;
    else
        rseed = _funcId;

    rrseed = rseed + 10000 * _trialId;
    gauss(gval, 1, rrseed);
    gauss(gval2, 1, rrseed + 1);
    return fmin1(1000., fmax1(-1000., (round1(100.*100.*gval[0]/gval2[0])/100.)));
}

void setGlobalVariables(ParamStruct params) {
    DIM = params.DIM;
    trialid = params.instanceId;
    isInitDone = 0;
    return;
}

void initbenchmarkshelper(void *data) {
    experiment_total *exp = (experiment_total *) data;
    gval = (double*)malloc(sizeof(double) * 1);
    gval2 = (double*)malloc(sizeof(double) * 1);
    gvect = (double*)malloc(sizeof(double) * DIM * DIM);
    uniftmp = (double*)malloc(sizeof(double) * 2 * DIM * DIM);
    tmpvect = (double*)malloc(sizeof(double) * DIM);
    Xopt = (double*)malloc(sizeof(double) * DIM);
    TOL = (*exp).test.VTR;
    return;
}

void finibenchmarkshelper(void) {
    free(gval);
    //free(gval2);
    free(gvect);
    free(uniftmp);
    free(tmpvect);
    free(Xopt);
    return;
}

/* error handling routines - same arguments as printf, i.e. format first, then list of things to print */
/* this one exits after printing - severe error, not recoverable */
void ERROR(const char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "ERROR: ");
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");
/* and EXIT */
  exit(1);
}

/* same, but returns to the caller, mild error */
void WARNING(const char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "WARNING: ");
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");
/* and RETURN */
  return;
}

/* create complete pathName from filename and dirname 
   is SYSTEM dependent (should be some #ifdef WINDOWS etc ...)
   fullFileName should already be allocated, at least 1024 bytes long
*/
void createFullFileName(char *fullFileName, char *dirName, char *fileName)
{
char sLoc[1024];
if ( (strlen(fileName) + strlen(dirName)) > 1022 )
   ERROR("FileName will be too long for %s + %s", dirName, fileName);

sprintf(sLoc, "%s/%s", dirName, fileName);
strcpy(fullFileName, sLoc);
/* What if fullFileName is not long enough? */
return;
}

/* Checks if sDir exists, 
   creates it if not
   checks if is writable thereafter
   Fatal ERROR if anything fails
*/
void dirOK(char *sDir);

/* checks the existence of a file */
int existFile(char * fileName)
{
    FILE * fLoc = fopen(fileName, "r");
    if (fLoc == NULL) /* does not exist */
        return 0;
    fclose(fLoc);
    return 1; /* does exist */
}

/* opens a file after checking it is there */
FILE * bbobOpenFile(char * fileName)
{
    FILE * fileId;

    fileId = fopen(fileName,"a");
    if (fileId == NULL)
        ERROR("Could not open %s", fileName);
    return fileId;
}

