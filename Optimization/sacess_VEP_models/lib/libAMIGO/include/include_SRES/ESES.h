/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** (miu,lambda)-Evolution Strategy                                 **
 **                                                                 **
 ** For ACADEMIC RESEARCH, this is licensed with GPL license        **
 ** For COMMERCIAL ACTIVITIES, please contact the authors           **
 **                                                                 **
 ** Copyright (C) 2005 Xinglai Ji (jix1@ornl.gov)                   **
 **                                                                 **
 ** This program is distributed in the hope that it will be useful, **
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of  **
 ** MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the    **
 ** GNU General Public License for more details.                    **
 **                                                                 **
 ** You should have received a copy of the GNU General Public       **
 ** License along with is program; if not, write to the Free        **
 ** Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, **
 ** MA 02111-1307, USA.                                             **
 **                                                                 **
 ** Author: Xinglai Ji (jix1@ornl.gov)                              **
 ** Date:   Mar 2, 2005; Mar 3, 2005; Mar 4, 2005; Mar 7, 2005;     **
 **         Mar 8, 2005; Mar 10, 2005; Mar 21, 2005; Mar 22, 2005;  **
 ** Organization: Oak Ridge National Laboratory                     **
 ** Reference:                                                      **
 **   1. Thomas P. Runarsson and Xin Yao. 2000. Stochastic Ranking  **
 **      for Constrained Evolutionary Optimization. 4(3):284-294.   **
 **      http://cerium.raunvis.hi.is/~tpr/software/sres/            **
 **   2. Thomas Philip Runarsson and Xin Yao. 2005. Search Biases   **
 **      in Constrained Evolutionary Optimization. IEEE             **
 **      Transactions on Systems, Man and Cybernetics -- Part C:    **
 **      Applications and Reviews. 35(2):233-243.                   **
 *********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
    
#if defined win32 || defined _win32 || defined __win32 || defined WIN32 || defined _WIN32 || defined _WIN64 || defined WIN64
#else
#include <time.h>
#endif

#ifndef _es_es_h
#define _es_es_h

#define esDefPopsize 300
#define esDefGeneration 500
#define esDefGamma 0.85
#define esDefAlpha 0.2
#define esDefVarphi 1
#define esDefRetry 10
#define esDefESPlus 0
#define esDefESSlash 1

/*********************************************************************
 ** function of fitness and constraints                             **
 ** to calculate fitness and constraints and assign to ESIndividual **
 ** fg(x,dim, f, g)                                                 **
 *********************************************************************/
typedef void(*ESfcnFG) (double *, double *, double *, void *);

/*********************************************************************
 ** function to transform x(op) and sp                              **
 ** double f(double)                                                **
 *********************************************************************/
typedef double(*ESfcnTrsfm) (double );

/*********************************************************************
 ** ESParameter: struct for ES-parameter                            **
 ** fg: functions of fitness and constraints                        **
 ** trsfm: to transform sp/op                                       **
 ** es: ES process, esDefESPlus/esDefESSlash                        **
 ** eslambda: lambda+miu or lambda according to ES process          **
 ** seed: random seed                                               **
 ** constraint: number of constraints                               **
 ** dim: dimension/number of genes in genome                        **
 ** ub[dim]: up bounds                                              **
 ** lb[dim]: low bounds                                             **
 ** spb[dim]: bounds on sp , spb = (ub-lb)/sqrt(dim)                **
 ** miu: parent/population size                                     **
 ** lambda: offsping/population size                                **
 ** gen: number of generations                                      **
 ** gamma: usually esDefGamma=0.85                                  **
 ** alpha: usually esDefAlpha=0.2                                   **
 ** chi: chi = 1/2n +1/2sqrt(n)                                     **
 ** varphi: = sqrt((2/chi)*log((1/alpha)*(exp(varphi^2*chi/2)       **
 **                  -(1-alpha))))                                  **
 **         expected rate of convergence                            **
 ** retry: retry times to check bounds                              **
 ** tau: learning rates: tau = varphi/(sqrt(2*sqrt(dim)))           **
 ** tar_: learning rates: tau_ = varphi((sqrt(2*dim)                **
 *********************************************************************/
typedef struct
  {
    ESfcnFG fg;
    ESfcnTrsfm *trsfm;
    int seed;
    int constraint;
    int dim;
    double *ub;
    double* lb;
    double* spb;
      int miu;
    int lambda;
    int gen;
    double gamma;
    double alpha;
    double varphi;
    int retry;
    double chi;
      double tau;
      double tau_;
      int es;
      int eslambda;
  } ESParameter;

/*********************************************************************
 ** ESIndividual: struct for each individual/genome                 **
 ** op[dim]: genes/objective parameters                             **
 ** sp[dim]: strategy parameters                                    **
 ** f: fitness                                                      **
 ** g[constraint]: constraint value                                 **
 ** phi: phi = sum( max(0,g)^2 )                                    **
 *********************************************************************/
typedef struct
  {
    double *op;
    double *sp;
      double f;
      double phi;
    double *g;
  } ESIndividual;

/*********************************************************************
 ** ESPopulation: struct for population                             **
 ** member[lambda]: each individual in this population              **
 ** f[lambda]: fitness                                              **
 ** phi[lambda]: constraints                                        **
 ** index[lambda]: ranking index                                    **
 *********************************************************************/
typedef struct
  {
    ESIndividual** member;
    double *f;
    double *phi;
    int *index;
  } ESPopulation;

/*********************************************************************
 ** ESStatistics: struct for ES-statistics                          **
 ** begintime: begin time when intializing                          **
 ** nowtime: time when do statistics                                **
 ** dt: nowtime - begintime                                         **
 ** bestindvdl: best individual                                     **
 ** thisbestindvdl: best individual in this generation              **
 ** bestgen: generation of the bestindividual                       **
 ** curgen: current generation                                      **
 *********************************************************************/
typedef struct
  {
    time_t begintime;
    time_t nowtime;
    int dt;
    int bestgen;
    int curgen;
    ESIndividual* bestindvdl;
    ESIndividual* thisbestindvdl;
  } ESStatistics;

/*********************************************************************
 ** initialize: parameters,populations and random seed              **
 ** ESInitial(seed, param,trsfm, fg,es,constraint,dim,ub,lb,miu,    **
 **            lambda,gen, gamma, alpha, varphi, retry,             **
 **             population, stats)                                  **
 ** seed: random seed, usually esDefSeed=0 (pid*time)               **
 ** outseed: seed value assigned , for next use                     **
 ** param: point to parameter                                       **
 ** fg: functions of fitness and constraints                        **
 ** trsfm: to transform sp/op                                       **
 ** es: ES process, esDefESPlus/esDefESSlash                        **
 ** constraint: number of constraints                               **
 ** dim: dimension/number of genes in genome                        **
 ** ub[dim]: up bounds                                              **
 ** lb[dim]: low bounds                                             **
 ** miu: parent/population size                                     **
 ** lambda: offsping/population size                                **
 ** gen: number of generations                                      **
 ** gamma: usually esDefGamma=0.85                                  **
 ** alpha: usually esDefAlpha=0.2                                   **
 ** chi: chi = 1/2n +1/2sqrt(n)                                     **
 ** varphi: = sqrt((2/chi)*log((1/alpha)*(exp(varphi^2*chi/2)       **
 **                  -(1-alpha))))                                  **
 **         expected rate of convergence                            **
 ** retry: retry times to check bounds                              **
 ** tau: learning rates: tau = varphi/(sqrt(2*sqrt(dim)))           **
 ** tar_: learning rates: tau_ = varphi((sqrt(2*dim)                **
 ** population: point to this population                            **
 ** stats: point to statistics                                      **
 **                                                                 **
 ** ESDeInitial(param,populationi,stats)                            **
 ** free param and population                                       **
 *********************************************************************/
void ESInitial(unsigned int, ESParameter**, ESfcnTrsfm *,   \
               ESfcnFG,int, int,int,double*,double*,int,int,int,  \
               double, double, double, int,  \
               ESPopulation**, ESStatistics**,void* data);
void ESDeInitial(ESParameter*, ESPopulation*, ESStatistics*,void* data);
/*********************************************************************
 ** initialize parameters                                           **
 ** ESInitialParam(param,trsfm,fg,es,constraint,                    **
 **                dim,ub,lb,miu,lambda,gen)                        **
 ** param: point to parameter                                       **
 ** fg: functions of fitness and constraints                        **
 ** trsfm: to transform sp/op                                       **
 ** es: ES process, esDefESPlus/esDefESSlash                        **
 ** seed: reserve seed for next use                                 **
 ** constraint: number of constraints                               **
 ** dim: dimension/number of genes in genome                        **
 ** ub[dim]: up bounds                                              **
 ** lb[dim]: low bounds                                             **
 ** spb[dim]: bounds on sp , spb = (ub-lb)/sqrt(dim)                **
 ** miu: parent/population size                                     **
 ** lambda: offsping/population size                                **
 ** gen: number of generations                                      **
 ** gamma: usually esDefGamma=0.85                                  **
 ** alpha: usually esDefAlpha=0.2                                   **
 ** chi: chi = 1/2n +1/2sqrt(n)                                     **
 ** varphi: = sqrt((2/chi)*log((1/alpha)*(exp(varphi^2*chi/2)       **
 **                  -(1-alpha))))                                  **
 **         expected rate of convergence                            **
 ** retry: retry times to check bounds                              **
 ** tau: learning rates: tau = varphi/(sqrt(2*sqrt(dim)))           **
 ** tar_: learning rates: tau_ = varphi((sqrt(2*dim)                **
 **                                                                 **
 ** ESDeInitialParam(param)                                         **
 ** free param                                                      **
 *********************************************************************/
void ESInitialParam(ESParameter **, ESfcnTrsfm *, ESfcnFG, int,   \
                    unsigned int,  \
                    int,int,double*,double*,int,int,int,  \
                    double, double, double, int,void*);
void ESDeInitialParam(ESParameter *);
/*********************************************************************
 ** initialize population                                           **
 ** ESInitialPopulation(population,param)                           **
 ** population: point to this population                            **
 ** param: point to this parameter                                  **
 **   -> index: 0->lambda-1                                         **
 **   -> individual[lambda]                                         **
 **   -> fg(individual)                                             **
 **   -> f,phi                                                      **
 ** the initialization is looked as first generation                **
 **                                                                 **
 ** ESDeInitialPopulation(population, param)                        **
 ** free population                                                 **
 *********************************************************************/
void ESInitialPopulation(ESPopulation **, ESParameter *,void *data);
void ESDeInitialPopulation(ESPopulation *, ESParameter *,void* data);
/*********************************************************************
 ** initialize individual                                           **
 ** ESInitialIndividual(indvdl, param)                              **
 ** to calculate f,g,and phi                                        **
 ** to initialize op and sp                                         **
 ** phi=sum{(g>0)^2}                                                **
 ** op = rand(lb, ub)                                               **
 ** sp = (ub - lb)/sqrt(dim)                                        **
 **                                                                 **
 **                                                                 **
 ** ESDeInitialIndividual(indvdl, param)                            **
 ** free individual                                                 **
 **                                                                 **
 ** ESPrintOp(indvdl, param)                                        **
 ** print individual information, indvdl->op                        **
 ** ESPrintSp(indvdl, param)                                        **
 ** print individual information, indvdl->sp                        **
 *********************************************************************/
void ESInitialIndividual(ESIndividual **, ESParameter *,void* data);
void ESDeInitialIndividual(ESIndividual *);
void ESPrintIndividual(ESIndividual *, ESParameter *);
void ESPrintOp(ESIndividual *, ESParameter *);
void ESPrintSp(ESIndividual *, ESParameter *);
/*********************************************************************
 ** copy a individual                                               **
 ** ESCopyIndividual(from, to, param)                               **
 *********************************************************************/
void ESCopyIndividual(ESIndividual *, ESIndividual *, ESParameter *);
/*********************************************************************
 ** initialize statistics                                           **
 ** ESInitialStat(stats, population, param)                         **
 ** to intialize time, curgen, bestindvdl,thisbestindvdl            **
 ** not to do the first statistics                                  **
 ** to set dt, bestgen                                              **
 **                                                                 **
 ** ESDeInitialStat(stats)                                          **
 ** free statistics                                                 **
 *********************************************************************/
void ESInitialStat(ESStatistics **, ESPopulation *, ESParameter *,void *data);
void ESDeInitialStat(ESStatistics *, void*data);
/*********************************************************************
 ** do statistics                                                   **
 ** ESDoStat(stats, population, param)                              **
 ** to set nowtime, dt, curgen, bestgen, (this)bestindvdl           **
 ** to do statistics                                                **
 ** if there's no feasible best, do nothing                         **
 ** the initialization is looked as zero generation                 **
 **                                                                 **
 ** ESPrintStat(stats, param)                                       **
 ** print statistics information                                    **
 ** gen=,time=,dt=,bestgen=,bestfitness=,bestindividual=,           **
 *********************************************************************/
void ESDoStat(ESStatistics *, ESPopulation *, ESParameter *);
void ESPrintStat(ESStatistics *, ESParameter *);

/*********************************************************************
 ** stepwise evolution                                              **
 ** ESStep(population, param, stats, pf)                            **
 **                                                                 **
 ** -> Stochastic ranking -> sort population based on ranking index **
 ** -> Mutate (recalculate f/g/phi) -> do statistics analysis on    **
 ** this generation -> print statistics information                 **
 *********************************************************************/
void ESStep(ESPopulation *, ESParameter *, ESStatistics *, double,void* data);

/*********************************************************************
 ** sort population based on Index by ESSRSort                      **
 ** ESSortPopulation(population, param)                             **
 *********************************************************************/
void ESSortPopulation(ESPopulation *, ESParameter *);

/*********************************************************************
 ** select the next generation                                      **
 ** ESSelectPopulation(population, param)                           **
 ** select first miu offsprings to fill up the next generation      **
 ** miu -> lambda : 1..miu,1..miu,..,lambda                         **
 *********************************************************************/
void ESSelectPopulation(ESPopulation *, ESParameter *);

/*********************************************************************
 ** mutate                                                          **
 ** ESMutate(population, param)                                     **
 **                                                                 **
 ** sp_ : copy of sp                                                **
 ** op_ : copy of op                                                **
 ** update sp                                                       **
 ** traditional technique using exponential smoothing               **
 ** sp(1->miu-1) : unchanged                                        **
 ** sp(miu->lambda): sp = sp_*exp(tau_*N(0,1) + tau*Nj(0,1))        **
 **                  Nj : random number generated for each j        **
 ** check sp bound                                                  **
 ** if(sp > bound) then sp = bound                                  **
 ** differential variation                                          **
 ** op(1->miu-1) = op_ + gamma*(op_[1] - op_[i+1])                  **
 ** mutation                                                        **
 ** op(miu->lambda): op = op_ +sp * N(0,1)                          **
 ** check op bound                                                  **
 ** if(op > ub || op < lb) then try retry times                     **
 **                        op = op_ + sp * N(0,1)                   **
 ** if still not in bound then op = op_                             **
 ** exponential smoothing                                           **
 ** sp(miu->lambda): sp = sp_ + alpha * (sp - sp_)                  **
 **                                                                 **
 ** re-calculate f/g/phi                                            **
 *********************************************************************/
void ESMutate(ESPopulation *, ESParameter *, void*data);

#endif

#ifdef __cplusplus
}
#endif

