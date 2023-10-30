/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** deal with shared functions                                      **
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
 ** Date:   Mar 2, 2005; Mar 3, 2005; Mar 7, 2005;                  **
 ** Organization: Oak Ridge National Laboratory                     **
 *********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _share_func_h
#define _share_func_h

#define shareDefSeed 0
#define shareDefTrue 1
#define shareDefFalse 0
#define shareDefMinZero 1e-20
#define shareDefMaxLine 4096
#define shareDefNullYes 0
#define shareDefNullNo 1

/*********************************************************************
 ** uniform random                                                  **
 ** to output a random value between min and max                    **
 ** double ShareRand(min,max)                                       **
 ** min: min value                                                  **
 ** max: max value                                                  **
 **                                                                 **
 ** return value = min + (max-min)*rand()/(RAND_MAX)                **
 **                                                                 **
 ** void ShareRandVec(s, n, min, max)                               **
 ** return s=vec(n)                                                 **
 *********************************************************************/
double ShareRand(double , double );
void ShareRandVec(double *, int, double, double);

/*********************************************************************
 ** gaussian random: normal distribution                            **
 ** N(mean, dev)                                                    **
 **                                                                 **
 ** void ShareNormalRandVec(s, n, mean, dev)                        **
 ** return s=vec(n)                                                 **
 **                                                                 **
 ** source code from                                                **
 **        http://remus.rutgers.edu/~rhoads/Code/code.html          **
 *********************************************************************/
double ShareNormalRand(double , double );
void ShareNormalRandVec(double *, int, double, double);

/*********************************************************************
 ** to set random seed                                              **
 ** void ShareSeed(inseed, outseed)                                 **
 ** if inseed==0, then use pid*time as seed                         **
 ** if inseed!=0, use this seed users set                           **
 ** outseed is to be used in the next step                          **
 *********************************************************************/
void ShareSeed(unsigned int, unsigned int *);

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1c(size): size*char                                 **
 ** to realloc memories                                             **
 ** ShareReallocM1c(s, size)                                        **
 ** to free memories                                                **
 ** ShareFreeM1c(s)                                                 **
 *********************************************************************/
char * ShareMallocM1c(int);
char * ShareReallocM1c(char *, int);
void ShareFreeM1c(char *);
/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM2c(size1, size2): size1*(char*), size2*char         **
 ** to realloc memories                                             **
 ** ShareReallocM2c(s, size1, size2)                                **
 ** to free memories                                                **
 ** ShareFreeM2c(s, size)                                           **
 *********************************************************************/
char ** ShareMallocM2c(int,int);
char ** ShareReallocM2c(char **, int, int);
void ShareFreeM2c(char **, int);
/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1i(size): size*int                                  **
 *********************************************************************/
int * ShareMallocM1i(int);
/*********************************************************************
 ** to free memories                                                **
 ** ShareFreeM1i(s)                                                 **
 *********************************************************************/
void ShareFreeM1i(int *);
/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1d(size): size*double                               **
 *********************************************************************/
double * ShareMallocM1d(int);
/*********************************************************************
 ** to free memories                                                **
 ** ShareFreeM1d(s)                                                 **
 *********************************************************************/
void ShareFreeM1d(double *);
/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM2d(size1, size2)                                    **
 **                                                                 **
 ** to free memories                                                **
 ** ShareFreeM2d(s, size)                                           **
 *********************************************************************/
double **ShareMallocM2d(int , int);
void ShareFreeM2d(double **, int);
/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM3d(size1, size2, siez3)                             **
 **                                                                 **
 ** to free memories                                                **
 ** ShareFreeM3d(s, size1, size2)                                   **
 *********************************************************************/
double ***ShareMallocM3d(int , int , int );
void ShareFreeM3d(double ***, int , int );
/*********************************************************************
 ** to check if it's equal to zero                                  **
 ** if(x<min) return true                                           **
 ** true=shareDefTrue=1, false=shareDefFalse=0                      **
 *********************************************************************/
int ShareIsZero(double);

/*********************************************************************
 ** ShareSplitStr(buf0, sepa, len, flag)                            **
 **   parse string based on separator string                        **
 **   similar to strtok                                             **
 **   return sub string array                                       **
 **   sub string array in strarr                                    **
 **   flag = 0 allow str=0                                          **
 **   flag = 1 only allow non-0 string                              **
 **   buf0: string to be parsed                                     **
 **   sepa: mark string                                             **
 **   len: number of sub strings                                    **
 *********************************************************************/
char ** ShareSplitStr(const char *, const char *, int *, int );
/*********************************************************************
 ** ShareChop(s)                                                    **
 ** cut the '\n' from the tail of a string                          **
 *********************************************************************/
void ShareChop(char *);

#endif

#ifdef __cplusplus
}
#endif
