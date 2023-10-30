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
 ** Date:   Mar 2, 2005; Mar 3, 2005                                **
 ** Organization: Oak Ridge National Laboratory                     **
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "sharefunc.h"

/*********************************************************************
 ** uniform random                                                  **
 ** double ShareRand(min,max)                                       **
 ** min: min value                                                  **
 ** max: max value                                                  **
 **                                                                 **
 ** return value = min + (max-min)*rand()/(RAND_MAX)                **
 **                                                                 **
 ** double ShareRandVec(n, min, max)                                **
 ** return s=vec(n)                                                 **
 *********************************************************************/
double ShareRand(double min, double max)
{
  double delta;
  double value;

  delta = max - min;

  value = (double)rand()/(RAND_MAX);
  value = min + delta*value;

  return value;
}

void ShareRandVec(double *s, int n, double min, double max)
{
  int i;

  for(i=0; i<n; i++)
    s[i] = ShareRand(min, max);

  return;
}

/*********************************************************************
 ** void ShareSeed(inseed, outseed)                                 **
 ** if inseed==0, then use pid*time as seed                         **
 ** if inseed!=0, use this seed users set                           **
 ** outseed is to be used in the next step                          **
 *********************************************************************/
void ShareSeed(unsigned int inseed, unsigned int *outseed)
{
  time_t nowtime;

/*********************************************************************
 ** shareDefSeed = 0                                                **
 *********************************************************************/
  if(inseed != shareDefSeed)
  {
    srand(inseed);
    *outseed = inseed;
    return;
  }

  time(&nowtime);
  inseed = nowtime;
  srand(inseed);
  *outseed = inseed;

  return;
}

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
double ShareNormalRand(double mean, double dev)
{
  static double V2, fac;
  static int phase = 0;
  double S, Z, U1, U2, V1;

  if (phase)
    Z = V2 * fac;
  else
  {
    do 
    {
      U1 = (double)rand() / RAND_MAX;
      U2 = (double)rand() / RAND_MAX;
      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;
    } while(S >= 1 || S ==0.0);

    fac = sqrt (-2 * log(S) / S);
    Z = V1 * fac;
  }

  phase = 1 - phase;

  Z = mean + dev*Z;

  return Z;
}

void ShareNormalRandVec(double *s, int n, double mean, double dev)
{
  int i;

  for(i=0; i<n; i++)
    s[i] = ShareNormalRand(mean, dev);

  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1c(size): size*char                                 **
 ** to realloc memories                                             **
 ** ShareReallocM1c(s, size)                                           **
 ** to free memories                                                **
 ** SharefreeM1c(s)                                                 **
 *********************************************************************/
char * ShareMallocM1c(int size)
{
  char *s = NULL;
  if( (s=(char*)calloc(size, sizeof(char))) == NULL )
  {
    printf("char * malloc error!\n");
    exit(1);
  }
  return s;
}
char * ShareReallocM1c(char *s, int size)
{
  if( !s ) 
  {
    s = ShareMallocM1c(size);
    return s;
  }
  if( (s = (char *)realloc(s, size*sizeof(char))) == NULL )
  {
    printf("char * realloc error!");
    exit(1);
  }
  return s;
}
void ShareFreeM1c(char *s)
{
  if(s)
    free((void *)s);
  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM2c(size1, size2): size1*(char*), size2*char         **
 **                                                                 **
 ** to realloc memories                                             **
 ** ShareReallocM2c(s, size1, size2)                                **
 ** to free memories                                                **
 ** ShareFreeM2c(s, size)                                           **
 *********************************************************************/
char ** ShareMallocM2c(int size1, int size2)
{
  int i;
  char **s = NULL;

  if( (s = (char **)calloc(size1,sizeof(char *))) == NULL )
  {
    printf("char ** malloc error!\n");
    exit(1);
  }
  if( size2 > 0 ) 
  {
    for( i=0; i<size1; i++ )
      s[i] = ShareMallocM1c(size2);
  }
  return s;
}

char ** ShareReallocM2c(char **s, int size1, int size2)
{
  int i;
                                                                                
  if( (s = (char **)realloc(s, size1*sizeof(char *))) == NULL )
  {
    printf("char ** realloc error!");
    exit(1);
  }
  if( size2 > 0 ) 
  {
    for( i=0; i<size1; i++ )
      s[i] = ShareReallocM1c(s[i], size2);
  }
  return s;
}

void ShareFreeM2c(char **s, int size)
{
  int i;

  if(!s)
    return;
  for(i=0; i<size; i++)
    ShareFreeM1c(s[i]);
  free((void *)s);

  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1i(size): size*int                                  **
 *********************************************************************/
int * ShareMallocM1i(int size)
{
  int *s = NULL;
  if( (s=(int*)calloc(size, sizeof(int))) == NULL )
  {
    printf("int * malloc error!\n");
    exit(1);
  }
  return s;
}

/*********************************************************************
 ** to free memories                                                **
 ** SharefreeM1i(s)                                                 **
 *********************************************************************/
void ShareFreeM1i(int *s)
{
  if(s)
    free((void *)s);
  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM1d(size): size*double                               **
 *********************************************************************/
double * ShareMallocM1d(int size)
{
  double *s = NULL;
  if( (s=(double*)calloc(size, sizeof(double))) == NULL )
  {
    printf("double * malloc error!\n");
    exit(1);
  }
  return s;
}

/*********************************************************************
 ** to free memories                                                **
 ** SharefreeM1d(s)                                                 **
 *********************************************************************/
void ShareFreeM1d(double *s)
{
  if(s)
    free((void *)s);
  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM2d(size1, size2)                                    **
 **                                                                 **
 ** to free memories                                                **
 ** ShareFreeM2d(s, size)                                           **
 *********************************************************************/
double **ShareMallocM2d(int size1, int size2)
{
  int i;
  double **s = NULL;

  if( (s=(double **)calloc(size1, sizeof(double *))) == NULL)
  {
    printf("double ** malloc error!\n");
    exit(1);
  }
  if( size2 > 0 )
  {
    for(i=0; i<size1; i++)
      s[i] = ShareMallocM1d(size2);
  }
 
  return s;
}

void ShareFreeM2d(double **s, int size)
{
  int i;

  if(!s)
    return;
  for(i=0; i<size; i++)
    ShareFreeM1d(s[i]);
  free((void *)s); 

  return;
}

/*********************************************************************
 ** to malloc memories                                              **
 ** ShareMallocM3d(size1, size2, siez3)                             **
 **                                                                 **
 ** to free memories                                                **
 ** ShareFreeM3d(s, size1, size2)                                   **
 *********************************************************************/
double ***ShareMallocM3d(int size1, int size2, int size3)
{
  int i;
  double ***s = NULL;

  if( (s=(double ***)calloc(size1, sizeof(double **))) == NULL)
  {
    printf("double *** malloc error!\n");
    exit(1);
  }
  if( size2 > 0 )
  {
    for(i=0; i<size1; i++)
      s[i] = ShareMallocM2d(size2, size3);
  }
 
  return s;
}

void ShareFreeM3d(double ***s, int size1, int size2)
{
  int i;

  if(!s)
    return;
  for(i=0; i<size1; i++)
    ShareFreeM2d(s[i], size2);
  free((void *)s); 

  return;
}

/*********************************************************************
 ** to check if it's equal to zero                                  **
 ** if(x<shareDefMinZero) return true                               **
 ** true=shareDefTrue=1, false=shareDefFalse=0                      **
 *********************************************************************/
int ShareIsZero(double x)
{
  if(fabs(x) <= shareDefMinZero)
    return shareDefTrue;
  else
    return shareDefFalse;
}

/*********************************************************************
 ** ShareSplitStr(buf0, sepa, len, flag)                            **
 **   parse string based on separator string                        **
 **   similar to strtok                                             **
 **   return sub string array                                       **
 **   sub string array in strarr                                    **
 **   flag = shareDefNullYes=0 allow str=0                          **
 **   flag = shareDefNullNo=1 only allow non-0 string               **
 **   buf0: string to be parsed                                     **
 **   sepa: mark string                                             **
 **   len: number of sub strings                                    **
 *********************************************************************/
char ** ShareSplitStr(const char *buf0, const char *sepa, int *len, int flag)
{
  char buf[shareDefMaxLine];
  char *str0, *str1, *str2;
  int num=0;
  char *s=NULL;
  char **strarr=NULL;
  int n;

  n=(strlen(buf0)+1>shareDefMaxLine)?shareDefMaxLine-1:strlen(buf0);
  strncpy(buf, buf0, n);
  buf[n] = 0;
  str0 = str1 = buf;
  str2 = buf + strlen(buf);
    
  strarr = ShareMallocM2c(1, 0);
  
  while( str0 <= str2 ) {
    str1 = strstr(str0, sepa);
    if(str1 == NULL)
      str1 = str2;
    
    if( ((str1>=str0)&&(flag==shareDefNullYes))   \
           || ((str1>str0)&&(flag==shareDefNullNo)) ) {
      strarr = ShareReallocM2c(strarr,num+1,0);
      s = ShareMallocM1c(str1-str0+1);
      
      strncpy(s, str0, str1-str0);
      s[str1-str0] = 0;
      strarr[num] = s;
      num++;
    }
    str0 = str1 + strlen(sepa);
  }

  *len = num;
  return strarr;
}

/*********************************************************************
 ** ShareChop(s)                                                    **
 ** cut the '\n' from the tail of a string                          **
 *********************************************************************/
void ShareChop(char *source)
{
  int nLen;

  nLen = strlen(source);
  
  if ( source[nLen-1] == 10 || source[nLen-1] == 13 )
    source[nLen-1] = 0;
  if ( source[nLen-2] == 10 || source[nLen-2] == 13 )
    source[nLen-2] = 0;
}
