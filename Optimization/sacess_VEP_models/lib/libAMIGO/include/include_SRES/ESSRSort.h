/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** Stochastic Bubble Sort                                          **
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
 ** Date:   Mar 1, 2005; Mar 2, 2005; Mar 3, 2005; Mar 21, 2005;    **
 ** Organization: Oak Ridge National Laboratory                     **
 ** Reference:                                                      **
 **   Thomas P. Runarsson and Xin Yao. 2000. Stochastic Ranking     **
 **   for Constrained Evolutionary Optimization. 4(3):284-294.      **
 **   http://cerium.raunvis.hi.is/~tpr/software/sres/               **
 *********************************************************************/

#ifdef __cplusplus 
extern "C" {
#endif

#ifndef _es_srsort_h
#define _es_srsort_h

#define essrDefPf 0.45

/*********************************************************************
 ** Stochastic Bubble Sort                                          **
 ** void ESSRSort(f,phi,pf,eslambda,N,I)                            **
 ** f[eslambda]: fitness                                            **
 ** phi[eslambda]: constraints                                      **
 ** pf: stochastic ranking, in (0,1), generally pf<0.5, pf=0.45     **
 ** eslambda: population size -- offspring or parent+offspring      **
 ** N: usually N=eslambda                                           **
 ** I[eslambda]: sort index                                         **
 **                                                                 **
 ** for i=1 to N do                                                 **
 **   for j=1 to eslambda-1 do                                      **
 **     rand u in (0,1)                                             **
 **     if( phi(I(j)) == phi(I(j+1)) == 0 || u<pf)                  **
 **       if( f(I(j)) > f(I(j+1)) )                                 **
 **         swap( I(j), I(j+1) )                                    **
 **     else                                                        **
 **       if( phi(I(j)) > phi(I(j+1)) )                             **
 **         swap( I(j), I(j+1) )                                    **
 **   if(numberOFswap == 0)                                         **
 **     break                                                       **
 *********************************************************************/
void ESSRSort(double *, double *, double , int , int , int *);

#endif

#ifdef __cplusplus 
}
#endif
