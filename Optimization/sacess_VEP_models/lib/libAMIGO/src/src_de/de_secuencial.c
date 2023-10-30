/*************************************************************************
**									**
**		D I F F E R E N T I A L    E V O L U T I O N		**
**									**
** This C-code implements Differential Evolution (DE) algorithm		**
** (DE/rand/1/bin version) described in:				**
**									**
** Rainer Storn and Kenneth V. Price, Differential evolution - a Simple **
** and Efficient Adaptive Scheme for Global Optimization Over		**
** Continuous paces, Technical Report, TR-95-012, ICSI, March, 1995,	**
** Available from							**
** www.icsi.berkeley.edu/ftp/global/pub/techreports/1995/tr-95-012.pdf	**
**									**
** Rainer Storn and Kenneth V. Price, Differential evolution - a Simple **
** and Efficient Heuristic for Global Optimization Over Continuous	**
** paces, Journal of Global Optimization, 11 (4), pp. 341-359, Dec,	**
** 1997, Kluwer Academic Publisher					**
**									**
** Kenneth V. Price, Rainer Storn and Jouni Lampinen, Differential 	**
** Evolution: A Practical Approach to  Global Optimization, 		**
** Springer-Verlag, Berlin, ISBN: 3-540-20950-6, 2005			**
**									**
** Code is free for scientific and academic use. Use for other purpose	**
** is not allowed without a permission of the author. There is no	**	
** warranty of any kind about correctness of the code and if you find	**
** a bug, please, inform the author.					**
**									**
**									**
** Author: Saku Kukkonen						**
**	   Lappeenranta University of Technology			**
**	   Department of Information Technology				**
**	   P.O.Box 20, FIN-53851 LAPPEENRANTA, Finland			**
**	   E-mail: saku.kukkonen@lut.fi					**
**									**
** Date: 14.6.2005							**
** Modified: 3.1.2007	(replacement/removal of				**
**			non-standard C-functions)			**
**									**
**									**
** Program: de.c (needs function func.c)				**
**									**
** Compile: gcc -Wall -pedantic -ansi -O -o de de.c problem.c -lm	**
**									**
** Run: ./de -h		for options 					**
**									**
** If an output file name is defined, then final population is written	**
** to this file so that the file has one individual on each row listing **
** first decision variable values and finally the objective function	**
** value, i.e.,								**
**									**
**	x1_1 x1_2 x1_3  ...     x1_D f1					**
**	x2_1 x2_2 x2_3  ...     x2_D f2					**
**	x3_1 x3_2 x3_3  ...     x3_D f3					**
**	 :    :    :             :					**
**	xN_1 xN_2 xN_3  ...     xN_D fN					**
**									**
**									**
** Please, acknowledge and inform the author if you use this code.	**
**									**
*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <de.h>

/* Function definitions		*/



/* Random number generator defined by URAND should return
   double-precision floating-point values uniformly distributed
   over the interval [0.0, 1.0)					*/

#define URAND	((double)rand()/((double)RAND_MAX + 1.0))

/* Definition for random number generator initialization	*/

#define INITRAND srand(time(0))

/* Usage for the program	*/


double func(double *, int D, double* Xl, double* Xu);
int usage(char *);

int usage(char *str)
{
   fprintf(stderr,"Usage: %s [-h] [-u] [-s] [-N NP (20*D)] ", str);
   fprintf(stderr,"[-G Gmax (1000)]\n");
   fprintf(stderr,"\t[-C crossover constant, CR (0.9)]\n");
   fprintf(stderr,"\t[-F mutation scaling factor, F (0.9)]\n");
   fprintf(stderr,"\t[-o <outputfile>]\n\n");
   fprintf(stderr,"\t-s does not initialize random number generator\n");
   exit(-1);
}

int de(int D,DE_opts* de_opts,double *Xl,double* Xu, double(*func)(double*,int,double*,double*,void*),void* data){

   register int i, j, k, r1, r2, r3, jrand, numofFE = 0;
   int c, index = -1, s = 1;
   double **popul, **next, **ptr, *iptr, *U, CR = 0.9, F = 0.9,
	min_value = DBL_MAX, totaltime = 0.0;
   char *ofile = NULL;
   FILE *fid;
   clock_t starttime, endtime;
   
   int NP, Gmax;
   NP=de_opts->NP;
   Gmax=de_opts->Gmax;


   /* Parse command line arguments given by user	*/
   s = 1;	/* different runs				*/
     

   if (s) INITRAND;

   /* Printing out information about optimization process for the user	*/

   printf("Program parameters: ");
   printf("NP = %d, Gmax = %d, CR = %.2f, F = %.2f\n",
	NP, Gmax, CR, F);

   printf("Dimension of the problem: %d\n", D);

   /* Starting timer    */

   starttime = clock();


   /* Allocating memory for current and next populations, intializing
      current population with uniformly distributed random values and
      calculating value for the objective function	*/

   popul = (double **)malloc(NP*sizeof(double *));
   if (popul == NULL) perror("malloc");

   next = (double **)malloc(NP*sizeof(double *));
   if (next == NULL) perror("malloc");

   for (i=0; i < NP; i++)
   {
      popul[i] = (double *)malloc((D+1)*sizeof(double));
      if (popul[i] == NULL) perror("malloc");

      for (j=0; j < D; j++)
         popul[i][j] = Xl[j] + (Xu[j] - Xl[j])*URAND;

      popul[i][D] = func(popul[i],D, Xl,Xu,data);
      numofFE++;

      next[i] = (double *)malloc((D+1)*sizeof(double));
      if (next[i] == NULL) perror("malloc");
   }

   /* Allocating memory for a trial vector U	*/

   U = (double *)malloc((D+1)*sizeof(double));
   if (U == NULL) perror("malloc");


   /* The main loop of the algorithm	*/
   for (k=0; k < Gmax; k++){

      for (i=0; i < NP; i++){	/* Going through whole population	*/

         /* Selecting random indeces r1, r2, and r3 to individuls of
            the population such that i != r1 != r2 != r3	*/
         do{
            r1 = (int)(NP*URAND);
         } while( r1==i );
         
		 do{
            r2 = (int)(NP*URAND);
         } while( r2==i || r2==r1);
         
		 do{
            r3 = (int)(NP*URAND);
         } while( r3==i || r3==r1 || r3==r2 );

         jrand = (int)(D*URAND);

         /* Mutation and crossover	*/

         for (j=0; j < D; j++)
         {
            if (URAND < CR || j == jrand)
            {
               U[j] = popul[r3][j] + F*(popul[r1][j] - popul[r2][j]);
            }
            else
               U[j] = popul[i][j];
         }

         U[D] = func(U,D,Xl,Xu,data);
         numofFE++;

         /* Comparing the trial vector 'U' and the old individual
            'next[i]' and selecting better one to continue in the
            next population.	*/

         if (U[D] <= popul[i][D]){
            iptr = U;
            U = next[i];
            next[i] = iptr;
         }
         else{
            for (j=0; j <= D; j++)
            next[i][j] = popul[i][j];
         }

      }	/* End of the going through whole population	*/
	 // printf("numofFE=%d best_found=%f\n",numofFE,popul[index][D]);

      /* Pointers of old and new populations are swapped	*/

      ptr = popul;
      popul = next;
      next = ptr;

   }	/* End of the main loop		*/


   /* Stopping timer	*/

   endtime = clock();
   totaltime = (double)(endtime - starttime);


   /* If user has defined output file, the whole final population is
      saved to the file						*/

   /* Finding best individual	*/

   for (i=0; i < NP; i++){

      if (popul[i][D] < min_value){
         min_value = popul[i][D];
         index = i;
      }
   }
   /* Printing out information about optimization process for the user	*/

   printf("Execution time: %.3f s\n", totaltime / (double)CLOCKS_PER_SEC);
   printf("Number of objective function evaluations: %d\n", numofFE);

   printf("Solution:\nValues of variables: ");
   for (i=0; i < D; i++)
      printf("%.15f ", popul[index][i]);

   printf("\nObjective function value: ");
   printf("%.15f\n", popul[index][D]);


   /* Freeing dynamically allocated memory	*/

   for (i=0; i < NP; i++){
      free(popul[i]);
      free(next[i]);
   }
   free(popul);
   free(next);
   free(U);

   return(0);
}

