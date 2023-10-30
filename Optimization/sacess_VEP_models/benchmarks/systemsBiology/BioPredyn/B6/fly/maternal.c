/**
 * @file maternal.c                                            
 * @author JR, modified by Yoginho
 *                                                               
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Defines a lot of constants, bias-related stuff and functions to deal
 * with the data structures
 *                                                               
 * 1. structs and constants for things used throughout the fly   
 *    model specific part of the code (since this file will al-  
 *    ways get included when we run the fly model); therefore    
 *    all fly-specific things that used to be in global.h are    
 *    now in here                                                
 * 2. stuff that deals with that part of the blastoderm which is 
 *    fixed by the maternal genotype, i.e. division schedules    
 *    (including the rules for zygotic.c, the number of nuclei   
 *    at each cleavage cycle and the lineage number of the most  
 *    anterior nucleus for each cleavage cycle), bicoid gra-     
 *    dients and diffusion schedules                             
 * 3. bias-related stuff is also in here since it's needed to    
 *    run the model                                              
 */

#include <float.h>              /* for DBL_EPSILON */
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
//#include <error.h>              /* for error and linked list functions  */

#include "maternal.h"
#include "integrate.h"          /* for HALF_EPSILON */
#include "zygotic.h"            /* for derivative function rules */
#include "fly_io.h"             /* i/o of parameters and data */

#include "ioTools.h"

double custom_gast;
double maxconc;
int olddivstyle;
/*** MITOSIS SCHEDULE: hard-wired cell division tables *********************
 *                                                                         *
 * From Foe & Alberts, '83 (Table 1, TOTAL ELAPSED added by JR and JJ):    *
 *                                                                         *
 * The authors observed anterior tips of growing embryos during cycles 10  *
 * to 14 and measured the time of somatic bud cycles and recorded the      *
 * visibility of the nuclei, which indicate the duration of mitoses.       * 
 *                                                                         *
 * TOTAL ELAPSED is the total time which elapsed since the start of a sim- *
 * ulation. I.e. we assume                                                 *
 *                                                                         *
 * t=0 at end of cell cycle 10 for NDIVS=3                                 *
 * t=0 1 min into interphase of cell cycle 10 for NDIVS=4 (see (*) below)  *
 *                                                                         *
 * Times are given in minutes, standard devs are between parentheses)      *
 *                                                                         *
 * CYCLE  TOTAL DURATION      NO NUCLEUS            TOTAL ELAPSED          *
 *                                              NDIVS=3       NDIVS=4      *
 *                                                                         *
 *  10 (*)  7.8 (0.6)          3.3 (0.9)                        0.0        *
 *  11      9.5 (0.7)          3.0 (0.9)          0.0           7.8        *
 *  12     12.4 (0.9)          3.3 (0.9)          9.5          17.3        * 
 *  13     21.1 (1.5)          5.1 (0.9)         21.9          29.7        *
 *     	             	       				                   *
 *  14        50+          mitosis post gast     43.0          50.8        *
 *                                                                         *
 *  gastrulation time:                           93.0         100.8        *
 *                                                                         *
 * (*) NOTE: cell cycle 10 actually takes 8.8 minutes, but the migration   *
 *           of nuclei to the cell surface is only completed 1 min into    *
 *           interphase of cell cycle 10. This is our simulation starting  *
 *           point for NDIVS=4.                                            * 
 *                                                                         *
 * The following arrays are the hard-wired stuff about times and durations *
 * of cell divisions and time of gastrulation; each function in maternal.c *
 * that returns information about timing in the blastoderm needs to choose *
 * the appropriate tables depending on the problem, i.e. the mitosis sche- *
 * dule used.                                                              *
 *                                                                         *
 * Oldstyle division times were a coarse approximation based on the 4 min  *
 * time units intially used by JR and DS. They were refined into 4/3 min   *
 * units in pre 9.2 code.                                                  *
 *                                                                         *
 * IMPORTANT: Post 9.2 code does NOT require rounded division times any-   *
 * more. These old division times are only included for backward compati-  *
 * bility and should not be used for annealing to new data!                *
 *                                                                         *
 * NOTE: The times must be in reverse oredr!!!                             *
 * ALSO NOTE: divtimes are at the END of the cell division!!!              *
 * LAST NOTE: there are 3 arrays each, which are (in order of appearance): *
 *            - oldstyle (see above)                                       *
 *            - 3-celldivision schedule (starting at cycle 11)             *
 *            - 4-celldivision schedule (starting at cycle 10)             *
 *                                                                         *
 * JJ (04/04/02):                                                          *
 *            I've added a 'zero'-division schedule which we will use to   *
 *            check if we can get patterns without cell divisions; since   *
 *            this division schedule doesn't have any divisions, it also   *
 *            doesn't need any division tables below                       *
 *                                                                         *
 ***************************************************************************/


/* division times: at end of mitosis! */

const double old_divtimes[3] = { 52.0, 28.0, 12.0 };

const double divtimes1[1] = { 21.1 };                   //NDIVS = 1
const double divtimes2[2] = { 33.5, 12.4 };             //NDIVS = 2   
const double divtimes3[3] = { 43.0, 21.9, 9.5 };        //NDIVS = 3
const double divtimes4[4] = { 50.8, 29.7, 17.3, 7.8 };  //NDIVS = 4   

/* division durations */

static const double old_div_duration[3] = { 4.0, 4.0, 4.0 };

const double div_duration1[1] = { 5.1 };
const double div_duration2[2] = { 5.1, 3.3 };
const double div_duration3[3] = { 5.1, 3.3, 3.0 };
const double div_duration4[4] = { 5.1, 3.3, 3.0, 3.3 };

/* gastrulation times */

const double old_gast_time = 88.;

const double gast_time0 = 50.;
const double gast_time1 = 71.1;
const double gast_time2 = 83.5;
const double gast_time3 = 93.;
const double gast_time4 = 100.8;

/* full division times: including t<0 */

static const double full_divtimes0[TOTAL_DIVS] = { 0.0, -21.1, -33.5, -43.0, -51.8, -57.8 };
static const double full_divtimes1[TOTAL_DIVS] = { 21.1, 0.0, -12.4, -21.9, -30.7, -36.7 };
static const double full_divtimes2[TOTAL_DIVS] = { 33.5, 12.4, 0.0, -9.5, -18.3, -24.3 };
static const double full_divtimes3[TOTAL_DIVS] = { 43.0, 21.9, 9.5, 0.0, -8.8, -14.8 };
static const double full_divtimes4[TOTAL_DIVS] = { 50.8, 29.7, 17.3, 7.8, -1.0, -7.0 };

/* division durations */

static const double full_div_durations[TOTAL_DIVS] = { 5.1, 3.3, 3.0, 3.3, 3.0, 3.0 };





/*** STATIC VARIABLES ******************************************************/


/* following contains lineage numbers at which nuclei at each ccycle start */

//static int         *lin_start;
//int              *full_lin_start=NULL;


//static int         bt_init_flag = 0;                   /* flag for BTtable */
//static int         d_flag       = 0;        /* flag for first call to GetD */
//static int         rule_flag    = 0;     /* flag for first call to GetRule */
static int theta_flag = 0;      /* flag for first call to
                                   Theta */




/*** CONSTANTS *************************************************************/

/* The following is the maximum stepsize for solvers (now set to the whole */
/* time it takes to model to run (100.8 minutes))                          */

const double MAX_STEPSIZE = 100.8;

/* masks for lineage numbers */

const int CYCLE1 = 1;
const int CYCLE2 = 2;
const int CYCLE3 = 4;
const int CYCLE4 = 8;
const int CYCLE5 = 16;
const int CYCLE6 = 32;
const int CYCLE7 = 64;
const int CYCLE8 = 128;
const int CYCLE9 = 256;
const int CYCLE10 = 512;
const int CYCLE11 = 1024;
const int CYCLE12 = 2048;
const int CYCLE13 = 4096;
const int CYCLE14 = 8192;





/*** INITIALIZATION FUNCTIONS **********************************************/

/**  InitTimes: Reading hardcoded times in case that we don't find the
 *              times section in the input file                       
 */
Times
InitTimes( TheProblem defs ) {
    Times times;

    times.total_divs = TOTAL_DIVS;
    times.gast_time = GetGastTime( defs.ndivs );
    times.div_times = GetDivtable( defs.ndivs );
    times.div_duration = GetDurations( defs.ndivs );
    
    times.full_div_durations = ( double * ) full_div_durations;
    times.full_div_times = NULL;

    if( defs.ndivs == 0 )
        times.full_div_times = ( double * ) full_divtimes0;
    else if( defs.ndivs == 1 )
        times.full_div_times = ( double * ) full_divtimes1;
    else if( defs.ndivs == 2 )
        times.full_div_times = ( double * ) full_divtimes2;
    else if( defs.ndivs == 3 )
        times.full_div_times = ( double * ) full_divtimes3;
    else if( defs.ndivs == 4 )
        times.full_div_times = ( double * ) full_divtimes4;
    else
        error( "InitTimes: can't handle %d cell divisions!", defs.ndivs );

    return times;
}

/**  InitBicoid: Copies the Blist read by ReadBicoid into the DArrPtr 
 *               structure; the bicoid DArrPtr contains pointers to bcd    
 *               arrays for each cell division cycle                       
 */
GenoType *
InitBicoid( FILE * fp, Zygote * zyg ) {
    int i;                      /* local loop counter */

    Blist *inlist;              /* temporary linked list for bcd */

    Slist *genotypes;           /* temporary linked list for geno- */
    Slist *current;             /* types from data file */

    GenoType *bcdtype;

    genotypes = ReadGenotypes( fp, zyg->defs.ngenes );
    if( zyg->nalleles == 0 ) {
        zyg->nalleles = count_Slist( genotypes );
    }

    if( !( bcdtype = ( GenoType * ) calloc( zyg->nalleles, sizeof( GenoType ) ) ) )
        error( "InitBicoid: Could not allocate bcdtype struct" );

    /*** for loop: read bicoid for each genotype *******************************/

    for( current = genotypes, i = 0; current; current = current->next, i++ ) {

        if( !( inlist = ReadBicoid( fp, current->bcd_section ) ) )      /* read bicoid */
            error( "InitBicoid: error reading %s", current->bcd_section );
        else {

            if( !( bcdtype[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
                error( "InitBicoid: could not allocate bcd genotype string" );
            bcdtype[i].genotype = strcpy( bcdtype[i].genotype, current->genotype );
            bcdtype[i].ptr.bicoid = List2Bicoid( inlist, zyg );
            free_Blist( inlist );
        }
    }
    free_Slist( genotypes );
    return bcdtype;
}

/**  InitBias:  puts bias records in a form where get_bias can use them; 
 *              it expects times in increasing order; it expects a non-    
 *              sparse entry, with no genes or nuclei missing              
 */
Bias
InitBias( FILE * fp, Zygote * zyg ) {
    /* The following two store bicoid gradients and bias static to maternal.c  */
    /* bias is found here because it contains the maternal contributions to    */
    /* some zygotic genes (such as hunchback), but it can also be used to add  */
    /* heatshocks during the simulation                                        */

    Bias bias;

    /* *bt is the bias times for each genotype */

    /* these two static structs are used for storing the times for which       */
    /* there is bias; these times can be retrieved by using GetBTimes          */




    int i;                      /* loop counters */

    Dlist *inlist;              /* temporary linked list for bias */

    Slist *genotypes;           /* temporary linked list for geno- */
    Slist *current;             /* types from data file */

    int ndp = 0;                /* dummy for ReadData, no need to count datapts here */

    genotypes = ReadGenotypes( fp, zyg->defs.ngenes );
    if( zyg->nalleles == 0 )
        zyg->nalleles = count_Slist( genotypes );

    if( !( bias.biastype = ( GenoType * ) calloc( zyg->nalleles, sizeof( GenoType ) ) ) )
        error( "InitBias: Could not allocate biastype struct" );

    /*** for loop: read bicoid for each genotype *******************************/

    for( current = genotypes, i = 0; current; current = current->next, i++ ) {

        if( !( inlist = ReadData( fp, current->bias_section, &ndp, &( zyg->defs ) ) ) ) /* read bias */
            error( "InitBias: error reading %s", current->bias_section );
        else {

            if( !( bias.biastype[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
                error( "InitBias: could not allocate bias genotype string" );
            bias.biastype[i].genotype = strcpy( bias.biastype[i].genotype, current->genotype );
            bias.biastype[i].ptr.bias = List2Bias( inlist, zyg->defs );

            free_Dlist( inlist );
        }
    }
    bias.bt = InitBTs( bias.biastype, zyg->nalleles );  /* initialize the bias time struct */
    free_Slist( genotypes );
    return bias;
}

/**  InitBTs: initializes the BT struct that holds all times for 
 *            which we have bias.                                          
 */
GenoType *
InitBTs( GenoType * biastype, int nalleles ) {
    int i, j;
    GenoType *bt;


    if( !( bt = ( GenoType * ) calloc( nalleles, sizeof( GenoType ) ) ) )
        error( "InitBTs: could not allocate bt struct" );

    for( i = 0; i < nalleles; i++ ) {

        if( !( bt[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "InitBTs: could not allocate BT genotype string" );
        bt[i].genotype = strcpy( bt[i].genotype, biastype[i].genotype );
        bt[i].ptr.times.size = biastype[i].ptr.bias.size;
        if( !( bt[i].ptr.times.array = ( double * ) calloc( bt[i].ptr.times.size, sizeof( double ) ) ) )
            error( "InitBTs: could not allocate bt array" );

        for( j = 0; j < bt[i].ptr.times.size; j++ )
            bt[i].ptr.times.array[j] = biastype[i].ptr.bias.array[j].time;
    }

    /*  bt_init_flag = 1; */

    return bt;
}

/**  InitNNucs: takes the global defs.nnucs and calculates number of nucs 
 *              for each cleavage cycle which are then stored in reverse   
 *              order in the static nnucs[] array                          
 *   CAUTION:   defs struct and lin_start need to be initialized before!   
 */
int *
InitNNucs( Zygote * zyg ) {
    /* following is number of nucs in each cleavage cycle, reverse order */

    int *nnucs;
    int i;                      /* loop counter */
    int n;                      /* used to calculate number of nucs for each cell cycle */

    if( !( nnucs = ( int * ) calloc( zyg->defs.ndivs + 1, sizeof( int ) ) ) )
        error( "InitNNucs: could not allocate nnucs array" );

    n = zyg->defs.nnucs;

    /* below we have to take into account two cases: a) most anterior lineage  *
     * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
     * b) most anterior lineage number is even-numbered: just add an additio-  *
     * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
     * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

    for( i = 0; i <= zyg->defs.ndivs; i++ ) {
        nnucs[i] = n;
        if( zyg->lin_start[i] % 2 )
            n = n / 2 + 1;
        else
            n = ( n % 2 ) ? n / 2 + 1 : n / 2;
    }
    return nnucs;
}

/**  InitFullNNucs: takes the global defs.nnucs and calculates number of nucs 
 *              for each cleavage cycle which are then stored in reverse   
 *              order in the static nnucs[] array                          
 *   CAUTION:   defs struct and lin_start need to be initialized before!   
 */
void
InitFullNNucs( Zygote * zyg, int *full_lin_start ) {

    if( zyg->full_nnucs ) {
        return;
    }
    int *full_nnucs = NULL;
    int i;                      /* loop counter */
    int n;                      /* used to calculate number of nucs for each cell cycle */


    if( !( full_nnucs = ( int * ) calloc( zyg->defs.full_ccycles, sizeof( int ) ) ) )
        error( "InitFullNNucs: could not allocate full_nnucs array" );

    n = zyg->defs.nnucs;

    /* below we have to take into account two cases: a) most anterior lineage  *
     * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
     * b) most anterior lineage number is even-numbered: just add an additio-  *
     * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
     * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

    for( i = 0; i < zyg->defs.full_ccycles; i++ ) {
        full_nnucs[i] = n;
        if( full_lin_start[i] % 2 )
            n = n / 2 + 1;
        else
            n = ( n % 2 ) ? n / 2 + 1 : n / 2;
    }

    /*  for (i=0; i < zyg->defs.full_ccycles; i++)
       printf("History lineages %d, nnucs %d\n", full_lin_start[i],full_nnucs[i]); */

    zyg->full_nnucs = full_nnucs;
}




/*** FUNCTIONS THAT RETURN INFO ABOUT THE EMBRYO **************************/

/**  GetBicoid: returns a bicoid gradients (in form of a DArrPtr) for a 
 *              specific time and genotype.                                
 */
DArrPtr
GetBicoid( double time, int genindex, GenoType * bcdtype, Zygote * zyg ) {
    int i;
    unsigned int ccycle;

    if( genindex < 0 || genindex >= zyg->nalleles )
        error( "GetBicoid: invalid genotype index %d", genindex );

    ccycle = GetCCycle( time, zyg->defs.ndivs, &( zyg->times ) );

    for( i = 0; i <= bcdtype[genindex].ptr.bicoid.size; i++ )
        if( ccycle == bcdtype[genindex].ptr.bicoid.array[i].ccycle )
            return bcdtype[genindex].ptr.bicoid.array[i].gradient;

    error( "GetBicoid: no bicoid gradient for ccycle %d", ccycle );

    return bcdtype[genindex].ptr.bicoid.array[i].gradient;
    /* just to make the compiler happy! */
}

/** GetBias: This function returns bias values for a given time and 
 *            genotype.                                                    
 */
DArrPtr
GetBias( double time, int genindex, Zygote * zyg ) {
    int i;

    if( genindex < 0 || genindex >= zyg->nalleles )
        error( "GetBias: invalid genotype index %d", genindex );

    for( i = 0; i < zyg->bias.biastype[genindex].ptr.bias.size; i++ ) {
        if( zyg->bias.biastype[genindex].ptr.bias.array[i].time == time )
            break;
    }
    return zyg->bias.biastype[genindex].ptr.bias.array[i].state;
}

/** GetBTimes: returns a sized array of times for which there is bias */
DArrPtr
GetBTimes( char *genotype, Zygote * zyg ) {
    int index;                  /* loop counter */
    for( index = 0; index < zyg->nalleles; index++ ) {
        if( !( strcmp( zyg->bias.biastype[index].genotype, genotype ) ) ) {
            break;
        }
    }
    /* if no explicit bias times for this genotype -> use wt bias times */

    if( index == zyg->nalleles )
        index = 0;

    /* check if we actually have biastimes at all or otherwise -> error! */
    return zyg->bias.bt[index].ptr.times;
}

/** GetNNucs: returns number of nuclei for a given time */
int
GetNNucs( TheProblem * defs, int *nnucs, double t, Times * times ) {
    int i;                      /* loop counter */
    double *table;              /* local copy of divtimes table */

    /* assign 'table' to the appropriate division schedule */

    if( olddivstyle ) {
        if( defs->ndivs != 3 )
            error( "GetNNucs: only 3 cell divisions allowed for oldstyle (-o)" );
        table = ( double * ) old_divtimes;
    } else if( defs->ndivs == 0 )
        return nnucs[0];
    else
        table = ( double * ) times->div_times;

    /*if ( defs->ndivs == 1 )
       table = (double *)divtimes1;
       else if ( defs->ndivs == 2 )
       table = (double *)divtimes2;
       else if ( defs->ndivs == 3 )
       table = (double *)divtimes3;
       else if ( defs->ndivs == 4 )
       table = (double *)divtimes4;
       else
       error("GetNNucs: can't handle %d cell divisions!", defs->ndivs); */

    /* evaluate nnucs for current time; note that for the *exact* time of cell *
     * division, we'll return the number of nuclei before the division has ac- *
     * tually occurred                                                         */

    for( i = 0; i < defs->ndivs; i++ )
        if( t > table[i] )
            return nnucs[i];

    return nnucs[i];
}

/** GetStartLin: returns the lineage number of the most anterior nucleus 
 *                for a given time                                         
 */
int
GetStartLin( double t, TheProblem defs, int *lin_start, Times * times ) {
    int i;                      /* loop counter */
    double *table;              /* local copy of divtimes table */

    /* assign 'table' to the appropriate division schedule */

    if( olddivstyle ) {
        if( defs.ndivs != 3 )
            error( "GetStartLin: only 3 cell divisions allowed for oldstyle (-o)" );
        table = ( double * ) old_divtimes;
    } else if( defs.ndivs == 0 )
        return lin_start[0];
    else
        table = ( double * ) times->div_times;

    /*if ( defs.ndivs == 1 )
       table = (double *)divtimes1;
       else if ( defs.ndivs == 2 )
       table = (double *)divtimes2; 
       else if ( defs.ndivs == 3 )
       table = (double *)divtimes3;
       else if ( defs.ndivs == 4 ) 
       table = (double *)divtimes4;
       else
       error("GetStartLin: can't handle %d cell divisions!", defs.ndivs); */

    /* evaluate lineage number of most anterior nucleus for current time; note *
     * that for the *exact* time of cell division, we'll return the lineage    *
     * number of the most anterior nucleus of the previous cell cycle          */

    for( i = 0; i < defs.ndivs; i++ )
        if( t > table[i] )
            return lin_start[i];

    return lin_start[i];
}

/** Index2StartLin: get starting lineage from index */
int
Index2StartLin( int index, int *full_lin_start ) {
    return full_lin_start[index];
}

/** Index2NNuc: get number of nuclei from index */
int
Index2NNuc( int index, int *full_nnucs ) {
    return full_nnucs[index];
}

/**  GetStartLinIndex: returns the index of the lineage array 
 *                 for a given time                                         
 */
int
GetStartLinIndex( double t, TheProblem * defs, Times * times ) {
    int i;                      /* loop counter */
    double *table;              /* local copy of divtimes table */

    /* assign 'table' to the appropriate division schedule */

    if( olddivstyle ) {
        if( defs->ndivs != 3 )
            error( "GetStartLin: only 3 cell divisions allowed for oldstyle (-o)" );
        table = ( double * ) old_divtimes;
    } else {
        /*if ( defs->ndivs == 0 )
           table = (double *)full_divtimes0;
           else if ( defs->ndivs == 1 )
           table = (double *)full_divtimes1;
           else if ( defs->ndivs == 2 )
           table = (double *)full_divtimes2; 
           else if ( defs->ndivs == 3 )
           table = (double *)full_divtimes3;
           else if ( defs->ndivs == 4 )
           table = (double *)full_divtimes4;
           else
           error("GetStartLin: can't handle %d cell divisions!", defs->ndivs); */

        table = ( double * ) times->full_div_times;

        //test
        //for (i=0; i<defs->full_ccycles - 1; i++)
        //printf("TIMES %d = %lg\n", i, times->full_div_times[i]); 
        //test

    }
    /* evaluate lineage number of most anterior nucleus for current time; note *
     * that for the *exact* time of cell division, we'll return the lineage    *
     * number of the most anterior nucleus of the previous cell cycle          */

    for( i = 0; i < defs->full_ccycles - 1; i++ )
        //printf("TABLE%d=%lg !!!!!\n", i, table[i]);
        if( t > table[i] + HALF_EPSILON ) {
            //printf("%lg > %lg, LININDEX=%d\n", t, table[i], i);
            return i;
        }
    //printf("FINALTABLE%d=%lg !!!!!\n", i, table[i]);
    //printf("FINALLININDEX=%d !!!!!\n", i);
    return i;
}

/** GetCCycle: returns cleavage cycle number for a given time */
unsigned int
GetCCycle( double time, int ndivs, Times * times ) {
    int i;                      /* loop counter */
    double *table;              /* local copy of divtimes table */

    /* assign 'table' to the appropriate division schedule */

    if( olddivstyle ) {
        if( ndivs != 3 )
            error( "GetCCycle: only 3 cell divisions allowed for oldstyle (-o)" );
        table = ( double * ) old_divtimes;
    } else if( ndivs == 0 )
        return 14;              /* if defs.ndivs == 0: we're always in cycle 14 */
    else
        table = ( double * ) times->div_times;

    /* evaluate number of cell cycle for current time; note that for the exact *
     * time of cell division, we'll return the number of the previous cell cy- *
     * cyle                                                                    */

    for( i = 0; i < ndivs; i++ )
        if( time > table[i] )
            return 14 - i;

    return 14 - ( i++ );
}

/** ParseLineage: takes lineage number as input and returns the cleavage 
 *                 cycle the nucleus belongs to.                           
 */
unsigned int
ParseLineage( unsigned int lin ) {
    if( lin & CYCLE14 )
        return 14;
    else if( lin & CYCLE13 )
        return 13;
    else if( lin & CYCLE12 )
        return 12;
    else if( lin & CYCLE11 )
        return 11;
    else if( lin & CYCLE10 )
        return 10;
    else
        error( "ParseLineage: illegal lineage number %d", lin );

    return 0;                   /* just to make the compiler happy! */
}

/** GetDivtable: returns times of cell divisions depending on ndivs and 
 *                olddivstyle; returns NULL in case of an error            
 */
double *
GetDivtable( int ndivs ) {
    /*if ( olddivstyle ) {
       if ( ndivs != 3 ) 
       error("GetDivtable: only 3 cell divisions allowed for oldstyle (-o)");
       return (double *)old_divtimes;
       } else */
    if( ndivs == 0 )
        return NULL;            /* return error if ndivs == 0 */
    else if( ndivs == 1 )
        return ( double * ) divtimes1;
    else if( ndivs == 2 )
        return ( double * ) divtimes2;
    else if( ndivs == 3 )
        return ( double * ) divtimes3;
    else if( ndivs == 4 )
        return ( double * ) divtimes4;
    else
        error( "GetDivtable: can't handle %d cell divisions!", ndivs );

    return NULL;
}

/** GetDurations: returns pointer to durations of cell divisions de- 
 *                 pending on ndivs and olddivstyle; returns NULL in case  
 *                 of an error                                             
 */
double *
GetDurations( int ndivs ) {
    /*if ( olddivstyle ) {
       if ( ndivs != 3 ) 
       error("GetDurations: only 3 cell divisions allowed for oldstyle (-o)");
       return (double *)old_div_duration;
       } else */
    if( ndivs == 0 )
        return NULL;            /* return error if ndivs == 0 */
    else if( ndivs == 1 )
        return ( double * ) div_duration1;
    else if( ndivs == 2 )
        return ( double * ) div_duration2;
    else if( ndivs == 3 )
        return ( double * ) div_duration3;
    else if( ndivs == 4 )
        return ( double * ) div_duration4;
    else
        error( "GetDurations: can't handle %d cell divisions!", ndivs );

    return NULL;
}

/** GetGastTime: returns time of gastrulation depending on ndivs and 
 *                olddivstyle; returns 0 in case of an error; if a custom  
 *                gastrulation time is chosen with -S, it will be returned 
 *                only if it's bigger than the normal gastrulation time    
 */
double
GetGastTime( int ndivs ) {
    if( olddivstyle ) {
        if( ndivs != 3 )
            error( "GetDurations: only 3 cell divisions allowed for oldstyle (-o)" );

        if( custom_gast > old_gast_time )
            return custom_gast;
        else
            return old_gast_time;

    } else
     if( ndivs == 0 )

        if( custom_gast > gast_time0 )
            return custom_gast;
        else
            return gast_time0;

    else if( ndivs == 1 )

        if( custom_gast > gast_time1 )
            return custom_gast;
        else
            return gast_time1;

    else if( ndivs == 2 )

        if( custom_gast > gast_time2 )
            return custom_gast;
        else
            return gast_time2;

    else if( ndivs == 3 )

        if( custom_gast > gast_time3 )
            return custom_gast;
        else
            return gast_time3;

    else if( ndivs == 4 )

        if( custom_gast > gast_time4 )
            return custom_gast;
        else
            return gast_time4;

    else
        error( "GetGastTime: can't handle %d cell divisions!", ndivs );

    return 0;
}

/** GetD: returns diffusion parameters D according to the diff. params. 
 *         in the data file and the diffusion schedule used.               
 *   NOTE: Caller must allocate D_tab                                      
 */
void
GetD( double t, double *d, double *D_tab, Zygote * zyg ) {
    double *table;              /* local copy of divtimes */

    int i;                      /* loop counter */
    double gast;                /* gastrulation time */
    double cutoff;              /* cutoff time */
    double lscale = 1;          /* scaling factor */

    /* first time GetD is called: set pointer to the right division table */
    if( olddivstyle ) {
        if( zyg->defs.ndivs != 3 )
            error( "GetD: only 3 cell divisions allowed for oldstyle (-o)" );
        table = ( double * ) old_divtimes;
    } else if( zyg->defs.ndivs == 0 )
        table = NULL;
    /* no need for table, if there are no cell divs */
    else {
        /*
           }if ( defs.ndivs == 1 ) {
           table = (double *)divtimes1;
           }
           else if ( defs.ndivs == 2 ) 
           table = (double *)divtimes2;
           else if ( defs.ndivs == 3 ) 
           table = (double *)divtimes3;
           else if ( defs.ndivs == 4 ) 
           table = (double *)divtimes4;
           else 
           error("GetD: can't handle %d cell divisions!", defs.ndivs); */
        table = ( double * ) ( zyg->times.div_times );
    }


    /* this loop takes lscale square for each cell division, i.e. the earlier  *
     * we are the bigger lscale (and the smaller the D's that we return        */
    for( i = 0; i < zyg->defs.ndivs; i++ )
        if( t <= table[i] )
            lscale *= 2;
    /* diffusion schedule A: all Ds always the same */

    if( zyg->defs.diff_schedule == 'A' )
        for( i = 0; i < zyg->defs.ngenes; i++ )
            D_tab[i] = d[0];

    /* diffusion schedule B: all genes have different D's that depend on in-   *
     * verse l-square                                                          */

    else if( zyg->defs.diff_schedule == 'B' ) {

        for( i = 0; i < zyg->defs.ngenes; i++ ) {
            D_tab[i] = d[i] / ( lscale * lscale );
        }
    }
    /* diffusion schedule C: all genes have the same D that depends on inverse *
     * l-square                                                                */
    else if( zyg->defs.diff_schedule == 'C' )
        for( i = 0; i < zyg->defs.ngenes; i++ )
            D_tab[i] = d[0] / ( lscale * lscale );

    /* diffusion schedule D: all genes have different D's which don't change   *
     * over time                                                               */

    else if( zyg->defs.diff_schedule == 'D' )
        for( i = 0; i < zyg->defs.ngenes; i++ )
            D_tab[i] = d[i];

    /* diffusion schedule E: used cutoff at gast-12 otherwise just like B */

    else if( zyg->defs.diff_schedule == 'E' ) {
        if( olddivstyle ) {
            if( zyg->defs.ndivs != 3 )
                error( "GetD: only 3 cell divisions allowed for oldstyle (-o)" );
            gast = old_gast_time;
        } else {

            gast = zyg->times.gast_time;

        }

        cutoff = gast - 12.0;   /* This value probably wrong; see Merrill88 */
        for( i = 0; i < zyg->defs.ngenes; i++ )
            D_tab[i] = ( t < cutoff ) ? d[i] / ( lscale * lscale ) : 0.;

        /* any other diffusion schedule: error! */
    } else
        error( "GetD: no code for schedule %c!", zyg->defs.diff_schedule );
}

/** Theta: Returns the value of theta(t) in the autonomous 
 *            version of the equations                              
 */
int
Theta( double time, Zygote * zyg ) {
    int i, ndivs;

    static double *dt;          /* pointer to division time table */
    static double *dd;          /* pointer to division duration table */

    ndivs = zyg->defs.ndivs;

    if( !theta_flag ) {         /* only do this once */
        if( olddivstyle ) {     /* get pointers to division time table */
            dt = ( double * ) old_divtimes;     /* and division duration table */
            dd = ( double * ) old_div_duration;
        } else {
            dt = ( double * ) zyg->times.full_div_times;
            dd = ( double * ) zyg->times.full_div_durations;

        }
        theta_flag = 1;
    }

    /* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
     * on Linux which can't handle truncation errors very well                 */
    //printf("totaldivs=%d HALF_EPSILON=%lg\n", zyg->times.total_divs, HALF_EPSILON);
    for( i = 0; i < zyg->times.total_divs; i++ ) {
        //printf("Thetatime=%.16f,[%.16f,%.16f]\n", time, (*(dt+i) - *(dd+i) + HALF_EPSILON), (*(dt+i) + HALF_EPSILON)); 
        if( ( time <= ( *( dt + i ) + HALF_EPSILON ) )
            && ( time >= ( *( dt + i ) - *( dd + i ) + HALF_EPSILON ) ) )
            return MITOSIS;
    }

    return INTERPHASE;
}

/** GetRule: returns the appropriate rule for a given time; used by the 
 *            derivative function                                          
 */
int
GetRule( double time, Zygote * zyg ) {
    int i;

    double *dt;                 /* pointer to division time table */
    double *dd;                 /* pointer to division duration table */

    if( olddivstyle ) {         /* get pointers to division time table */
        dt = ( double * ) old_divtimes; /* and division duration table */
        dd = ( double * ) old_div_duration;
    } else {
        if( zyg->defs.ndivs == 0 ) {
            return INTERPHASE;  /* no cell division? no MITOSIS! */
        } else {
            dt = ( double * ) zyg->times.div_times;
            dd = ( double * ) zyg->times.div_duration;

        }
    }
    /* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
     * on Linux which can't handle truncation errors very well                 */
    for( i = 0; i < zyg->defs.ndivs; i++ ) {
        if( ( time <= ( *( dt + i ) + HALF_EPSILON ) ) && ( time >= ( *( dt + i ) - *( dd + i ) + HALF_EPSILON ) ) ) {
            return MITOSIS;
        }
    }
    return INTERPHASE;
}

/** GetIndex: this functions returns the genotype index for a given 
 *             genotype number for reading the GenoType struct.            
 */
int
GetIndex( char *genotype, GenoType * bcdtype, int nalleles ) {
    int i;                      /* local loop counter */

    for( i = 0; i < nalleles; i++ )     /* nalleles static to score.c */
        if( !( strcmp( bcdtype[i].genotype, genotype ) ) )
            return i;

    error( "GetIndex: could not find index for genotype %s", genotype );
    return -1000;
}

/** MakeTable: this function constructs a time table for Blastoderm 
 *             based on a comand line option (print_stepsize) and the      
 *             cell division tables here in maternal.c                     
 * DISCLAIMER: this thing is written in a very bad way; but hey it does    
 *             its job!!!!!!!                                              
 */
DArrPtr
MakeTable( double p_stepsize, Zygote * zyg ) {
    int i;                      // local loop counter
    int t;                      // time counter

    double time;                // double counter
    double gast;                // gastrulation time

    DArrPtr table;              // time table to be returned



    gast = zyg->times.gast_time;        // get gastrulation time

    if( p_stepsize > gast )
        error( "GetTable: output stepsize can't be larger than gast time" );

    t = 0;
    for( time = 0; time < gast; time += p_stepsize ) {
        t++;
    }

    table.size = t + 1;         // add one for gastrulation time itself
    if( !( table.array = ( double * ) calloc( table.size, sizeof( double ) ) ) )
        error( "GetTable: error allocating times array" );

    time = 0;
    for( i = 0; i < table.size - 1; i++ ) {
        table.array[i] = time;
        time += p_stepsize;
    }
    table.array[i] = gast;

    return table;
}

/** List2Bicoid: takes a Blist and returns the corresponding BArrPtr 
 *                structure; also initializes the lin_start array which    
 *                contains the lineage number of the most anterior nucleus 
 *                for each cell cycle (we need this for cell division and  
 *                printing model output)                                   
 */
BArrPtr
List2Bicoid( Blist * inlist, Zygote * zyg ) {
    int i = 0;                  /* local loop counter */
    int n;                      /* used to evaluate # of nuclei */
    int lin_count = 0;          /* counter for initializing lin_start */

    unsigned int ccycle;        /* what cleavage cycle are we in? */
    unsigned int samecycle = 0; /* same cleavage cycle as before? */

    Blist *current;             /* current element of Blist */
    Blist *start;               /* used to evaluate # of nuclei */

    BArrPtr bicoid;             /* BArrPtr to be returned */

    bicoid.size = 0;
    bicoid.array = NULL;

    if( !( zyg->lin_start = ( int * ) calloc( zyg->defs.ndivs + 1, sizeof( int ) ) ) )
        error( "List2Bicoid: could not allocate lin_start array" );

    /*** for loop: step through linked list and copies values into an array    *
     *             of BcdGrads; there's one BcdGrad for each cleavage cycle;   *
     *             each BcdGrad has: - ccycle (the cleavage cycle number)      *
     *                               - gradient.size (# of nuclei for grad.)   *
     *                               - gradient.array (pointer to array)       *
     ***************************************************************************/

    for( current = inlist; current; current = current->next ) {

        ccycle = ParseLineage( current->lineage );      /* which cycle is it? */

        /* allocate new gradient struct for new cell cycle */

        if( ccycle != samecycle ) {

            samecycle = ccycle;

            bicoid.size++;
            bicoid.array = ( BcdGrad * ) realloc( bicoid.array, bicoid.size * sizeof( BcdGrad ) );

            /* this loop determines the number of nuclei in each cycle */

            n = 0;
            for( start = current; ParseLineage( current->lineage ) == samecycle; current = current->next ) {
                n++;
                if( !( current->next ) )        /* don't count garbage -> core dumps... */
                    break;
            }
            current = start;    /* reset pointer to where we were */

            /* initialize lin_start: this array is later used by Blastoderm and such   */

            zyg->lin_start[zyg->defs.ndivs - lin_count] = current->lineage;
            lin_count++;

            /* allocate array for gradient here and reset the counter */

            bicoid.array[bicoid.size - 1].ccycle = samecycle;   /* next three */
            /* lines define BcdGrad for each cleavage cycle */
            bicoid.array[bicoid.size - 1].gradient.array = ( double * ) calloc( n, sizeof( double ) );
            bicoid.array[bicoid.size - 1].gradient.size = n;
            i = 0;

        }

        /* in any case: read concentration into gradient array */

        bicoid.array[bicoid.size - 1].gradient.array[i] = current->conc;
        i++;
    }

    return bicoid;
}

/** List2Bias: takes a Dlist and returns the corresponding DArrPtr 
 *              structure.                                                 
 */
NArrPtr
List2Bias( Dlist * inlist, TheProblem defs ) {
    int i = 0;
    int j;                      /* local loop counters */

    int n;                      /* used to evaluate # of nuclei */
    double now = -999999999.0;  /* variable for time */

    Dlist *current;             /* holds current element of Dlist */
    Dlist *start;               /* pointer used for evaluating # of */

    NArrPtr bias;               /* NArrPtr to be returned */

    bias.size = 0;
    bias.array = NULL;

    /*** for loop: steps through linked list and copies values into an array   *
     *             of NucStates. There's one NucState for each time step       *
     *             Each NucState has: - time (the time)                        *
     *                                - state.size (# of genes * # of nucs)    *
     *                                - state.array (pointer to array)         *
     ***************************************************************************/

    for( current = inlist; current; current = current->next ) {

        if( current->d[0] != now ) {    /* a new time point: allocate */
            now = current->d[0];        /* the time is now! */
            bias.size++;        /* one NucState for each time step */
            bias.array =        /* allocate for one more NucState */
                ( NucState * ) realloc( bias.array, bias.size * sizeof( NucState ) );

            /* determine number of nuclei per time step */

            n = 0;
            for( start = current; current->d[0] == now; current = current->next ) {
                n++;
                if( !( current->next ) )        /* don't count garbage and cause dis- */
                    //wow! that was a very 'braintraining' commentary...
                    break;      /* may and core dumps... */
            }
            current = start;    /* reset list ptr to where we were */

            /* allocate a bias array for each biastime */

            bias.array[bias.size - 1].time = now;
            bias.array[bias.size - 1].state.array = ( double * ) calloc( n * defs.ngenes, sizeof( double ) );
            bias.array[bias.size - 1].state.size = n * defs.ngenes;

            i = 0;
        }

        /* always: read concs into array */

        for( j = 1; j <= defs.ngenes; j++ ) {
            bias.array[bias.size - 1].state.array[i]
                = current->d[j];
            i++;
        }
    }

    return bias;
}

/** List2Interp: takes a Dlist and returns the corresponding  
 *              interpolation DataTable structure.                        
 */
DataTable *
List2Interp( Dlist * inlist, Input * inp, int ngenes ) {


    int i = 0;
    int j;                      /* local loop counters */

    double now = -999999999.;   /* assigns data to specific time */

    Dlist *current;             /* holds current element of Dlist */

    DataTable *D;               /* local copy of DataTable */



    if( !inp->zyg.full_lin_start ) {
        if( !( inp->zyg.full_lin_start = ( int * ) calloc( inp->zyg.defs.ndivs + 1, sizeof( int ) ) ) )
            error( "List2Interp: could not allocate full_lin_start array" );
        inp->zyg.defs.full_ccycles = inp->zyg.defs.ndivs + 1;
        for( i = 0; i < inp->zyg.defs.full_ccycles; i++ )
            inp->zyg.full_lin_start[i] = inp->zyg.lin_start[i];
    }

    D = ( DataTable * ) malloc( sizeof( DataTable ) );
    /* Initialize DataTable structure */
    D->size = 0;
    D->record = NULL;

    /*** for loop: steps through linked list and transfers facts into Data-    *
     *             Records, one for each time step                             *
     ***************************************************************************/
    for( current = inlist; current; current = current->next ) {

        if( current->d[0] != now ) {    /* a new time point: allocate */
            now = current->d[0];        /* the time is now! */
            //printf("NOW = %lg\n", now);
            D->size++;          /* one DataRecord for each time */
            D->record =         /* allocate DataRecord */
                ( DataRecord * ) realloc( D->record, D->size * sizeof( DataRecord ) );

            D->record[D->size - 1].time = now;  /* next three lines define */
            D->record[D->size - 1].size = 0;    /* DataRecord for each */
            D->record[D->size - 1].array = NULL;        /* time step */
            i = 0;
        }
        for( j = 1; j <= ngenes; j++ ) {        /* always: read concs into array */
            D->record[D->size - 1].size++;      /* one more in this record */
            D->record[D->size - 1].array =      /* reallocate memory for array! */
                realloc( D->record[D->size - 1].array, D->record[D->size - 1].size * sizeof( DataPoint ) );

            /* the following two lines assign concentration value and index *********** */

            D->record[D->size - 1].array[D->record[D->size - 1].size - 1].conc = current->d[j];

            D->record[D->size - 1].array[D->record[D->size - 1].size - 1].index = i;
            i++;

        }
        /* initialize lin_start: this array is later used by Blastoderm and such   */
        if( ParseLineage( inp->zyg.full_lin_start[inp->zyg.defs.full_ccycles - 1] ) != ParseLineage( current->lineage ) ) {
            inp->zyg.defs.full_ccycles++;
            if( !( inp->zyg.full_lin_start = ( int * ) realloc( inp->zyg.full_lin_start, inp->zyg.defs.full_ccycles * sizeof( int ) ) ) )
                error( "List2Interp: could not allocate full_lin_start array" );
            inp->zyg.full_lin_start[inp->zyg.defs.full_ccycles - 1] = current->lineage;
        }

    }
    qsort( ( void * ) inp->zyg.full_lin_start, inp->zyg.defs.full_ccycles, sizeof( int ), ( int ( * )( const void *, const void * ) ) descend );


    /* Now lets remove duplicates */
    i = 0;
    while( i < inp->zyg.defs.full_ccycles - 1 ) {
        if( inp->zyg.full_lin_start[i] == inp->zyg.full_lin_start[i + 1] ) {
            memmove( ( inp->zyg.full_lin_start + i ), ( inp->zyg.full_lin_start + i + 1 ), ( inp->zyg.defs.full_ccycles - i - 1 ) * sizeof( int ) );
            /*          printf("Shifted %d elements to %d\n",inp->zyg.defs.full_ccycles-i-1,i); */
            inp->zyg.defs.full_ccycles--;
            i--;
        }

        i++;
        inp->zyg.full_lin_start = ( int * ) realloc( inp->zyg.full_lin_start, inp->zyg.defs.full_ccycles * sizeof( int ) );
    }
    return D;
}

/** Comparison function for qsort */
int
descend( int *x, int *y ) {

    if( *x < *y )
        return 1;
    else if( *x > *y )
        return -1;
    else
        return 0;

}

/* FOLLOWING FUNCTIONS ARE UTILITY FUNCTIONS FOR DIFFERENT LINKED LISTS ****
 * which are used to read in data of unknown size. All utility functions   *
 * follow the same scheme (X stands for the different letters below):      *
 *                                                                         *
 * - init_Xlist:         allocates first element and returns pointer to it *
 * - addto_Xlist:        adds the adduct to an already existing linkd list *
 * - free_Xlist:         frees memory of linked list                       *
 * - count_Xlist:        counts elements in the linked list                *
 *                                                                         *
 ***************************************************************************/

/*** Utility functions for Blist *******************************************/

/** init_Blist:         allocates first element and returns pointer to it */
Blist *
init_Blist( void ) {
    Blist *p;

    if( ( p = ( Blist * ) malloc( sizeof( Blist ) ) ) ) {

        p->lineage = 0;
        p->conc = 0;
        p->next = NULL;

    } else
        error( "init_Blist: couldn't allocate!" );

    return p;
}

/** addto_Blist:       adds the adduct to an already existing linkd list */
Blist *
addto_Blist( Blist * start, Blist * adduct ) {
    Blist *current;

    if( !start )
        return adduct;

    current = start;

    while( current->next ) {
        current = current->next;
    }

    current->next = adduct;
    return start;
}

/** free_Blist:         frees memory of linked list */
void
free_Blist( Blist * start ) {
    if( start->next )
        free_Blist( start->next );

    free( start );
}

/** count_Blist:        counts elements in the linked list */
int
count_Blist( Blist * start ) {
    int n = 0;

    while( start != NULL ) {
        n++;
        start = start->next;
    }

    return n;
}

/*** Utility functions for Dlist *******************************************/

/** init_Dlist:         allocates first element and returns pointer to it */
Dlist *
init_Dlist( int size ) {
    Dlist *p;

    if( ( p = ( Dlist * ) malloc( sizeof( Dlist ) ) ) ) {

        if( !( p->d = ( double * ) calloc( size, sizeof( double ) ) ) )
            error( "initDlist: couldn't allocate!" );

        p->lineage = 0;
        p->next = NULL;

    } else
        error( "initDlist: couldn't allocate!" );

    return p;
}

/** addto_Dlist:       adds the adduct to an already existing linkd list */
Dlist *
addto_Dlist( Dlist * start, Dlist * adduct ) {
    Dlist *current;

    if( !start )
        return adduct;

    current = start;

    while( current->next ) {
        current = current->next;
    }

    current->next = adduct;
    return start;
}

/** free_Dlist:         frees memory of linked list */
void
free_Dlist( Dlist * start ) {
    if( start->next )
        free_Dlist( start->next );

    free( start->d );
    free( start );
}

/** count_Dlist:        counts elements in the linked list */
int
count_Dlist( Dlist * start ) {
    int n = 0;
    while( start != NULL ) {
        n++;
        start = start->next;
    }

    return n;
}

/*** Utility functions for Slist *******************************************/

/** init_Slist:         allocates first element and returns pointer to it */
Slist *
init_Slist( void ) {
    Slist *p;

    if( ( p = ( Slist * ) malloc( sizeof( Slist ) ) ) ) {

        p->bias_section = NULL;
        p->fact_section = NULL;
        p->bcd_section = NULL;
        p->hist_section = NULL;
        p->ext_section = NULL;
        p->genotype = NULL;
        p->next = NULL;

    } else
        error( "initSlist: couldn't allocate" );

    return p;
}

/** addto_Slist:       adds the adduct to an already existing linkd list */
Slist *
addto_Slist( Slist * start, Slist * adduct ) {
    Slist *current;

    if( !start )
        return adduct;

    current = start;

    while( current->next ) {
        current = current->next;
    }

    current->next = adduct;
    return start;
}

/** free_Slist:         frees memory of linked list */
void
free_Slist( Slist * start ) {
    if( start->next )
        free_Slist( start->next );

    free( start->bias_section );
    free( start->fact_section );
    free( start->bcd_section );
    free( start->hist_section );
    free( start->ext_section );
    free( start->weights_section );
    free( start->genotype );
    free( start );
}

/** count_Slist:        counts elements in the linked list */
int
count_Slist( Slist * start ) {

    int n = 0;
    while( start != NULL ) {
        n++;
        start = start->next;
    }
    return n;
}

double *
Get_Theta_Discons( int *theta_discon_size, Zygote * zyg ) {

    int s, i, j, ndivs;
    double *dt, *dd, *disc_array;

    ndivs = zyg->defs.ndivs;

    s = ( zyg->times.total_divs - ndivs );
    disc_array = ( double * ) calloc( 2 * s, sizeof( double ) );

    dt = ( double * ) zyg->times.full_div_times;
    dd = ( double * ) zyg->times.full_div_durations;

    for( i = 0, j = 0; i < s; i++, j += 2 ) {
        *( disc_array + j ) = dt[ndivs + i];
        *( disc_array + j + 1 ) = dt[ndivs + i] - dd[ndivs + i];
    }

    *theta_discon_size = 2 * s;
    return disc_array;
}
