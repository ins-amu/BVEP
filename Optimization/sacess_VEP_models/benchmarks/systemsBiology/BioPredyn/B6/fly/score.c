/**
 * @file score.c                                               
 * @author JR, modified by Yoginho,
 * scoregut functions by Yousong Wang in June 2002 
 *                                                                  
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Function implementation for reading and defining facts/limits, 
 * and for scoring.
 * 
 * The functions defined here initialize or manipulate facts or 
 * data time tables, read and initialize limits and penalty (if needed) 
 * and do the actual scoring by least squares.                              
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "score.h"              /* obviously */
#include "global.h"

//extern int proc_id;
//extern int debug;

/* some major structs */


FILE *fp;

GutInfo gutparms;               /* Yousong's gutinfo structure used to print guts below                      */

/* ... plus an init flag */

static int tt_init_flag = 0;    /* flag: InitTT called yet or not? */

double best_score = 2000000000; //a very large number

//
//  added by YF for the  residuals computation
//
static NArrPtr resComp2;

static int resC;                /* do we compute the residuals? */
static int nbScore;             /* number of times we ran score */

//the different possible type of objective functions YF
static const int LSE = 0;
static const int MAD = 1;

const int SLEEP_LGTH = 2;
const int NPOINTS = 50;

/*** INITIALIZATION FUNCTIONS **********************************************/

/** InitScoring: intializes a) facts-related structs and TTimes and 
 *                           b) parameter ranges for the Score function.   
 */
Scoring
InitScoring( FILE * fp, int method, Input * inp ) {

    Scoring scoring = ( const struct Scoring ) { 0 };
    scoring.method = method;
    scoring.facts = InitFacts( fp, inp );       /* initializes facts */
    if( method == 0 ) {
        scoring.weights = InitWeights( fp, inp );       /* initializes weights */
    }
    scoring.searchspace = InitLimits( fp, inp );        /* installs search space into score.c */

    return scoring;
}

/*** FUNCTION DEFINITIONS **************************************************/

/** InitTweak: installs tweak as a static variable in translate.c; tweak 
 *              is read from the $tweak section in newstyle data files     
 */
Tweak
InitTweak( FILE * fp, int *mask, TheProblem defs ) {
    Tweak tweak;                /* which parameters to tweak */
    tweak = ReadTweak( fp, mask, defs );
    return tweak;
}

/** InitFacts: puts facts records into the appropriate DataTable. 
 *              Returns a pointer to a Facts struture, which contains a    
 *              Datatable, which in turn points to a sized array of        
 *              DataRecords, one for each time step.                       
 */
Facts
InitFacts( FILE * fp, Input * inp ) {
    int i;                      /* local loop counter */

    Facts facts;

    Dlist *inlist;              /* temporary linked list for facts */

    Slist *genotypes;           /* temporary linked list for geno- */
    Slist *current;             /* types from data file */

    genotypes = ( Slist * ) ReadGenotypes( fp, inp->zyg.defs.ngenes );  /* get the genotypes into an SList */
    if( inp->zyg.nalleles == 0 )
        inp->zyg.nalleles = ( int ) count_Slist( genotypes );

    if( !( facts.facttype = ( GenoType * ) calloc( inp->zyg.nalleles, sizeof( GenoType ) ) ) )
        error( "InitFacts: could not allocate facttype struct" );

    /*** for loop: read the data for each genotype *****************************/

    for( current = genotypes, i = 0; current; current = current->next, i++ ) {

        if( !( inlist = ( Dlist * ) ReadData( fp, current->fact_section, &( inp->zyg.ndp ), &( inp->zyg.defs ) ) ) )
            error( "InitFacts: no Dlist to initialize facts" );
        else {
            if( !( facts.facttype[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
                error( "InitFacts: could not allocate facts genotype string" );
            facts.facttype[i].genotype = strcpy( facts.facttype[i].genotype, current->genotype );
            facts.facttype[i].ptr.facts = List2Facts( inlist, inp->zyg.defs.ngenes );
            free_Dlist( inlist );
        }
    }

    facts.tt = InitTTs( facts.facttype, inp->zyg.nalleles );
    free_Slist( genotypes );
    return facts;
}

/** getFacts: gets facts records from the appropriate DataTable. 
 *              Returns a pointer to a Facts structure, which contains a   
 *              DataTable, which in turn points to a sized array of        
 *              DataRecords, one for each time step.                       
 */
DataTable
getFact( int i, GenoType * facttype ) {
    return *( facttype[i].ptr.facts );
}

/** InitWeights: puts facts records into the appropriate DataTable. 
 *              Returns a pointer to a DataTable, which in turn points to  
 *              a sized array of DataRecords                               
 */
Weights
InitWeights( FILE * fp, Input * inp ) {
    int i;                      //, j; /* local loop counters */
    int ndp = 0;                // dummy
    //double* medians_wt;
    //double* medians;
    //double* multipliers;
    Weights weights = ( const struct Weights ){ 0 };

    Dlist *inlist;              /* temporary linked list for facts */

    Slist *genotypes;           /* temporary linked list for geno- */
    Slist *current;             /* types from data file */

    genotypes = ( Slist * ) ReadGenotypes( fp, inp->zyg.defs.ngenes );  /* get the genotypes into an SList */
    if( inp->zyg.nalleles == 0 )
        inp->zyg.nalleles = ( int ) count_Slist( genotypes );

    if( !( weights.weighttype = ( GenoType * ) calloc( inp->zyg.nalleles, sizeof( GenoType ) ) ) )
        error( "InitWeights: could not allocate struct" );

    /*** for loop: read the data for each genotype *****************************/

    for( current = genotypes, i = 0; current; current = current->next, i++ ) {
        if( !*( current->weights_section ) ) {
            weights.weighttype[i].genotype = NULL;
            weights.weighttype[i].ptr.facts = NULL;
        } else {
            if( !( inlist = ( Dlist * ) ReadData( fp, current->weights_section, &ndp, &( inp->zyg.defs ) ) ) )
                error( "InitWeights: no Dlist to initialize weights" );
            else {
                if( !( weights.weighttype[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
                    error( "InitWeights: could not allocate weights genotype string" );
                weights.weighttype[i].genotype = strcpy( weights.weighttype[i].genotype, current->genotype );
                weights.weighttype[i].ptr.facts = List2Facts( inlist, inp->zyg.defs.ngenes );

                //Here we rescale weights of mutants to their median, in order to get data more uniform
                /*if (i == 0) {
                   medians_wt = getMediansFromDataTable(weights.weighttype[i].ptr.facts, inp->zyg.defs.ngenes);                    
                   } else {                                        
                   medians = getMediansFromDataTable(weights.weighttype[i].ptr.facts, inp->zyg.defs.ngenes);
                   multipliers = (double *) calloc(inp->zyg.defs.ngenes, sizeof(double));
                   for (j=0; j<inp->zyg.defs.ngenes; j++) {
                   multipliers[j] = medians_wt[j]/medians[j];
                   }
                   Rescale(weights.weighttype[i].ptr.facts, multipliers, inp->zyg.defs.ngenes);
                   free(medians);
                   free(multipliers);
                   } */
                free_Dlist( inlist );
            }
        }
    }
    free_Slist( genotypes );
    //free(medians_wt);
    return weights;
}

/** getMediansFromDataTable: Function to calculate medians per gene from data */
double *
getMediansFromDataTable( DataTable * data, int ngenes ) {
    int g, i, j, ind, elements = 0, tempint, count;
    double *array;
    double *medians = ( double * ) malloc( ngenes * sizeof( double ) );

    for( i = 0; i < data->size; i++ ) { //get number of elements per gene in the weight matrix (total number of weights divided by number of genes)
        tempint = ( int ) data->record[i].size / ( double ) ngenes;
        elements += tempint;
    }
    array = ( double * ) malloc( elements * sizeof( double ) );
    for( g = 0; g < ngenes; g++ ) {
        count = 0;
        for( i = 0; i < data->size; i++ ) {
            for( j = 0; j < data->record[i].size; j++ ) {
                if( ( j - g ) % ngenes == 0 ) {
                    array[count] = data->record[i].array[j].conc;
                    count++;
                }
            }
        }
        qsort( ( void * ) array, elements, sizeof( double ), ( int ( * )( const void *, const void * ) ) compare );
        if( elements % 2 == 0 ) {
            ind = ( int ) ( elements / 2 );
            medians[g] = ( array[ind] + array[ind + 1] ) / 2;
        } else {
            ind = ( int ) ( ( elements - 1 ) / 2 ) + 1;
            medians[g] = array[ind];
        }
    }

    free( array );

    return medians;
}

/** Rescale: Function to rescale data - for every gene multiply data per previously calculated multiplier */
void
Rescale( DataTable * data, double *multipliers, int ngenes ) {
    int g, i, j;
    for( i = 0; i < data->size; i++ ) {
        for( j = 0; j < data->record[i].size; j = j + 4 ) {
            for( g = 0; g < ngenes; g++ ) {
                data->record[i].array[j + g].conc *= multipliers[g];
            }
        }
    }
}


/** InitLimits: reads limits section from the data file into the struct 
 *               limits, which is static to score.c. Then, it initializes  
 *               the penalty function if necessary.                        
 *   NOTE:       lambda limits are stored as protein half lives in the     
 *               data file and therefore get converted upon reading        
 */
SearchSpace *
InitLimits( FILE * fp, Input * inp ) {
    SearchSpace *limits;        /* structure that holds the limits for   */
    /* the search space                      */
    //printf( "# init limits\n" );
    limits = ReadLimits( fp, inp->zyg.defs );
    if( limits->pen_vec != NULL ) {
        //printf( "# init penalty\n" );
        InitPenalty( fp, inp->zyg.defs, limits );       /* installs mmax and vmax into penalty vector */
    }
    return limits;
}

/** InitPenalty: initializes vmax[] and mmax static to score.c; these 
 *                variables are used to calculate the penalty function     
 *                for scoring.                                             
 *         NOTE: Should only get called, if penalty function is used       
 */
void
InitPenalty( FILE * fp, TheProblem defs, SearchSpace * limits ) {
    int i, j;                   /* local loop counters */

    char *g_type;               /* genotype string of penalty data to be read */

    char *factpen_section;      /* title of penalty data section */
    char *extpen_section;       /* title of penalty data
                                   section for external inputs */
    char *bcdpen_section;       /* title of maternal penalty data section */

    Slist *s_inlist;            /* linked list for genotypes */
    Slist *s_current;

    Dlist *d_inlist;            /* linked list for reading penalty data */
    Dlist *d_current;

    Blist *b_inlist;            /* linked list for reading maternal */
    Blist *b_current;           /* penalty data */

    double *mmax;               /* max. bcd protein conc. in maternal penalty data */
    double *vmax;               /* array for max. prot concs. in penalty data */

    int ndp = 0;                /* dummy for ReadData, no need to count datapts here */

    int genot_id = -1;

    factpen_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    bcdpen_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    extpen_section = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    g_type = ( char * ) calloc( MAX_RECORD, sizeof( char ) );

    mmax = ( limits->pen_vec ) + 1;     /* mmax stored in limits->pen_vec[1] */
    vmax = ( limits->pen_vec ) + 2;     /* vmax stored in penalty vector array */

    *mmax = -1.;
    for( i = 0; i < defs.ngenes + defs.egenes; i++ )
        vmax[i] = -1.;

    s_inlist = ( Slist * ) ReadGenotypes( fp, defs.ngenes );

    /* loop to read penalty data for all genotypes */

    for( s_current = s_inlist; s_current; s_current = s_current->next ) {
        genot_id++;
        g_type = strcpy( g_type, s_current->genotype );

        /* this part reads the facts data */

        if( !( d_inlist = ( Dlist * ) ReadData( fp, s_current->fact_section, &ndp, &defs ) ) )
            error( "InitPenalty: no Dlist to initialize penalty" );
        else {
            for( d_current = d_inlist; d_current; d_current = d_current->next )
                for( j = 0; j < defs.ngenes; j++ )
                    if( d_current->d[j + 1] > vmax[j] )
                        vmax[j] = d_current->d[j + 1];

            free_Dlist( d_inlist );
        }

        //for every gene take maximum of all facts, and all penalties

        /* read additional penalty sections if necessary */
//printf("vmax1 = %lg %lg %lg %lg %lg %lg %lg %lg\n", vmax[0], vmax[1], vmax[2], vmax[3], vmax[4], vmax[5], vmax[6], vmax[7]);
        sprintf( factpen_section, "penalty_data.%d", genot_id );
        if( FindSection( fp, factpen_section ) ) {
            if( ( d_inlist = ( Dlist * ) ReadData( fp, factpen_section, &ndp, &defs ) ) ) {
                for( d_current = d_inlist; d_current; d_current = d_current->next )
                    for( j = 0; j < defs.ngenes; j++ )
                        if( d_current->d[j + 1] > vmax[j] )
                            vmax[j] = d_current->d[j + 1];
            } else {
                error( "InitPenalty: error reading penalty section for genotype %s", genot_id );
            }

            free_Dlist( d_inlist );
        }

        /* just in case if all values are -1 */

        for( i = 0; i < defs.ngenes; i++ )
            if( vmax[i] == -1. )
                vmax[i] = 255.; /* set max value to 255 */


        //for every gene take maximum of all external imputs, and all external penalties

        /* this part reads the external inputs data */
//printf("vmax2 = %lg %lg %lg %lg %lg %lg %lg %lg\n", vmax[0], vmax[1], vmax[2], vmax[3], vmax[4], vmax[5], vmax[6], vmax[7]);
        if( !( d_inlist = ( Dlist * ) ReadInterpData( fp, s_current->ext_section, defs.egenes, &ndp ) ) )
            error( "InitPenalty: no Dlist to initialize penalty" );
        else {
            for( d_current = d_inlist; d_current; d_current = d_current->next )
                for( j = 0; j < defs.egenes; j++ )
                    if( d_current->d[j + 1] > vmax[j + defs.ngenes] )
                        vmax[j + defs.ngenes] = d_current->d[j + 1];

            free_Dlist( d_inlist );
        }

        /* read additional penalty section if necessary */
//printf("vmax3 = %lg %lg %lg %lg %lg %lg %lg %lg\n", vmax[0], vmax[1], vmax[2], vmax[3], vmax[4], vmax[5], vmax[6], vmax[7]);
        sprintf( extpen_section, "external_penalty_data.%d", genot_id );
        if( FindSection( fp, extpen_section ) ) {
            if( ( d_inlist = ( Dlist * ) ReadInterpData( fp, extpen_section, defs.egenes, &ndp ) ) ) {
                for( d_current = d_inlist; d_current; d_current = d_current->next )
                    for( j = 0; j < defs.egenes; j++ )
                        if( d_current->d[j + 1] > vmax[j + defs.ngenes] )
                            vmax[j + defs.ngenes] = d_current->d[j + 1];
            } else {
                error( "InitPenalty: error reading external penalty section for genotype %d", genot_id );
            }

            free_Dlist( d_inlist );
        }

        /* just in case if all values are -1 */
        for( i = 0; i < defs.egenes; i++ )
            if( vmax[i + defs.ngenes] == -1. )
                vmax[i] = 255.; /* set max value to 255 */

        /* this part reads the bicoid data */
        //for bcde take maximum of all bcd concentrations, and all maternal penalties

        if( !( b_inlist = ( Blist * ) ReadBicoid( fp, s_current->bcd_section ) ) )
            error( "InitPenalty: no Blist to initialize penalty" );
        else {
            for( b_current = b_inlist; b_current; b_current = b_current->next )
                if( b_current->conc > *mmax )
                    *mmax = b_current->conc;

            free_Blist( b_inlist );
        }

        /* this part read maternal penalty data */

        sprintf( bcdpen_section, "maternal_penalty_data.%d", genot_id );
        if( FindSection( fp, bcdpen_section ) ) {
            if( ( b_inlist = ( Blist * ) ReadBicoid( fp, bcdpen_section ) ) ) {
                for( b_current = b_inlist; b_current; b_current = b_current->next )
                    if( b_current->conc > *mmax )
                        *mmax = b_current->conc;
            } else {
                error( "InitPenalty: error reading maternal penalty section for genotype %d", genot_id );
            }

            free_Blist( b_inlist );
        }

        /* just in case there is no bicoid gradient in the data file */

        if( *mmax == -1. )
            *mmax = 255.;       /* set max value to 255 */

    }

    //printf("vmax = %lg %lg %lg %lg %lg %lg %lg %lg\n", vmax[0], vmax[1], vmax[2], vmax[3], vmax[4], vmax[5], vmax[6], vmax[7]);

    free( factpen_section );
    free( bcdpen_section );
    free( extpen_section );
    free( g_type );

    free_Slist( s_inlist );
}

/** InitTTs: Initializes the time points for which we need model output 
 *            i.e. the time points for which we have data.                 
 *      NOTE: we do not allow data at t=0 but it still gets included into  
 *            TTs                                                          
 */
GenoType *
InitTTs( GenoType * facttype, int nalleles ) {
    int i, j;                   /* local loop counters */

    GenoType *tt;

    /* tt.ptr.times is a DArrPtr static to score.c.                            */
    /* Each DArrPtr in the array points to an array of times for which we have */
    /* data for each genotype (if you know what I mean...)                     */

    if( !( tt = ( GenoType * ) calloc( nalleles, sizeof( GenoType ) ) ) )
        error( "InitTTs: could not allocate tt struct" );

    /* the following loop copies the times from the Facts section of GenoTab   */
    /* to the tt array of DArrPtrs                                             */

    for( i = 0; i < nalleles; i++ ) {

        if( !( tt[i].genotype = ( char * ) calloc( MAX_RECORD, sizeof( char ) ) ) )
            error( "InitTTs: could not allocate tt genotype string" );
        tt[i].genotype = strcpy( tt[i].genotype, facttype[i].genotype );
        tt[i].ptr.times.size = 1 + facttype[i].ptr.facts->size;
        if( !( tt[i].ptr.times.array = ( double * ) calloc( tt[i].ptr.times.size, sizeof( double ) ) ) )
            error( "InitTTs: could not allocate bt array" );
        ;
        tt[i].ptr.times.array[0] = 0.;

        for( j = 1; j < tt[i].ptr.times.size; j++ )
            tt[i].ptr.times.array[j] = facttype[i].ptr.facts->record[j - 1].time;
    }
    tt_init_flag = 1;
    return tt;
}

/** InitStepsize: the only thing this function does is putting stepsize 
 *                 and accuracy in a structure to be passed to the solver 
 *                 by the Score() function                             
 */
Step_Acc
InitStepsize( double step, double acc, FILE * slog, char *infile ) {
    Step_Acc step_acc;
    step_acc.stepsize = step;
    step_acc.accuracy = acc;
    step_acc.slogptr = slog;
    step_acc.filename = infile;
    return step_acc;
}



/*** REAL SCORING CODE HERE ************************************************/

/** Score: as the name says, score runs the simulation, gets a solution 
 *          and then compares it to the data using the Eval least squares  
 *          function                                                       
 *   NOTE:  both InitZygote and InitScoring have to be called first!       
 */

void
Score( Input * inp, ScoreOutput * out, int jacobian ) {
    //name of the output dir
    //extern char *outname;
    ScoreEval eval;
    int i, j, ii;
    double totalscore = 0;
    // file pointer
    //FILE *fp;
    // name of debug full filename (with path)
    char *debugfile = NULL;
    // stores the Solution from Blastoderm
    NArrPtr answer;

    // summed squared differences
    double chisq = 0;

    // variable for penalty
    double penalty = 0;

    /* debugging mode: need debugging file name */
    if( debug ) {
        debugfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    }

    /*if( !tt_init_flag ) {          // tt_init_flag is a flag static to score.c        //REMOVE? - InitTTs is done in getFacts() function
       InitTTs();                       // initializes tabulated times (tt)
       } */

    /* The following will be called after parms are tweaked, hence it must    *
     * check signs. If it appears cleaner, sign checking could be done by the *
     * tweaker                                                                */
    //printf("Score: checking limits... \n");
    for( ii = 0; ii < inp->zyg.defs.ngenes; ii++ ) {
        if( inp->zyg.parm.R[ii] < 0 )
            inp->zyg.parm.R[ii] = -inp->zyg.parm.R[ii];
	
        if( inp->zyg.parm.lambda[ii] < 0 )
            inp->zyg.parm.lambda[ii] = -inp->zyg.parm.lambda[ii];
    }
    /* The following fors and ifs check the searchspace in such a way as to   *
     * return after as few calculations as possible                           */
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {
        if( inp->zyg.parm.R[i] > inp->sco.searchspace->Rlim[i]->upper ) {
            //printf("OUT_OF_BOUND_R: %.10lf > %.10lf | %d\n", inp->zyg.parm.R[i], inp->sco.searchspace->Rlim[i]->upper, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.R[i] < inp->sco.searchspace->Rlim[i]->lower ) {
            //printf("OUT_OF_BOUND_R: %.10lf < %.10lf | %d\n", inp->zyg.parm.R[i], inp->sco.searchspace->Rlim[i]->lower, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.lambda[i] > inp->sco.searchspace->lambdalim[i]->upper ) {
            //printf("OUT_OF_BOUND_lambda: %.10lf > %.10lf | %d\n", inp->zyg.parm.lambda[i], inp->sco.searchspace->lambdalim[i]->upper, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.lambda[i] < inp->sco.searchspace->lambdalim[i]->lower ) {
            //printf("OUT_OF_BOUND_lambda: %.10lf < %.10lf | %d\n", inp->zyg.parm.lambda[i], inp->sco.searchspace->lambdalim[i]->lower, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.tau[i] > inp->sco.searchspace->taulim[i]->upper ) {
            //printf("OUT_OF_BOUND_tau>\n");
            out->score = FORBIDDEN_MOVE;
            return;
       }
        if( inp->zyg.parm.tau[i] < inp->sco.searchspace->taulim[i]->lower ) {
            //printf("OUT_OF_BOUND_tau<\n");
            out->score = FORBIDDEN_MOVE;
            return;
        }
    }
   // printf("SCR STEP2\n");
    if( ( inp->zyg.defs.diff_schedule == 'A' ) || ( inp->zyg.defs.diff_schedule == 'C' ) ) {
        if( inp->zyg.parm.d[0] > inp->sco.searchspace->dlim[0]->upper ) {
            //printf("OUT_OF_BOUND_d: %.10lf > %.10lf\n", inp->zyg.parm.d[0], inp->sco.searchspace->dlim[0]->upper);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.d[0] < inp->sco.searchspace->dlim[0]->lower ) {
            //printf("OUT_OF_BOUND_d: %.10lf < %.10lf\n", inp->zyg.parm.d[0], inp->sco.searchspace->dlim[0]->lower);
            out->score = FORBIDDEN_MOVE;
            return;
        }
    } else {
        for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {
            if( inp->zyg.parm.d[i] > inp->sco.searchspace->dlim[i]->upper ) {
                //printf("OUT_OF_BOUND_d: %.10lf > %.10lf for i = %d\n", inp->zyg.parm.d[i], inp->sco.searchspace->dlim[i]->upper, i);
                out->score = FORBIDDEN_MOVE;
                return;
            }
            if( inp->zyg.parm.d[i] < inp->sco.searchspace->dlim[i]->lower ) {
                //printf("OUT_OF_BOUND_d: %.10lf < %.10lf for i = %d\n", inp->zyg.parm.d[i], inp->sco.searchspace->dlim[i]->lower, i);
                out->score = FORBIDDEN_MOVE;
                return;
            }
        }
    }
    
    /* Penalty stuff below:
     * With asym limits, it is easier now to always check for limits. In case
     * we are doing partly limits, partly penalty, we still need to check for
     * the limits... so the test if pen_vec == NULL is not sufficient.
     *
     * Perhaps in the future (if this part of the code appears to be a 
     * bottleneck) we can make a new test for only limits or only penalty, but
     * as said above, now it is simply easier to just check for penalties
     * always.
     */
    /* If you're using limits on contributors to u, check'em here */
    
    /*if( inp->sco.searchspace->pen_vec == NULL ) {*/
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {
        for( j = 0; j < inp->zyg.defs.ngenes; j++ ) {

            if( inp->zyg.parm.T[( i * inp->zyg.defs.ngenes ) + j] > inp->sco.searchspace->Tlim[( i * inp->zyg.defs.ngenes ) + j]->upper ) {
                //printf("OUT_OF_BOUND_T: %.10lf > %.10lf | %d | %d\n", inp->zyg.parm.T[(i * inp->zyg.defs.ngenes) + j], inp->sco.searchspace->Tlim[(i * inp->zyg.defs.ngenes) + j]->upper, i, j);
                out->score = FORBIDDEN_MOVE;
                return;
            }
            if( inp->zyg.parm.T[( i * inp->zyg.defs.ngenes ) + j] < inp->sco.searchspace->Tlim[( i * inp->zyg.defs.ngenes ) + j]->lower ) {
                //printf("OUT_OF_BOUND_T: %.10lf < %.10lf | %d | %d\n", inp->zyg.parm.T[(i * inp->zyg.defs.ngenes) + j], inp->sco.searchspace->Tlim[(i * inp->zyg.defs.ngenes) + j]->lower, i, j);
                out->score = FORBIDDEN_MOVE;
                return;
            }
        }
        for( j = 0; j < inp->zyg.defs.egenes; j++ ) {
            if( inp->zyg.parm.E[( i * inp->zyg.defs.egenes ) + j] > inp->sco.searchspace->Elim[( i * inp->zyg.defs.egenes ) + j]->upper ) {
                //printf("OUT_OF_BOUND_E: %.10lf > %.10lf | %d | %d\n", inp->zyg.parm.E[(i * inp->zyg.defs.egenes) + j], inp->sco.searchspace->Elim[(i * inp->zyg.defs.egenes) + j]->upper, i, j);
                out->score = FORBIDDEN_MOVE;
                return;
            }
            if( inp->zyg.parm.E[( i * inp->zyg.defs.egenes ) + j] < inp->sco.searchspace->Elim[( i * inp->zyg.defs.egenes ) + j]->lower ) {
                //printf("OUT_OF_BOUND_E: %.10lf < %.10lf | %d | %d\n", inp->zyg.parm.E[(i * inp->zyg.defs.egenes) + j], inp->sco.searchspace->Elim[(i * inp->zyg.defs.egenes) + j]->lower, i, j);
                out->score = FORBIDDEN_MOVE;
                return;
            }
        }
        if( inp->zyg.parm.m[i] > inp->sco.searchspace->mlim[i]->upper ) {
            //printf("OUT_OF_BOUND_m: %.10lf > %.10lf for i = %d\n", inp->zyg.parm.m[i], inp->sco.searchspace->mlim[i]->upper, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.m[i] < inp->sco.searchspace->mlim[i]->lower ) {
            //printf("OUT_OF_BOUND_m: %.10lf < %.10lf for i = %d\n", inp->zyg.parm.m[i], inp->sco.searchspace->mlim[i]->lower, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.h[i] > inp->sco.searchspace->hlim[i]->upper ) {
            //printf("OUT_OF_BOUND_h: %.10lf > %.10lf for i = %d\n", inp->zyg.parm.h[i], inp->sco.searchspace->hlim[i]->upper, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
        if( inp->zyg.parm.h[i] < inp->sco.searchspace->hlim[i]->lower ) {
            //printf("OUT_OF_BOUND_h: %.10lf < %.10lf for i = %d\n", inp->zyg.parm.h[i], inp->sco.searchspace->hlim[i]->lower, i);
            out->score = FORBIDDEN_MOVE;
            return;
        }
    }
    out->penalty = 0;

    /* if you're going to calculate penalty on u, do it here */
    /* following lines calculate exp of sum of squares penalty function */

    /*} else {*/
    if( inp->sco.searchspace->pen_vec != NULL ) {    
        penalty = GetPenalty( inp, inp->sco.searchspace );
        if( penalty == FORBIDDEN_MOVE ) {
    //        printf("FORBIDDEN_MOVE_Penalty\n");
            out->penalty = FORBIDDEN_MOVE;
    //        return;
        }
    /*if (penalty > 0) {
       printf( "PENALTY = %lg\n", penalty);
    }*/
        out->penalty = penalty;
    }
    
    
    /* runs the model and sums squared differences for all genotypes */
    for( i = 0; i < inp->zyg.nalleles; i++ ) {
        answer = Blastoderm( i, inp->sco.facts.facttype[i].genotype, inp, inp->ste.slogptr );
        if( debug ) {
            sprintf( debugfile, "%s.%s.pout", inp->ste.filename, inp->sco.facts.facttype[i].genotype );
            fp = fopen( debugfile, "w" );
            if( !fp ) {
                perror( "printscore" );
                exit( 1 );
            }
            PrintBlastoderm( fp, answer, "debug_output", MAX_PRECISION, &( inp->zyg ) );
            fclose( fp );
        }
        if( gutparms.flag )     //change this to the new Eval() format, if we want to use it
            GutEval( &eval, &answer, i, inp );
        else {
            Eval( &eval, &answer, i, inp );
        }
        chisq += eval.chisq;
        if( i == 0 ) {
            out->residuals = ( double * ) realloc( out->residuals, eval.residuals_size * sizeof( double ) );
            for( j = 0; j < eval.residuals_size; j++ ) {
                out->residuals[j] = 0;
            }
        }
        for( j = 0; j < eval.residuals_size; j++ ) {
            out->residuals[j] += eval.residuals[j];
        }
        free( eval.residuals );
        for( j = 0; j < answer.size; j++ ) {
            free( answer.array[j].state.array );
        }
        free( answer.array );
    }
    if ( debug ) {
        free( debugfile );
    }
    nbScore++;

    out->score = chisq;
//    printf("score=%lg\n", out->score);
//    totalscore = out->score + out->penalty;
//    printf("penalty=%lg\n", out->penalty);
//    printf("totalscore=%lg debug=%d\n", totalscore, debug);
//    printf("================================================\n");     

    if( ( debug ) && ( totalscore < best_score ) ) {
        //if (totalscore < best_score) {
        //    printf("%d ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BEST RMS SCORE %lg -> rms %lg\n", proc_id, totalscore, sqrt(totalscore/inp->zyg.ndp));
        printf( "%d ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BEST SCORE %lg \n", proc_id, totalscore );
        best_score = totalscore;
    }
    out->size_resid_arr = eval.residuals_size;
}

/** Eval: scores the summed squared differences between equation solution 
 *         and data. Because the times for states written to the Solution  
 *         structure are read out of the data file itself, we do not check 
 *         for consistency of times in this function---all times with data 
 *         will be in the table, but the table may also contain additional 
 *         times.                                                          
 */
void
Eval( ScoreEval * eval, NArrPtr * Solution, int gindex, Input * inp ) {
    const double big_epsilon = BIG_EPSILON;     /* used to recognize new time */

    DataTable fact_tab;         /* stores a copy of the Facts (from GenoTab) */
    DataTable weight_tab = ( const struct DataTable ){ 0 };       
                                /* stores a copy of the Weights (from GenoTab) */
    DataPoint point;            /* used to extract an element of DataTable */
    DataPoint weight;           /* used to extract an element of DataTable */

    GenoType *facttype = inp->sco.facts.facttype;
    GenoType *weighttype = NULL;

    int tindex;                 /* index for facts timepoints */
    int sindex;                 /* index for Solution timepoints */
    int vindex;                 /* index for facts datapoint */

    double time;                /* time for each facts timepoint */
    double difference;          /* diff btw data and model (per datapoint) */
    double *v;                  /* ptr to solution for each timepoint */
    double chisq = 0;           /* the final score to be returned */

    //int count=0;
    //double testsum = 0.0;

    int currsize = 0;
    eval->residuals = NULL;     /* the function returns two distances - chi-square and */
    /* the classic difference between the sum of solution  */
    /* and the sum of facts */
    /* extablish a pointer to the appropriate facts section in GenoTab */

    fact_tab = *( facttype[gindex].ptr.facts );

    if( inp->sco.method == 0 ) {
        weighttype = inp->sco.weights.weighttype;
        if( weighttype[gindex].ptr.facts != NULL ) {
            weight_tab = *( weighttype[gindex].ptr.facts );
        }
    }

    sindex = 0;

    /* the following loop adds all squared differences for the whole Solution */
    for( tindex = 0; tindex < fact_tab.size; tindex++ ) {
        time = fact_tab.record[tindex].time;
        /* new time step to be compared? -> get the Solution for this time */

        while( fabs( time - Solution->array[sindex].time ) >= big_epsilon ) {
            sindex++;
        }

        v = Solution->array[sindex].state.array;
        /* this loop steps through the Solution for a specific timepoint and       */
        /* evaluates the squared diffs                                             */
        eval->residuals = ( double * ) realloc( eval->residuals, ( currsize + fact_tab.record[tindex].size ) * sizeof( double ) );
        for( vindex = 0; vindex < fact_tab.record[tindex].size; vindex++ ) {
            point = fact_tab.record[tindex].array[vindex];
            difference = point.conc - v[point.index];
            if( inp->sco.method == 0 ) {
                if( weighttype[gindex].ptr.facts != NULL ) {    //no weights in input file
                    weight = weight_tab.record[tindex].array[vindex];
                    difference *= weight.conc;
                    //testsum += weight.conc;
                    //count++;
                } else {
                    printf( "WARNING: Error reading weights from input file - using OLS\n" );
                    inp->sco.method = 1;
                }
            }
            chisq += difference * difference;
            /* residuals are weighted; if you want residuals before adding weights, 
             *  put these two lines three lines up :) 
             */
            currsize++;
            eval->residuals[currsize - 1] = fabs( difference );
        }
    }
    //testsum = testsum / count;
    //printf("TESTSUM GI %d EVAL = %lg\n", gindex, testsum);
    eval->chisq = chisq;
    eval->residuals_size = currsize;
}

/*** SCOREGUT FUNCTIONS ****************************************************/

/** SetGuts: sets the gut info in score.c for printing out guts */
void
SetGuts( int gutflag, int ndigits ) {
    gutparms.flag = gutflag;
    gutparms.ndigits = ndigits;
}

/** GutEval: this is the same as Eval, i.e it calculates the summed squa- 
 *            red differences between equation solution and data, with the 
 *            addition that individual squared differences between data-   
 *            points are written to STDOUT in the unfold output format     
 */
void
GutEval( ScoreEval * eval, NArrPtr * Solution, int gindex, Input * inp ) {

    const double big_epsilon = BIG_EPSILON;     // used to recognize new time 

    int i, j;                   // loop counters 

    DataTable fact_tab;         // stores a copy of the Facts (from GenoTab) 
    DataTable weight_tab = ( const struct DataTable ){ 0 };
                                // stores a copy of the Weights (from GenoTab) 
    DataPoint point;            // used to extract an element of DataTable 
    DataPoint weight;           // used to extract an element of DataTable 

    GenoType *facttype = inp->sco.facts.facttype;
    GenoType *weighttype = NULL;

    NArrPtr gut;                // individual square root diff for a datapoint 
    NArrPtr outgut;             // output gut structure 

    int tindex;                 // index for facts timepoints 
    int sindex;                 // index for Solution timepoints 
    int vindex;                 // index for facts datapoint 

    double time;                // time for each facts timepoint 
    double difference;          // diff btw data and model (per datapoint) 
    double *v;                  // ptr to solution for each timepoint 
    double chisq = 0;           // the final score to be returned 
    char gen_print[MAX_RECORD]; // for PrintBlastoderm 

    int currsize = 0;
    eval->residuals = NULL;     // the function returns two distances - chi-square and 
    // the classic difference between the sum of solution 
    // and the sum of facts 
    // extablish a pointer to the appropriate facts section in GenoTab 

    fact_tab = *( facttype[gindex].ptr.facts );


    if( inp->sco.method == 0 ) {
        weighttype = inp->sco.weights.weighttype;
        if( weighttype[gindex].ptr.facts != NULL ) {
            weight_tab = *( weighttype[gindex].ptr.facts );
        }
    }

    sindex = 0;
    gen_print[0] = ( char ) 48 + gindex;
    gen_print[1] = '\0';

    // initialize the gut structure 
    gut.array = ( NucState * ) calloc( Solution->size, sizeof( NucState ) );
    gut.size = Solution->size;
    for( i = 0; i < gut.size; i++ ) {
        gut.array[i].time = Solution->array[i].time;
        gut.array[i].state.size = Solution->array[i].state.size;
        gut.array[i].state.array = ( double * ) calloc( gut.array[i].state.size, sizeof( double ) );
        for( j = 0; j < gut.array[i].state.size; j++ )
            gut.array[i].state.array[j] = -1.0;
    }

    // the following loop adds all squared differences for the whole Solution 
    for( tindex = 0; tindex < fact_tab.size; tindex++ ) {
        time = fact_tab.record[tindex].time;
        // new time step to be compared? -> get the Solution for this time 

        while( fabs( time - Solution->array[sindex].time ) >= big_epsilon ) {
            sindex++;
        }
        v = Solution->array[sindex].state.array;
        // this loop steps through the Solution for a specific timepoint and       
        // evaluates the squared diffs                                             
        eval->residuals = ( double * ) realloc( eval->residuals, ( currsize + fact_tab.record[tindex].size ) * sizeof( double ) );
        for( vindex = 0; vindex < fact_tab.record[tindex].size; vindex++ ) {
            point = fact_tab.record[tindex].array[vindex];
            difference = point.conc - v[point.index];
            if( inp->sco.method == 0 ) {
                if( weighttype[gindex].ptr.facts != NULL ) {
                    weight = weight_tab.record[tindex].array[vindex];
                    difference *= weight.conc;
                } else {
                    printf( "WARNING: Error reading weights from input file - using OLS\n" );
                    inp->sco.method = 1;
                }
            }
            chisq += difference * difference;
            gut.array[sindex].state.array[point.index] = difference * difference;
            currsize++;
            eval->residuals[currsize - 1] = fabs( difference );
        }

    }


    // strip gut struct of cell division times and print it to stdout 

    outgut = ConvertAnswer( gut, inp->sco.facts.tt[gindex].ptr.times );
    PrintBlastoderm( stdout, outgut, strcat( gen_print, " genotype\n" ), gutparms.ndigits, &( inp->zyg ) );

    eval->chisq = chisq;
    eval->residuals_size = currsize;

    FreeSolution( &outgut );
    FreeSolution( &gut );
}

/*** FUNCTIONS THAT RETURN SCORE.C-SPECIFIC STUFF **************************/


/*** GetTTimes: this function returns the times for which there's data *****
 *              for a given genotype                                       *
 *     CAUTION: InitTTs has to be called first!                            *
 ***************************************************************************/

/*
DArrPtr GetTTimes(char *genotype)
{
  int         i;

  if ( tt_init_flag ) {
    i = GetIndex(genotype);
    return tt[i].ptr.times;
  }
  else
    error("GetTTimes: called without initialized TTs");

  return tt[0].ptr.times;               // just to make the compiler happy
 }
 */

/** GetLimits:  returns a pointer to the static limits struct in score.c 
 *     CAUTION:  InitScoring must be called first!                         
 */
SearchSpace *
GetLimits( Input * inp ) {
    return inp->sco.searchspace;
}

/**
Experimental feature: how to calculate penalties in the new model
@author Anton Crombach 
@date 2012, August
*/
double
DoCalculatePenalty( double param ) {
    // how do I want to calculate the penalty?!
    return param * param;
}

/** Calculate penalties for a simple parameter, like R, m, h, d, lambda */
double
CalculateSinglePenalty( double *param, Range **lim, int ncols, double max ) {
    double penalty = 0.0;
    int i; 
    for( i = 0; i < ncols; ++i ) {
        /* calculate the penalty if at least one of the two limits is set to
           DBL_MAX */
        if( fabs(lim[ i ]->lower + DBL_MAX) < EPSILON || 
            fabs(lim[ i ]->upper - DBL_MAX) < EPSILON ) {
            penalty += DoCalculatePenalty( param[ i ] * max );
        }
    }
    return penalty;
}

/** Calculate penalties for a matrix parameter, like T_h, T_b, E_h, E_b */
double
CalculateCompoundPenalty( double *param, Range **lim, int nrows, int ncols, double *vmax ) {
    double penalty = 0.0;
    int i,j;
    for( i = 0; i < nrows; i++ ) {
        for( j = 0; j < ncols; j++ ) {
            /* calculate the penalty if at least one of the two limits is set 
               to DBL_MAX */
            if( fabs(lim[ i*ncols + j ]->lower + DBL_MAX) < EPSILON || 
                fabs(lim[ i*ncols + j ]->upper - DBL_MAX) < EPSILON ) {
                penalty += DoCalculatePenalty( param[ i*ncols + j ] * vmax[j]);
                //printf("penalty(%d,%d)=%lg vmax=%lg, acc=%lg PARAM=%lg\n", j, i, DoCalculatePenalty( param[ i*ncols + j ] * vmax[j]), vmax[j], penalty, param[ i*ncols + j ]);
            }
        }
    }
    return penalty;
}

/** GetPenalty: calculates penalty from static limits, vmax and mmax 
*      CAUTION: InitPenalty must be called first!                         
*/
double
GetPenalty( Input * inp, SearchSpace * limits ) {
    double penalty = 0.0;

    double Lambda;          /* penalty Lambda */
    double mmax;            /* mmax in penalty vector (deprecated!) */
    double *vmax;           /* pointer to vmax array in penalty vector */
    
    EqParms *parm;          /* local pointer to parameters */
    
    static int donethis = 0;
    int i;
    
    parm = &( inp->lparm );
    if( limits->pen_vec == NULL ) 
        return -1;
    
    /* Lambda: first entry in penalty vector */
    Lambda = *( ( limits->pen_vec ) );
    /* locate mmax and vmax */
    mmax = *( ( limits->pen_vec ) + 1 );
    vmax = ( limits->pen_vec ) + 2;
    
    if( debug && !donethis ) {
        printf( "Penalty:\n" );
        printf( "Lambda:    %10.8f\n", Lambda );
        printf( "mmax:    %6.2f\n", mmax );
        for( i = 0; i < inp->zyg.defs.ngenes; i++ )
            printf( "vmax[%d]: %6.2f\n", i, vmax[i] );
        for( i = 0; i < inp->zyg.defs.egenes; i++ )
            printf( "vmax[%d]: %6.2f\n", i, vmax[i + inp->zyg.defs.ngenes] );
        printf( "\n" );
    }

    /* AC: calculate penalty, I need limits to check which params I want. 
       Also, do I want to normalize for the number of params used?
    */
    
    //penalty = CalculateSinglePenalty( parm->R, limits->Rlim, inp->zyg.defs.ncols );
    penalty += CalculateCompoundPenalty( parm->T, limits->Tlim, inp->zyg.defs.ngenes, inp->zyg.defs.ngenes, vmax );
    penalty += CalculateCompoundPenalty( parm->E, limits->Elim, inp->zyg.defs.ngenes, inp->zyg.defs.egenes, ( vmax ) + inp->zyg.defs.ngenes );
    penalty += CalculateSinglePenalty( parm->m, limits->mlim, inp->zyg.defs.ngenes, mmax );
    penalty += CalculateSinglePenalty( parm->h, limits->hlim, inp->zyg.defs.ngenes, 1 );
    
    //penalty += CalculateSinglePenalty( parm->d, limits->dlim, inp->zyg.defs.ncols );
    //penalty += CalculateSinglePenalty( parm->lambda, limits->lambdalim, inp->zyg.defs.ncols );
    //penalty += CalculateSinglePenalty( parm->tau, limits->taulim, inp->zyg.defs.ncols );
    
    
    /* with 64-bit machines, the maximum safe value is about 709, but we stuck to  */
    /* the old value 88.7228391 because we don't want circuits with big penalties anyway */
    
    
    
    if( Lambda * penalty > 88.7228391 ) {
    //if( Lambda * penalty > 709 ) {
        //printf( "# Argument too big: %lf > %lf\n", penalty, ( 85.19565 / Lambda ) );
        return FORBIDDEN_MOVE;
    } else {
        penalty = exp( Lambda * penalty ) - 2.718281828459045;
    }

    donethis = 1;

    return ( penalty < 0 ) ? 0 : penalty;
}

//------------------------------------------------------------------------------
// End of experimental feature
//------------------------------------------------------------------------------

/**  GetPenalty_OLD: calculates penalty from static limits, vmax and mmax 
 *   CAUTION:    InitPenalty must be called first!                         
 */
double
GetPenalty_OLD( Input * inp, SearchSpace * limits ) {

    int i, j;                   /* local loop counter */
    double argument = 0;        /* variable for penalty function argument */
    double penalty = 0;         /* holds penalty result */

    EqParms *parm;              /* local copy of eqparms */

    double Lambda;              /* penalty Lambda */
    double mmax;                /* mmax in penalty vector */
    double *vmax;               /* pointer to vmax array in penalty vector */

    static int donethis = 0;    /* KLUDGE: only print this info once */

    if( limits->pen_vec == NULL )
        return -1;

    Lambda = *( ( limits->pen_vec ) );  /* Lambda: first entry in penalty vector */
    mmax = *( ( limits->pen_vec ) + 1 );        /* locate mmax and vmax */
    vmax = ( limits->pen_vec ) + 2;

    /* print debugging info */
    if( debug && !donethis ) {
        printf( "Penalty:\n" );
        printf( "Lambda:    %10.8f\n", Lambda );
        printf( "mmax:    %6.2f\n", mmax );
        for( i = 0; i < inp->zyg.defs.ngenes; i++ )
            printf( "vmax[%d]: %6.2f\n", i, vmax[i] );
        for( i = 0; i < inp->zyg.defs.egenes; i++ )
            printf( "vmax[%d]: %6.2f\n", i, vmax[i + inp->zyg.defs.ngenes] );
        printf( "\n" );
    }

    /* calculate penalty */
    parm = &( inp->lparm );
    //printf("Argument00 = %lf\n", parm->T[0]);
    for( i = 0; i < inp->zyg.defs.ngenes; i++ ) {
        for( j = 0; j < inp->zyg.defs.ngenes; j++ ) {
            argument += ( parm->T[( i * inp->zyg.defs.ngenes ) + j] * vmax[j] ) * ( parm->T[( i * inp->zyg.defs.ngenes ) + j] * vmax[j] );
        }
        for( j = 0; j < inp->zyg.defs.egenes; j++ ) {
            argument +=
                ( parm->E[( i * inp->zyg.defs.egenes ) + j] * vmax[j + inp->zyg.defs.ngenes] ) * ( parm->E[( i * inp->zyg.defs.egenes ) + j] *
                                                                                                   vmax[j + inp->zyg.defs.ngenes] );
        }
        argument += ( parm->m[i] * mmax ) * ( parm->m[i] * mmax );
        argument += parm->h[i] * parm->h[i];
    }
    /* e was 2.718281828 in old code---doubt it matters!! */
    /* 88.7228391 is the maximum safe value to use with exp (see man exp) */
    // 
    // According to cplusplus.com exp return HUGE_VAL if the value is too big..
    if( ( Lambda * argument ) > 88.7228391 ) {
        //printf( "argument too big: %lf > %lf\n", argument, (88.7228391*(1.0/Lambda)));
        return FORBIDDEN_MOVE;
    } else {
        penalty = exp( Lambda * argument ) - 2.718281828459045;
    }
    if( penalty <= 0 )
        penalty = 0;
    donethis = 1;
    return penalty;
}

/*double GetCurPenalty() { return  currPenalty; }*/


/*** A FUNCTION TO CONVERT PENALTY INTO EXPLICIT LIMITS ********************/

/** Penalty2Limits: uses the inverse function of g(u) to calculate upper 
 *                   and lower limits for T, m and h using the penalty     
 *                   lambda parameter; these limits can only be an appro-  
 *                   ximation, since all T, m and h are added up to yield  
 *                   u for g(u); we try to compensate for this summation   
 *                   dividing the limits by sqrt(n); this function then    
 *                   sets the penalty vector to NULL and supplies explicit 
 *                   limits for T, m and h, which can be used for scram-   
 *                   bling parameters and such                             
 *          CAUTION: this func DOES NOT RESET pen_vec, caller must do this 
 ***************************************************************************
 *                                                                         
 * A short comment on penalty limits (JJ, Aug 7, 2001):                    
 *                                                                         
 * In order to produce random values of T, m and h that are                
 * within reasonable limits when no explicit limits are used for           
 * them, we do the following approximation:                                
 *                                                                          
 * First, we determine the limits of u in g(u) within which                
 * there is no penalty assigned. This is true for all g(u) be-             
 * tween Lambda and (1-Lambda). Therefore, we calculate the in-            
 * verse function of g(u) for those to values to get upper and             
 * lower limits for u (based on ideas by Eric Mjolsness).                  
 * Since we need to sum up all T * vmax, m * mmax and h to get             
 * u, we'll compensate by dividing the limits by the sqrt of the           
 * number of genes in the problem. This way we think, we'll get            
 * reasonable limits for single parameters (idea by JR).                   
 * All the above happens in the Penalty2Limits function in sco-            
 * re.c. When comparing parameters to limits, don't forget to              
 * multiply Ts with vmax and ms with mmax. This happens in main            
 * of scramble below. hs are compared as they are.                         
 * See JJs lab notes for further detail on g(u)-inverse and such.          
 *                                                                         
 */
void
Penalty2Limits( SearchSpace * limits, TheProblem defs ) {
    int i, j;                   /* local loop counters */

    double Lambda;              /* penalty Lambda */

    Range u;                    /* explicit range for u (as in g(u)) */
    Range gu;                   /* range of g(u), in which there's no penalty */

    Range x;                    /* these two are used to store intermediate results below */
    Range y;                    /* x = 2gu - 1 ; y = sqrt( 1 - x * x ) */

    Lambda = limits->pen_vec[0];        /* use explicit variable name for Lambda */

    gu.lower = Lambda;          /* range within which there's no penalty */
    gu.upper = 1 - Lambda;

    x.lower = ( 2 * gu.lower - 1 );     /* the following calculates the inverse */
    x.upper = ( 2 * gu.upper - 1 );     /* function of g(u) for gu limits above */
    /* (see JJs lab notes for details) */
    y.lower = sqrt( 1 - x.lower * x.lower );
    y.upper = sqrt( 1 - x.upper * x.upper );

    u.lower = x.lower / y.lower;
    u.upper = x.upper / y.upper;

    u.lower = u.lower / sqrt( defs.ngenes );    /* this is to compensate for the */
    u.upper = u.upper / sqrt( defs.ngenes );    /* summing up of parameters */

    /*
    limits->Tlim = ( Range ** ) calloc( defs.ngenes * defs.ngenes, sizeof( Range * ) );
    limits->Elim = ( Range ** ) calloc( defs.ngenes * defs.egenes, sizeof( Range * ) );
    limits->mlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );
    limits->hlim = ( Range ** ) calloc( defs.ngenes, sizeof( Range * ) );
    */

    /* Only change limits where the value was set to (-)DBLMAX, and 
     * _importantly_ we do not need to claim memory anymore, since this has 
     * been done during the reading in of limits.
     */
    for( i = 0; i < defs.ngenes; i++ ) {
        for( j = 0; j < defs.ngenes; j++ ) {
            //limits->Tlim[( i * defs.ngenes ) + j] = ( Range * ) malloc( sizeof( Range ) );
            if( fabs( limits->Tlim[( i * defs.ngenes ) + j]->lower + DBL_MAX ) < EPSILON ) {
                limits->Tlim[( i * defs.ngenes ) + j]->lower = u.lower;
            }
            if( fabs( limits->Tlim[( i * defs.ngenes ) + j]->upper - DBL_MAX ) < EPSILON ) {
                limits->Tlim[( i * defs.ngenes ) + j]->upper = u.upper;
            }
        }

        for( j = 0; j < defs.egenes; j++ ) {
            //limits->Elim[( i * defs.egenes ) + j] = ( Range * ) malloc( sizeof( Range ) );
            if( fabs( limits->Elim[( i * defs.egenes ) + j]->lower + DBL_MAX ) < EPSILON ) {
                limits->Elim[( i * defs.egenes ) + j]->lower = u.lower;
            }
            if( fabs( limits->Elim[( i * defs.egenes ) + j]->upper - DBL_MAX ) < EPSILON ) {
                limits->Elim[( i * defs.egenes ) + j]->upper = u.upper;
            }
        }

        //limits->mlim[ i ] = ( Range * ) malloc( sizeof( Range ) );
        if( fabs( limits->mlim[ i ]->lower + DBL_MAX ) < EPSILON ) {
            limits->mlim[ i ]->lower = u.lower;
        }
        if( fabs( limits->mlim[ i ]->upper - DBL_MAX ) < EPSILON ) {
            limits->mlim[ i ]->upper = u.upper;
        }

        //limits->hlim[ i ] = ( Range * ) malloc( sizeof( Range ) );
        if( fabs( limits->hlim[ i ]->lower + DBL_MAX ) < EPSILON ) {
            limits->hlim[ i ]->lower = u.lower;
        }
        if( fabs( limits->hlim[ i ]->upper - DBL_MAX ) < EPSILON ) {
            limits->hlim[ i ]->upper = u.upper;
        }
    }
}

/** List2Facts: takes a Dlist and returns the corresponding DataTable 
 *               structure we use for facts data.                          
 *                                                                         
 * An extensive comment about indices: (by JJ) 
 *                                                                         
 *               All data points with -1 as a value are NOT read from the  
 *               data file. The difference between such ignored and zero   
 *               values is crucial: -1 data points WILL NOT BE COMPARED TO 
 *               simulation data, whereas 0 means NO PROTEIN AT THAT TIME  
 *               IN THAT NUCLEUS.                                          
 *               Index numbers help maintain the integrity of the data.    
 *               An index number is defined as the array index at which a  
 *               protein concentration would be if the data was complete,  
 *               i.e. available for all nuclei at all times. In this way   
 *               a sparse set of data can be compared to a complete set of 
 *               simulation output.                                        
 *               Thus, indices are defined as starting from 1 for each     
 *               DataRecord (each time step) and increase by one for each  
 *               gene in each nucleus in the order predefined by JR.       
 *                                                                         
 */
DataTable *
List2Facts( Dlist * inlist, int ngenes ) {
    int i = 0;
    int j;                      /* local loop counters */

    double now = -999999999.;   /* assigns data to specific time */

    Dlist *current;             /* holds current element of Dlist */

    DataTable *D;               /* local copy of DataTable */


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
            D->size++;          /* one DataRecord for each time */
            D->record =         /* allocate DataRecord */
                ( DataRecord * ) realloc( D->record, D->size * sizeof( DataRecord ) );

            D->record[D->size - 1].time = now;  /* next three lines define */
            D->record[D->size - 1].size = 0;    /* DataRecord for each */
            D->record[D->size - 1].array = NULL;        /* time step */
            i = 0;
        }

        for( j = 1; j <= ngenes; j++ ) {        /* always: read concs into array */
            if( current->d[j] != IGNORE ) {     /* valid conc? if IGNORE -> ignore! */
                D->record[D->size - 1].size++;  /* one more in this record */
                /* reallocate memory for array! */
                D->record[D->size - 1].array = realloc( D->record[D->size - 1].array, D->record[D->size - 1].size * sizeof( DataPoint ) );

                /* the following two lines assign concentration value and index *********** */

                D->record[D->size - 1].array[D->record[D->size - 1].size - 1].conc = current->d[j];
                D->record[D->size - 1].array[D->record[D->size - 1].size - 1].index = i;
            }

            i++;

        }
    }

    return D;
}

/** FreeFacts: Function to fre tha DataTable object */
void
FreeFacts( DataTable * D ) {
    int i;
    for( i = 0; i < D->size; i++ )
        free( D->record[i].array );
    free( D->record );
    free( D );
}

/** InitHistory: Initializing the full set of nuclei based on the lineages of the history,
       please make sure that all of the alleles' lineages are the same */
InterpObject *
InitHistory( FILE * fp, Input * inp ) {

    InterpObject *polations;    /* array of interpobjects for the interpolating functions for history for each genotype */
    int i, ii;
    DataTable **temp_table;
    double *theta_discons;
    int theta_discons_size;
    double *temp_divtable;
    double *temp_durations;
    Slist *geno, *curr;         /* We will read in the genotypes to set the interp_dat, bias_dat tables */
    geno = ( Slist * ) ReadGenotypes( fp, inp->zyg.defs.ngenes );
    if( inp->zyg.nalleles == 0 )
        inp->zyg.nalleles = count_Slist( geno );
    if( !( temp_table = ( DataTable ** ) calloc( 1, sizeof( DataTable * ) ) ) )
        error( "InitHistory: could not allocate temp_table struct" );
    if( !( polations = ( InterpObject * ) calloc( inp->zyg.nalleles, sizeof( InterpObject ) ) ) )
        error( "InitHistory: could not allocate interp_tab struct" );

    for( curr = geno, i = 0; curr; curr = curr->next, i++ ) {
        GetInterp( fp, curr->hist_section, inp, inp->zyg.defs.ngenes, temp_table );
        DoInterp( *temp_table, polations + i, inp->zyg.defs.ngenes, &( inp->zyg ) );
        FreeFacts( *temp_table );
        theta_discons = Get_Theta_Discons( &theta_discons_size, &( inp->zyg ) );
        ( polations + i )->fact_discons =
            ( double * ) realloc( ( polations + i )->fact_discons, ( ( polations + i )->fact_discons_size + theta_discons_size ) * sizeof( double ) );
        for( ii = 0; ii < theta_discons_size; ii++ )
            ( polations + i )->fact_discons[( polations + i )->fact_discons_size + ii] = theta_discons[ii];
        ( polations + i )->fact_discons_size += theta_discons_size;
        free( theta_discons );
        if( inp->zyg.defs.ndivs > 0 ) {
            if( !( temp_divtable = inp->zyg.times.div_times ) )
                error( "InitHistory: error getting temp_divtable" );
            if( !( temp_durations = inp->zyg.times.div_duration ) )
                error( "Inithistory: error getting division temp_durations" );

            for( ii = 0; ii < inp->zyg.defs.ndivs; ii++ ) {
                ( polations + i )->fact_discons = ( double * ) realloc( ( polations + i )->fact_discons,
                                                                        ( ( polations + i )->fact_discons_size + 4 ) * sizeof( double ) );

                ( polations + i )->fact_discons[( polations + i )->fact_discons_size] = temp_divtable[ii];
                ( polations + i )->fact_discons[( polations + i )->fact_discons_size + 1] = temp_divtable[ii] + EPSILON;
                ( polations + i )->fact_discons[( polations + i )->fact_discons_size + 2] = temp_divtable[ii] - temp_durations[ii];
                ( polations + i )->fact_discons[( polations + i )->fact_discons_size + 3] = temp_divtable[ii] - temp_durations[ii] + EPSILON;

                ( polations + i )->fact_discons_size += 4;
            }
        }
        qsort( ( void * ) ( polations + i )->fact_discons, ( polations + i )->fact_discons_size, sizeof( double ),
               ( int ( * )( const void *, const void * ) ) compare );
    }

    free_Slist( geno );
    free( temp_table );
    return polations;
}

void
GetInterp( FILE * fp, char *title, Input * inp, int num_genes, DataTable ** interp_tables ) {
    int ndatapoints = 0;        //dummy?
    Dlist *list_dat;
    if( !( list_dat = ( Dlist * ) ReadInterpData( fp, title, num_genes, &ndatapoints ) ) )
        error( "SetFacts: no facts section for history for delay solver" );
    *interp_tables = List2Interp( list_dat, inp, num_genes );

    InitFullNNucs( &( inp->zyg ), inp->zyg.full_lin_start );
    free_Dlist( list_dat );
}

void
FreeHistory( int nalleles, InterpObject * polations ) {
    int i;
    if( nalleles == 0 )
        error( "No alleles to free history for!\n" );
    for( i = 0; i < nalleles; i++ ) {
        FreeInterpObject( polations + i );
    }
}

void
FreeExternalInputs( int nalleles, InterpObject * extinp_polations ) {
    int i;

    if( nalleles == 0 )
        error( "No alleles to free history for!\n" );

    for( i = 0; i < nalleles; i++ ) {
        FreeInterpObject( extinp_polations + i );
    }
}

InterpObject *
InitExternalInputs( FILE * fp, Input * inp ) {
    InterpObject *extinp_polations;     /* array of interpobjects for the interpolating functions for external inputs for each genotype */

    int i;
    DataTable **temp_table;
    Slist *geno, *curr;         /* We will read in the genotypes to set the interp_dat, bias_dat tables */
    geno = ( Slist * ) ReadGenotypes( fp, inp->zyg.defs.ngenes );


    if( inp->zyg.nalleles == 0 )
        inp->zyg.nalleles = count_Slist( geno );

    if( !( temp_table = ( DataTable ** ) calloc( 1, sizeof( DataTable * ) ) ) )
        error( "InitExternalInputs: could not allocate temp_table struct" );


    if( !( extinp_polations = ( InterpObject * ) calloc( inp->zyg.nalleles, sizeof( InterpObject ) ) ) )
        error( "InitExternalInputs: could not allocate extinp_polations struct" );

    for( curr = geno, i = 0; curr; curr = curr->next, i++ ) {

        GetInterp( fp, curr->ext_section, inp, inp->zyg.defs.egenes, temp_table );
        DoInterp( *temp_table, extinp_polations + i, inp->zyg.defs.egenes, &( inp->zyg ) );
        FreeFacts( *temp_table );

    }

    /* Initializing the full set of nuclei based on the lineages of the history,
       please make sure that all of the alleles' lineages are the same */


    free_Slist( geno );
    free( temp_table );
    return extinp_polations;
}

/*
  YF added
 */

/***************************************************************************/

void
initResComp( NArrPtr inData ) {

    int i, j;
    resComp2.array = ( NucState * ) calloc( inData.size, sizeof( NucState ) );
    resComp2.size = inData.size;
    for( i = 0; i < resComp2.size; i++ ) {
        resComp2.array[i].time = inData.array[i].time;
        resComp2.array[i].state.size = inData.array[i].state.size;
        resComp2.array[i].state.array = ( double * ) calloc( resComp2.array[i].state.size, sizeof( double ) );
        for( j = 0; j < resComp2.array[i].state.size; j++ )
            resComp2.array[i].state.array[j] = -1.0;
    }
}

/*void initResCompPar(char *outname){
   resC=1;
   nbScore=0;   
   ResCompFileName = (char *)calloc(MAX_RECORD, sizeof(char)); 
   ResCompFileName = strcpy(ResCompFileName, outname);   
   ResCompFileName = strcat(ResCompFileName, ".resComp");    
 }*/

int
isResComp(  ) {
    return ( resC == 1 ) ? 1 : 0;
}
