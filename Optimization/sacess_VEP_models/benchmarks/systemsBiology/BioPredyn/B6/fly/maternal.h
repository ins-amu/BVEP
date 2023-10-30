/**
 * @file maternal.h
 * @author JR, modified by Yoginho
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 *                                                               
 * @brief Contains various structs, constants and bias-related stuff.
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

#ifndef MATERNAL_INCLUDED
#define MATERNAL_INCLUDED

#include <stdio.h>

//#include "global.h"
//#include "distributions.h"

extern const double old_divtimes[3];

/*** CONSTANTS *************************************************************/

/* The following is the maximum stepsize for solvers (now set to the whole */
/* time it takes to model to run (100.8 minutes))                          */

extern const double MAX_STEPSIZE;

/* masks for lineage numbers */

extern const int CYCLE1;
extern const int CYCLE2;
extern const int CYCLE3;
extern const int CYCLE4;
extern const int CYCLE5;
extern const int CYCLE6;
extern const int CYCLE7;
extern const int CYCLE8;
extern const int CYCLE9;
extern const int CYCLE10;
extern const int CYCLE11;
extern const int CYCLE12;
extern const int CYCLE13;
extern const int CYCLE14;

/* define total number of divisions including the ones before cycle 10
 */

#define TOTAL_DIVS 6
//#define TOTAL_DIVS 3

/* Data points that have a concentration of -1 will be ignored! */
/* this needs to be here rather than in score.h, since it's used to */
/* calculate the number of data points in ReadData */

#define     IGNORE          -1.


#define c13NucStart 17
#define c13NucEnd 46

#define c14NucStart  34
#define c14NucEnd 91


/*** STRUCTS ***************************************************************/
 
/** @brief This is the problem at hand */
typedef struct TheProblem {
    int ngenes;                 /* number of genes to consider */
    int egenes;                 /* number of external inputs */
    char *gene_ids;             /* pointer to char with gene IDs */
    char *egene_ids;            /* pointer to char with external gene IDs */
    int ndivs;                  /* number of cell divisions */
    int nnucs;                  /* number of nuclei at gastrulation */
    char diff_schedule;         /* diffusion schedule */
    int full_ccycles;
} TheProblem;

/** @brief Holds the equation parameters */
typedef struct EqParms {
    double *R;                  /* strength of each promoter--always >= 0. */
    double *T;                  /* the genetic interconnect matrix */
    double *E;                  /* the external input regulatory matrix */
    double *m;                  /* regulatory coefficient of bcd on gene */
    double *h;                  /* reg. coeff. for generic TFs on synthesis of gene */
    double *d;                  /* spatial interaction at gastrulation--always >= 0. */
    double *lambda;             /*protein half lives--always >= 0. */
    double *tau;                /* delay times for the proteins */
} EqParms;

/** @brief Tweaking individual parameters or not.
 *
 * Each pointer points to an array of ints which represent each parameter          
 * 1 means tweak, 0 means leave it alone */
typedef struct Tweak {
    int *Rtweak;                // which Rs to be tweaked 
    int *Ttweak;                // which Ts to be tweaked 
    int *Etweak;                // which Es to be tweaked 
    int *mtweak;                // which ms to be tweaked 
    int *htweak;                // which hs to be tweaked 
    int *dtweak;                // which ds to be tweaked 
    int *lambdatweak;           // which lambdas to be tweaked 
    int *tautweak;              // which taus to be tweaked 
} Tweak;

/** @brief Output of the (Gut)Eval function */
typedef struct ScoreEval {      
    double chisq;
    size_t residuals_size;
    double *residuals;
} ScoreEval;

/** @brief Information from History and ExternalInputs */
typedef struct FactDiscons {    
    double fact_discons_size;   
    double *fact_discons;
} FactDiscons;

/** @brief History and ExternalInputs to solvers */
typedef struct SolverInput {    
    double time;                
    int genindex;
    FactDiscons all_fact_discons;
} SolverInput;

/** @brief Valid param range */
typedef struct Range {
    double lower;
    double upper;
} Range;

/** @brief Pointers to parameters to be tweaked */
typedef struct ParamList {
    double *param;              
    Range *param_range;         /* pointers to corresponding range limits */
} ParamList;

typedef struct PArrPtr {
    int size;                   /* size of the ParamList array */
    ParamList *array;           /* points to 1st element of ParamList array */
    double *pen_vec;            /* penalty vector: see score.h, struct SearchSpace */
} PArrPtr;

/** @brief General struct used for sized array of doubles */
typedef struct DArrPtr {
    int size;
    double *array;
} DArrPtr;

/** @brief Bicoid gradients */
typedef struct BcdGrad {
    int ccycle;                 /* the cleavage cycle */
    DArrPtr gradient;           /* the array of concs for a gradient */
} BcdGrad;

/** @brief More bicoid gradients */
typedef struct BArrPtr {        /* points to an array of BcdStates of 'size' */
    int size;
    BcdGrad *array;
} BArrPtr;

/** @brief Bias data and Blastoderm() output (solution) */
typedef struct NucState {
    double time;                /* the time */
    DArrPtr state;              /* the array of v's, and how many of them */
} NucState;

/** @brief Bias data and Blastoderm() output (solution) */
typedef struct NArrPtr {
    int size;                   /* How many NucState elements in array? */
    NucState *array;            /* points to 1st element of array of NucState elems. */
} NArrPtr;

                                                   
/* NOTE: The following three structs needs to be in this generic header for union definition below     *
 ***************************************************************************/

/** @brief Store facts data.
 *
 * The idea is to encode sparse data against non-sparse v[index] state     
 * arrays. This will in most of the present cases use more memory, but     
 * that won't always be true. Also, it is more convenient for those        
 * 'don't care' points. 
 */
typedef struct DataPoint {
    int index;
    double conc;
} DataPoint;

/** @brief Store facts data.
 *
 * The idea is to encode sparse data against non-sparse v[index] state     
 * arrays. This will in most of the present cases use more memory, but     
 * that won't always be true. Also, it is more convenient for those        
 * 'don't care' points. 
 */
typedef struct DataRecord {
    int size;
    double time;
    DataPoint *array;
} DataRecord;

/** @brief Store facts data.
 *
 * The idea is to encode sparse data against non-sparse v[index] state     
 * arrays. This will in most of the present cases use more memory, but     
 * that won't always be true. Also, it is more convenient for those        
 * 'don't care' points. 
 */
typedef struct DataTable {
    int size;
    DataRecord *record;
} DataTable;

/** @brief Used for bicoid, bias, facts and time tables */
typedef union DataPtr {
    DArrPtr times;              /* used for bias and tab times */
    BArrPtr bicoid;             /* for bicoid stuff            */
    NArrPtr bias;               /* for the bias                */
    DataTable *facts;           /* for facts                   */
} DataPtr;

/** @brief Genotype number and pointer to data */
typedef struct GenoType {
    char *genotype;
    DataPtr ptr;
} GenoType;

/** @brief linked lists for reading bicoid (Blist) */
typedef struct Blist {          
    unsigned int lineage;
    double conc;
    struct Blist *next;
} Blist;

/** @brief linked list used to read bias & facts */
typedef struct Dlist {          
    unsigned int lineage;
    double *d;
    struct Dlist *next;
} Dlist;

/** @brief linked list used to read in genotype specific sections from datafile */
typedef struct Slist {
    char *genotype;             /* genotype string (see dataformatX.X) */
    char *bcd_section;          /* bcd section title */
    char *bias_section;         /* bias section title */
    char *fact_section;         /* fact section title */
    char *hist_section;         /* hist section title */
    char *ext_section;          /* ext section title */
    char *weights_section;      /* weights section title */
    struct Slist *next;
} Slist;

/**@brief Just groups the two bias arrays */
typedef struct Bias {
    GenoType *biastype;         /* array of structs that hold bias for each genotype  */
    GenoType *bt;
} Bias;

/** @brief Just groups the two facts arrays */
typedef struct Facts {
    GenoType *facttype;         /* array of structs that hold facts for each genotype */
    GenoType *tt;
} Facts;

/** @brief Just groups the two weights arrays */
typedef struct Weights {
    GenoType *weighttype;       /* array of structs that hold weights for each genotype */
    GenoType *tt;
} Weights;

/** @brief Parameter limits and the penalty vector */
typedef struct SearchSpace {
    double *pen_vec;            /* pointer to array that defines penalty function */
    Range **Rlim;               /* limits fore promoter strengths */
    Range **Tlim;               /* limits for T matrix, NULL if pen_vec != NULL */
    Range **Elim;               /* limits for E matrix, NULL if pen_vec != NULL */
    Range **mlim;               /* limits for m, NULL if pen_vec != NULL */
    Range **hlim;               /* limits for h, NULL if pen_vec != NULL */
    Range **dlim;               /* limit(s) for diffusion parameter(s) */
    Range **lambdalim;          /* limits for lambda (prot. decay rate) */
    Range **taulim;             /* limits for tau (delays) */
} SearchSpace;

/** @brief This is returned by InitScoring function */
typedef struct Scoring {        
    SearchSpace *searchspace;
    Weights weights;
    Facts facts;
    int method;
} Scoring;

/** @brief Interpolation object */
typedef struct InterpObject {
    double *fact_discons;
    int fact_discons_size;
    NArrPtr func;
    NArrPtr slope;
    int maxsize;
    double maxtime;
} InterpObject;

/** @brief Stepsize, accuracy, solver log file pointer and input file name */
typedef struct Step_Acc {
    double stepsize;            /* solver stepsize */
    double accuracy;            /* solver accuracy */
    FILE *slogptr;              /* solver log file pointer */
    char *filename;             /* infile name */
} Step_Acc;

/** @brief Times from the config file or the ones hardcoded in maternal.c  */
typedef struct Times {          
    int total_divs;             
    double gast_time;
    double *div_times;
    double *div_duration;
    double *full_div_times;
    double *full_div_durations;
} Times;

/** @brief Problem information (ngenes, ndivs etc...).
 * 
 * Contains parameters, times, number of nuclei, number of alleles, etc. 
 */
typedef struct Zygote {         
    TheProblem defs;
    EqParms parm;
    GenoType *bcdtype;
    Bias bias;
    Times times;
    int *nnucs;
    int *full_nnucs;
    int *lin_start;
    int *full_lin_start;

    int nalleles;
    int ndp;
} Zygote;

/** @brief The whole input, and nothing but the input.
 *
 * Contains pointers to all other relevant structures with data taken from 
 * the input file.
 */
typedef struct Input {          
    Zygote zyg;
    Scoring sco;
    InterpObject *his;
    InterpObject *ext;
    Step_Acc ste;
    Tweak twe;
    PArrPtr tra;
    //DistParms dis;
    EqParms lparm;
} Input;


/*** GLOBALS ** ************************************************************/

extern int olddivstyle;                /* flag: old or new division times? */
extern double maxconc;                 /* max prot conc: 12 (old-), 255 (newstyle) */
extern double custom_gast;             /* custom gastrulation time set by -S */


/*** FUNCTION PROTOTYPES ***************************************************/

/* Initialization Functions */

/**  InitBicoid: Copies the Blist read by ReadBicoid into the DArrPtr 
 *               structure; the bicoid DArrPtr contains pointers to bcd    
 *               arrays for each cell division cycle                       
 */
GenoType *InitBicoid( FILE * fp, Zygote * zyg );

/**  InitBias:  puts bias records in a form where get_bias can use them; 
 *              it expects times in increasing order; it expects a non-    
 *              sparse entry, with no genes or nuclei missing              
 */
Bias InitBias( FILE * fp, Zygote * zyg );

/**  InitBTs: initializes the BT struct that holds all times for 
 *            which we have bias.                                          
 */
GenoType *InitBTs( GenoType * biastype, int nalleles );

/**  InitNNucs: takes the global defs.nnucs and calculates number of nucs 
 *              for each cleavage cycle which are then stored in reverse   
 *              order in the static nnucs[] array                          
 *   CAUTION:   defs struct and lin_start need to be initialized before!   
 */
int *InitNNucs( Zygote * zyg );

int getStartNuc( double time, char *typedata, char *genderdata );
int getEndNuc( double time, char *typedata, char *genderdata );


/* Functions that return info about embryo */

/** GetBicoid: returns a bicoid gradients (in form of a DArrPtr) for a 
 *              specific time and genotype.                            
 */
DArrPtr GetBicoid( double time, int genindex, GenoType * bcdtype, Zygote * zyg );

/** GetBias: This function returns bias values for a given time and 
 *            genotype.                                                    
 */
DArrPtr GetBias( double time, int genindex, Zygote * zig );

/** GetBTimes: returns a sized array of times for which there is bias */
DArrPtr GetBTimes( char *genotype, Zygote * zyg );

/** GetNNucs: returns number of nuclei for a given time */
int GetNNucs( TheProblem * defs, int *nnucs, double t, Times * times );

/** GetStartLin: returns the lineage number of the most anterior nucleus 
 *                for a given time                                        
 */
int GetStartLin( double t, TheProblem defs, int *lin_start, Times * times );

/** GetCCycle: returns cleavage cycle number for a given time */
unsigned int GetCCycle( double time, int ndivs, Times * times );

/** ParseLineage: takes lineage number as input and returns the cleavage 
 *                 cycle the nucleus belongs to.                           
 */
unsigned int ParseLineage( unsigned int lin );

/** GetDivtable: returns times of cell divisions depending on ndivs and 
 *                olddivstyle; returns NULL in case of an error         
 */
double *GetDivtable( int ndivs );

/** GetDurations: returns pointer to durations of cell divisions de- 
 *                 pending on ndivs and olddivstyle; returns NULL in case  
 *                 of an error                                             
 */
double *GetDurations( int ndivs );

/** GetGastTime: returns time of gastrulation depending on ndivs and 
 *                olddivstyle; returns 0 in case of an error; if a custom  
 *                gastrulation time is chosen with -S, it will be returned 
 *                only if it's bigger than the normal gastrulation time    
 */
double GetGastTime( int ndivs );

/** GetD: returns diffusion parameters D according to the diff. params. 
 *         in the data file and the diffusion schedule used.               
 *   NOTE: Caller must allocate D_tab                                      
 */
void GetD( double t, double *d, double *D_tab, Zygote * zyg );

/** GetRule: returns the appropriate rule for a given time; used by the 
 *            derivative function                                          
 */
int GetRule( double time, Zygote * zyg );

/** Theta: Returns the value of theta(t) in the autonomous 
 *            version of the equations                              
 */
int Theta( double time, Zygote * zyg );

/** GetIndex: this functions returns the genotype index for a given 
 *             genotype number for reading the GenoType struct.            
 */
int GetIndex( char *genotype, GenoType * bcdindex, int nalleles );

/** MakeTable: this function constructs a time table for Blastoderm 
 *             based on a comand line option (print_stepsize) and the      
 *             cell division tables here in maternal.c                     
 * DISCLAIMER: this thing is written in a very bad way; but hey it does    
 *             its job!!!!!!!                                              
 */DArrPtr MakeTable( double p_stepsize, Zygote * zyg );


/** List2Bicoid: takes a Blist and returns the corresponding BArrPtr 
 *                structure; also initializes the lin_start array which    
 *                contains the lineage number of the most anterior nucleus 
 *                for each cell cycle (we need this for cell division and  
 *                printing model output)                                   
 */
BArrPtr List2Bicoid( Blist * inlist, Zygote * zyg );

/** List2Bias: takes a Dlist and returns the corresponding DArrPtr 
 *              structure.                                         
 */
NArrPtr List2Bias( Dlist * inlist, TheProblem defs );

/* Following functions are utility functions for different linked lists 
 * which are used to read in data of unknown size. All utility functions   
 * follow the same scheme (X stands for the different letters below):      
 *                                                                         
 * - init_Xlist:         allocates first element and returns pointer to it 
 * - addto_Xlist:        adds the adduct to an already existing linkd list 
 * - free_Xlist:         frees memory of linked list                       
 * - count_Xlist:        counts elements in the linked list                
 *                                                                         
 */

/* Utility functions for Blist (for reading bicoid) */

/** init_Blist:         allocates first element and returns pointer to it */
Blist *init_Blist( void );

/** addto_Blist:       adds the adduct to an already existing linkd list */
Blist *addto_Blist( Blist * start, Blist * adduct );

/** free_Blist:         frees memory of linked list */
void free_Blist( Blist * start );

/** count_Blist:        counts elements in the linked list */
int count_Blist( Blist * start );


/* Utility functions for Dlist (for reading bias and facts) */

/** init_Dlist:         allocates first element and returns pointer to it */
Dlist *init_Dlist( int size );

/** addto_Dlist:       adds the adduct to an already existing linkd list */
Dlist *addto_Dlist( Dlist * start, Dlist * adduct );

/** free_Dlist:         frees memory of linked list */
void free_Dlist( Dlist * start );

/** count_Dlist:        counts elements in the linked list */
int count_Dlist( Dlist * start );


/* Utility functions for Slist (for reading genotypes ) */

/** init_Slist:         allocates first element and returns pointer to it */
Slist *init_Slist( void );

/** addto_Slist:       adds the adduct to an already existing linkd list */
Slist *addto_Slist( Slist * start, Slist * adduct );

/** free_Slist:         frees memory of linked list */
void free_Slist( Slist * start );

/** count_Slist:        counts elements in the linked list */
int count_Slist( Slist * start );

/**  InitFullNNucs: takes the global defs.nnucs and calculates number of nucs 
 *              for each cleavage cycle which are then stored in reverse   
 *              order in the static nnucs[] array                          
 *   CAUTION:   defs struct and lin_start need to be initialized before!   
 */
void InitFullNNucs( Zygote * zyg, int *full_lin_start );

/** Index2StartLin: get starting lineage from index */
int Index2StartLin( int index, int *full_lin_start );

/** Index2NNuc: get number of nuclei from index */
int Index2NNuc( int index, int *full_nnucs );

/**  GetStartLinIndex: returns the index of the lineage array 
 *                 for a given time                                         
 */
int GetStartLinIndex( double t, TheProblem * defs, Times * times );

double *Get_Theta_Discons( int *theta_discon_size, Zygote * zyg );

/** List2Interp: takes a Dlist and returns the corresponding  
 *              interpolation DataTable structure.                         
 */
DataTable *List2Interp( Dlist * inlist, Input * inp, int ngenes );

/** Comparison function for qsort */
int descend( int *x, int *y );

/**  InitTimes: Reading hardcoded times in case that we don't find the
 *              times section in the input file                            
 */
Times InitTimes( TheProblem defs );

/**  Translate: creates an array of pointers that point to the parameters 
 *              to be tweaked; the PArrPtr that is returned also includes  
 *              pointers to the corresponding parameter ranges for each    
 *              parameter, although this feature is not yet used anywhere  
 *              in the annealing code                                      
 *     CAUTION: InitZygote and InitScoring have to be called first!        
 */
PArrPtr Translate( Input * inp );


#endif
