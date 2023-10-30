/** @file mathLib.h */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MATHLIB_INCLUDED
#define MATHLIB_INCLUDED


#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

#ifndef MAX_RECORD
#define MAX_RECORD 256
#endif
    
    static const double ZeroValue = 0;
    static const int BoolTrue = 1;
    static const int BoolFalse = 0;
    /*extern const int MAX_RECORD; */       /* max. length of lines read from file */

    /*
       //
       // Strange stuff here. AC
       //

       // The following defines the maximum float precision that is supported by
       // the code.

       #define MAX_PRECISION          16

       #define BIG_EPSILON     0.001
       #define DBL_EPSILON     2.2204460492503131e-016 
       #define FLT_EPSILON     1.192092896e-07F

       #define BIG_EQ(x,v) (((v - BIG_EPSILON) < x) && (x <( v + BIG_EPSILON)))
       #define DBL_EPSILON_EQ(x,v) (((v - DBL_EPSILON) < x) && (x <( v + DBL_EPSILON)))
       #define FLT_EPSILON_EQ(x,v) (((v - FLT_EPSILON) < x) && (x <( v + FLT_EPSILON)))
     */

    /** trunc x to n decimals */
    double Trunc( double x, int n );    

    double AbsDiff( double x, double y );
    double SqrDiff( double x, double y );

    /** to malloc memories AllocChar1D(size): size*char */
    char *AllocChar1D( int );
    
    /** to realloc memories ReallocChar1D(s, size) */
    char *ReallocChar1D( char *, int );
    
    /** to free memories FreeChar1D(s) */
    void FreeChar1D( char * );
    
    
     /** to malloc memories AllocChar2D(size1, size2): size1*(char*), size2*char */
    char **AllocChar2D( int, int );
    
    /** to realloc memories ReallocChar2D(s, size1, size2) */
    char **ReallocChar2D( char **, int, int );
    
    /** to free memories FreeChar2D(s, size) */
    void FreeChar2D( char **, int );
    
    /** to malloc memories AllocInt1D(size): size*int */
    int *AllocInt1D( int );
    
    /** to free memories FreeInt1D(s) */
    void FreeInt1D( int * );
    
    /** to malloc memories AllocDbl1D(size): size*double */
    double *AllocDbl1D( int );
    
    /** to free memories FreeDbl1D(s) */
    void FreeDbl1D( double * );
    
    /** to malloc memories AllocDbl2D(size1, size2) */
    double **AllocDbl2D( int, int );
    
    /** to free memories FreeDbl2D(s, size) */
    void FreeDbl2D( double **, int );

    /** to malloc memories AllocDbl3D(size1, size2, siez3) */
    double ***AllocDbl3D( int, int, int );
    
    /** to free memories FreeDbl3D(s, size1, size2) */
    void FreeDbl3D( double ***, int, int );
    
    /** to check if it's equal to zero
     * if(x<min) return true
     * true=BoolTrue=1, false=BoolFalse=0
     */
    int IsZero( double );

#endif

#ifdef __cplusplus
}                               /* closing brace for extern "C" */
#endif
