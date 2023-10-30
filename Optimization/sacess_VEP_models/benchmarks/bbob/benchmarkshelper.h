

double round1( double );
double fmin1(  double ,  double );
double fmax1(  double ,  double );
void unif(double* , int , int );
void gauss(double * , int , int );
void computeXopt(int , int );
void monotoneTFosc(double* );
void freeStarStar(double** , int );
double** reshape(double** , double*, int , int );
void computeRotation(double ** , int, int);
double myrand(void);
double randn(void);
double FGauss(double , double );
double FUniform(double , double , double );
double FCauchy(double , double , double );
int compare_doubles (const void *, const void *);
void initbenchmarkshelper(void *);
void finibenchmarkshelper(void);
double computeFopt(int , int );
void setNoiseSeed(unsigned int , unsigned int );

/* error handling routines - same arguments as printf, i.e. format first, then list of things to print */
/* this one exits after printing - severe error, not recoverable */
void ERROR(const char *, ...);
/* same, but returns to the caller, mild error */
void WARNING(const char *, ...);

/* Checks if sDir exists, 
   creates it if not
   checks if is writable thereafter
   Fatal ERROR if anything fails
*/
void dirOK(char *);

/* create complete pathName from filename and dirname 
   is SYSTEM dependent (should be some #ifdef WINDOWS etc ...)
   fullFileName should already be allocated, at least 1024 bytes long
*/
void createFullFileName(char *, char *, char *);

/* checks the existence of a file */
int existFile(char * );

/* opens a file after checking it is there */
FILE * bbobOpenFile(char * );
