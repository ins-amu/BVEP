
void printsolution_(void *, double *, double *);
int savehdf5solutions_(double *, int *, int *);
int savehdf5solutionspar_(double *, int *, int *, int *, const char *);
void initprintfile_(void *, double *, int *, int *, double *, long *, int *);
void printiterationcesslog_(void *,  double *, long *, double *, int *) ;
void printdescartsolution_(void *,  double *, int  *);
void printinititeration_(void *,  int *, double *);
void printlsinitlog_(void *, double * );
void printlocalsolverinsert_(void *, double *, double *, int *, int * );
void printinitlocalsolverinsert_(void *, double *, double *, int * );
void printlsendlog_(void *, double *, double *, long *);
#ifdef MPI2
void printmasterocurrencesend_(void *, int *, int  *, int *, int *, int *, double * );
#endif
void printcomparenewsolutionmasterlog_(void *, double *, double *, double *, double *, double * );
void printputmasterlog_(void *, double *, double *, int * );
void printputcheckmasterlog_(void *, double *new, double * );
void printreturnselectsolution_(void *,  double *, double *, double *, int * );
void printcheckvtr_(void *,  double *, long *,  double *, double * );
void printrecvmasterlog_(void *,  double *, double *, double *, int * ) ;
void printputmasterendlog_(void *,  double *, double *, int* ) ;
void printrefset_(void *, double *, int *, double *, int * ) ;
void printrestart_(void *,  double * ) ;
void printrestartslave_(void *,  double * );
void printreplaceslavelog_(void *, double *, double *, int * );
void printreceivedslave_(void *, double *, double *) ;
void printdiscardreceivedslave_(void *, double *, double *, double * );
void printdiscardreceivedslave2_(void *, double *, double *, double * );
void printcomenewsolutionslavelog_(void *, double *new, double *, double *, double * ) ;
void printfinalsendslavelog_(void *, double *new, double * ) ;
void printiteration_(void *, int *, double *, long *, double *);
void printgant_(void *, double *, int *);
void printiterationcess_(void *, double *, long *, double *,
        int *, int *);
void printverboselocaloutput_(void *, int *, double *, double *,int *) ;
void printverboselocaloutput2_(void *, double *, int *);
void improvelocalsolver_(void *, double *, double * ) ;
void print_verbose_local_success(experiment_total , int ) ;
void print_end_file_(experiment_total *, double *, double *, result_solver *) ;
void verboseiteration_(void *, int ,  double , double , double , long , int );
void verboseiterationfortran_(void *, int *, double *,  double *, double *, long *,
        int *, int *);
char * concat_char(int *, char *, char * );
void matlab_plot_file(experiment_total ,  char *,  char *, const char *, int , 
        const char *marca, int , int , int , int , int , int  );
void matlab_plot_file_gant(experiment_total ,  int , int , char* , int ) ;
int row_count(char* ) ;
double varianze(double *array, double , int ) ;
char select_color(int );
char select_marca(int ) ;
void plot_file(experiment_total , int , int , int );
void plot_file_cess(experiment_total , int , int , int );
void initoutputvars_(void *);
void updateresultsandprint_(void *, void *, double *, long *, double *, double * ) ;
void updateresultsandprintess_(void *, void *, double *, long *, double *, double *, double *, double *) ;
int updateresultsrandomsearch_(void *, void *, double *, double *, double * );
void init_message(int ,  experiment_total * , int) ;
void printresults_(void *, double *, long *, double *, double *, double *);
void printpercentage_(void *, double *, double *, double *);
void bechmark_message(experiment_total *, const char *);
void printresults_end(experiment_total *, result_solver );
