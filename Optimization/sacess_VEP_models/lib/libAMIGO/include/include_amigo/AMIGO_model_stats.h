typedef struct{

	int error_flag;

	int sens;

	long int nst;
	long int nfe;
	long int nsetups;
	long int nni;
	long int ncfn;
	long int netf;
	long int nfSe;
	long int nfeS;
	long int nsetupsS; 
	long int nniS; 
	long int ncfnS;
	long int netfS;
	long int nli;
	long int ncfl;
	long int npe;
	long int nps;

	int npars;
	double* pars;

}AMIGO_model_stats;

void free_AMIGO_model_stats(AMIGO_model_stats* amigo_model_stats);
