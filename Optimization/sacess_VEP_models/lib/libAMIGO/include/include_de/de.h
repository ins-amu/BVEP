
#pragma once
#include <AMIGO_problem.h>
typedef struct
{
	int NP;
	int Gmax;
	
}DE_opts;

DE_opts* get_default_DE_opts();

int DE_AMIGO_pe(AMIGO_problem* amigo_problem);

double DE_AMIGO_pe_objective(double *X,int D, double* Xl, double* Xu, void* data);

