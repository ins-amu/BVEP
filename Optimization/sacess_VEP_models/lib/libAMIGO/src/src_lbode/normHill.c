/*$Id: normHill.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double normHill(double x,double n,double k)
{
	return(pow(x,n)/(pow(x,n)+pow(k,n))*(1+pow(k,n)));
}
