/*$Id: hill_function.c 1372 2012-01-24 11:16:50Z cokelaer $ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double hill_function(double x,double n,double k)
{
	return(pow(x,n)/(pow(x,n)+pow(k,n)));
}
