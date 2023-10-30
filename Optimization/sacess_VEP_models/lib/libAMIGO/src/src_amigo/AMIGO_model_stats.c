#pragma once
#include <AMIGO_model_stats.h>
#include <stdlib.h>

void free_AMIGO_model_stats(AMIGO_model_stats* amigo_model_stats){
	
	if(amigo_model_stats->npars>0){
		free(amigo_model_stats->pars);
	}
	free(amigo_model_stats);
}