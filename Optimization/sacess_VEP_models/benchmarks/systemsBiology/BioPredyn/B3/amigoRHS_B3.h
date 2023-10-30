#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h> 

int amigoRHS_B3(realtype , N_Vector , N_Vector , void *);

void amigoRHS_get_OBS_B3(void* );

void amigoRHS_get_sens_OBS_B3(void* );
 
int amigoRHS_B3(realtype , N_Vector , N_Vector , void *);

void amigo_Y_at_tcon_B3(void* , realtype , N_Vector );

int amigoJAC_B3(int , realtype , N_Vector , N_Vector , DlsMat , void *, N_Vector , N_Vector , N_Vector );

int amigoSensRHS_B3(int , realtype , N_Vector , N_Vector , int , N_Vector , N_Vector , void *, N_Vector , N_Vector );

