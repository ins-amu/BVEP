#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h> 
/* *** Definition of the algebraic variables *** */

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_CIRCADIAN(realtype , N_Vector , N_Vector , void *);

/* Jacobian of the system (dfdx)*/
int amigoJAC_CIRCADIAN(int , realtype , N_Vector , N_Vector , DlsMat , void *, N_Vector , N_Vector , N_Vector );

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_CIRCADIAN(int , realtype , N_Vector , N_Vector , int , N_Vector , N_Vector , void *, N_Vector , N_Vector );

void amigoRHS_get_OBS_CIRCADIAN(void* );

void amigoRHS_get_sens_OBS_CIRCADIAN(void* );

void amigo_Y_at_tcon_CIRCADIAN(void* , realtype , N_Vector );
