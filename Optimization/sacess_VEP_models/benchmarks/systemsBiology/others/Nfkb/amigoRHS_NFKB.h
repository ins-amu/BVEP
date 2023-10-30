#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h> 

/* *** Definition of the algebraic variables *** */

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_NFKB(realtype , N_Vector , N_Vector , void *);

/* Jacobian of the system (dfdx)*/
int amigoJAC_NFKB(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector );

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_NFKB(int , realtype , N_Vector , N_Vector , int , N_Vector , N_Vector , void *, N_Vector, N_Vector);

void amigoRHS_get_OBS_NFKB(void* );

void amigoRHS_get_sens_OBS_NFKB(void* );

void amigo_Y_at_tcon_NFKB(void* , realtype , N_Vector );
