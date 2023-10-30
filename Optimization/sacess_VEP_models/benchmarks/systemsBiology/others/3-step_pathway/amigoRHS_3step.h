

/* *** Definition of the algebraic variables *** */

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_3step_pathway(realtype , N_Vector , N_Vector , void *);

/* Jacobian of the system (dfdx)*/
int amigoJAC_3step_pathway(int , realtype , N_Vector , N_Vector , DlsMat , void *, N_Vector , N_Vector , N_Vector);

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_3step_pathway(int , realtype , N_Vector , N_Vector , int , N_Vector , N_Vector , void *, N_Vector , N_Vector);

void amigoRHS_get_OBS_3step_pathway(void* );

void amigoRHS_get_sens_OBS_3step_pathway(void* );


void amigo_Y_at_tcon_3step_pathway(void* , realtype , N_Vector );
