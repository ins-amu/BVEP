#include <amigoRHS.h>
#include <math.h>
#include <AMIGO_model.h>

	/* *** Definition of the states *** */

#define	G1 Ith(y,0)
#define	G2 Ith(y,1)
#define	G3 Ith(y,2)
#define	E1 Ith(y,3)
#define	E2 Ith(y,4)
#define	E3 Ith(y,5)
#define	M1 Ith(y,6)
#define	M2 Ith(y,7)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dG1 Ith(ydot,0)
#define	dG2 Ith(ydot,1)
#define	dG3 Ith(ydot,2)
#define	dE1 Ith(ydot,3)
#define	dE2 Ith(ydot,4)
#define	dE3 Ith(ydot,5)
#define	dM1 Ith(ydot,6)
#define	dM2 Ith(ydot,7)

	/* *** Definition of the parameters *** */

#define	V1    (*amigo_model).pars[0]
#define	Ki1   (*amigo_model).pars[1]
#define	ni1   (*amigo_model).pars[2]
#define	Ka1   (*amigo_model).pars[3]
#define	na1   (*amigo_model).pars[4]
#define	k_1   (*amigo_model).pars[5]
#define	V2    (*amigo_model).pars[6]
#define	Ki2   (*amigo_model).pars[7]
#define	ni2   (*amigo_model).pars[8]
#define	Ka2   (*amigo_model).pars[9]
#define	na2   (*amigo_model).pars[10]
#define	k_2   (*amigo_model).pars[11]
#define	V3    (*amigo_model).pars[12]
#define	Ki3   (*amigo_model).pars[13]
#define	ni3   (*amigo_model).pars[14]
#define	Ka3   (*amigo_model).pars[15]
#define	na3   (*amigo_model).pars[16]
#define	k_3   (*amigo_model).pars[17]
#define	V4    (*amigo_model).pars[18]
#define	K4    (*amigo_model).pars[19]
#define	k_4   (*amigo_model).pars[20]
#define	V5    (*amigo_model).pars[21]
#define	K5    (*amigo_model).pars[22]
#define	k_5   (*amigo_model).pars[23]
#define	V6    (*amigo_model).pars[24]
#define	K6    (*amigo_model).pars[25]
#define	k_6   (*amigo_model).pars[26]
#define	kcat1 (*amigo_model).pars[27]
#define	Km1   (*amigo_model).pars[28]
#define	Km2   (*amigo_model).pars[29]
#define	kcat2 (*amigo_model).pars[30]
#define	Km3   (*amigo_model).pars[31]
#define	Km4   (*amigo_model).pars[32]
#define	kcat3 (*amigo_model).pars[33]
#define	Km5   (*amigo_model).pars[34]
#define	Km6   (*amigo_model).pars[35]
#define S	((*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim])
#define P	((*amigo_model).controls_v[1][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[1][(*amigo_model).index_t_stim])

	/* *** Definition of the algebraic variables *** */

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_MENDES(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	/* *** Equations *** */

	dG1=V1/(1+pow(P/Ki1,ni1)+pow(Ka1/S,na1))-k_1*G1;
	dG2=V2/(1+pow(P/Ki2,ni2)+pow(Ka2/M1,na2))-k_2*G2;
	dG3=V3/(1+pow(P/Ki3,ni3)+pow(Ka3/M2,na3))-k_3*G3;
	dE1=V4*G1/(K4+G1)-k_4*E1;
	dE2=V5*G2/(K5+G2)-k_5*E2;
	dE3=V6*G3/(K6+G3)-k_6*E3;
	dM1=kcat1*E1*(1/Km1)*(S-M1)/(1+S/Km1+M1/Km2)-kcat2*E2*(1/Km3)*(M1-M2)/(1+M1/Km3+M2/Km4);
	dM2=kcat2*E2*(1/Km3)*(M1-M2)/(1+M1/Km3+M2/Km4)-kcat3*E3*(1/Km5)*(M2-P)/(1+M2/Km5+P/Km6);

	return(0);

}


/* Jacobian of the system (dfdx)*/
int amigoJAC_MENDES(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_MENDES(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}

#define	 G1 (amigo_model->sim_results[0][j]) 
#define	 G2 (amigo_model->sim_results[1][j]) 
#define	 G3 (amigo_model->sim_results[2][j]) 
#define	 E1 (amigo_model->sim_results[3][j]) 
#define	 E2 (amigo_model->sim_results[4][j]) 
#define	 E3 (amigo_model->sim_results[5][j]) 
#define	 M1 (amigo_model->sim_results[6][j]) 
#define	 M2 (amigo_model->sim_results[7][j]) 



void amigoRHS_get_OBS_MENDES(void* data){

	int j;
	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 1:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 2:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 3:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 4:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 5:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 6:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 7:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 8:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 9:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 10:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 11:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 12:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 13:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 14:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;
		#define	 obsG1 amigo_model->obs_results[0][j] 
		#define	 obsG2 amigo_model->obs_results[1][j] 
		#define	 obsG3 amigo_model->obs_results[2][j] 
		#define	 obsE1 amigo_model->obs_results[3][j] 
		#define	 obsE2 amigo_model->obs_results[4][j] 
		#define	 obsE3 amigo_model->obs_results[5][j] 
		#define	 obsM1 amigo_model->obs_results[6][j] 
		#define	 obsM2 amigo_model->obs_results[7][j] 

		 case 15:


			 for (j = 0; j < amigo_model->n_times; ++j){
				obsG1=G1;
				obsG2=G2;
				obsG3=G3;
				obsE1=E1;
				obsE2=E2;
				obsE3=E3;
				obsM1=M1;
				obsM2=M2;

			}

		 break;

	}

	return(amigo_model);

}

#define	 G1 (amigo_model->sens_results[0][j][k]) 
#define	 G2 (amigo_model->sens_results[1][j][k]) 
#define	 G3 (amigo_model->sens_results[2][j][k]) 
#define	 E1 (amigo_model->sens_results[3][j][k]) 
#define	 E2 (amigo_model->sens_results[4][j][k]) 
#define	 E3 (amigo_model->sens_results[5][j][k]) 
#define	 M1 (amigo_model->sens_results[6][j][k]) 
#define	 M2 (amigo_model->sens_results[7][j][k]) 



void amigoRHS_get_sens_OBS_MENDES(void* data){
	int j,k;

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 1:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 2:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 3:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 4:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 5:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 6:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 7:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 8:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 9:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 10:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 11:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 12:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 13:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 14:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;

		 case 15:

		#define	 obsG1 amigo_model->sens_obs[0][j][k] 
		#define	 obsG2 amigo_model->sens_obs[1][j][k] 
		#define	 obsG3 amigo_model->sens_obs[2][j][k] 
		#define	 obsE1 amigo_model->sens_obs[3][j][k] 
		#define	 obsE2 amigo_model->sens_obs[4][j][k] 
		#define	 obsE3 amigo_model->sens_obs[5][j][k] 
		#define	 obsM1 amigo_model->sens_obs[6][j][k] 
		#define	 obsM2 amigo_model->sens_obs[7][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					obsG1=G1;
					obsG2=G2;
					obsG3=G3;
					obsE1=E1;
					obsE2=E2;
					obsE3=E3;
					obsM1=M1;
					obsM2=M2;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon_MENDES(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
}

