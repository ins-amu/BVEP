#include <amigoRHS.h>
#include <AMIGO_model.h>
#include <math.h>


	/* *** Definition of the states *** */

#define	IKKn         Ith(y,0)
#define	IKKa         Ith(y,1)
#define	IKKi         Ith(y,2)
#define	IKKaIkBa     Ith(y,3)
#define	IKKaIkBaNFkB Ith(y,4)
#define	NFkB         Ith(y,5)
#define	NFkBn        Ith(y,6)
#define	A20          Ith(y,7)
#define	A20t         Ith(y,8)
#define	IkBa         Ith(y,9)
#define	IkBan        Ith(y,10)
#define	IkBat        Ith(y,11)
#define	IkBaNFkB     Ith(y,12)
#define	IkBanNFkBn   Ith(y,13)
#define	cgent        Ith(y,14)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dIKKn         Ith(ydot,0)
#define	dIKKa         Ith(ydot,1)
#define	dIKKi         Ith(ydot,2)
#define	dIKKaIkBa     Ith(ydot,3)
#define	dIKKaIkBaNFkB Ith(ydot,4)
#define	dNFkB         Ith(ydot,5)
#define	dNFkBn        Ith(ydot,6)
#define	dA20          Ith(ydot,7)
#define	dA20t         Ith(ydot,8)
#define	dIkBa         Ith(ydot,9)
#define	dIkBan        Ith(ydot,10)
#define	dIkBat        Ith(ydot,11)
#define	dIkBaNFkB     Ith(ydot,12)
#define	dIkBanNFkBn   Ith(ydot,13)
#define	dcgent        Ith(ydot,14)

	/* *** Definition of the parameters *** */

#define	a1    (*amigo_model).pars[0]
#define	a2    (*amigo_model).pars[1]
#define	t1    (*amigo_model).pars[2]
#define	a3    (*amigo_model).pars[3]
#define	t2    (*amigo_model).pars[4]
#define	c1a   (*amigo_model).pars[5]
#define	c2a   (*amigo_model).pars[6]
#define	c3a   (*amigo_model).pars[7]
#define	c4a   (*amigo_model).pars[8]
#define	c5a   (*amigo_model).pars[9]
#define	c6a   (*amigo_model).pars[10]
#define	c1    (*amigo_model).pars[11]
#define	c2    (*amigo_model).pars[12]
#define	c3    (*amigo_model).pars[13]
#define	c4    (*amigo_model).pars[14]
#define	c5    (*amigo_model).pars[15]
#define	k1    (*amigo_model).pars[16]
#define	k2    (*amigo_model).pars[17]
#define	k3    (*amigo_model).pars[18]
#define	kprod (*amigo_model).pars[19]
#define	kdeg  (*amigo_model).pars[20]
#define	kv    (*amigo_model).pars[21]
#define	i1    (*amigo_model).pars[22]
#define	e2a   (*amigo_model).pars[23]
#define	i1a   (*amigo_model).pars[24]
#define	e1a   (*amigo_model).pars[25]
#define	c1c   (*amigo_model).pars[26]
#define	c2c   (*amigo_model).pars[27]
#define	c3c   (*amigo_model).pars[28]
#define Tr	((*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim])

	/* *** Definition of the algebraic variables *** */

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_NFKB(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	/* *** Equations *** */

	dIKKn=kprod-kdeg*IKKn-Tr*k1*IKKn;
	dIKKa=Tr*k1*IKKn-k3*IKKa-Tr*k2*IKKa*A20-kdeg*IKKa-a2*IKKa*IkBa+t1*IKKaIkBa-a3*IKKa*IkBaNFkB+t2*IKKaIkBaNFkB;
	dIKKi=k3*IKKa+Tr*k2*IKKa*A20-kdeg*IKKi;
	dIKKaIkBa=a2*IKKa*IkBa-t1*IKKaIkBa;
	dIKKaIkBaNFkB=a3*IKKa*IkBaNFkB-t2*IKKaIkBaNFkB;
	dNFkB=c6a*IkBaNFkB-a1*NFkB*IkBa+t2*IKKaIkBaNFkB-i1*NFkB;
	dNFkBn=i1*kv*NFkB-a1*IkBan*NFkBn;
	dA20=c4*A20t-c5*A20;
	dA20t=c2+c1*NFkBn-c3*A20t;
	dIkBa=-a2*IKKa*IkBa-a1*IkBa*NFkB+c4a*IkBat-c5a*IkBa-i1a*IkBa+e1a*IkBan;
	dIkBan=-a1*IkBan*NFkBn+i1a*kv*IkBa-e1a*kv*IkBan;
	dIkBat=c2a+c1a*NFkBn-c3a*IkBat;
	dIkBaNFkB=a1*IkBa*NFkB-c6a*IkBaNFkB-a3*IKKa*IkBaNFkB+e2a*IkBanNFkBn;
	dIkBanNFkBn=a1*IkBan*NFkBn-e2a*kv*IkBanNFkBn;
	dcgent=c2c+c1c*NFkBn-c3c*cgent;

	return(0);

}


/* Jacobian of the system (dfdx)*/
int amigoJAC_NFKB(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_NFKB(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}

#define	 IKKn         (amigo_model->sim_results[0][j]) 
#define	 IKKa         (amigo_model->sim_results[1][j]) 
#define	 IKKi         (amigo_model->sim_results[2][j]) 
#define	 IKKaIkBa     (amigo_model->sim_results[3][j]) 
#define	 IKKaIkBaNFkB (amigo_model->sim_results[4][j]) 
#define	 NFkB         (amigo_model->sim_results[5][j]) 
#define	 NFkBn        (amigo_model->sim_results[6][j]) 
#define	 A20          (amigo_model->sim_results[7][j]) 
#define	 A20t         (amigo_model->sim_results[8][j]) 
#define	 IkBa         (amigo_model->sim_results[9][j]) 
#define	 IkBan        (amigo_model->sim_results[10][j]) 
#define	 IkBat        (amigo_model->sim_results[11][j]) 
#define	 IkBaNFkB     (amigo_model->sim_results[12][j]) 
#define	 IkBanNFkBn   (amigo_model->sim_results[13][j]) 
#define	 cgent        (amigo_model->sim_results[14][j]) 



void amigoRHS_get_OBS_NFKB(void* data){

	int j;
	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 NFkB_n  amigo_model->obs_results[0][j] 
		#define	 TIkBa_c amigo_model->obs_results[1][j] 
		#define	 A20mRNA amigo_model->obs_results[2][j] 
		#define	 TIKK    amigo_model->obs_results[3][j] 
		#define	 IKK_a   amigo_model->obs_results[4][j] 
		#define	 IkBa_t  amigo_model->obs_results[5][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){
				NFkB_n=NFkBn;
				TIkBa_c=IkBa+IkBaNFkB;
				A20mRNA=A20t;
				TIKK=IKKn+IKKa+IKKi;
				IKK_a=IKKa;
				IkBa_t=IkBat;

			}

		 break;

	}

	return(amigo_model);

}

#define	 IKKn         (amigo_model->sens_results[0][j][k]) 
#define	 IKKa         (amigo_model->sens_results[1][j][k]) 
#define	 IKKi         (amigo_model->sens_results[2][j][k]) 
#define	 IKKaIkBa     (amigo_model->sens_results[3][j][k]) 
#define	 IKKaIkBaNFkB (amigo_model->sens_results[4][j][k]) 
#define	 NFkB         (amigo_model->sens_results[5][j][k]) 
#define	 NFkBn        (amigo_model->sens_results[6][j][k]) 
#define	 A20          (amigo_model->sens_results[7][j][k]) 
#define	 A20t         (amigo_model->sens_results[8][j][k]) 
#define	 IkBa         (amigo_model->sens_results[9][j][k]) 
#define	 IkBan        (amigo_model->sens_results[10][j][k]) 
#define	 IkBat        (amigo_model->sens_results[11][j][k]) 
#define	 IkBaNFkB     (amigo_model->sens_results[12][j][k]) 
#define	 IkBanNFkBn   (amigo_model->sens_results[13][j][k]) 
#define	 cgent        (amigo_model->sens_results[14][j][k]) 



void amigoRHS_get_sens_OBS_NFKB(void* data){
	int j,k;

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 NFkB_n  amigo_model->sens_obs[0][j][k] 
		#define	 TIkBa_c amigo_model->sens_obs[1][j][k] 
		#define	 A20mRNA amigo_model->sens_obs[2][j][k] 
		#define	 TIKK    amigo_model->sens_obs[3][j][k] 
		#define	 IKK_a   amigo_model->sens_obs[4][j][k] 
		#define	 IkBa_t  amigo_model->sens_obs[5][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					NFkB_n=NFkBn;
					TIkBa_c=IkBa+IkBaNFkB;
					A20mRNA=A20t;
					TIKK=IKKn+IKKa+IKKi;
					IKK_a=IKKa;
					IkBa_t=IkBat;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon_NFKB(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
}