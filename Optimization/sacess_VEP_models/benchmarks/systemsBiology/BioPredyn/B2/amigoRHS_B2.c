/*
Check README.txt
*/

#include <AMIGO_problem.h>

int amigoRHS_B2(realtype t, N_Vector y, N_Vector ydot, void *data);

void amigoRHS_get_OBS_B2(void* data);

void amigoRHS_get_sens_OBS_B2(void* data);
 
int amigoRHS_B2(realtype t, N_Vector y, N_Vector ydot, void *data);

void amigo_Y_at_tcon_B2(void* amigo_model, realtype t, N_Vector y);


//Define the objective function
//This is an useful wrapper since most solvers in C accept a void pointer
//The problem data is passed as a void pointer
double objective_function_B2(double* x, void* data){
    
    //Convert void pointer to AMIGO_problem pointer
    AMIGO_problem* amigo_problem=(AMIGO_problem*)data;
    
	//Copy x to AMIGO_problem->x
    set_AMIGO_problem_pars(x,amigo_problem);
    
	//Increment the number of evaluations counter
    amigo_problem->nevals++;
	
    return(eval_AMIGO_problem_LSQ(amigo_problem));
	
}


/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_B2(realtype t, N_Vector y, N_Vector ydot, void *data){
	/* *** Definition of the states *** */

#define	cdhap   Ith(y,0)
#define	ce4p    Ith(y,1)
#define	cf6p    Ith(y,2)
#define	cfdp    Ith(y,3)
#define	cg1p    Ith(y,4)
#define	cg6p    Ith(y,5)
#define	cgap    Ith(y,6)
#define	cpep    Ith(y,7)
#define	cpg     Ith(y,8)
#define	cpg2    Ith(y,9)
#define	cpg3    Ith(y,10)
#define	cpgp    Ith(y,11)
#define	cpyr    Ith(y,12)
#define	crib5p  Ith(y,13)
#define	cribu5p Ith(y,14)
#define	csed7p  Ith(y,15)
#define	cxyl5p  Ith(y,16)
#define	cglcex  Ith(y,17)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dcdhap   Ith(ydot,0)
#define	dce4p    Ith(ydot,1)
#define	dcf6p    Ith(ydot,2)
#define	dcfdp    Ith(ydot,3)
#define	dcg1p    Ith(ydot,4)
#define	dcg6p    Ith(ydot,5)
#define	dcgap    Ith(ydot,6)
#define	dcpep    Ith(ydot,7)
#define	dcpg     Ith(ydot,8)
#define	dcpg2    Ith(ydot,9)
#define	dcpg3    Ith(ydot,10)
#define	dcpgp    Ith(ydot,11)
#define	dcpyr    Ith(ydot,12)
#define	dcrib5p  Ith(ydot,13)
#define	dcribu5p Ith(ydot,14)
#define	dcsed7p  Ith(ydot,15)
#define	dcxyl5p  Ith(ydot,16)
#define	dcglcex  Ith(ydot,17)

	/* *** Definition of the parameters *** */

#define	kALDOdhap          (*amigo_model).pars[0]
#define	kALDOeq            (*amigo_model).pars[1]
#define	kALDOfdp           (*amigo_model).pars[2]
#define	kALDOgap           (*amigo_model).pars[3]
#define	kALDOgapinh        (*amigo_model).pars[4]
#define	KDAHPSe4p          (*amigo_model).pars[5]
#define	KDAHPSpep          (*amigo_model).pars[6]
#define	KENOeq             (*amigo_model).pars[7]
#define	KENOpep            (*amigo_model).pars[8]
#define	KENOpg2            (*amigo_model).pars[9]
#define	KG1PATatp          (*amigo_model).pars[10]
#define	KG1PATfdp          (*amigo_model).pars[11]
#define	KG1PATg1p          (*amigo_model).pars[12]
#define	KG3PDHdhap         (*amigo_model).pars[13]
#define	KG6PDHg6p          (*amigo_model).pars[14]
#define	KG6PDHnadp         (*amigo_model).pars[15]
#define	KG6PDHnadphg6pinh  (*amigo_model).pars[16]
#define	KG6PDHnadphnadpinh (*amigo_model).pars[17]
#define	KGAPDHeq           (*amigo_model).pars[18]
#define	KGAPDHgap          (*amigo_model).pars[19]
#define	KGAPDHnad          (*amigo_model).pars[20]
#define	KGAPDHnadh         (*amigo_model).pars[21]
#define	KGAPDHpgp          (*amigo_model).pars[22]
#define	KPDHpyr            (*amigo_model).pars[23]
#define	KpepCxylasefdp     (*amigo_model).pars[24]
#define	KpepCxylasepep     (*amigo_model).pars[25]
#define	KPFKadpa           (*amigo_model).pars[26]
#define	KPFKadpb           (*amigo_model).pars[27]
#define	KPFKadpc           (*amigo_model).pars[28]
#define	KPFKampa           (*amigo_model).pars[29]
#define	KPFKampb           (*amigo_model).pars[30]
#define	KPFKatps           (*amigo_model).pars[31]
#define	KPFKf6ps           (*amigo_model).pars[32]
#define	KPFKpep            (*amigo_model).pars[33]
#define	KPGDHatpinh        (*amigo_model).pars[34]
#define	KPGDHnadp          (*amigo_model).pars[35]
#define	KPGDHnadphinh      (*amigo_model).pars[36]
#define	KPGDHpg            (*amigo_model).pars[37]
#define	KPGIeq             (*amigo_model).pars[38]
#define	KPGIf6p            (*amigo_model).pars[39]
#define	KPGIf6ppginh       (*amigo_model).pars[40]
#define	KPGIg6p            (*amigo_model).pars[41]
#define	KPGIg6ppginh       (*amigo_model).pars[42]
#define	KPGKadp            (*amigo_model).pars[43]
#define	KPGKatp            (*amigo_model).pars[44]
#define	KPGKeq             (*amigo_model).pars[45]
#define	KPGKpg3            (*amigo_model).pars[46]
#define	KPGKpgp            (*amigo_model).pars[47]
#define	KPGluMueq          (*amigo_model).pars[48]
#define	KPGluMupg2         (*amigo_model).pars[49]
#define	KPGluMupg3         (*amigo_model).pars[50]
#define	KPGMeq             (*amigo_model).pars[51]
#define	KPGMg1p            (*amigo_model).pars[52]
#define	KPGMg6p            (*amigo_model).pars[53]
#define	KPKadp             (*amigo_model).pars[54]
#define	KPKamp             (*amigo_model).pars[55]
#define	KPKatp             (*amigo_model).pars[56]
#define	KPKfdp             (*amigo_model).pars[57]
#define	KPKpep             (*amigo_model).pars[58]
#define	KPTSa1             (*amigo_model).pars[59]
#define	KPTSa2             (*amigo_model).pars[60]
#define	KPTSa3             (*amigo_model).pars[61]
#define	KPTSg6p            (*amigo_model).pars[62]
#define	KR5PIeq            (*amigo_model).pars[63]
#define	KRPPKrib5p         (*amigo_model).pars[64]
#define	KRu5Peq            (*amigo_model).pars[65]
#define	KSerSynthpg3       (*amigo_model).pars[66]
#define	KSynth1pep         (*amigo_model).pars[67]
#define	KSynth2pyr         (*amigo_model).pars[68]
#define	KTAeq              (*amigo_model).pars[69]
#define	kTISdhap           (*amigo_model).pars[70]
#define	kTISeq             (*amigo_model).pars[71]
#define	kTISgap            (*amigo_model).pars[72]
#define	KTKaeq             (*amigo_model).pars[73]
#define	KTKbeq             (*amigo_model).pars[74]
#define	LPFK               (*amigo_model).pars[75]
#define	LPK                (*amigo_model).pars[76]
#define	nDAHPSe4p          (*amigo_model).pars[77]
#define	nDAHPSpep          (*amigo_model).pars[78]
#define	nG1PATfdp          (*amigo_model).pars[79]
#define	nPDH               (*amigo_model).pars[80]
#define	npepCxylasefdp     (*amigo_model).pars[81]
#define	nPFK               (*amigo_model).pars[82]
#define	nPK                (*amigo_model).pars[83]
#define	nPTSg6p            (*amigo_model).pars[84]
#define	rmaxALDO           (*amigo_model).pars[85]
#define	rmaxDAHPS          (*amigo_model).pars[86]
#define	rmaxENO            (*amigo_model).pars[87]
#define	rmaxG1PAT          (*amigo_model).pars[88]
#define	rmaxG3PDH          (*amigo_model).pars[89]
#define	rmaxG6PDH          (*amigo_model).pars[90]
#define	rmaxGAPDH          (*amigo_model).pars[91]
#define	rmaxMetSynth       (*amigo_model).pars[92]
#define	rmaxMurSynth       (*amigo_model).pars[93]
#define	rmaxPDH            (*amigo_model).pars[94]
#define	rmaxpepCxylase     (*amigo_model).pars[95]
#define	rmaxPFK            (*amigo_model).pars[96]
#define	rmaxPGDH           (*amigo_model).pars[97]
#define	rmaxPGI            (*amigo_model).pars[98]
#define	rmaxPGK            (*amigo_model).pars[99]
#define	rmaxPGluMu         (*amigo_model).pars[100]
#define	rmaxPGM            (*amigo_model).pars[101]
#define	rmaxPK             (*amigo_model).pars[102]
#define	rmaxPTS            (*amigo_model).pars[103]
#define	rmaxR5PI           (*amigo_model).pars[104]
#define	rmaxRPPK           (*amigo_model).pars[105]
#define	rmaxRu5P           (*amigo_model).pars[106]
#define	rmaxSerSynth       (*amigo_model).pars[107]
#define	rmaxSynth1         (*amigo_model).pars[108]
#define	rmaxSynth2         (*amigo_model).pars[109]
#define	rmaxTA             (*amigo_model).pars[110]
#define	rmaxTIS            (*amigo_model).pars[111]
#define	rmaxTKa            (*amigo_model).pars[112]
#define	rmaxTKb            (*amigo_model).pars[113]
#define	rmaxTrpSynth       (*amigo_model).pars[114]
#define	VALDOblf           (*amigo_model).pars[115]
#define	cfeed              (*amigo_model).pars[116]
#define	Dil                (*amigo_model).pars[117]
#define	mu                 (*amigo_model).pars[118]
#define	cytosol            (*amigo_model).pars[119]
#define	extracellular      (*amigo_model).pars[120]

	/* *** Definition of the algebraic variables *** */

	double	cadp;
	double	camp;
	double	catp;
	double	cnad;
	double	cnadh;
	double	cnadp;
	double	cnadph;
	double	vALDO;
	double	vDAHPS;
	double	vDHAP;
	double	vE4P;
	double	vENO;
	double	vEXTER;
	double	vG1PAT;
	double	vG3PDH;
	double	vG6P;
	double	vG6PDH;
	double	vGAP;
	double	vGAPDH;
	double	vGLP;
	double	vMURSyNTH;
	double	vMethSynth;
	double	vPDH;
	double	vPEP;
	double	vPFK;
	double	vPG;
	double	vPG3;
	double	vPGDH;
	double	vPGI;
	double	vPGK;
	double	vPGM;
	double	vPGP;
	double	vPK;
	double	vPPK;
	double	vPTS;
	double	vR5PI;
	double	vRIB5P;
	double	vRibu5p;
	double	vRu5P;
	double	vSED7P;
	double	vSynth1;
	double	vSynth2;
	double	vTA;
	double	vTIS;
	double	vTKA;
	double	vTKB;
	double	vTRPSYNTH;
	double	vXYL5P;
	double	vf6P;
	double	vfdP;
	double	vpepCxylase;
	double	vpg2;
	double	vpyr;
	double	vrpGluMu;
	double	vsersynth;

	AMIGO_model* amigo_model=(AMIGO_model*)data;

	/* *** Equations *** */

	cadp=0.582+1.73*pow(2.731,-0.15*t)*(0.12*t+0.000214*pow(t,3));
	camp=0.123+7.25*(t/(7.25+1.47*t+0.17*pow(t,2)))+1.073/(1.29+8.05*t);
	catp=4.27-4.163*(t/(0.657+1.43*t+0.0364*pow(t,2)));
	cnad=1.314+1.314*pow(2.73,-0.0435*t-0.342)-(t+7.871)*(pow(2.73,-0.0218*t-0.171)/(8.481+t));
	cnadh=0.0934+0.00111*pow(2.371,-0.123*t)*(0.844*t+0.104*pow(t,3));
	cnadp=0.159-0.00554*(t/(2.8-0.271*t+0.01*pow(t,2)))+0.182/(4.82+0.526*t);
	cnadph=0.062+0.332*pow(2.718,-0.464*t)*(0.0166*pow(t,1.58)+0.000166*pow(t,4.73)+0.1312*pow(10,-9)*pow(t,7.89)+0.1362*pow(10,-12)*pow(t,11)+0.1233*pow(10,-15)*pow(t,14.2));
	vALDO=cytosol*rmaxALDO*(cfdp-cgap*cdhap/kALDOeq)/(kALDOfdp+cfdp+kALDOgap*cdhap/(kALDOeq*VALDOblf)+kALDOdhap*cgap/(kALDOeq*VALDOblf)+cfdp*cgap/kALDOgapinh+cgap*cdhap/(VALDOblf*kALDOeq));
	vDAHPS=cytosol*rmaxDAHPS*pow(ce4p,nDAHPSe4p)*pow(cpep,nDAHPSpep)/((KDAHPSe4p+pow(ce4p,nDAHPSe4p))*(KDAHPSpep+pow(cpep,nDAHPSpep)));
	vDHAP=cytosol*mu*cdhap;
	vE4P=cytosol*mu*ce4p;
	vENO=cytosol*rmaxENO*(cpg2-cpep/KENOeq)/(KENOpg2*(1.0+cpep/KENOpep)+cpg2);
	vEXTER=extracellular*Dil*(cfeed-cglcex);
	vG1PAT=cytosol*rmaxG1PAT*cg1p*catp*(1.0+pow(cfdp/KG1PATfdp,nG1PATfdp))/((KG1PATatp+catp)*(KG1PATg1p+cg1p));
	vG3PDH=cytosol*rmaxG3PDH*cdhap/(KG3PDHdhap+cdhap);
	vG6P=cytosol*mu*cg6p;
	vG6PDH=cytosol*rmaxG6PDH*cg6p*cnadp/((cg6p+KG6PDHg6p)*(1.0+cnadph/KG6PDHnadphg6pinh)*(KG6PDHnadp*(1.0+cnadph/KG6PDHnadphnadpinh)+cnadp));
	vGAP=cytosol*mu*cgap;
	vGAPDH=cytosol*rmaxGAPDH*(cgap*cnad-cpgp*cnadh/KGAPDHeq)/((KGAPDHgap*(1.0+cpgp/KGAPDHpgp)+cgap)*(KGAPDHnad*(1.0+cnadh/KGAPDHnadh)+cnad));
	vGLP=cytosol*mu*cg1p;
	vMURSyNTH=cytosol*rmaxMurSynth;
	vMethSynth=cytosol*rmaxMetSynth;
	vPDH=cytosol*rmaxPDH*pow(cpyr,nPDH)/(KPDHpyr+pow(cpyr,nPDH));
	vPEP=cytosol*mu*cpep;
	vPFK=cytosol*rmaxPFK*catp*cf6p/((catp+KPFKatps*(1.0+cadp/KPFKadpc))*(cf6p+KPFKf6ps*(1.0+cpep/KPFKpep+cadp/KPFKadpb+camp/KPFKampb)/(1.0+cadp/KPFKadpa+camp/KPFKampa))*(1.0+LPFK/pow(1.0+cf6p*(1.0+cadp/KPFKadpa+camp/KPFKampa)/(KPFKf6ps*(1.0+cpep/KPFKpep+cadp/KPFKadpb+camp/KPFKampb)),nPFK)));
	vPG=cytosol*mu*cpg;
	vPG3=cytosol*mu*cpg3;
	vPGDH=cytosol*rmaxPGDH*cpg*cnadp/((cpg+KPGDHpg)*(cnadp+KPGDHnadp*(1.0+cnadph/KPGDHnadphinh)*(1.0+catp/KPGDHatpinh)));
	vPGI=cytosol*rmaxPGI*(cg6p-cf6p/KPGIeq)/(KPGIg6p*(1.0+cf6p/(KPGIf6p*(1.0+cpg/KPGIf6ppginh))+cpg/KPGIg6ppginh)+cg6p);
	vPGK=cytosol*rmaxPGK*(cadp*cpgp-catp*cpg3/KPGKeq)/((KPGKadp*(1.0+catp/KPGKatp)+cadp)*(KPGKpgp*(1.0+cpg3/KPGKpg3)+cpgp));
	vPGM=cytosol*rmaxPGM*(cg6p-cg1p/KPGMeq)/(KPGMg6p*(1.0+cg1p/KPGMg1p)+cg6p);
	vPGP=cytosol*mu*cpgp;
	vPK=cytosol*rmaxPK*cpep*pow(cpep/KPKpep+1.0,nPK-1.0)*cadp/(KPKpep*(LPK*pow((1.0+catp/KPKatp)/(cfdp/KPKfdp+camp/KPKamp+1.0),nPK)+pow(cpep/KPKpep+1.0,nPK))*(cadp+KPKadp));
	vPPK=cytosol*rmaxRPPK*crib5p/(KRPPKrib5p+crib5p);
	vPTS=extracellular*rmaxPTS*cglcex*(cpep/cpyr)/((KPTSa1+KPTSa2*(cpep/cpyr)+KPTSa3*cglcex+cglcex*(cpep/cpyr))*(1.0+pow(cg6p,nPTSg6p)/KPTSg6p));
	vR5PI=cytosol*rmaxR5PI*(cribu5p-crib5p/KR5PIeq);
	vRIB5P=cytosol*mu*crib5p;
	vRibu5p=cytosol*mu*cribu5p;
	vRu5P=cytosol*rmaxRu5P*(cribu5p-cxyl5p/KRu5Peq);
	vSED7P=cytosol*mu*csed7p;
	vSynth1=cytosol*rmaxSynth1*cpep/(KSynth1pep+cpep);
	vSynth2=cytosol*rmaxSynth2*cpyr/(KSynth2pyr+cpyr);
	vTA=cytosol*rmaxTA*(cgap*csed7p-ce4p*cf6p/KTAeq);
	vTIS=cytosol*rmaxTIS*(cdhap-cgap/kTISeq)/(kTISdhap*(1.0+cgap/kTISgap)+cdhap);
	vTKA=cytosol*rmaxTKa*(crib5p*cxyl5p-csed7p*cgap/KTKaeq);
	vTKB=cytosol*rmaxTKb*(cxyl5p*ce4p-cf6p*cgap/KTKbeq);
	vTRPSYNTH=cytosol*rmaxTrpSynth;
	vXYL5P=cytosol*mu*cxyl5p;
	vf6P=cytosol*mu*cf6p;
	vfdP=cytosol*mu*cfdp;
	vpepCxylase=cytosol*rmaxpepCxylase*cpep*(1.0+pow(cfdp/KpepCxylasefdp,npepCxylasefdp))/(KpepCxylasepep+cpep);
	vpg2=cytosol*mu*cpg2;
	vpyr=cytosol*mu*cpyr;
	vrpGluMu=cytosol*rmaxPGluMu*(cpg3-cpg2/KPGluMueq)/(KPGluMupg3*(1.0+cpg2/KPGluMupg2)+cpg3);
	vsersynth=cytosol*rmaxSerSynth*cpg3/(KSerSynthpg3+cpg3);
	dcdhap=(vALDO-vDHAP-vG3PDH-vTIS)/cytosol;
	dce4p=(-vDAHPS-vE4P+vTA-vTKB)/cytosol;
	dcf6p=(-2.0*vMURSyNTH-vPFK+vPGI+vTA+vTKB-vf6P)/cytosol;
	dcfdp=(-vALDO+vPFK-vfdP)/cytosol;
	dcg1p=(-vG1PAT-vGLP+vPGM)/cytosol;
	dcg6p=(-vG6P-vG6PDH-vPGI-vPGM+65.0*vPTS)/cytosol;
	dcgap=(vALDO-vGAP-vGAPDH-vTA+vTIS+vTKA+vTKB+vTRPSYNTH)/cytosol;
	dcpep=(-vDAHPS+vENO-vPEP-vPK-65.0*vPTS-vSynth1-vpepCxylase)/cytosol;
	dcpg=(vG6PDH-vPG-vPGDH)/cytosol;
	dcpg2=(-vENO-vpg2+vrpGluMu)/cytosol;
	dcpg3=(-vPG3+vPGK-vrpGluMu-vsersynth)/cytosol;
	dcpgp=(vGAPDH-vPGK-vPGP)/cytosol;
	dcpyr=(vMethSynth-vPDH+vPK+65.0*vPTS-vSynth2+vTRPSYNTH-vpyr)/cytosol;
	dcrib5p=(-vPPK+vR5PI-vRIB5P-vTKA)/cytosol;
	dcribu5p=(vPGDH-vR5PI-vRibu5p-vRu5P)/cytosol;
	dcsed7p=(-vSED7P-vTA+vTKA)/cytosol;
	dcxyl5p=(vRu5P-vTKA-vTKB-vXYL5P)/cytosol;
	dcglcex=(vEXTER-vPTS)/extracellular;

	return(0);

}


/* Jacobian of the system (dfdx)*/
int amigoJAC_B2(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_B2(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}




void amigoRHS_get_OBS_B2(void* data){

#define	 cdhap   (amigo_model->sim_results[0][j]) 
#define	 ce4p    (amigo_model->sim_results[1][j]) 
#define	 cf6p    (amigo_model->sim_results[2][j]) 
#define	 cfdp    (amigo_model->sim_results[3][j]) 
#define	 cg1p    (amigo_model->sim_results[4][j]) 
#define	 cg6p    (amigo_model->sim_results[5][j]) 
#define	 cgap    (amigo_model->sim_results[6][j]) 
#define	 cpep    (amigo_model->sim_results[7][j]) 
#define	 cpg     (amigo_model->sim_results[8][j]) 
#define	 cpg2    (amigo_model->sim_results[9][j]) 
#define	 cpg3    (amigo_model->sim_results[10][j]) 
#define	 cpgp    (amigo_model->sim_results[11][j]) 
#define	 cpyr    (amigo_model->sim_results[12][j]) 
#define	 crib5p  (amigo_model->sim_results[13][j]) 
#define	 cribu5p (amigo_model->sim_results[14][j]) 
#define	 csed7p  (amigo_model->sim_results[15][j]) 
#define	 cxyl5p  (amigo_model->sim_results[16][j]) 
#define	 cglcex  (amigo_model->sim_results[17][j]) 

	int j;
	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 cpep_obs amigo_model->obs_results[0][j] 
		#define	 cg6p_obs amigo_model->obs_results[1][j] 
		#define	 cpyr_obs amigo_model->obs_results[2][j] 
		#define	 cf6p_obs amigo_model->obs_results[3][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){
				cpep_obs=cpep;
				cg6p_obs=cg6p;
				cpyr_obs=cpyr;
				cf6p_obs=cf6p;

			}

		 break;
		#define	 cglcex_obs amigo_model->obs_results[0][j] 

		 case 1:


			 for (j = 0; j < amigo_model->n_times; ++j){
				cglcex_obs=cglcex;

			}

		 break;
		#define	 cg1p_obs amigo_model->obs_results[0][j] 

		 case 2:


			 for (j = 0; j < amigo_model->n_times; ++j){
				cg1p_obs=cg1p;

			}

		 break;
		#define	 cpg_obs amigo_model->obs_results[0][j] 

		 case 3:


			 for (j = 0; j < amigo_model->n_times; ++j){
				cpg_obs=cpg;

			}

		 break;
		#define	 cfdp_obs amigo_model->obs_results[0][j] 
		#define	 cgap_obs amigo_model->obs_results[1][j] 

		 case 4:


			 for (j = 0; j < amigo_model->n_times; ++j){
				cfdp_obs=cfdp;
				cgap_obs=cgap;

			}

		 break;

	}

	return(amigo_model);

}



void amigoRHS_get_sens_OBS_B2(void* data){
	int j,k;

#define	 cdhap   (amigo_model->sens_results[0][j][k]) 
#define	 ce4p    (amigo_model->sens_results[1][j][k]) 
#define	 cf6p    (amigo_model->sens_results[2][j][k]) 
#define	 cfdp    (amigo_model->sens_results[3][j][k]) 
#define	 cg1p    (amigo_model->sens_results[4][j][k]) 
#define	 cg6p    (amigo_model->sens_results[5][j][k]) 
#define	 cgap    (amigo_model->sens_results[6][j][k]) 
#define	 cpep    (amigo_model->sens_results[7][j][k]) 
#define	 cpg     (amigo_model->sens_results[8][j][k]) 
#define	 cpg2    (amigo_model->sens_results[9][j][k]) 
#define	 cpg3    (amigo_model->sens_results[10][j][k]) 
#define	 cpgp    (amigo_model->sens_results[11][j][k]) 
#define	 cpyr    (amigo_model->sens_results[12][j][k]) 
#define	 crib5p  (amigo_model->sens_results[13][j][k]) 
#define	 cribu5p (amigo_model->sens_results[14][j][k]) 
#define	 csed7p  (amigo_model->sens_results[15][j][k]) 
#define	 cxyl5p  (amigo_model->sens_results[16][j][k]) 
#define	 cglcex  (amigo_model->sens_results[17][j][k]) 

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 cpep_obs amigo_model->sens_obs[0][j][k] 
		#define	 cg6p_obs amigo_model->sens_obs[1][j][k] 
		#define	 cpyr_obs amigo_model->sens_obs[2][j][k] 
		#define	 cf6p_obs amigo_model->sens_obs[3][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					cpep_obs=cpep;
					cg6p_obs=cg6p;
					cpyr_obs=cpyr;
					cf6p_obs=cf6p;
				}
			}
		 break;

		 case 1:

		#define	 cglcex_obs amigo_model->sens_obs[0][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					cglcex_obs=cglcex;
				}
			}
		 break;

		 case 2:

		#define	 cg1p_obs amigo_model->sens_obs[0][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					cg1p_obs=cg1p;
				}
			}
		 break;

		 case 3:

		#define	 cpg_obs amigo_model->sens_obs[0][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					cpg_obs=cpg;
				}
			}
		 break;

		 case 4:

		#define	 cfdp_obs amigo_model->sens_obs[0][j][k] 
		#define	 cgap_obs amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					cfdp_obs=cfdp;
					cgap_obs=cgap;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon_B2(void* data,realtype t, N_Vector y){
	AMIGO_model* amigo_model=(AMIGO_model*)data;


}
