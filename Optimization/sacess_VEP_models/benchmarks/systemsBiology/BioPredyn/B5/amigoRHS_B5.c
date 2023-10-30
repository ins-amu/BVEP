
/*
Check README.txt
*/

#include <AMIGO_problem.h>

int amigoRHS_B5(realtype t, N_Vector y, N_Vector ydot, void *data);

void amigo_Y_at_tcon_B5(void* amigo_model, realtype t, N_Vector y);

void amigoRHS_get_OBS_B5(void* data);

void amigoRHS_get_sens_OBS_B5(void* data);
 
int amigoRHS_B5(realtype t, N_Vector y, N_Vector ydot, void *data);


#define nik	Ith(y,0)
#define mkk4	Ith(y,1)
#define ask1	Ith(y,2)
#define map3k7	Ith(y,3)
#define mkk7	Ith(y,4)
#define tnfr	Ith(y,5)
#define egfr	Ith(y,6)
#define ph	Ith(y,7)
#define ex	Ith(y,8)
#define mek	Ith(y,9)
#define ras	Ith(y,10)
#define traf2	Ith(y,11)
#define ikk	Ith(y,12)
#define akt	Ith(y,13)
#define pi3k	Ith(y,14)
#define ikb	Ith(y,15)
#define nfkb	Ith(y,16)
#define cjun	Ith(y,17)
#define jnk	Ith(y,18)
#define map3k1	Ith(y,19)
#define erk	Ith(y,20)
#define raf1	Ith(y,21)
#define sos	Ith(y,22)
#define p38	Ith(y,23)
#define gsk3	Ith(y,24)
#define ap1	Ith(y,25)

#define dnik	Ith(ydot,0)
#define dmkk4	Ith(ydot,1)
#define dask1	Ith(ydot,2)
#define dmap3k7	Ith(ydot,3)
#define dmkk7	Ith(ydot,4)
#define dtnfr	Ith(ydot,5)
#define degfr	Ith(ydot,6)
#define dph	Ith(ydot,7)
#define dex	Ith(ydot,8)
#define dmek	Ith(ydot,9)
#define dras	Ith(ydot,10)
#define dtraf2	Ith(ydot,11)
#define dikk	Ith(ydot,12)
#define dakt	Ith(ydot,13)
#define dpi3k	Ith(ydot,14)
#define dikb	Ith(ydot,15)
#define dnfkb	Ith(ydot,16)
#define dcjun	Ith(ydot,17)
#define djnk	Ith(ydot,18)
#define dmap3k1	Ith(ydot,19)
#define derk	Ith(ydot,20)
#define draf1	Ith(ydot,21)
#define dsos	Ith(ydot,22)
#define dp38	Ith(ydot,23)
#define dgsk3	Ith(ydot,24)
#define dap1	Ith(ydot,25)

//"raf1" "erk"  "ap1"  "gsk3" "p38"  "nfkb"
//22     21      26      25      24      17


#define map3k7_n_nik	(*amigo_model).pars[0]
#define map3k7_k_nik	(*amigo_model).pars[1]
#define tau_nik         (*amigo_model).pars[2]
#define map3k7_n_mkk4	(*amigo_model).pars[3]
#define map3k7_k_mkk4	(*amigo_model).pars[4]
#define map3k1_n_mkk4	(*amigo_model).pars[5]
#define map3k1_k_mkk4	(*amigo_model).pars[6]
#define tau_mkk4        (*amigo_model).pars[7]
#define traf2_n_ask1	(*amigo_model).pars[8]
#define traf2_k_ask1	(*amigo_model).pars[9]
#define tau_ask1        (*amigo_model).pars[10]
#define traf2_n_map3k7	(*amigo_model).pars[11]
#define traf2_k_map3k7	(*amigo_model).pars[12]
#define tau_map3k7      (*amigo_model).pars[13]
#define ask1_n_mkk7     (*amigo_model).pars[14]
#define ask1_k_mkk7     (*amigo_model).pars[15]
#define map3k1_n_mkk7	(*amigo_model).pars[16]
#define map3k1_k_mkk7	(*amigo_model).pars[17]
#define tau_mkk7        (*amigo_model).pars[18]
#define tnfa_n_tnfr     (*amigo_model).pars[19]
#define tnfa_k_tnfr     (*amigo_model).pars[20]
#define tau_tnfr        (*amigo_model).pars[21]
#define egf_n_egfr      (*amigo_model).pars[22]
#define egf_k_egfr      (*amigo_model).pars[23]
#define tau_egfr        (*amigo_model).pars[24]
#define erk_n_ph        (*amigo_model).pars[25]
#define erk_k_ph        (*amigo_model).pars[26]
#define tau_ph          (*amigo_model).pars[27]
#define nfkb_n_ex       (*amigo_model).pars[28]
#define nfkb_k_ex       (*amigo_model).pars[29]
#define tau_ex          (*amigo_model).pars[30]
#define raf1_n_mek      (*amigo_model).pars[31]
#define raf1_k_mek      (*amigo_model).pars[32]
#define tau_mek         (*amigo_model).pars[33]
#define sos_n_ras       (*amigo_model).pars[34]
#define sos_k_ras       (*amigo_model).pars[35]
#define tau_ras         (*amigo_model).pars[36]
#define tnfr_n_traf2	(*amigo_model).pars[37]
#define tnfr_k_traf2	(*amigo_model).pars[38]
#define tau_traf2       (*amigo_model).pars[39]
#define nik_n_ikk       (*amigo_model).pars[40]
#define nik_k_ikk       (*amigo_model).pars[41]
#define tau_ikk         (*amigo_model).pars[42]
#define pi3k_n_akt      (*amigo_model).pars[43]
#define pi3k_k_akt      (*amigo_model).pars[44]
#define tau_akt         (*amigo_model).pars[45]
#define egfr_n_pi3k     (*amigo_model).pars[46]
#define egfr_k_pi3k     (*amigo_model).pars[47]
#define tau_pi3k        (*amigo_model).pars[48]
#define ex_n_ikb        (*amigo_model).pars[49]
#define ex_k_ikb        (*amigo_model).pars[50]
#define ikk_n_ikb       (*amigo_model).pars[51]
#define ikk_k_ikb       (*amigo_model).pars[52]
#define tau_ikb         (*amigo_model).pars[53]
#define ikb_n_nfkb      (*amigo_model).pars[54]
#define ikb_k_nfkb      (*amigo_model).pars[55]
#define tau_nfkb        (*amigo_model).pars[56]
#define jnk_n_cjun      (*amigo_model).pars[57]
#define jnk_k_cjun      (*amigo_model).pars[58]
#define tau_cjun        (*amigo_model).pars[59]
#define mkk7_n_jnk      (*amigo_model).pars[60]
#define mkk7_k_jnk      (*amigo_model).pars[61]
#define tau_jnk         (*amigo_model).pars[62]
#define ras_n_map3k1	(*amigo_model).pars[63]
#define ras_k_map3k1	(*amigo_model).pars[64]
#define tau_map3k1      (*amigo_model).pars[65]
#define mek_n_erk       (*amigo_model).pars[66]
#define mek_k_erk       (*amigo_model).pars[67]
#define tau_erk         (*amigo_model).pars[68]
#define ras_n_raf1      (*amigo_model).pars[69]
#define ras_k_raf1      (*amigo_model).pars[70]
#define tau_raf1        (*amigo_model).pars[71]
#define egfr_n_sos      (*amigo_model).pars[72]
#define egfr_k_sos      (*amigo_model).pars[73]
#define ph_n_sos        (*amigo_model).pars[74]
#define ph_k_sos        (*amigo_model).pars[75]
#define tau_sos         (*amigo_model).pars[76]
#define mkk4_n_p38      (*amigo_model).pars[77]
#define mkk4_k_p38      (*amigo_model).pars[78]
#define tau_p38         (*amigo_model).pars[79]
#define akt_n_gsk3      (*amigo_model).pars[80]
#define akt_k_gsk3      (*amigo_model).pars[81]
#define tau_gsk3        (*amigo_model).pars[82]
#define cjun_n_ap1      (*amigo_model).pars[83]
#define cjun_k_ap1      (*amigo_model).pars[84]
#define tau_ap1         (*amigo_model).pars[85]


#define map3k7_nik_w                (*amigo_model).pars[86]     //1
//Original is map3k7  AND mkk4
#define map3k7_mkk4_w               (*amigo_model).pars[87]     //0
#define map3k1_mkk4_w               (*amigo_model).pars[88]     //0
#define map3k7_map3k1_mkk4_AND_w    (*amigo_model).pars[89]     //1
#define traf2_ask1_w                (*amigo_model).pars[90]     //1
#define traf2_map3k7_w              (*amigo_model).pars[91]     //1
#define ask1_mkk7_w                 (*amigo_model).pars[92]     //1 
#define map3k1_mkk7_w               (*amigo_model).pars[93]     //1  
#define ask1_map3k1_mkk7_AND_w      (*amigo_model).pars[94]     //0    
#define tnfa_tnfr_w                 (*amigo_model).pars[95]     //1
#define egf_egfr_w                  (*amigo_model).pars[96]     //1
#define erk_ph_w                    (*amigo_model).pars[97]     //1 
#define nfkb_ex_w                   (*amigo_model).pars[98]     //1
#define raf1_mek_w                  (*amigo_model).pars[99]     //1 
#define sos_ras_w                   (*amigo_model).pars[100]    //1
#define tnfr_traf2_w                (*amigo_model).pars[101]    //1
#define nik_ikk_w                   (*amigo_model).pars[102]    //1 
#define pi3k_akt_w                  (*amigo_model).pars[103]    //1
#define egfr_pi3k_w                 (*amigo_model).pars[104]    //1
// Original model was ex OR NOT ikk  
#define ex_ikb_w                    (*amigo_model).pars[105]    //1  
#define ikk_ikb_w                   (*amigo_model).pars[106]    //1   
#define ex_ikk_ikb_AND_w            (*amigo_model).pars[107]    //0
#define ikb_nfkb_w                  (*amigo_model).pars[108]    //1
#define jnk_cjun_w                  (*amigo_model).pars[109]    //1
#define mkk7_jnk_w                  (*amigo_model).pars[110]    //1
#define ras_map3k1_w                (*amigo_model).pars[111]    //1
#define mek_erk_w                   (*amigo_model).pars[112]    //1 
#define ras_raf1_w                  (*amigo_model).pars[113]    //1
//define Original was egfr AND NOT ph    
#define egfr_sos_w                  (*amigo_model).pars[114]    //0    
#define ph_sos_w                    (*amigo_model).pars[115]    //0
#define egfr_ph_sos_AND_w           (*amigo_model).pars[116]    //1   
#define mkk4_p38_w                  (*amigo_model).pars[117]    //1
#define akt_gsk3_w                  (*amigo_model).pars[118]    //1
#define cjun_ap1_w                  (*amigo_model).pars[119]    //1    


//'map3k7_nik_w','map3k7_mkk4_w','map3k1_mkk4_w','map3k7_map3k1_mkk4_AND_w','traf2_ask1_w','traf2_map3k7_w','ask1_mkk7_w','map3k1_mkk7_w','ask1_map3k1_mkk7_AND_w','tnfa_tnfr_w','egf_egfr_w,'erk_ph_w','nfkb_ex_w','raf1_mek_w','sos_ras_w','tnfr_traf2_w','nik_ikk_w','pi3k_akt_w','egfr_pi3k_w','ex_ikb_w','ikk_ikb_w','ex_ikk_ikb_AND_w','ikb_nfkb_w','jnk_cjun_w','mkk7_jnk_w','ras_map3k1_w','mek_erk_w','ras_raf1_w','egfr_sos_w','ph_sos_w','egfr_ph_sos_AND_w','mkk4_p38_w','akt_gsk3_w','cjun_ap1_w'
//1,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1 
        
        



///////////////////////////////////////////////////////////////////////////////////////////

#define egf (*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+ \
(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim]
        
#define tnfa (*amigo_model).controls_v[1][(*amigo_model).index_t_stim]+ \
        (t-(*amigo_model).tlast)*(*amigo_model).slope[1][(*amigo_model).index_t_stim]
        
#define pi3k_inh (*amigo_model).controls_v[2][(*amigo_model).index_t_stim]+ \
        (t-(*amigo_model).tlast)*(*amigo_model).slope[2][(*amigo_model).index_t_stim]
        
#define raf1_inh (*amigo_model).controls_v[3][(*amigo_model).index_t_stim]+ \
        (t-(*amigo_model).tlast)*(*amigo_model).slope[3][(*amigo_model).index_t_stim]


double OR(double x1,double x2,double x3){
    
    double f=
	/*0 0 0*/0*(1-x1)*(1-x2)*(1-x3)+ 
	/*0 0 1*/(1-x1)*(1-x2)*(x3)+
	/*0 1 0*/(1-x1)*(x2)*(1-x3)+
	/*0 1 1*/(1-x1)*(x2)*(x3)+
	/*1 0 0*/(x1)*(1-x2)*(1-x3)+
	/*1 0 1*/(x1)*(1-x2)*(x3)+
	/*1 1 0*/(x1)*(x2)*(1-x3)+
	/*1 1 1*/(x1)*(x2)*(x3);
    
    return(f);
}

        
        
    int amigoRHS_B5(realtype t, N_Vector y, N_Vector ydot, void *data){
    
    double DNAfree,Dp,D,E,k_trans_f,k_deg_f,psi,R1,R2,fabc_mini,fabc_mini_neg,mrna_mini, r3_mini;
    
    int i;
    
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
   
       double map3k7_fhill_nik;
double map3k7_fhill_mkk4;
double map3k1_fhill_mkk4;
double traf2_fhill_ask1;
double traf2_fhill_map3k7;
double ask1_fhill_mkk7;
double map3k1_fhill_mkk7;
double tnfa_fhill_tnfr;
double egf_fhill_egfr;
double erk_fhill_ph;
double nfkb_fhill_ex;
double raf1_fhill_mek;
double sos_fhill_ras;
double tnfr_fhill_traf2;
double nik_fhill_ikk;
double pi3k_fhill_akt;
double egfr_fhill_pi3k;
double ex_fhill_ikb;
double ikk_fhill_ikb;
double ikb_fhill_nfkb;
double jnk_fhill_cjun;
double mkk7_fhill_jnk;
double ras_fhill_map3k1;
double mek_fhill_erk;
double ras_fhill_raf1;
double egfr_fhill_sos;
double ph_fhill_sos;
double mkk4_fhill_p38;
double akt_fhill_gsk3;
double cjun_fhill_ap1; 
        
        
    map3k7_fhill_nik=(pow(map3k7,map3k7_n_nik)/(pow(map3k7_k_nik,map3k7_n_nik)+pow(map3k7,map3k7_n_nik)))*(1+pow(map3k7_k_nik,map3k7_n_nik));
    
    //map3k7_nik_w 1
        
    dnik=(
            0*(1-map3k7_fhill_nik)+
            map3k7_nik_w * map3k7_fhill_nik
            -nik)*tau_nik;
    
    map3k7_fhill_mkk4=(pow(map3k7,map3k7_n_mkk4)/(pow(map3k7_k_mkk4,map3k7_n_mkk4)+pow(map3k7,map3k7_n_mkk4)))*(1+pow(map3k7_k_mkk4,map3k7_n_mkk4));
    map3k1_fhill_mkk4=(pow(map3k1,map3k1_n_mkk4)/(pow(map3k1_k_mkk4,map3k1_n_mkk4)+pow(map3k1,map3k1_n_mkk4)))*(1+pow(map3k1_k_mkk4,map3k1_n_mkk4));
    
        //Original is map3k7  AND mkk4
        //map3k7_mkk4_w 0
        //map3k1_mkk4_w 0  
        //map3k7_map3k1_mkk4_AND_w 1
    
    dmkk4=(
            0*(1-map3k7_fhill_mkk4)*(1-map3k1_fhill_mkk4)+
            map3k1_mkk4_w *(1-map3k7_fhill_mkk4)*map3k1_fhill_mkk4+
            map3k7_mkk4_w *map3k7_fhill_mkk4*(1-map3k1_fhill_mkk4)+
            OR(map3k7_mkk4_w,map3k1_mkk4_w,map3k7_map3k1_mkk4_AND_w)* map3k7_fhill_mkk4*map3k1_fhill_mkk4
            -mkk4)*tau_mkk4;
    
    
    
    traf2_fhill_ask1=(pow(traf2,traf2_n_ask1)/(pow(traf2_k_ask1,traf2_n_ask1)+pow(traf2,traf2_n_ask1)))*(1+pow(traf2_k_ask1,traf2_n_ask1));
    
    //traf2_ask1_w 1
        
    dask1=(
            0 * (1-traf2_fhill_ask1)+
            traf2_ask1_w * traf2_fhill_ask1
            -ask1)*tau_ask1;
    
    
    traf2_fhill_map3k7=(pow(traf2,traf2_n_map3k7)/(pow(traf2_k_map3k7,traf2_n_map3k7)+pow(traf2,traf2_n_map3k7)))*(1+pow(traf2_k_map3k7,traf2_n_map3k7));
    
    //traf2_map3k7_w 1   
    
    dmap3k7=(
            0*(1-traf2_fhill_map3k7)+
            traf2_map3k7_w *traf2_fhill_map3k7
            -map3k7)*tau_map3k7;
    
    
    ask1_fhill_mkk7=(pow(ask1,ask1_n_mkk7)/(pow(ask1_k_mkk7,ask1_n_mkk7)+pow(ask1,ask1_n_mkk7)))*(1+pow(ask1_k_mkk7,ask1_n_mkk7));
    map3k1_fhill_mkk7=(pow(map3k1,map3k1_n_mkk7)/(pow(map3k1_k_mkk7,map3k1_n_mkk7)+pow(map3k1,map3k1_n_mkk7)))*(1+pow(map3k1_k_mkk7,map3k1_n_mkk7));
    
        
        
    //ask1_mkk7_w             1 
    //map3k1_mkk7_w           1  
    //ask1_map3k1_mkk7_AND_w  0    
        
    dmkk7=(
            0*(1-ask1_fhill_mkk7)*(1-map3k1_fhill_mkk7)+
            map3k1_mkk7_w * (1-ask1_fhill_mkk7)*map3k1_fhill_mkk7+
            ask1_mkk7_w * ask1_fhill_mkk7*(1-map3k1_fhill_mkk7)+
            OR(ask1_mkk7_w,map3k1_mkk7_w,ask1_map3k1_mkk7_AND_w)*ask1_fhill_mkk7*map3k1_fhill_mkk7
            -mkk7)*tau_mkk7;
    
    
    tnfa_fhill_tnfr=(pow(tnfa,tnfa_n_tnfr)/(pow(tnfa_k_tnfr,tnfa_n_tnfr)+pow(tnfa,tnfa_n_tnfr)))*(1+pow(tnfa_k_tnfr,tnfa_n_tnfr));
    
    //tnfa_tnfr_w 1    
        
    dtnfr=(
            0*(1-tnfa_fhill_tnfr)+
            tnfa_tnfr_w * tnfa_fhill_tnfr
            -tnfr)*tau_tnfr;
    
   //egf_egfr_w 1 
        
    egf_fhill_egfr=(pow(egf,egf_n_egfr)/(pow(egf_k_egfr,egf_n_egfr)+pow(egf,egf_n_egfr)))*(1+pow(egf_k_egfr,egf_n_egfr));
    
    degfr=(
            0*(1-egf_fhill_egfr)+
            egf_egfr_w * egf_fhill_egfr
            -egfr)*tau_egfr;
    
    
    erk_fhill_ph=(pow(erk,erk_n_ph)/(pow(erk_k_ph,erk_n_ph)+pow(erk,erk_n_ph)))*(1+pow(erk_k_ph,erk_n_ph));
    
   
    //erk_ph_w 1       
        
    dph=(
            0 * (1-erk_fhill_ph)+
            erk_ph_w * erk_fhill_ph
            -ph)*tau_ph;
    
    nfkb_fhill_ex=(pow(nfkb,nfkb_n_ex)/(pow(nfkb_k_ex,nfkb_n_ex)+pow(nfkb,nfkb_n_ex)))*(1+pow(nfkb_k_ex,nfkb_n_ex));
    
    //nfkb_ex_w 1    
    
    dex=(
            0 * (1-nfkb_fhill_ex)+
            nfkb_ex_w * nfkb_fhill_ex
            -ex)*tau_ex;
    
    raf1_fhill_mek=(pow(raf1,raf1_n_mek)/(pow(raf1_k_mek,raf1_n_mek)+pow(raf1,raf1_n_mek)))*(1+pow(raf1_k_mek,raf1_n_mek));
    
    //raf1_mek_w    
    
    dmek=(
            0*(1-raf1_fhill_mek)+
            raf1_mek_w*raf1_fhill_mek
            -mek)*tau_mek;
    
    
    sos_fhill_ras=(pow(sos,sos_n_ras)/(pow(sos_k_ras,sos_n_ras)+pow(sos,sos_n_ras)))*(1+pow(sos_k_ras,sos_n_ras));
    
    //sos_ras_w 1        
        
    dras=(
            0*(1-sos_fhill_ras)+
            sos_ras_w * sos_fhill_ras
            -ras)*tau_ras;
    
    
    tnfr_fhill_traf2=(pow(tnfr,tnfr_n_traf2)/(pow(tnfr_k_traf2,tnfr_n_traf2)+pow(tnfr,tnfr_n_traf2)))*(1+pow(tnfr_k_traf2,tnfr_n_traf2));
    
    //tnfr_traf2_w 1   
        
    dtraf2=(
            0 * (1-tnfr_fhill_traf2)+
            tnfr_traf2_w  * tnfr_fhill_traf2
            -traf2)*tau_traf2;
    
    
    nik_fhill_ikk=(pow(nik,nik_n_ikk)/(pow(nik_k_ikk,nik_n_ikk)+pow(nik,nik_n_ikk)))*(1+pow(nik_k_ikk,nik_n_ikk));
    
    //nik_ikk_w 1   
        
    dikk=(
            0*(1-nik_fhill_ikk)+
            nik_ikk_w*nik_fhill_ikk
            -ikk)*tau_ikk;
    
    
    pi3k_fhill_akt=(pow(pi3k,pi3k_n_akt)/(pow(pi3k_k_akt,pi3k_n_akt)+pow(pi3k,pi3k_n_akt)))*(1+pow(pi3k_k_akt,pi3k_n_akt));
    
    //pi3k_akt_w 1
        
    dakt=(
            0 * (1-pi3k_fhill_akt)+
            pi3k_akt_w * pi3k_fhill_akt
            -akt)*tau_akt;
    
    
    egfr_fhill_pi3k=(pow(egfr,egfr_n_pi3k)/(pow(egfr_k_pi3k,egfr_n_pi3k)+pow(egfr,egfr_n_pi3k)))*(1+pow(egfr_k_pi3k,egfr_n_pi3k));
    
    //egfr_pi3k_w 1    
        
    dpi3k=(
            0 * (1-egfr_fhill_pi3k)+
            egfr_pi3k_w * egfr_fhill_pi3k
            -pi3k)*tau_pi3k*(1-pi3k_inh);
    
    
    ex_fhill_ikb=(pow(ex,ex_n_ikb)/(pow(ex_k_ikb,ex_n_ikb)+pow(ex,ex_n_ikb)))*(1+pow(ex_k_ikb,ex_n_ikb));
    ikk_fhill_ikb=(pow(ikk,ikk_n_ikb)/(pow(ikk_k_ikb,ikk_n_ikb)+pow(ikk,ikk_n_ikb)))*(1+pow(ikk_k_ikb,ikk_n_ikb));
    

    // Original model was ex OR NOT ikk  
    //ex_ikb_w  1  
    //ikk_ikb_w 1   
    //ex_ikk_ikb_AND_w 0    
        
    dikb=(
            ikk_ikb_w * (1-ex_fhill_ikb) * (1-ikk_fhill_ikb)+
            0 * (1-ex_fhill_ikb)*ikk_fhill_ikb+
            OR(ex_ikb_w,ikk_ikb_w,ex_ikk_ikb_AND_w) * ex_fhill_ikb * (1-ikk_fhill_ikb)+
            ex_ikb_w * ex_fhill_ikb * ikk_fhill_ikb
            -ikb)*tau_ikb;
    
    
    ikb_fhill_nfkb=(pow(ikb,ikb_n_nfkb)/(pow(ikb_k_nfkb,ikb_n_nfkb)+pow(ikb,ikb_n_nfkb)))*(1+pow(ikb_k_nfkb,ikb_n_nfkb));
    
    //ikb_nfkb_w 1        
        
    dnfkb=(
            ikb_nfkb_w * (1-ikb_fhill_nfkb)+
            0*ikb_fhill_nfkb
            -nfkb)*tau_nfkb;
    
    
        
    jnk_fhill_cjun=(pow(jnk,jnk_n_cjun)/(pow(jnk_k_cjun,jnk_n_cjun)+pow(jnk,jnk_n_cjun)))*(1+pow(jnk_k_cjun,jnk_n_cjun));
        
    //jnk_cjun_w 1    
    
    dcjun=(
            0 * (1-jnk_fhill_cjun)+
            jnk_cjun_w * jnk_fhill_cjun
            -cjun) * tau_cjun;
    
    
    mkk7_fhill_jnk=(pow(mkk7,mkk7_n_jnk)/(pow(mkk7_k_jnk,mkk7_n_jnk)+pow(mkk7,mkk7_n_jnk)))*(1+pow(mkk7_k_jnk,mkk7_n_jnk));
    
    //mkk7_jnk_w 1    
        
    djnk=(
            0 * (1-mkk7_fhill_jnk)+
            mkk7_jnk_w * mkk7_fhill_jnk
            -jnk) * tau_jnk;
    
    ras_fhill_map3k1=(pow(ras,ras_n_map3k1)/(pow(ras_k_map3k1,ras_n_map3k1)+pow(ras,ras_n_map3k1)))*(1+pow(ras_k_map3k1,ras_n_map3k1));
    
    //ras_map3k1_w 1
        
    dmap3k1=(
            0*(1-ras_fhill_map3k1)+
            ras_map3k1_w*ras_fhill_map3k1
            -map3k1)*tau_map3k1;
    
    
    mek_fhill_erk=(pow(mek,mek_n_erk)/(pow(mek_k_erk,mek_n_erk)+pow(mek,mek_n_erk)))*(1+pow(mek_k_erk,mek_n_erk));
    
    //mek_erk_w 1    
    
    derk=(
            0 * (1-mek_fhill_erk)+
            mek_erk_w * mek_fhill_erk
            -erk) * tau_erk;
    
    
    ras_fhill_raf1=(pow(ras,ras_n_raf1)/(pow(ras_k_raf1,ras_n_raf1)+pow(ras,ras_n_raf1)))*(1+pow(ras_k_raf1,ras_n_raf1));
    
    //ras_raf1_w 1    
        
    draf1=(
            0 * (1-ras_fhill_raf1)+
            ras_raf1_w * ras_fhill_raf1
            -raf1) * tau_raf1 * (1-raf1_inh);
    
    
    egfr_fhill_sos=(pow(egfr,egfr_n_sos)/(pow(egfr_k_sos,egfr_n_sos)+pow(egfr,egfr_n_sos)))*(1+pow(egfr_k_sos,egfr_n_sos));
    ph_fhill_sos=(pow(ph,ph_n_sos)/(pow(ph_k_sos,ph_n_sos)+pow(ph,ph_n_sos)))*(1+pow(ph_k_sos,ph_n_sos));
    
        
        
    //Original was egfr AND NOT ph    
    //egfr_sos_w 0    
    //ph_sos_w   0
    //egfr_ph_sos_AND_w 1   
        
    dsos=(
            ph_sos_w * (1-egfr_fhill_sos) * (1-ph_fhill_sos)+
            0 * (1-egfr_fhill_sos) * ph_fhill_sos+
            OR(egfr_sos_w,ph_sos_w,egfr_ph_sos_AND_w) * egfr_fhill_sos * (1-ph_fhill_sos)+
            egfr_sos_w * egfr_fhill_sos * ph_fhill_sos
            -sos) * tau_sos;
    
    
    mkk4_fhill_p38=(pow(mkk4,mkk4_n_p38)/(pow(mkk4_k_p38,mkk4_n_p38)+pow(mkk4,mkk4_n_p38)))*(1+pow(mkk4_k_p38,mkk4_n_p38));
    
    //mkk4_p38_w 1    
    
    dp38=(
            0 * (1-mkk4_fhill_p38)+
            mkk4_p38_w * mkk4_fhill_p38
            -p38) * tau_p38;
    
    
    akt_fhill_gsk3=(pow(akt,akt_n_gsk3)/(pow(akt_k_gsk3,akt_n_gsk3)+pow(akt,akt_n_gsk3)))*(1+pow(akt_k_gsk3,akt_n_gsk3));
    
    //akt_gsk3_w 1    
        
    dgsk3=(
            akt_gsk3_w * (1-akt_fhill_gsk3)+
            0 * akt_fhill_gsk3
            -gsk3) * tau_gsk3;
    
    
    cjun_fhill_ap1=(pow(cjun,cjun_n_ap1)/(pow(cjun_k_ap1,cjun_n_ap1)+pow(cjun,cjun_n_ap1)))*(1+pow(cjun_k_ap1,cjun_n_ap1));
    
    //cjun_ap1_w 1    
        
    dap1=(
            0 * (1-cjun_fhill_ap1)+
            cjun_ap1_w * cjun_fhill_ap1
            -ap1) * tau_ap1;
           
    
    
    return(0);
}
    
    int amigoJAC_B5(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        return(0);
    }
    
    int amigoSensRHS_B5(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
        return(0);
    };
        
void amigoRHS_get_OBS_B5(void* data){

}

void amigoRHS_get_sens_OBS_B5(void* data){

}

void amigo_Y_at_tcon_B5(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
}
    