#include <AMIGO_model.h>


#define axOD 1
#define axACT 2
#define axGLC 3
#define axACoA 4
#define axAKG 5
#define axcAMP 6
#define axFBP 7
#define axG6P 8
#define axGLX 9
#define axICT 10
#define axMAL 11
#define axOAA 12
#define axPEP 13
#define axPG3 14
#define axPYR 15
#define axAceA 16
#define axAceB 17
#define axAceK 18
#define axAcoa2act 19
#define axAcs 20
#define axAkg2mal 21
#define axCAMPdegr 22
#define axCya 23
#define axEmp 24
#define axEno 25
#define axFdp 26
#define axGltA 27
#define axIcd 28
#define axIcdP 29
#define axMdh 30
#define axMaeAB 31
#define axPckA 32
#define axPdh 33
#define axPfkA 34
#define axPpc 35
#define axPpsA 36
#define axPykF 37
#define axEIIA 38
#define axEIIAP 39
#define axEIICB 40
#define axCra 41
#define axCraFBP 42
#define axCrp 43
#define axCrpcAMP 44
#define axIclR 45
#define axPdhR 46
#define axPdhRPYR 47

#define afenvgrowth 1
#define afenvGLCup 2
#define afenvACTup 3
#define afenvACTex 4
#define afeAceA 5
#define afeAceB 6
#define afeAceKKi 7
#define afeAceKPh 8
#define afeAcoa2act 9
#define afeAcs 10
#define afeAkg2mal 11
#define afeCAMPdegr 12
#define afeCya 13
#define afeEmp 14
#define afeEno 15
#define afeFdp 16
#define afeGltA 17
#define afeIcd 18
#define afeMaeAB 19
#define afeMdh 20
#define afePckA 21
#define afePdh 22
#define afePfkA 23
#define afePpc 24
#define afePpsA 25
#define afePykF 26
#define afptsr1 27
#define afptsr4 28
#define aftfCra 29
#define aftfCrp 30
#define aftfPdhR 31
#define afgaceA 32
#define afgaceB 33
#define afgaceK 34
#define afgacoa2act 35
#define afgacs 36
#define afgakg2mal 37
#define afgcampdegr 38
#define afgcra 39
#define afgcrp 40
#define afgcya 41
#define afgEIIA 42
#define afgEIICB 43
#define afgemp 44
#define afgeno 45
#define afgfdp 46
#define afggltA 47
#define afgicd 48
#define afgiclr 49
#define afgmdh 50
#define afgmaeAB 51
#define afgpckA 52
#define afgpdh 53
#define afgpdhr 54
#define afgpfkA 55
#define afgppc 56
#define afgppsA 57
#define afgpykF 58
#define afdACoA 59
#define afdAKG 60
#define afdcAMP 61
#define afdFBP 62
#define afdG6P 63
#define afdGLX 64
#define afdICT 65
#define afdMAL 66
#define afdOAA 67
#define afdPEP 68
#define afdPG3 69
#define afdPYR 70
#define afdAceA 71
#define afdAceB 72
#define afdAceK 73
#define afdAcoa2act 74
#define afdAcs 75
#define afdAkg2mal 76
#define afdCAMPdegr 77
#define afdCra 78
#define afdCraFBP 79
#define afdCrp 80
#define afdCrpcAMP 81
#define afdCya 82
#define afdEIIA 83
#define afdEIIAP 84
#define afdEIICB 85
#define afdEmp 86
#define afdEno 87
#define afdFdp 88
#define afdGltA 89
#define afdIcd 90
#define afdIcdP 91
#define afdIclR 92
#define afdMdh 93
#define afdMaeAB 94
#define afdPckA 95
#define afdPdh 96
#define afdPdhR 97
#define afdPdhRPYR 98
#define afdPfkA 99
#define afdPpc 100
#define afdPpsA 101
#define afdPykF 102
#define afbmACoA 103
#define afbmAKG 104
#define afbmG6P 105
#define afbmOAA 106
#define afbmPEP 107
#define afbmPG3 108
#define afbmPYR 109

// Definition of parameter aliases

#define apenvMACT 1
#define apenvMGLC  2
#define apenvuc  3
#define apeAceAkcat  4
#define apeAceAn  5
#define apeAceAL  6
#define apeAceAKict  7
#define apeAceAKpep  8
#define apeAceAKpg3  9
#define apeAceAKakg  10
#define apeAceBkcat  11
#define apeAceBKglx  12
#define apeAceBKacoa  13
#define apeAceBKglxacoa  14
#define apeAceKkcatki  15
#define apeAceKkcatph  16
#define apeAceKn  17
#define apeAceKL  18
#define apeAceKKicd  19
#define apeAceKKicdP  20
#define apeAceKKpep  21
#define apeAceKKpyr  22
#define apeAceKKoaa  23
#define apeAceKKglx  24
#define apeAceKKakg  25
#define apeAceKKpg3  26
#define apeAceKKict  27
#define apeAcoa2actkcat  28
#define apeAcoa2actn  29
#define apeAcoa2actL  30
#define apeAcoa2actKacoa  31
#define apeAcoa2actKpyr  32
#define apeAcskcat  33
#define apeAcsKact  34
#define apeAkg2malkcat  35
#define apeAkg2malKakg  36
#define apeCAMPdegrkcat  37
#define apeCAMPdegrKcAMP  38
#define apeCyakcat  39
#define apeCyaKEIIA  40
#define apeEmpkcatf  41
#define apeEmpkcatr  42
#define apeEmpKfbp  43
#define apeEmpKpg3  44
#define apeEnokcatf  45
#define apeEnokcatr  46
#define apeEnoKpg3  47
#define apeEnoKpep  48
#define apeFdpkcat 49
#define apeFdpn  50
#define apeFdpL  51
#define apeFdpKfbp  52
#define apeFdpKpep  53
#define apeGltAkcat  54
#define apeGltAKoaa  55
#define apeGltAKacoa  56
#define apeGltAKoaaacoa  57
#define apeGltAKakg  58
#define apeIcdkcat  59
#define apeIcdn  60
#define apeIcdL  61
#define apeIcdKict  62
#define apeIcdKpep  63
#define apeMdhkcat  64
#define apeMdhn  65
#define apeMdhKmal  66
#define apeMaeABkcat  67
#define apeMaeABn  68
#define apeMaeABL  69
#define apeMaeABKmal  70
#define apeMaeABKacoa  71
#define apeMaeABKcamp  72
#define apePckAkcat  73
#define apePckAKoaa  74
#define apePckAKpep  75
#define apePdhkcat  76
#define apePdhn  77
#define apePdhL  78
#define apePdhKpyr  79
#define apePdhKpyrI  80
#define apePdhKglx  81
#define apePfkAkcat  82
#define apePfkAn  83
#define apePfkAL  84
#define apePfkAKg6p  85
#define apePfkAKpep  86
#define apePpckcat  87
#define apePpcn  88
#define apePpcL  89
#define apePpcKpep  90
#define apePpcKfbp  91
#define apePpsAkcat  92
#define apePpsAn  93
#define apePpsAL  94
#define apePpsAKpyr  95
#define apePpsAKpep  96
#define apePykFkcat  97
#define apePykFn  98
#define apePykFL  99
#define apePykFKpep  100
#define apePykFKfbp  101
#define apptsk1  102
#define apptskm1  103
#define apptsk4  104
#define apptsKEIIA  105
#define apptsKglc  106
#define aptfCrascale  107
#define aptfCrakfbp  108
#define aptfCran  109
#define aptfCrpscale  110
#define aptfCrpkcamp  111
#define aptfCrpn  112
#define aptfPdhRscale  113
#define aptfPdhRkpyr  114
#define aptfPdhRn  115
#define apgaceBAKvcraunbound  116
#define apgaceBAKvcrabound  117
#define apgaceBAKKcra  118
#define apgaceBAKaceBfactor  119
#define apgaceBAKaceKfactor  120
#define apgaceBAKKDNA  121
#define apgaceBAKKP  122
#define apgaceBAKKPprime  123
#define apgaceBAKKG  124
#define apgaceBAKL  125
#define apgaceBAKkcaticlr  126
#define apgaceBAKDNA  127
#define apgaceBAKvcrpbound  128
#define apgaceBAKvcrpunbound  129
#define apgaceBAKKcrp  130
#define apgacsvcrpunbound  131
#define apgacsvcrpbound  132
#define apgacsn  133
#define apgacsKcrp  134
#define apgakg2malvcrpunbound  135
#define apgakg2malvcrpbound  136
#define apgakg2malKcrp  137
#define apgakg2maln  138
#define apgempvcraunbound  139
#define apgempvcrabound  140
#define apgempKcra  141
#define apgempvcrpunbound  142
#define apgempvcrpbound  143
#define apgempKcrp  144
#define apgenovcraunbound  145
#define apgenovcrabound  146
#define apgenoKcra  147
#define apgfdpvcraunbound  148
#define apgfdpvcrabound  149
#define apgfdpKcra  150
#define apggltAvcrpunbound  151
#define apggltAvcrpbound  152
#define apggltAKcrp  153
#define apggltAn  154
#define apgicdvcraunbound  155
#define apgicdvcrabound  156
#define apgicdKcra  157
#define apgmdhvcrpunbound  158
#define apgmdhvcrpbound  159
#define apgmdhKcrp  160
#define apgpckAvcraunbound  161
#define apgpckAvcrabound  162
#define apgpckAKcra  163
#define apgpdhvpdhrunbound  164
#define apgpdhvpdhrbound  165
#define apgpdhKpdhr  166
#define apgpfkAvcraunbound  167
#define apgpfkAvcrabound  168
#define apgpfkAKcra  169
#define apgppsAvcraunbound  170
#define apgppsAvcrabound  171
#define apgppsAKcra  172
#define apgpykFvcraunbound  173
#define apgpykFvcrabound  174
#define apgpykFKcra  175
#define apdkdegr  176
#define apbmkexpr  177
#define apbmmuACT  178
#define apbmmuGLC  179
#define apbmGLCACoA  180
#define apbmGLCAKG  181
#define apbmGLCG6P  182
#define apbmGLCOAA  183
#define apbmGLCPEP  184
#define apbmGLCPG3  185
#define apbmGLCPYR  186
#define apbmACTACoA  187
#define apbmACTAKG  188
#define apbmACTG6P  189
#define apbmACTOAA  190
#define apbmACTPEP  191
#define apbmACTPG3  192
#define apbmACTPYR  193

#define NSTATES 47
#define NREACS 109
#define NPARS 193

double f[NREACS];
double ssACT[NSTATES];
double ssGLC[NSTATES];

#define x(i) ( NV_DATA_S(y)[i-1] )
#define xdot ( NV_DATA_S(ydot) )
#define p(i) ( amigo_model->pars[i-1] )
#define f(i) ( f[i-1] )
#define ssGLC(i) ( ssGLC[i-1] )
#define ssACT(i) ( ssACT[i-1] )

double S[NSTATES*NREACS];

#define S(i,j) (S[(j-1)*NSTATES+(i-1)])

int ode_step(realtype t, N_Vector y, N_Vector ydot, void *data);

int amigoRHS_B3(realtype t, N_Vector y, N_Vector ydot, void *data){
    
    double  alpha = 1, beta =0;
    int l=NSTATES, n=1 ,m=NREACS;
    int  lda=NSTATES, ldb=1, ldc=NSTATES;

    char transa='N';
    char transb='T';
    
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
    double alphaGLC;
    
    double alphaACT;
    double mu;
    double kbmACoA;
    double kbmAKG;
    double kbmG6P;
    double kbmOAA;
    double kbmPEP;
    double kbmPG3;
    double kbmPYR;
    
    double SSMaeAB;
    double SSPpc;
    


    ssACT(axACoA) = 2.045970006;
    ssACT(axAKG) = 1.11161563;
    ssACT(axcAMP) = 4.021170806;
    ssACT(axFBP) = 0.272276143;
    ssACT(axG6P) = 1.145931991;
    ssACT(axGLX) = 1.322800153;
    ssACT(axICT) = 1.480239304;
    ssACT(axMAL) = 6.396156035;
    ssACT(axOAA) = 0.064573174;
    ssACT(axPEP) = 0.557008135;
    ssACT(axPG3) = 1.293579693;
    ssACT(axPYR) = 0.037723095;
    ssACT(axAceA) = 0.101787757;
    ssACT(axAceB) = 0.030536327;
    ssACT(axAceK) = 0.003053633;
    ssACT(axAcoa2act) = 0.001;
    ssACT(axAcs) = 0.010144567*0.000036201/0.001096222;
    ssACT(axAkg2mal) = 0.002192398;
    ssACT(axCAMPdegr) = 0.001;
    ssACT(axCya) = 0.001;
    ssACT(axEmp) = 0.009751833*0.011389032/0.011515593;
    ssACT(axEno) = 0.006304314*0.011389032/0.011552813;
    ssACT(axFdp) = 0.000513512*0.000074810/0.000157492;
    ssACT(axGltA) = 0.003539257*0.000292771/0.001029612;
    ssACT(axIcd) = 0.002404566;
    ssACT(axIcdP) = 0.007516755;
    ssACT(axMdh) = 0.010969029*0.000491491/0.00345727;
    ssACT(axMaeAB) = 0.003399346;
    ssACT(axPckA) = 0.018902966*0.000336947/0.002290892;
    ssACT(axPdh) = 0.001760705*0.001/0.004647401;
    ssACT(axPfkA) = 8.89703E-05*0.000242131/0.000143816;
    ssACT(axPpc) = 0.000279893*0.000377962/0.000999714;
    ssACT(axPpsA) = 0.012844496;
    ssACT(axPykF) = 0.001305745*0.002501893/0.005977168;
    ssACT(axCra) = 0.007009039;
    ssACT(axCraFBP) = 0.000280931;
    ssACT(axCrp) = 0.001327161;
    ssACT(axCrpcAMP) = 0.005962839;
    ssACT(axIclR) = 0.00729;
    ssACT(axPdhR) = 0.005926738;
    ssACT(axPdhRPYR) = 0.001363262;
    ssACT(axEIIA) = 0.002631995;
    ssACT(axEIIAP) = 0.097368005;
    ssACT(axEIICB) = 0.003;
    
    ssGLC(axACoA) = 0.351972298;
    ssGLC(axAKG) = 0.191190619;
    ssGLC(axcAMP) = 0.202804098;
    ssGLC(axFBP) = 6.57504207;
    ssGLC(axG6P) = 1.908140784;
    ssGLC(axGLX) = 5.70593E-09;
    ssGLC(axICT) = 0.001408116;
    ssGLC(axMAL) = 3.278779135;
    ssGLC(axOAA) = 0.050535354;
    ssGLC(axPEP) = 0.210455879;
    ssGLC(axPG3) = 5.720977255;
    ssGLC(axPYR) = 0.863278018;
    ssGLC(axAceA) = 0.00472323;
    ssGLC(axAceB) = 0.001416969;
    ssGLC(axAceK) = 0.000141697;
    ssGLC(axAcoa2act) = 0.001;
    ssGLC(axAcs) = 0.000036201;//old: 0.001096222;
    ssGLC(axAkg2mal) = 0.001026848;
    ssGLC(axCAMPdegr) = 0.001;
    ssGLC(axCya) = 0.001;
    ssGLC(axEmp) = 0.011389032;//old: 0.011515593;
    ssGLC(axEno) = 0.011389032;//old: 0.011552813;
    ssGLC(axFdp) = 0.000074810;//old: 0.000157492;
    ssGLC(axGltA) = 0.000292771;//old: 0.001029612;
    ssGLC(axIcd) = 0.004290789;
    ssGLC(axIcdP) = 0.000220477;
    ssGLC(axMdh) = 0.000491491;// old: 0.00345727;
    ssGLC(axMaeAB) = 0.000999714;
    ssGLC(axPckA) = 0.000336947;//old: 0.002290892;
    ssGLC(axPdh) = 0.001;// old: 0.004647401;
    ssGLC(axPfkA) = 0.000242131;//old: 0.000143816;
    ssGLC(axPpc) = 0.000377962;//old: 0.000999714;
    ssGLC(axPpsA) = 0.000987493;
    ssGLC(axPykF) = 0.002501893;//old: 0.005977168;
    ssGLC(axCra) = 0.000299098;
    ssGLC(axCraFBP) = 0.006990902;
    ssGLC(axCrp) = 0.005943273;
    ssGLC(axCrpcAMP) = 0.001346727;
    ssGLC(axIclR) = 0.00729;
    ssGLC(axPdhR) = 0.001163813;
    ssGLC(axPdhRPYR) = 0.006126187;
    ssGLC(axEIIA) = 0.09647707;
    ssGLC(axEIIAP) = 0.00352292;
    ssGLC(axEIICB) = 0.003;
    
    // Write steady state values into auxiliary variables
    // The two auxiliary variables ssACT and ssGLC contain the
    // steady states on acetate and glucose, respectively.
    // Remember that these values, which have been obtained with
    // simulations on acetate and glucose as sole carbon sources,
    // depend on the chosen parameter set.
    // The sole purpose of these two auxiliary variables is to allow
    // for the initial conditions to be set to either steady state.
    // An initial condition equal to a steady state is interesting
    // when a transition to the other steady state is investigated.
    
    alphaGLC = x(axGLC)/(x(axGLC)+p(apptsKglc));
    alphaACT = x(axACT)/(x(axACT)+p(apeAcsKact))*(1-x(axGLC)/(x(axGLC)+p(apptsKglc)));
    
    // Calculate the growth rate 'mu'
    mu = alphaGLC*p(apbmmuGLC) + alphaACT*p(apbmmuACT);
    
    // Calculate the first order rate constants of the seven biomass reactions
    kbmACoA = alphaGLC*p(apbmGLCACoA) +  alphaACT*p(apbmACTACoA);
    kbmAKG = alphaGLC*p(apbmGLCAKG)   +  alphaACT*p(apbmACTAKG);
    kbmG6P = alphaGLC*p(apbmGLCG6P)   +  alphaACT*p(apbmACTG6P);
    kbmOAA = alphaGLC*p(apbmGLCOAA)   +  alphaACT*p(apbmACTOAA);
    kbmPEP = alphaGLC*p(apbmGLCPEP)   +  alphaACT*p(apbmACTPEP);
    kbmPG3 = alphaGLC*p(apbmGLCPG3)   +  alphaACT*p(apbmACTPG3);
    kbmPYR = alphaGLC*p(apbmGLCPYR)   +  alphaACT*p(apbmACTPYR);
    
    f(afbmACoA) = kbmACoA * x(axACoA);
    f(afbmAKG)  = kbmAKG  * x(axAKG);
    f(afbmG6P)  = kbmG6P  * x(axG6P);
    f(afbmOAA)  = kbmOAA  * x(axOAA);
    f(afbmPEP)  = kbmPEP  * x(axPEP);
    f(afbmPG3)  = kbmPG3  * x(axPG3);
    f(afbmPYR)  = kbmPYR  * x(axPYR);
    
    //Calculate the actual steady state levels of MaeAB and Ppc
    SSMaeAB = alphaGLC*ssGLC(axMaeAB) + alphaACT*ssACT(axMaeAB);
    SSPpc = alphaGLC*ssGLC(axPpc) + alphaACT*ssACT(axPpc);
    
    S(axOD,afenvgrowth) = 1;
    S(axGLC,afenvGLCup) = -1;
    S(axACT,afenvACTup) = -1;
    S(axACT,afenvACTex) = 1;
    S(axAceA,afgaceA) = 1;
    S(axAceA,afdAceA) = -1;
    S(axAceB,afgaceB) = 1;
    S(axAceB,afdAceB) = -1;
    S(axAceK,afgaceK) = 1;
    S(axAceK,afdAceK) = -1;
    S(axAcoa2act,afgacoa2act) = 1;
    S(axAcoa2act,afdAcoa2act) = -1;
    S(axAcs,afgacs) = 1;
    S(axAcs,afdAcs) = -1;
    S(axAkg2mal,afgakg2mal) = 1;
    S(axAkg2mal,afdAkg2mal) = -1;
    S(axCAMPdegr,afgcampdegr) = 1;
    S(axCAMPdegr,afdCAMPdegr) = -1;
    S(axCra,afgcra) = 1;
    S(axCra,afdCra) = -1;
    S(axCraFBP,afdCraFBP) = -1;
    S(axCrp,afgcrp) = 1;
    S(axCrp,afdCrp) = -1;
    S(axCrpcAMP,afdCrpcAMP) = -1;
    S(axCya,afgcya) = 1;
    S(axCya,afdCya) = -1;
    S(axEmp,afgemp) = 1;
    S(axEmp,afdEmp) = -1;
    S(axEno,afgeno) = 1;
    S(axEno,afdEno) = -1;
    S(axFdp,afgfdp) = 1;
    S(axFdp,afdFdp) = -1;
    S(axGltA,afggltA) = 1;
    S(axGltA,afdGltA) = -1;
    S(axIcd,afgicd) = 1;
    S(axIcd,afdIcd) = -1;
    S(axIcdP,afdIcdP) = -1;
    S(axIclR,afgiclr) = 1;
    S(axIclR,afdIclR) = -1;
    S(axMdh,afgmdh) = 1;
    S(axMdh,afdMdh) = -1;
    S(axMaeAB,afgmaeAB) = 1;
    S(axMaeAB,afdMaeAB) = -1;
    S(axPckA,afgpckA) = 1;
    S(axPckA,afdPckA) = -1;
    S(axPdh,afgpdh) = 1;
    S(axPdh,afdPdh) = -1;
    S(axPdhR,afgpdhr) = 1;
    S(axPdhR,afdPdhR) = -1;
    S(axPdhRPYR,afdPdhRPYR) = -1;
    S(axPfkA,afgpfkA) = 1;
    S(axPfkA,afdPfkA) = -1;
    S(axPpc,afgppc) = 1;
    S(axPpc,afdPpc) = -1;
    S(axPpsA,afgppsA) = 1;
    S(axPpsA,afdPpsA) = -1;
    S(axPykF,afgpykF) = 1;
    S(axPykF,afdPykF) = -1;
    S(axIcd,afeAceKKi) = -1;
    S(axIcd,afeAceKPh) = 1;
    S(axIcdP,afeAceKKi) = 1;
    S(axIcdP,afeAceKPh) = -1;
    S(axEIIA,afptsr1) = -1;
    S(axEIIA,afptsr4) = 1;
    S(axEIIA,afgEIIA) = 1;
    S(axEIIA,afdEIIA) = -1;
    S(axEIIAP,afptsr1) = 1;
    S(axEIIAP,afptsr4) = -1;
    S(axEIIAP,afdEIIAP) = -1;
    S(axEIICB,afgEIICB) = 1;
    S(axEIICB,afdEIICB) = -1;
    S(axCraFBP,aftfCra) = 1;
    S(axCra,aftfCra) = -1;
    S(axCrpcAMP,aftfCrp) = 1;
    S(axCrp,aftfCrp) = -1;
    S(axPdhRPYR,aftfPdhR) = 1;
    S(axPdhR,aftfPdhR) = -1;
    S(axACoA,afeAcs) = 1;
    S(axACoA,afePdh) = 1;
    S(axACoA,afeAcoa2act) = -1;
    S(axACoA,afeGltA) = -1;
    S(axACoA,afeAceB) = -1;
    S(axACoA,afdACoA) = -1;
    S(axACoA,afbmACoA) = -1;
    S(axAKG,afeIcd) = 1;
    S(axAKG,afeAceA) = 1;
    S(axAKG,afeAkg2mal) = -1;
    S(axAKG,afdAKG) = -1;
    S(axAKG,afbmAKG) = -1;
    S(axcAMP,afeCya) = 1;
    S(axcAMP,afeCAMPdegr) = -1;
    S(axcAMP,afdcAMP) = -1;
    S(axFBP,afePfkA) = 1;
    S(axFBP,afeEmp) = -0.5;
    S(axFBP,afeFdp) = -1;
    S(axFBP,afdFBP) = -1;
    S(axG6P,afptsr4) = 1;
    S(axG6P,afeFdp) = 1;
    S(axG6P,afePfkA) = -1;
    S(axG6P,afdG6P) = -1;
    S(axG6P,afbmG6P) = -1;
    S(axGLX,afeAceA) = 1;
    S(axGLX,afeAceB) = -1;
    S(axGLX,afdGLX) = -1;
    S(axICT,afeGltA) = 1;
    S(axICT,afeAceA) = -1;
    S(axICT,afeIcd) = -1;
    S(axICT,afdICT) = -1;
    S(axMAL,afeAceB) = 1;
    S(axMAL,afeAkg2mal) = 1;
    S(axMAL,afeMaeAB) = -1;
    S(axMAL,afeMdh) = -1;
    S(axMAL,afdMAL) = -1;
    S(axOAA,afePckA) = -1;
    S(axOAA,afeGltA) = -1;
    S(axOAA,afePpc) = 1;
    S(axOAA,afeMdh) = 1;
    S(axOAA,afdOAA) = -1;
    S(axOAA,afbmOAA) = -1;
    S(axPEP,afePckA) = 1;
    S(axPEP,afePpsA) = 1;
    S(axPEP,afePykF) = -1;
    S(axPEP,afePpc) = -1;
    S(axPEP,afeEno) = 1;
    S(axPEP,afptsr1) = -1;
    S(axPEP,afdPEP) = -1;
    S(axPEP,afbmPEP) = -1;
    S(axPG3,afeEmp) = 1;
    S(axPG3,afeEno) = -1;
    S(axPG3,afdPG3) = -1;
    S(axPG3,afbmPG3) = -1;
    S(axPYR,afeMaeAB) = 1;
    S(axPYR,afePykF) = 1;
    S(axPYR,afePpsA) = -1;
    S(axPYR,afePdh) = -1;
    S(axPYR,afptsr1) = 1;
    S(axPYR,afdPYR) = -1;
    S(axPYR,afbmPYR) = -1;
    
    
    
    
    // Protein phosphorylation rates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PTS phosphorylation kinetics
    f(afptsr1) = p(apptsk1)*x(axEIIA)*x(axPEP)-p(apptskm1)*x(axEIIAP)*x(axPYR);
    
    f(afptsr4) = (p(apptsk4)*x(axEIICB)*x(axEIIAP)*x(axGLC))/((p(apptsKEIIA)+x(axEIIAP))*(p(apptsKglc)+x(axGLC)));
    
    // AceKki kinetics: MWC, substrate: Icd, inhibitors: GLX, ICT, OAA,
    // PEP, PG3, PYR, AKG
    f(afeAceKKi) = (p(apeAceKkcatki)*x(axAceK)*x(axIcd)*pow(x(axIcd)/p(apeAceKKicd)+1.0,p(apeAceKn)-1.0))/(p(apeAceKKicd)*(p(apeAceKL)*pow(x(axAKG)/p(apeAceKKakg)+x(axGLX)/p(apeAceKKglx)+x(axICT)/p(apeAceKKict)+x(axOAA)/p(apeAceKKoaa)+x(axPG3)/p(apeAceKKpg3)+x(axPEP)/p(apeAceKKpep)+x(axPYR)/p(apeAceKKpyr)+1.0,p(apeAceKn))+pow(x(axIcd)/p(apeAceKKicd)+1.0,p(apeAceKn))));
    
    // AceKph kinetics: MWC, substrate: IcdP, activators: OAA, PEP,
    // PG3, PYR, AKG
    f(afeAceKPh) = (p(apeAceKkcatph)*x(axAceK)*x(axIcdP)*pow(x(axIcdP)/p(apeAceKKicdP)+1.0,p(apeAceKn)-1.0))/(p(apeAceKKicdP)*(p(apeAceKL)*pow(x(axAKG)/p(apeAceKKakg)+x(axOAA)/p(apeAceKKoaa)+x(axPG3)/p(apeAceKKpg3)+x(axPEP)/p(apeAceKKpep)+x(axPYR)/p(apeAceKKpyr)+1.0,-p(apeAceKn))+pow(x(axIcdP)/p(apeAceKKicdP)+1.0,p(apeAceKn))));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    //// Metabolite- transcription factor binding rates
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // The IclR-GLX-PYR binding state is incorporated into the gene
    // expression kinetics of the aceBAK operon
    
    // Cra-FBP binding kinetics: Hill
    f(aftfCra) = -p(aptfCrascale)*(x(axCraFBP)-(pow(x(axFBP),p(aptfCran))*(x(axCra)+x(axCraFBP)))/(pow(p(aptfCrakfbp),p(aptfCran))+pow(x(axFBP),p(aptfCran))));
    
    // Crp-cAMP binding kinetics: Hill
    f(aftfCrp) = -p(aptfCrpscale)*(x(axCrpcAMP)-(pow(x(axcAMP),p(aptfCrpn))*(x(axCrp)+x(axCrpcAMP)))/(pow(p(aptfCrpkcamp),p(aptfCrpn))+pow(x(axcAMP),p(aptfCrpn))));
    
    // PdhR-PYR binding kinetics: Hill
    f(aftfPdhR) = -p(aptfPdhRscale)*(x(axPdhRPYR)-(pow(x(axPYR),p(aptfPdhRn))*(x(axPdhR)+x(axPdhRPYR)))/(pow(p(aptfPdhRkpyr),p(aptfPdhRn))+pow(x(axPYR),p(aptfPdhRn))));
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    //// Metabolic reaction rates
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // AceA kinetics: MWC, substrate: ICT, inhibitors: PG3, PEP, AKG
    f(afeAceA) = (p(apeAceAkcat)*x(axAceA)*x(axICT)*pow(x(axICT)/p(apeAceAKict)+1.0,p(apeAceAn)-1.0))/(p(apeAceAKict)*(pow(x(axICT)/p(apeAceAKict)+1.0,p(apeAceAn))+p(apeAceAL)*pow(x(axAKG)/p(apeAceAKakg)+x(axPG3)/p(apeAceAKpg3)+x(axPEP)/p(apeAceAKpep)+1.0,p(apeAceAn))));
    
    // AceB kinetics: Two-substrate MM, substrates: GLX, ACoA
    f(afeAceB) = (p(apeAceBkcat)*x(axACoA)*x(axAceB)*x(axGLX))/(p(apeAceBKacoa)*p(apeAceBKglxacoa)+p(apeAceBKglx)*x(axACoA)+p(apeAceBKacoa)*x(axGLX)+x(axACoA)*x(axGLX));
    
    // Acoa2act kinetics: MWC, substrate: ACoA, activator: PYR
    f(afeAcoa2act) = (p(apeAcoa2actkcat)*x(axACoA)*x(axAcoa2act)*pow(x(axACoA)/p(apeAcoa2actKacoa)+1.0,p(apeAcoa2actn)-1.0))/(p(apeAcoa2actKacoa)*(p(apeAcoa2actL)*pow(x(axPYR)/p(apeAcoa2actKpyr)+1.0,-p(apeAcoa2actn))+pow(x(axACoA)/p(apeAcoa2actKacoa)+1.0,p(apeAcoa2actn))));
    
    // Acs kinetics: MM, substrate: ACT
    f(afeAcs) = (p(apeAcskcat)*x(axACT)*x(axAcs))/(p(apeAcsKact)+x(axACT));
    
    // Akg2mal kinetics: MM, substrate: AKG
    f(afeAkg2mal) = (p(apeAkg2malkcat)*x(axAKG)*x(axAkg2mal))/(p(apeAkg2malKakg)+x(axAKG));
    
    // CAMPdegr kinetics: MM, substrate: cAMP
    f(afeCAMPdegr) = (p(apeCAMPdegrkcat)*x(axCAMPdegr)*x(axcAMP))/(p(apeCAMPdegrKcAMP)+x(axcAMP));
    
    // Cya kinetics: MM, substrate: Cya
    f(afeCya) = (p(apeCyakcat)*x(axCya)*x(axEIIAP))/(p(apeCyaKEIIA)+x(axEIIAP));
    
    // Emp kinetics: reversible MM, substrates: FBP, PG3
    f(afeEmp) = ((p(apeEmpkcatf)*x(axEmp)*x(axFBP))/p(apeEmpKfbp)-(p(apeEmpkcatr)*x(axEmp)*x(axPG3))/p(apeEmpKpg3))/(x(axFBP)/p(apeEmpKfbp)+x(axPG3)/p(apeEmpKpg3)+1.0);
    
    // Eno kinetics: reversible MM, substrates: PG3, PEP
    f(afeEno) = ((p(apeEnokcatf)*x(axEno)*x(axPG3))/p(apeEnoKpg3)-(p(apeEnokcatr)*x(axEno)*x(axPEP))/p(apeEnoKpep))/(x(axPG3)/p(apeEnoKpg3)+x(axPEP)/p(apeEnoKpep)+1.0);
    
    // Fdp kinetics: MWC, substrate: FBP, activator: PEP
    f(afeFdp) = (p(apeFdpkcat)*x(axFBP)*x(axFdp)*pow(x(axFBP)/p(apeFdpKfbp)+1.0,p(apeFdpn)-1.0))/(p(apeFdpKfbp)*(p(apeFdpL)*pow(x(axPEP)/p(apeFdpKpep)+1.0,-p(apeFdpn))+pow(x(axFBP)/p(apeFdpKfbp)+1.0,p(apeFdpn))));
    
    // GltA kinetics: Two-substrate MM, substrates: OAA, ACoA,
    // competitive inhibitor: AKG
    f(afeGltA) = (p(apeGltAkcat)*x(axACoA)*x(axGltA)*x(axOAA))/(p(apeGltAKacoa)*x(axOAA)+x(axACoA)*x(axOAA)+p(apeGltAKacoa)*p(apeGltAKoaaacoa)*(x(axAKG)/p(apeGltAKakg)+1.0)+p(apeGltAKoaa)*x(axACoA)*(x(axAKG)/p(apeGltAKakg)+1.0));
    
    // Icd kinetics: MWC, substrate: ICT, inhibitor: PEP
    f(afeIcd) = (p(apeIcdkcat)*x(axICT)*x(axIcd)*pow(x(axICT)/p(apeIcdKict)+1.0,p(apeIcdn)-1.0))/(p(apeIcdKict)*(p(apeIcdL)*pow(x(axPEP)/p(apeIcdKpep)+1.0,p(apeIcdn))+pow(x(axICT)/p(apeIcdKict)+1.0,p(apeIcdn))));
    
    // MaeAB kinetics: MWC, substrate: MAL, inhibitors: AcoA, cAMP
    f(afeMaeAB) = (p(apeMaeABkcat)*x(axMAL)*x(axMaeAB)*pow(x(axMAL)/p(apeMaeABKmal)+1.0,p(apeMaeABn)-1.0))/(p(apeMaeABKmal)*(p(apeMaeABL)*pow(x(axACoA)/p(apeMaeABKacoa)+x(axcAMP)/p(apeMaeABKcamp)+1.0,p(apeMaeABn))+pow(x(axMAL)/p(apeMaeABKmal)+1.0,p(apeMaeABn))));
    
    // Mdh kinetics: Hill, substrate: MAL
    f(afeMdh) = (p(apeMdhkcat)*pow(x(axMAL),p(apeMdhn))*x(axMdh))/(pow(p(apeMdhKmal),p(apeMdhn))+pow(x(axMAL),p(apeMdhn)));
    
    // PckA kinetics: MM, substrate: OAA, competitive inhibitor: PEP
    f(afePckA) = (p(apePckAkcat)*x(axOAA)*x(axPckA))/(x(axOAA)+p(apePckAKoaa)*(x(axPEP)/p(apePckAKpep)+1.0));
    
    // Pdh kinetics: MWC, substrate: PYR, inhibitors: GLX, PYR
    f(afePdh) = (p(apePdhkcat)*x(axPYR)*x(axPdh)*pow(x(axPYR)/p(apePdhKpyr)+1.0,p(apePdhn)-1.0))/(p(apePdhKpyr)*(p(apePdhL)*pow(x(axGLX)/p(apePdhKglx)+x(axPYR)/p(apePdhKpyrI)+1.0,p(apePdhn))+pow(x(axPYR)/p(apePdhKpyr)+1.0,p(apePdhn))));
    
    // PfkA kinetics: MWC, substrate: G6P, inhibitor: PEP
    f(afePfkA) = (p(apePfkAkcat)*x(axG6P)*x(axPfkA)*pow(x(axG6P)/p(apePfkAKg6p)+1.0,p(apePfkAn)-1.0))/(p(apePfkAKg6p)*(p(apePfkAL)*pow(x(axPEP)/p(apePfkAKpep)+1.0,p(apePfkAn))+pow(x(axG6P)/p(apePfkAKg6p)+1.0,p(apePfkAn))));
    
    // Ppc kinetics: MWC, substrate: PEP, activator: FBP
    f(afePpc) = (p(apePpckcat)*x(axPEP)*x(axPpc)*pow(x(axPEP)/p(apePpcKpep)+1.0,p(apePpcn)-1.0))/(p(apePpcKpep)*(p(apePpcL)*pow(x(axFBP)/p(apePpcKfbp)+1.0,-p(apePpcn))+pow(x(axPEP)/p(apePpcKpep)+1.0,p(apePpcn))));
    
    // PpsA kinetics: MWC, substrate: PYR, inhibitor: PEP
    f(afePpsA) = (p(apePpsAkcat)*x(axPYR)*x(axPpsA)*pow(x(axPYR)/p(apePpsAKpyr)+1.0,p(apePpsAn)-1.0))/(p(apePpsAKpyr)*(p(apePpsAL)*pow(x(axPEP)/p(apePpsAKpep)+1.0,p(apePpsAn))+pow(x(axPYR)/p(apePpsAKpyr)+1.0,p(apePpsAn))));
    
    // PykF kinetics: MWC, substrate: PEP, activator: FBP
    f(afePykF) = (p(apePykFkcat)*x(axPEP)*x(axPykF)*pow(x(axPEP)/p(apePykFKpep)+1.0,p(apePykFn)-1.0))/(p(apePykFKpep)*(p(apePykFL)*pow(x(axFBP)/p(apePykFKfbp)+1.0,-p(apePykFn))+pow(x(axPEP)/p(apePykFKpep)+1.0,p(apePykFn))));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    //// Gene expression rates
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // aceBAK expression: sum of the following three kinetics
    // MM plus basal expression, substrate: Cra
    // MM plus basal expression, substrate: Crpcamp
    // MWC-like, substrate: IclR, activator: GLX, inhibitor: PYR
    // aceB and aceK expression are coupled to aceA expression
    // with constant factors
    f(afgaceA) = p(apbmkexpr)*mu*((1-x(axCra)/(x(axCra)+p(apgaceBAKKcra)))*p(apgaceBAKvcraunbound)
    +x(axCra)/(x(axCra)+p(apgaceBAKKcra))*p(apgaceBAKvcrabound)
    +(1-x(axCrpcAMP)/(x(axCrpcAMP)+p(apgaceBAKKcrp)))*p(apgaceBAKvcrpunbound)
    +x(axCrpcAMP)/(x(axCrpcAMP)+p(apgaceBAKKcrp))*p(apgaceBAKvcrpbound)
    +p(apgaceBAKkcaticlr)*x(axIclR)*(1-(p(apgaceBAKDNA)/p(apgaceBAKKDNA))
    *(1+x(axPYR)/p(apgaceBAKKPprime))/(1+(x(axGLX)/p(apgaceBAKKG)
    *(1+x(axGLX)/p(apgaceBAKKG)))/p(apgaceBAKL)
    +p(apgaceBAKDNA)/p(apgaceBAKKDNA)+x(axPYR)/p(apgaceBAKKP)
    +p(apgaceBAKDNA)*x(axPYR)/p(apgaceBAKKDNA)/p(apgaceBAKKPprime))));
    f(afgaceB) = f(afgaceA)*p(apgaceBAKaceBfactor);
    f(afgaceK) = p(apgaceBAKaceKfactor)*f(afgaceA);
    
    // acoa2act kinetics: constitutive expression
    f(afgacoa2act) = 0;
    
    // acs kinetics: Hill plus basal expression, substrate: Crpcamp
    f(afgacs) = -mu*p(apbmkexpr)*(p(apgacsvcrpunbound)*(pow(x(axCrpcAMP),p(apgacsn))/(pow(p(apgacsKcrp),p(apgacsn))+pow(x(axCrpcAMP),p(apgacsn)))-1.0)-(p(apgacsvcrpbound)*pow(x(axCrpcAMP),p(apgacsn)))/(pow(p(apgacsKcrp),p(apgacsn))+pow(x(axCrpcAMP),p(apgacsn))));
    
    // akg2mal kinetics: Hill plus basal expression, substrate: Crpcamp
    f(afgakg2mal) = -mu*p(apbmkexpr)*(p(apgakg2malvcrpunbound)*(pow(x(axCrpcAMP),p(apgakg2maln))/(pow(p(apgakg2malKcrp),p(apgakg2maln))+pow(x(axCrpcAMP),p(apgakg2maln)))-1.0)-(p(apgakg2malvcrpbound)*pow(x(axCrpcAMP),p(apgakg2maln)))/(pow(p(apgakg2malKcrp),p(apgakg2maln))+pow(x(axCrpcAMP),p(apgakg2maln))));
    
    // campdegr kinetics: constitutive expression
    f(afgcampdegr) = 0.0;
    
    // cra kinetics: constitutive expression
    f(afgcra) = 0.0;
    
    // crp kinetics: constitutive expression
    f(afgcrp) = 0.0;
    
    // cya kinetics: constitutive expression
    f(afgcya) = 0.0;
    
    // emp kinetics: MM plus basal expression, substrate: Cra
    f(afgemp) = p(apbmkexpr)*mu*((1-x(axCra)/(x(axCra)+p(apgempKcra)))*p(apgempvcraunbound)
    +x(axCra)/(x(axCra)+p(apgempKcra))*p(apgempvcrabound)
    +(1-x(axCrpcAMP)/(x(axCrpcAMP)+p(apgempKcrp)))*p(apgempvcrpunbound)
    +x(axCrpcAMP)/(x(axCrpcAMP)+p(apgempKcrp))*p(apgempvcrpbound));
    
    // eno kinetics: MM plus basal expression, substrate: Cra
    f(afgeno) = -mu*p(apbmkexpr)*(p(apgenovcraunbound)*(x(axCra)/(p(apgenoKcra)+x(axCra))-1.0)-(p(apgenovcrabound)*x(axCra))/(p(apgenoKcra)+x(axCra)));
    
    // fdp kinetics: MM plus basal expression, substrate: Cra
    f(afgfdp) = -mu*p(apbmkexpr)*(p(apgfdpvcraunbound)*(x(axCra)/(p(apgfdpKcra)+x(axCra))-1.0)-(p(apgfdpvcrabound)*x(axCra))/(p(apgfdpKcra)+x(axCra)));
    
    // gltA kinetics: Hill plus basal expression, substrate: Crpcamp
    f(afggltA) = -mu*p(apbmkexpr)*(p(apggltAvcrpunbound)*(pow(x(axCrpcAMP),p(apggltAn))/(pow(p(apggltAKcrp),p(apggltAn))+pow(x(axCrpcAMP),p(apggltAn)))-1.0)-(p(apggltAvcrpbound)*pow(x(axCrpcAMP),p(apggltAn)))/(pow(p(apggltAKcrp),p(apggltAn))+pow(x(axCrpcAMP),p(apggltAn))));
    
    // icd kinetics: MM plus basal expression, substrate: Cra
    f(afgicd) = -mu*p(apbmkexpr)*(p(apgicdvcraunbound)*(x(axCra)/(p(apgicdKcra)+x(axCra))-1.0)-(p(apgicdvcrabound)*x(axCra))/(p(apgicdKcra)+x(axCra)));
    
    // iclr expression: constitutive
    f(afgiclr) = 0.0;
    
    // mdh kinetics: MM plus basal expression, substrate: Crpcamp
    f(afgmdh) = -mu*p(apbmkexpr)*(p(apgmdhvcrpunbound)*(x(axCrpcAMP)/(p(apgmdhKcrp)+x(axCrpcAMP))-1.0)-(p(apgmdhvcrpbound)*x(axCrpcAMP))/(p(apgmdhKcrp)+x(axCrpcAMP)));
    
    // me kinetics: growth rate- dependent constitutive expression
    f(afgmaeAB) = SSMaeAB*(mu+p(apdkdegr));
    
    // pckA kinetics: MM plus basal expression, substrate: Cra
    f(afgpckA) = -mu*p(apbmkexpr)*(p(apgpckAvcraunbound)*(x(axCra)/(p(apgpckAKcra)+x(axCra))-1.0)-(p(apgpckAvcrabound)*x(axCra))/(p(apgpckAKcra)+x(axCra)));
    
    // pdh kinetics: MM plus basal expression, substrate: PdhR
    f(afgpdh) = -mu*p(apbmkexpr)*(p(apgpdhvpdhrunbound)*(x(axPdhR)/(p(apgpdhKpdhr)+x(axPdhR))-1.0)-(p(apgpdhvpdhrbound)*x(axPdhR))/(p(apgpdhKpdhr)+x(axPdhR)));
    
    // pdhr expression: constitutive
    f(afgpdhr) = 0.0;
    
    // pfkA kinetics: MM plus basal expression, substrate: Cra
    f(afgpfkA) = -mu*p(apbmkexpr)*(p(apgpfkAvcraunbound)*(x(axCra)/(p(apgpfkAKcra)+x(axCra))-1.0)-(p(apgpfkAvcrabound)*x(axCra))/(p(apgpfkAKcra)+x(axCra)));
    
    // ppc kinetics: growth rate- dependent constitutive expression
    f(afgppc) = SSPpc*(mu+p(apdkdegr));
    
    // ppsA kinetics: MM plus basal expression, substrate: Cra
    f(afgppsA) = -mu*p(apbmkexpr)*(p(apgppsAvcraunbound)*(x(axCra)/(p(apgppsAKcra)+x(axCra))-1.0)-(p(apgppsAvcrabound)*x(axCra))/(p(apgppsAKcra)+x(axCra)));
    
    // pykF kinetics: MM plus basal expression, substrate: Cra
    f(afgpykF) = -mu*p(apbmkexpr)*(p(apgpykFvcraunbound)*(x(axCra)/(p(apgpykFKcra)+x(axCra))-1.0)-(p(apgpykFvcrabound)*x(axCra))/(p(apgpykFKcra)+x(axCra)));
    
    // EIIA kinetics: constitutive expression
    f(afgEIIA) = 0.0;
    
    // EIICB kinetics: constitutive expression
    f(afgEIICB) = 0.0;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    //// Protein degradation and dilution rates
    // Constitutively produced proteins are neither produced,
    // nor degraded, nor diluted, to keep their levels constant
    // All other proteins degrade with a constant rate,
    // and dilute with the growth rate
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    f(afdAceA) = (mu+p(apdkdegr))*x(axAceA);
    f(afdAceB) = (mu+p(apdkdegr))*x(axAceB);
    f(afdAceK) = (mu+p(apdkdegr))*x(axAceK);
    f(afdAcoa2act) = 0;
    f(afdAcs) = (mu+p(apdkdegr))*x(axAcs);
    f(afdAkg2mal) = (mu+p(apdkdegr))*x(axAkg2mal);
    f(afdCAMPdegr) = 0;
    f(afdCra) = 0;
    f(afdCraFBP) = 0;
    f(afdCrp) = 0;
    f(afdCrpcAMP) = 0;
    f(afdCya)  = 0;
    f(afdEmp)  = (mu+p(apdkdegr))*x(axEmp);
    f(afdEno)  = (mu+p(apdkdegr))*x(axEno);
    f(afdFdp)  = (mu+p(apdkdegr))*x(axFdp);
    f(afdGltA) = (mu+p(apdkdegr))*x(axGltA);
    f(afdIcd) = (mu+p(apdkdegr))*x(axIcd);
    f(afdIcdP) = (mu+p(apdkdegr))*x(axIcdP);
    f(afdIclR) = 0;
    f(afdMdh)  = (mu+p(apdkdegr))*x(axMdh);
    f(afdMaeAB) = (mu+p(apdkdegr))*x(axMaeAB);
    f(afdPckA) = (mu+p(apdkdegr))*x(axPckA);
    f(afdPdh)  = (mu+p(apdkdegr))*x(axPdh);
    f(afdPdhR) = 0;
    f(afdPdhRPYR) = 0;
    f(afdPfkA) = (mu+p(apdkdegr))*x(axPfkA);
    f(afdPpc)  = (mu+p(apdkdegr))*x(axPpc);
    f(afdPpsA) = (mu+p(apdkdegr))*x(axPpsA);
    f(afdPykF) = (mu+p(apdkdegr))*x(axPykF);
    f(afdEIIA) = 0;
    f(afdEIIAP) = 0;
    f(afdEIICB) = 0;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //// Metabolite dilution rates
    // Intracellular metabolites do not degrade,
    // only dilute with the growth rate
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    f(afdACoA) = mu*x(axACoA);
    f(afdAKG)  = mu*x(axAKG);
    f(afdcAMP) = mu*x(axcAMP);
    f(afdFBP)  = mu*x(axFBP);
    f(afdG6P)  = mu*x(axG6P);
    f(afdGLX)  = mu*x(axGLX);
    f(afdICT)  = mu*x(axICT);
    f(afdMAL)  = mu*x(axMAL);
    f(afdOAA)  = mu*x(axOAA);
    f(afdPEP)  = mu*x(axPEP);
    f(afdPG3)  = mu*x(axPG3);
    f(afdPYR)  = mu*x(axPYR);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Let the cell population grow with the actual growth rate
    f(afenvgrowth) = x(axOD)*mu;
    
    // Scale glucose uptake from the environment with the
    // actual population size
    f(afenvGLCup) = p(apenvuc)*p(apenvMGLC)*x(axOD)*f(afptsr4);
    
    // Scale acetate uptake from the environment with the
    // actual population size
    f(afenvACTup) = p(apenvuc)*p(apenvMACT)*x(axOD)*f(afeAcs);
    
    // Scale acetate excretion to the environment with the
    // actual population size for all scenarios except 1 and 4,
    // for which glucose should remain the sole carbon source -
    // in these cases, the excreted acetate is directed to nowhere
    
    if(amigo_model->controls_v[0][0]==1 || amigo_model->controls_v[0][0]==4){
        f(afenvACTex) = 0;
    }
    else{
        f(afenvACTex) = p(apenvuc)*p(apenvMACT)*x(axOD)*f(afeAcoa2act);
    }
    //S*f'
    dgemm_(&transa, &transb, &l, &n, &m, &alpha, S, &lda, f, &ldb, &beta, xdot, &ldc);

    return(0);
}



int amigoJAC_B3(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
    return(0);
    
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_B3(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}

void amigoRHS_get_OBS_B3(void* data){

}

void amigoRHS_get_sens_OBS_B3(void* data){

}
    
#define y(i) ( NV_DATA_S(y)[i-1] )
void amigo_Y_at_tcon_B3(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
    if(t==29000.0001){
        x(1)=0.03;
        x(2)=5;
        x(3)=0;
    }
    
    if(t==99000.0001){
        x(1)=5e-4;
        x(2)=3;
        x(3)=3;
    }
    
}
