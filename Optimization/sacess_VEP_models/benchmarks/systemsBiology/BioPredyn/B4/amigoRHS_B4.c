/*
Check README.txt
*/

#include <AMIGO_problem.h>
 
	/* *** Definition of the states *** */

#define	x1  Ith(y,0)
#define	x2  Ith(y,1)
#define	x3  Ith(y,2)
#define	x4  Ith(y,3)
#define	x5  Ith(y,4)
#define	x6  Ith(y,5)
#define	x7  Ith(y,6)
#define	x8  Ith(y,7)
#define	x9  Ith(y,8)
#define	x10 Ith(y,9)
#define	x11 Ith(y,10)
#define	x12 Ith(y,11)
#define	x13 Ith(y,12)
#define	x14 Ith(y,13)
#define	x15 Ith(y,14)
#define	x16 Ith(y,15)
#define	x17 Ith(y,16)
#define	x18 Ith(y,17)
#define	x19 Ith(y,18)
#define	x20 Ith(y,19)
#define	x21 Ith(y,20)
#define	x22 Ith(y,21)
#define	x23 Ith(y,22)
#define	x24 Ith(y,23)
#define	x25 Ith(y,24)
#define	x26 Ith(y,25)
#define	x27 Ith(y,26)
#define	x28 Ith(y,27)
#define	x29 Ith(y,28)
#define	x30 Ith(y,29)
#define	x31 Ith(y,30)
#define	x32 Ith(y,31)
#define	x33 Ith(y,32)
#define	x34 Ith(y,33)
#define	x35 Ith(y,34)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dx1  Ith(ydot,0)
#define	dx2  Ith(ydot,1)
#define	dx3  Ith(ydot,2)
#define	dx4  Ith(ydot,3)
#define	dx5  Ith(ydot,4)
#define	dx6  Ith(ydot,5)
#define	dx7  Ith(ydot,6)
#define	dx8  Ith(ydot,7)
#define	dx9  Ith(ydot,8)
#define	dx10 Ith(ydot,9)
#define	dx11 Ith(ydot,10)
#define	dx12 Ith(ydot,11)
#define	dx13 Ith(ydot,12)
#define	dx14 Ith(ydot,13)
#define	dx15 Ith(ydot,14)
#define	dx16 Ith(ydot,15)
#define	dx17 Ith(ydot,16)
#define	dx18 Ith(ydot,17)
#define	dx19 Ith(ydot,18)
#define	dx20 Ith(ydot,19)
#define	dx21 Ith(ydot,20)
#define	dx22 Ith(ydot,21)
#define	dx23 Ith(ydot,22)
#define	dx24 Ith(ydot,23)
#define	dx25 Ith(ydot,24)
#define	dx26 Ith(ydot,25)
#define	dx27 Ith(ydot,26)
#define	dx28 Ith(ydot,27)
#define	dx29 Ith(ydot,28)
#define	dx30 Ith(ydot,29)
#define	dx31 Ith(ydot,30)
#define	dx32 Ith(ydot,31)
#define	dx33 Ith(ydot,32)
#define	dx34 Ith(ydot,33)
#define	dx35 Ith(ydot,34)

	/* *** Definition of the parameters *** */

#define	par1   (*amigo_model).pars[0]
#define	par2   (*amigo_model).pars[1]
#define	par3   (*amigo_model).pars[2]
#define	par4   (*amigo_model).pars[3]
#define	par5   (*amigo_model).pars[4]
#define	par6   (*amigo_model).pars[5]
#define	par7   (*amigo_model).pars[6]
#define	par8   (*amigo_model).pars[7]
#define	par9   (*amigo_model).pars[8]
#define	par10  (*amigo_model).pars[9]
#define	par11  (*amigo_model).pars[10]
#define	par12  (*amigo_model).pars[11]
#define	par13  (*amigo_model).pars[12]
#define	par14  (*amigo_model).pars[13]
#define	par15  (*amigo_model).pars[14]
#define	par16  (*amigo_model).pars[15]
#define	par17  (*amigo_model).pars[16]
#define	par18  (*amigo_model).pars[17]
#define	par19  (*amigo_model).pars[18]
#define	par20  (*amigo_model).pars[19]
#define	par21  (*amigo_model).pars[20]
#define	par22  (*amigo_model).pars[21]
#define	par23  (*amigo_model).pars[22]
#define	par24  (*amigo_model).pars[23]
#define	par25  (*amigo_model).pars[24]
#define	par26  (*amigo_model).pars[25]
#define	par27  (*amigo_model).pars[26]
#define	par28  (*amigo_model).pars[27]
#define	par29  (*amigo_model).pars[28]
#define	par30  (*amigo_model).pars[29]
#define	par31  (*amigo_model).pars[30]
#define	par32  (*amigo_model).pars[31]
#define	par33  (*amigo_model).pars[32]
#define	par34  (*amigo_model).pars[33]
#define	par35  (*amigo_model).pars[34]
#define	par36  (*amigo_model).pars[35]
#define	par37  (*amigo_model).pars[36]
#define	par38  (*amigo_model).pars[37]
#define	par39  (*amigo_model).pars[38]
#define	par40  (*amigo_model).pars[39]
#define	par41  (*amigo_model).pars[40]
#define	par42  (*amigo_model).pars[41]
#define	par43  (*amigo_model).pars[42]
#define	par44  (*amigo_model).pars[43]
#define	par45  (*amigo_model).pars[44]
#define	par46  (*amigo_model).pars[45]
#define	par47  (*amigo_model).pars[46]
#define	par48  (*amigo_model).pars[47]
#define	par49  (*amigo_model).pars[48]
#define	par50  (*amigo_model).pars[49]
#define	par51  (*amigo_model).pars[50]
#define	par52  (*amigo_model).pars[51]
#define	par53  (*amigo_model).pars[52]
#define	par54  (*amigo_model).pars[53]
#define	par55  (*amigo_model).pars[54]
#define	par56  (*amigo_model).pars[55]
#define	par57  (*amigo_model).pars[56]
#define	par58  (*amigo_model).pars[57]
#define	par59  (*amigo_model).pars[58]
#define	par60  (*amigo_model).pars[59]
#define	par61  (*amigo_model).pars[60]
#define	par62  (*amigo_model).pars[61]
#define	par63  (*amigo_model).pars[62]
#define	par64  (*amigo_model).pars[63]
#define	par65  (*amigo_model).pars[64]
#define	par66  (*amigo_model).pars[65]
#define	par67  (*amigo_model).pars[66]
#define	par68  (*amigo_model).pars[67]
#define	par69  (*amigo_model).pars[68]
#define	par70  (*amigo_model).pars[69]
#define	par71  (*amigo_model).pars[70]
#define	par72  (*amigo_model).pars[71]
#define	par73  (*amigo_model).pars[72]
#define	par74  (*amigo_model).pars[73]
#define	par75  (*amigo_model).pars[74]
#define	par76  (*amigo_model).pars[75]
#define	par77  (*amigo_model).pars[76]
#define	par78  (*amigo_model).pars[77]
#define	par79  (*amigo_model).pars[78]
#define	par80  (*amigo_model).pars[79]
#define	par81  (*amigo_model).pars[80]
#define	par82  (*amigo_model).pars[81]
#define	par83  (*amigo_model).pars[82]
#define	par84  (*amigo_model).pars[83]
#define	par85  (*amigo_model).pars[84]
#define	par86  (*amigo_model).pars[85]
#define	par87  (*amigo_model).pars[86]
#define	par88  (*amigo_model).pars[87]
#define	par89  (*amigo_model).pars[88]
#define	par90  (*amigo_model).pars[89]
#define	par91  (*amigo_model).pars[90]
#define	par92  (*amigo_model).pars[91]
#define	par93  (*amigo_model).pars[92]
#define	par94  (*amigo_model).pars[93]
#define	par95  (*amigo_model).pars[94]
#define	par96  (*amigo_model).pars[95]
#define	par97  (*amigo_model).pars[96]
#define	par98  (*amigo_model).pars[97]
#define	par99  (*amigo_model).pars[98]
#define	par100 (*amigo_model).pars[99]
#define	par101 (*amigo_model).pars[100]
#define	par102 (*amigo_model).pars[101]
#define	par103 (*amigo_model).pars[102]
#define	par104 (*amigo_model).pars[103]
#define	par105 (*amigo_model).pars[104]
#define	par106 (*amigo_model).pars[105]
#define	par107 (*amigo_model).pars[106]
#define	par108 (*amigo_model).pars[107]
#define	par109 (*amigo_model).pars[108]
#define	par110 (*amigo_model).pars[109]
#define	par111 (*amigo_model).pars[110]
#define	par112 (*amigo_model).pars[111]
#define	par113 (*amigo_model).pars[112]
#define	par114 (*amigo_model).pars[113]
#define	par115 (*amigo_model).pars[114]
#define	par116 (*amigo_model).pars[115]
#define	par117 (*amigo_model).pars[116]
#define s1	((*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim])

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_B4(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;
	/* *** Definition of the algebraic variables *** */

	double	rho;
	double	vol_1;
	double	vol_2;
	double	rss6;
	double	rss7;
	double	rss8;
	double	rss9;
	double	rss10;
	double	rss11;
	double	rss12;
	double	rss13;
	double	rss14;
	double	rss15;
	double	rss16;
	double	rss17;
	double	rss18;
	double	rss19;
	double	rss20;
	double	rss21;
	double	rss22;
	double	rss23;
	double	rss24;
	double	rss25;
	double	rss26;
	double	rss27;
	double	rss28;
	double	rss29;
	double	rss30;
	double	rss31;
	double	rss32;
	double	css1;
	double	css2;
	double	css3;
	double	css4;
	double	css5;
	double	css6;
	double	css7;
	double	css8;
	double	css9;
	double	css10;
	double	css11;
	double	css12;
	double	css13;
	double	css14;
	double	css15;
	double	css16;
	double	css17;
	double	css18;
	double	css19;
	double	css20;
	double	css21;
	double	css22;
	double	css23;
	double	css24;
	double	css25;
	double	css26;
	double	css27;
	double	css28;
	double	css29;
	double	css30;
	double	css31;
	double	css32;
	double	css33;
	double	css34;
	double	e1;
	double	e2;
	double	e3;
	double	e4;
	double	e5;
	double	e6;
	double	e7;
	double	e8;
	double	e14;
	double	e15;
	double	e16;
	double	e17;
	double	e18;
	double	e19;
	double	e20;
	double	e21;
	double	e22;
	double	e23;
	double	e24;
	double	e25;
	double	e26;
	double	e27;
	double	e28;
	double	e29;
	double	e30;
	double	e31;
	double	e32;
	double	e33;
	double	e34;
	double	e35;
	double	e36;
	double	e37;
	double	e38;
	double	e39;
	double	e40;
	double	e41;
	double	e42;
	double	e43;
	double	e44;
	double	e45;
	double	e46;
	double	e47;
	double	e48;
	double	e49;
	double	e50;
	double	e51;
	double	e52;
	double	e53;
	double	e54;
	double	e55;
	double	e56;
	double	e57;
	double	e58;
	double	e59;
	double	e60;
	double	e61;
	double	e62;
	double	e63;
	double	e64;
	double	e65;
	double	e66;
	double	e67;
	double	e68;
	double	e69;
	double	e70;
	double	e71;
	double	e72;
	double	e73;
	double	e74;
	double	e75;
	double	e76;
	double	e77;
	double	e78;
	double	e79;
	double	e80;
	double	e81;
	double	e82;
	double	e83;
	double	e84;
	double	e85;
	double	e86;
	double	e87;
	double	e88;
	double	e89;
	double	e90;
	double	e91;
	double	e92;
	double	e93;
	double	e94;
	double	e95;
	double	e96;
	double	e97;
	double	e98;
	double	e99;
	double	e100;
	double	e101;
	double	e102;
	double	e103;
	double	e104;
	double	e105;
	double	e106;
	double	e107;
	double	e108;
	double	e109;
	double	e110;
	double	e111;
	double	e112;
	double	e113;
	double	e114;
	double	e115;
	double	e116;
	double	e117;
	double	e118;
	double	e119;
	double	e120;
	double	e121;
	double	e122;
	double	e123;
	double	e124;
	double	e125;
	double	e126;
	double	e127;
	double	e128;
	double	clogk1;
	double	clogk2;
	double	clogk3;
	double	clogk4;
	double	clogk5;
	double	clogk6;
	double	clogk7;
	double	clogk8;
	double	clogk9;
	double	clogk10;
	double	clogk11;
	double	clogk12;
	double	clogk13;
	double	clogk14;
	double	clogk15;
	double	clogk16;
	double	clogk17;
	double	clogk18;
	double	clogk19;
	double	clogk20;
	double	clogk21;
	double	clogk22;
	double	clogk23;
	double	clogk24;
	double	clogk25;
	double	clogk26;
	double	clogk27;
	double	clogk28;
	double	clogk29;
	double	clogk30;
	double	clogk31;
	double	clogk32;
	double	clogk33;
	double	clogk34;
	double	r1;
	double	r2;
	double	r3;
	double	r4;
	double	r5;
	double	r6;
	double	r7;
	double	r8;
	double	r9;
	double	r10;
	double	r11;
	double	r12;
	double	r13;
	double	r14;
	double	r15;
	double	r16;
	double	r17;
	double	r18;
	double	r19;
	double	r20;
	double	r21;
	double	r22;
	double	r23;
	double	r24;
	double	r25;
	double	r26;
	double	r27;
	double	r28;
	double	r29;
	double	r30;
	double	r31;
	double	r32;



	/* *** Equations *** */

	rho=141.4710605;
	vol_1=0.9;
	vol_2=0.1;
	rss6=2413.6562849217303*vol_1*1.0/vol_2;
	rss7=603.414071180172;
	rss8=40133.38664791273*vol_1*1.0/vol_2;
	rss9=74836.04665510054*vol_1*1.0/vol_2;
	rss10=180754.26889005644*vol_1*1.0/vol_2;
	rss11=164390.2745505474;
	rss12=36512.902220891134;
	rss13=166803.93083553354;
	rss14=39529.972576829896;
	rss15=1810.242213701644;
	rss16=127877.37232994902;
	rss17=166200.51676427448;
	rss18=216060.34296802033*vol_1*1.0/vol_2;
	rss19=5.028450593613567;
	rss20=35306.0740779408*vol_1*1.0/vol_2;
	rss21=0.0*vol_1*1.0/vol_2;
	rss22=532531.1670427525;
	rss23=38926.55850561767;
	rss24=40133.38664808034;
	rss25=689982.1797742102;
	rss26=495414.85075147304*vol_1*1.0/vol_2;
	rss27=2413.6562849329107;
	rss28=83401.9654177705;
	rss29=603.4140712332287;
	rss30=603.4140712332287;
	rss31=2413.656284931604;
	rss32=533134.5811139438;
	css1=100000;
	css2=1;
	css3=1000;
	css4=1000;
	css5=30;
	css6=1000;
	css7=1000;
	css8=1000;
	css9=1000;
	css10=1000;
	css11=1000;
	css12=1000;
	css13=3000;
	css14=1000;
	css15=1000;
	css16=1000;
	css17=100;
	css18=100;
	css19=100;
	css20=100;
	css21=1000;
	css22=1000;
	css23=1000;
	css24=1000;
	css25=1000;
	css26=1000;
	css27=1000;
	css28=1000;
	css29=1000;
	css30=3000;
	css31=100000;
	css32=1000;
	css33=1000;
	css34=1000;
	e1=par1;
	e2=par2;
	e3=par3;
	e4=par4;
	e5=par5;
	e6=par6;
	e7=par7;
	e8=par8;
	e14=par9;
	e15=par10;
	e16=-par11;
	e17=-par12;
	e18=par13;
	e19=par14;
	e20=-par15;
	e21=-par16;
	e22=par17;
	e23=par18;
	e24=-par19;
	e25=-par20;
	e26=par21;
	e27=par22;
	e28=par23;
	e29=-par24;
	e30=-par25;
	e31=par26;
	e32=par27;
	e33=par28;
	e34=-par29;
	e35=-par30;
	e36=-par31;
	e37=par32;
	e38=par33;
	e39=-par34;
	e40=-par35;
	e41=0.0;
	e42=par36;
	e43=par37;
	e44=par38;
	e45=-par39;
	e46=-par40;
	e47=par41;
	e48=0.0;
	e49=-par42;
	e50=par43;
	e51=par44;
	e52=-par45;
	e53=-par46;
	e54=par47;
	e55=-par48;
	e56=-par49;
	e57=-par50;
	e58=par51;
	e59=par52;
	e60=par53;
	e61=-par54;
	e62=-par55;
	e63=-par56;
	e64=par57;
	e65=-par58;
	e66=par59;
	e67=par60;
	e68=par61;
	e69=-par62;
	e70=-par63;
	e71=-par64;
	e72=par65;
	e73=par66;
	e74=-par67;
	e75=-par68;
	e76=par69;
	e77=-par70;
	e78=par71;
	e79=par72;
	e80=-par73;
	e81=-par74;
	e82=par75;
	e83=par76;
	e84=par77;
	e85=par78;
	e86=par79;
	e87=-par80;
	e88=-par81;
	e89=-par82;
	e90=-par83;
	e91=0.0;
	e92=0.0;
	e93=-par84;
	e94=0.7;
	e95=-2.0e-1;
	e96=par85;
	e97=par86;
	e98=-par87;
	e99=-par88;
	e100=par89;
	e101=par90;
	e102=-par91;
	e103=-par92;
	e104=par93;
	e105=par94;
	e106=-par95;
	e107=-par96;
	e108=par97;
	e109=-par98;
	e110=par99;
	e111=par100;
	e112=par101;
	e113=-par102;
	e114=-par103;
	e115=par104;
	e116=par105;
	e117=-par106;
	e118=par107;
	e119=-par108;
	e120=par109;
	e121=-par110;
	e122=par111;
	e123=-par112;
	e124=par113;
	e125=-par114;
	e126=par115;
	e127=-par116;
	e128=-par117;
	clogk1=log(x1);
	clogk2=log(x2);
	clogk3=log(x3);
	clogk4=log(x4);
	clogk5=log(x5);
	clogk6=log(x6);
	clogk7=log(x7);
	clogk8=log(x8);
	clogk9=log(x9);
	clogk10=log(x10);
	clogk11=log(x11);
	clogk12=log(x12);
	clogk13=log(x13);
	clogk14=log(x14);
	clogk15=log(x15);
	clogk16=log(x16);
	clogk17=log(x17);
	clogk18=log(x18);
	clogk19=log(x19);
	clogk20=log(x20);
	clogk21=log(x21);
	clogk22=log(x22);
	clogk23=log(x23);
	clogk24=log(x24);
	clogk25=log(x25);
	clogk26=log(x26);
	clogk27=log(x27);
	clogk28=log(x28);
	clogk29=log(x29);
	clogk30=log(x30);
	clogk31=log(x31);
	clogk32=log(x32);
	clogk33=log(x33);
	clogk34=log(x34);
	r1=0.0;
	r2=0.0;
	r3=0.0;
	r4=0.0;
	r5=0.0;
	r6=((rss6*(1.0+(((e14*clogk8)+(e15*clogk9))+((e16*clogk6)+(e17*clogk7))))));
	r7=((rss7*(1.0+(((e18*clogk12)+(e19*clogk13))+((e20*clogk10)+(e21*clogk11))))));
	r8=((rss8*(1.0+((((e22*clogk15)+(e23*clogk6))+((e24*clogk8)+(e25*clogk14)))+(e26*clogk16)))));
	r9=((rss9*(1.0+(((e27*clogk16)+(e28*clogk7))+((e29*clogk15)+(e30*clogk9))))));
	r10=((rss10*(1.0+((((e31*clogk9)+(e32*clogk19))+(e33*clogk20))+(((e34*clogk7)+(e35*clogk17))+(e36*clogk18))))));
	r11=((rss11*(1.0+((((e37*clogk22)+(e38*clogk11))+((e39*clogk21)+(e40*clogk13)))+(e41*clogk13)))));
	r12=((rss12*(1.0+(((((((e42*clogk15)+(e43*clogk7))+(e44*clogk21))+((e45*clogk8)+(e46*clogk9)))+(e47*clogk32))+(e48*clogk9))+(e49*clogk30)))));
	r13=((rss13*(1.0+(((((((e50*clogk25)+(e51*clogk26))+((e52*clogk23)+(e53*clogk24)))+(e54*clogk11))+(e55*clogk11))+(e56*clogk13))+(e57*clogk22)))));
	r14=((rss14*(1.0+((((((e58*clogk28)+(e59*clogk29))+(e60*clogk23))+(((e61*clogk12)+(e62*clogk27))+(e63*clogk25)))+(e64*clogk27))+(e65*clogk10)))));
	r15=((rss15*(1.0+((((e66*clogk32)+(e67*clogk22))+(e68*clogk16))+(((e69*clogk15)+(e70*clogk30))+(e71*clogk27))))));
	r16=((rss16*(1.0+(((e72*clogk21)+(e73*clogk23))+((e74*clogk2)+(e75*clogk25))))));
	r17=((rss17*(1.0+((e76*clogk24)+(e77*clogk22)))));
	r18=((rss18*(1.0+(((e78*clogk20)+(e79*clogk17))+((e80*clogk18)+(e81*clogk19))))));
	r19=((rss19/(((((((((css24/e1)*(css25/e2))*(css12/e3))*(css33/e4))*(css34/e5))*(css29/e6))*(css10/e7))*(css13/e8))/((((((((1.0+(css24/e1))*(1.0+(css25/e2)))*(1.0+(css12/e3)))*(1.0+(css33/e4)))*(1.0+(css34/e5)))*(1.0+(css29/e6)))*(1.0+(css10/e7)))*(1.0+(css13/e8)))))*((((((((((x24*css24)/e1)*((x25*css25)/e2))*((x12*css12)/e3))*((x33*css33)/e4))*((x34*css34)/e5))*((x29*css29)/e6))*((x10*css10)/e7))*((x13*css13)/e8))/((((((((1.0+((x24*css24)/e1))*(1.0+((x25*css25)/e2)))*(1.0+((x12*css12)/e3)))*(1.0+((x33*css33)/e4)))*(1.0+((x34*css34)/e5)))*(1.0+((x29*css29)/e6)))*(1.0+((x10*css10)/e7)))*(1.0+((x13*css13)/e8)))));
	r20=((rss20*(1.0+(((((((((e82*clogk8)+(e83*clogk7))+(e84*clogk19))+(e85*clogk31))+(e86*clogk32))+((((e87*clogk16)+(e88*clogk9))+(e89*clogk17))+(e90*clogk30)))+(e91*clogk31))+(e92*clogk9))+(e93*clogk15)))));
	r21=((rss21*(1.0+((e94*clogk18)+(e95*clogk20)))));
	r22=((rss22*(1.0+(((e96*clogk11)+(e97*clogk30))+((e98*clogk32)+(e99*clogk13))))));
	r23=((rss23*(1.0+(((e100*clogk8)+(e101*clogk27))+((e102*clogk28)+(e103*clogk16))))));
	r24=((rss24*(1.0+(((e104*clogk14)+(e105*clogk12))+((e106*clogk29)+(e107*clogk6))))));
	r25=((rss25*(1.0+((e108*clogk13)+(e109*clogk11)))));
	r26=((rss26*(1.0+((((e110*clogk32)+(e111*clogk31))+(e112*clogk18))+((e113*clogk30)+(e114*clogk20))))));
	r27=((rss27*(1.0+(((e115*clogk27)+(e116*clogk31))+(e117*clogk16)))));
	r28=((rss28*(1.0+((e118*clogk1)+(e119*clogk26)))));
	r29=((rss29*(1.0+((e120*clogk3)+(e121*clogk33)))));
	r30=((rss30*(1.0+((e122*clogk4)+(e123*clogk34)))));
	r31=((rss31*(1.0+((e124*clogk6)+(e125*clogk12)))));
	r32=((rss32*(1.0+((e126*clogk18)+((e127*clogk31)+(e128*clogk20))))));
	dx35=0;
	dx1=(r1-r28*1.0/rho-x1*css1/x35*dx35)/css1;
	dx2=(-r2+r16*1.0/rho-x2*css2/x35*dx35)/css2;
	dx3=(r3-r29*1.0/rho-x3*css3/x35*dx35)/css3;
	dx4=(r4-r30*1.0/rho-x4*css4/x35*dx35)/css4;
	dx5=(-r5+r19*1.0/rho-x5*css5/x35*dx35)/css5;
	dx6=(r6-r8+r24*vol_1/vol_2-r31*vol_1/vol_2)/css6;
	dx7=(r6-r9+r10-2.0*r12*vol_1/vol_2-r20)/css7;
	dx8=(-r6+r8+r12*vol_1/vol_2-r20-r23*vol_1/vol_2)/css8;
	dx9=(-r6+r9-r10+2.0*r12*vol_1/vol_2+r20)/css9;
	dx10=(r7-120.0*r19)/css10;
	dx11=(r7-r11+1260.0*r19-r22+r25)/css11;
	dx12=(-r7+r14-240.0*r19-r24+r31)/css12;
	dx13=(-r7+r11-1260.0*r19+r22-r25)/css13;
	dx14=(r8-r24*vol_1/vol_2)/css14;
	dx15=(-r8+r9-r12*vol_1/vol_2+r15*vol_1/vol_2)/css15;
	dx16=(-r9-r15*vol_1/vol_2+r20+r23*vol_1/vol_2+r27*vol_1/vol_2)/css16;
	dx17=(2.0*r10-2.0*r18+2.0*r20)/css17;
	dx18=(4.0*r10+6.0*r18-r21-3.0*r26-r32*vol_1/vol_2)/css18;
	dx19=(-2.0*r10+2.0*r18-2.0*r20)/css19;
	dx20=(-4.0*r10-6.0*r18+r21+3.0*r26+r32*vol_1/vol_2)/css20;
	dx21=(r11-r12-r16)/css21;
	dx22=(-r11-r15+r17)/css22;
	dx23=(r13-r14-r16+120.0*r19)/css23;
	dx24=(r13-r17-120.0*r19)/css24;
	dx25=(-r13+r14+r16-120.0*r19)/css25;
	dx26=(-0.5*r13+r28)/css26;
	dx27=(r14+r15-r23-r27)/css27;
	dx28=(-r14+120.0*r19+r23)/css28;
	dx29=(-r14-120.0*r19+r24)/css29;
	dx30=(r15*vol_1/vol_2+r20-r22*vol_1/vol_2+r26)/css30;
	dx31=(-r20-r26-r27*vol_1/vol_2+r32*vol_1/vol_2)/css31;
	dx32=(-r15*vol_1/vol_2-r20+r22*vol_1/vol_2-r26)/css32;
	dx33=(-120.0*r19+r29)/css33;
	dx34=(-120.0*r19+r30)/css34;

	return(0);

}


/* Jacobian of the system (dfdx)*/
int amigoJAC_B4(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_B4(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}

#define	 x1  (amigo_model->sim_results[0][j]) 
#define	 x2  (amigo_model->sim_results[1][j]) 
#define	 x3  (amigo_model->sim_results[2][j]) 
#define	 x4  (amigo_model->sim_results[3][j]) 
#define	 x5  (amigo_model->sim_results[4][j]) 
#define	 x6  (amigo_model->sim_results[5][j]) 
#define	 x7  (amigo_model->sim_results[6][j]) 
#define	 x8  (amigo_model->sim_results[7][j]) 
#define	 x9  (amigo_model->sim_results[8][j]) 
#define	 x10 (amigo_model->sim_results[9][j]) 
#define	 x11 (amigo_model->sim_results[10][j]) 
#define	 x12 (amigo_model->sim_results[11][j]) 
#define	 x13 (amigo_model->sim_results[12][j]) 
#define	 x14 (amigo_model->sim_results[13][j]) 
#define	 x15 (amigo_model->sim_results[14][j]) 
#define	 x16 (amigo_model->sim_results[15][j]) 
#define	 x17 (amigo_model->sim_results[16][j]) 
#define	 x18 (amigo_model->sim_results[17][j]) 
#define	 x19 (amigo_model->sim_results[18][j]) 
#define	 x20 (amigo_model->sim_results[19][j]) 
#define	 x21 (amigo_model->sim_results[20][j]) 
#define	 x22 (amigo_model->sim_results[21][j]) 
#define	 x23 (amigo_model->sim_results[22][j]) 
#define	 x24 (amigo_model->sim_results[23][j]) 
#define	 x25 (amigo_model->sim_results[24][j]) 
#define	 x26 (amigo_model->sim_results[25][j]) 
#define	 x27 (amigo_model->sim_results[26][j]) 
#define	 x28 (amigo_model->sim_results[27][j]) 
#define	 x29 (amigo_model->sim_results[28][j]) 
#define	 x30 (amigo_model->sim_results[29][j]) 
#define	 x31 (amigo_model->sim_results[30][j]) 
#define	 x32 (amigo_model->sim_results[31][j]) 
#define	 x33 (amigo_model->sim_results[32][j]) 
#define	 x34 (amigo_model->sim_results[33][j]) 
#define	 x35 (amigo_model->sim_results[34][j]) 



void amigoRHS_get_OBS_B4(void* data){

	int j;
	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 y1  amigo_model->obs_results[0][j] 
		#define	 y2  amigo_model->obs_results[1][j] 
		#define	 y3  amigo_model->obs_results[2][j] 
		#define	 y4  amigo_model->obs_results[3][j] 
		#define	 y5  amigo_model->obs_results[4][j] 
		#define	 y6  amigo_model->obs_results[5][j] 
		#define	 y7  amigo_model->obs_results[6][j] 
		#define	 y8  amigo_model->obs_results[7][j] 
		#define	 y9  amigo_model->obs_results[8][j] 
		#define	 y10 amigo_model->obs_results[9][j] 
		#define	 y11 amigo_model->obs_results[10][j] 
		#define	 y12 amigo_model->obs_results[11][j] 
		#define	 y13 amigo_model->obs_results[12][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){
				y1=x5;
				y2=x4;
				y3=x3;
				y4=x2;
				y5=x1;
				y6=x29;
				y7=x27;
				y8=x21;
				y9=x15;
				y10=x13;
				y11=x30;
				y12=x32;
				y13=x11;

			}

		 break;

	}

	return(amigo_model);

}

#define	 x1  (amigo_model->sens_results[0][j][k]) 
#define	 x2  (amigo_model->sens_results[1][j][k]) 
#define	 x3  (amigo_model->sens_results[2][j][k]) 
#define	 x4  (amigo_model->sens_results[3][j][k]) 
#define	 x5  (amigo_model->sens_results[4][j][k]) 
#define	 x6  (amigo_model->sens_results[5][j][k]) 
#define	 x7  (amigo_model->sens_results[6][j][k]) 
#define	 x8  (amigo_model->sens_results[7][j][k]) 
#define	 x9  (amigo_model->sens_results[8][j][k]) 
#define	 x10 (amigo_model->sens_results[9][j][k]) 
#define	 x11 (amigo_model->sens_results[10][j][k]) 
#define	 x12 (amigo_model->sens_results[11][j][k]) 
#define	 x13 (amigo_model->sens_results[12][j][k]) 
#define	 x14 (amigo_model->sens_results[13][j][k]) 
#define	 x15 (amigo_model->sens_results[14][j][k]) 
#define	 x16 (amigo_model->sens_results[15][j][k]) 
#define	 x17 (amigo_model->sens_results[16][j][k]) 
#define	 x18 (amigo_model->sens_results[17][j][k]) 
#define	 x19 (amigo_model->sens_results[18][j][k]) 
#define	 x20 (amigo_model->sens_results[19][j][k]) 
#define	 x21 (amigo_model->sens_results[20][j][k]) 
#define	 x22 (amigo_model->sens_results[21][j][k]) 
#define	 x23 (amigo_model->sens_results[22][j][k]) 
#define	 x24 (amigo_model->sens_results[23][j][k]) 
#define	 x25 (amigo_model->sens_results[24][j][k]) 
#define	 x26 (amigo_model->sens_results[25][j][k]) 
#define	 x27 (amigo_model->sens_results[26][j][k]) 
#define	 x28 (amigo_model->sens_results[27][j][k]) 
#define	 x29 (amigo_model->sens_results[28][j][k]) 
#define	 x30 (amigo_model->sens_results[29][j][k]) 
#define	 x31 (amigo_model->sens_results[30][j][k]) 
#define	 x32 (amigo_model->sens_results[31][j][k]) 
#define	 x33 (amigo_model->sens_results[32][j][k]) 
#define	 x34 (amigo_model->sens_results[33][j][k]) 
#define	 x35 (amigo_model->sens_results[34][j][k]) 



void amigoRHS_get_sens_OBS_B4(void* data){
	int j,k;

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 y1  amigo_model->sens_obs[0][j][k] 
		#define	 y2  amigo_model->sens_obs[1][j][k] 
		#define	 y3  amigo_model->sens_obs[2][j][k] 
		#define	 y4  amigo_model->sens_obs[3][j][k] 
		#define	 y5  amigo_model->sens_obs[4][j][k] 
		#define	 y6  amigo_model->sens_obs[5][j][k] 
		#define	 y7  amigo_model->sens_obs[6][j][k] 
		#define	 y8  amigo_model->sens_obs[7][j][k] 
		#define	 y9  amigo_model->sens_obs[8][j][k] 
		#define	 y10 amigo_model->sens_obs[9][j][k] 
		#define	 y11 amigo_model->sens_obs[10][j][k] 
		#define	 y12 amigo_model->sens_obs[11][j][k] 
		#define	 y13 amigo_model->sens_obs[12][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){
					y1=x5;
					y2=x4;
					y3=x3;
					y4=x2;
					y5=x1;
					y6=x29;
					y7=x27;
					y8=x21;
					y9=x15;
					y10=x13;
					y11=x30;
					y12=x32;
					y13=x11;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon_B4(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
   
}