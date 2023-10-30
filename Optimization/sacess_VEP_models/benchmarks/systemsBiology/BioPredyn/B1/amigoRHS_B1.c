/*
Check README.txt
*/

#include <AMIGO_problem.h>

int amigoRHS_B1(realtype t, N_Vector y, N_Vector ydot, void *data);

void amigo_Y_at_tcon_B1(void* amigo_model, realtype t, N_Vector y);

void amigoRHS_get_OBS_B1(void* data);

void amigoRHS_get_sens_OBS_B1(void* data);
 
int amigoRHS_B1(realtype t, N_Vector y, N_Vector ydot, void *data);

double objective_function_B1(double* x, void* data); 
 
//Define the objective function
//This is an useful wrapper since most solvers in C accept a void pointer
//The problem data is passed as a void pointer
double objective_function_B1(double* x, void* data){
    
    //Convert void pointer to AMIGO_problem pointer
    AMIGO_problem* amigo_problem=(AMIGO_problem*)data;
    
	//Copy x to AMIGO_problem->x
    set_AMIGO_problem_pars(x,amigo_problem);
    
	//Increment the number of evaluations counter
    amigo_problem->nevals++;
	
    return(eval_AMIGO_problem_LSQ(amigo_problem));
	
}

/* *** Definition of the states *** */
#define	s_0002 Ith(y,0)
#define	s_0008 Ith(y,1)
#define	s_0009 Ith(y,2)
#define	s_0010 Ith(y,3)
#define	s_0015 Ith(y,4)
#define	s_0016 Ith(y,5)
#define	s_0018 Ith(y,6)
#define	s_0019 Ith(y,7)
#define	s_0025 Ith(y,8)
#define	s_0028 Ith(y,9)
#define	s_0033 Ith(y,10)
#define	s_0037 Ith(y,11)
#define	s_0039 Ith(y,12)
#define	s_0056 Ith(y,13)
#define	s_0061 Ith(y,14)
#define	s_0062 Ith(y,15)
#define	s_0063 Ith(y,16)
#define	s_0066 Ith(y,17)
#define	s_0075 Ith(y,18)
#define	s_0076 Ith(y,19)
#define	s_0077 Ith(y,20)
#define	s_0078 Ith(y,21)
#define	s_0082 Ith(y,22)
#define	s_0086 Ith(y,23)
#define	s_0089 Ith(y,24)
#define	s_0118 Ith(y,25)
#define	s_0120 Ith(y,26)
#define	s_0122 Ith(y,27)
#define	s_0126 Ith(y,28)
#define	s_0141 Ith(y,29)
#define	s_0142 Ith(y,30)
#define	s_0145 Ith(y,31)
#define	s_0146 Ith(y,32)
#define	s_0158 Ith(y,33)
#define	s_0162 Ith(y,34)
#define	s_0165 Ith(y,35)
#define	s_0176 Ith(y,36)
#define	s_0178 Ith(y,37)
#define	s_0180 Ith(y,38)
#define	s_0188 Ith(y,39)
#define	s_0190 Ith(y,40)
#define	s_0201 Ith(y,41)
#define	s_0204 Ith(y,42)
#define	s_0207 Ith(y,43)
#define	s_0209 Ith(y,44)
#define	s_0210 Ith(y,45)
#define	s_0211 Ith(y,46)
#define	s_0218 Ith(y,47)
#define	s_0231 Ith(y,48)
#define	s_0232 Ith(y,49)
#define	s_0258 Ith(y,50)
#define	s_0259 Ith(y,51)
#define	s_0260 Ith(y,52)
#define	s_0261 Ith(y,53)
#define	s_0262 Ith(y,54)
#define	s_0291 Ith(y,55)
#define	s_0295 Ith(y,56)
#define	s_0296 Ith(y,57)
#define	s_0297 Ith(y,58)
#define	s_0298 Ith(y,59)
#define	s_0299 Ith(y,60)
#define	s_0300 Ith(y,61)
#define	s_0301 Ith(y,62)
#define	s_0302 Ith(y,63)
#define	s_0304 Ith(y,64)
#define	s_0306 Ith(y,65)
#define	s_0312 Ith(y,66)
#define	s_0313 Ith(y,67)
#define	s_0314 Ith(y,68)
#define	s_0322 Ith(y,69)
#define	s_0324 Ith(y,70)
#define	s_0325 Ith(y,71)
#define	s_0326 Ith(y,72)
#define	s_0327 Ith(y,73)
#define	s_0328 Ith(y,74)
#define	s_0335 Ith(y,75)
#define	s_0340 Ith(y,76)
#define	s_0349 Ith(y,77)
#define	s_0359 Ith(y,78)
#define	s_0362 Ith(y,79)
#define	s_0367 Ith(y,80)
#define	s_0373 Ith(y,81)
#define	s_0380 Ith(y,82)
#define	s_0386 Ith(y,83)
#define	s_0390 Ith(y,84)
#define	s_0393 Ith(y,85)
#define	s_0394 Ith(y,86)
#define	s_0403 Ith(y,87)
#define	s_0409 Ith(y,88)
#define	s_0419 Ith(y,89)
#define	s_0423 Ith(y,90)
#define	s_0427 Ith(y,91)
#define	s_0434 Ith(y,92)
#define	s_0445 Ith(y,93)
#define	s_0454 Ith(y,94)
#define	s_0455 Ith(y,95)
#define	s_0456 Ith(y,96)
#define	s_0467 Ith(y,97)
#define	s_0471 Ith(y,98)
#define	s_0475 Ith(y,99)
#define	s_0481 Ith(y,100)
#define	s_0493 Ith(y,101)
#define	s_0499 Ith(y,102)
#define	s_0515 Ith(y,103)
#define	s_0516 Ith(y,104)
#define	s_0522 Ith(y,105)
#define	s_0526 Ith(y,106)
#define	s_0529 Ith(y,107)
#define	s_0539 Ith(y,108)
#define	s_0550 Ith(y,109)
#define	s_0551 Ith(y,110)
#define	s_0555 Ith(y,111)
#define	s_0557 Ith(y,112)
#define	s_0563 Ith(y,113)
#define	s_0567 Ith(y,114)
#define	s_0568 Ith(y,115)
#define	s_0573 Ith(y,116)
#define	s_0574 Ith(y,117)
#define	s_0577 Ith(y,118)
#define	s_0581 Ith(y,119)
#define	s_0582 Ith(y,120)
#define	s_0584 Ith(y,121)
#define	s_0586 Ith(y,122)
#define	s_0587 Ith(y,123)
#define	s_0589 Ith(y,124)
#define	s_0595 Ith(y,125)
#define	s_0602 Ith(y,126)
#define	s_0613 Ith(y,127)
#define	s_0615 Ith(y,128)
#define	s_0619 Ith(y,129)
#define	s_0625 Ith(y,130)
#define	s_0629 Ith(y,131)
#define	s_0633 Ith(y,132)
#define	s_0644 Ith(y,133)
#define	s_0645 Ith(y,134)
#define	s_0649 Ith(y,135)
#define	s_0654 Ith(y,136)
#define	s_0656 Ith(y,137)
#define	s_0657 Ith(y,138)
#define	s_0662 Ith(y,139)
#define	s_0664 Ith(y,140)
#define	s_0666 Ith(y,141)
#define	s_0672 Ith(y,142)
#define	s_0680 Ith(y,143)
#define	s_0700 Ith(y,144)
#define	s_0709 Ith(y,145)
#define	s_0710 Ith(y,146)
#define	s_0722 Ith(y,147)
#define	s_0725 Ith(y,148)
#define	s_0739 Ith(y,149)
#define	s_0743 Ith(y,150)
#define	s_0745 Ith(y,151)
#define	s_0750 Ith(y,152)
#define	s_0754 Ith(y,153)
#define	s_0764 Ith(y,154)
#define	s_0765 Ith(y,155)
#define	s_0767 Ith(y,156)
#define	s_0773 Ith(y,157)
#define	s_0782 Ith(y,158)
#define	s_0785 Ith(y,159)
#define	s_0835 Ith(y,160)
#define	s_0836 Ith(y,161)
#define	s_0837 Ith(y,162)
#define	s_0841 Ith(y,163)
#define	s_0849 Ith(y,164)
#define	s_0918 Ith(y,165)
#define	s_0940 Ith(y,166)
#define	s_0943 Ith(y,167)
#define	s_0951 Ith(y,168)
#define	s_0953 Ith(y,169)
#define	s_0955 Ith(y,170)
#define	s_0959 Ith(y,171)
#define	s_0965 Ith(y,172)
#define	s_0969 Ith(y,173)
#define	s_0973 Ith(y,174)
#define	s_0978 Ith(y,175)
#define	s_0979 Ith(y,176)
#define	s_0980 Ith(y,177)
#define	s_0981 Ith(y,178)
#define	s_0991 Ith(y,179)
#define	s_0999 Ith(y,180)
#define	s_1003 Ith(y,181)
#define	s_1006 Ith(y,182)
#define	s_1010 Ith(y,183)
#define	s_1011 Ith(y,184)
#define	s_1012 Ith(y,185)
#define	s_1014 Ith(y,186)
#define	s_1016 Ith(y,187)
#define	s_1020 Ith(y,188)
#define	s_1021 Ith(y,189)
#define	s_1025 Ith(y,190)
#define	s_1029 Ith(y,191)
#define	s_1032 Ith(y,192)
#define	s_1035 Ith(y,193)
#define	s_1038 Ith(y,194)
#define	s_1039 Ith(y,195)
#define	s_1045 Ith(y,196)
#define	s_1048 Ith(y,197)
#define	s_1051 Ith(y,198)
#define	s_1056 Ith(y,199)
#define	s_1059 Ith(y,200)
#define	s_1065 Ith(y,201)
#define	s_1073 Ith(y,202)
#define	s_1084 Ith(y,203)
#define	s_1101 Ith(y,204)
#define	s_1107 Ith(y,205)
#define	s_1151 Ith(y,206)
#define	s_1153 Ith(y,207)
#define	s_1161 Ith(y,208)
#define	s_1176 Ith(y,209)
#define	s_1182 Ith(y,210)
#define	s_1187 Ith(y,211)
#define	s_1191 Ith(y,212)
#define	s_1192 Ith(y,213)
#define	s_1194 Ith(y,214)
#define	s_1195 Ith(y,215)
#define	s_1198 Ith(y,216)
#define	s_1203 Ith(y,217)
#define	s_1207 Ith(y,218)
#define	s_1212 Ith(y,219)
#define	s_1233 Ith(y,220)
#define	s_1234 Ith(y,221)
#define	s_1238 Ith(y,222)
#define	s_1255 Ith(y,223)
#define	s_1266 Ith(y,224)
#define	s_1269 Ith(y,225)
#define	s_1270 Ith(y,226)
#define	s_1271 Ith(y,227)
#define	s_1275 Ith(y,228)
#define	s_1286 Ith(y,229)
#define	s_1302 Ith(y,230)
#define	s_1322 Ith(y,231)
#define	s_1331 Ith(y,232)
#define	s_1337 Ith(y,233)
#define	s_1342 Ith(y,234)
#define	s_1343 Ith(y,235)
#define	s_1346 Ith(y,236)
#define	s_1351 Ith(y,237)
#define	s_1360 Ith(y,238)
#define	s_1364 Ith(y,239)
#define	s_1365 Ith(y,240)
#define	s_1366 Ith(y,241)
#define	s_1376 Ith(y,242)
#define	s_1377 Ith(y,243)
#define	s_1386 Ith(y,244)
#define	s_1399 Ith(y,245)
#define	s_1405 Ith(y,246)
#define	s_1408 Ith(y,247)
#define	s_1413 Ith(y,248)
#define	s_1416 Ith(y,249)
#define	s_1426 Ith(y,250)
#define	s_1427 Ith(y,251)
#define	s_1429 Ith(y,252)
#define	s_1445 Ith(y,253)
#define	s_1447 Ith(y,254)
#define	s_1449 Ith(y,255)
#define	s_1454 Ith(y,256)
#define	s_1467 Ith(y,257)
#define	s_1469 Ith(y,258)
#define	s_1487 Ith(y,259)
#define	s_1520 Ith(y,260)
#define	s_1524 Ith(y,261)
#define	s_1535 Ith(y,262)
#define	s_1537 Ith(y,263)
#define	s_1538 Ith(y,264)
#define	s_1543 Ith(y,265)
#define	s_1545 Ith(y,266)
#define	s_1559 Ith(y,267)
#define	s_1565 Ith(y,268)
#define	s_1569 Ith(y,269)
#define	s_1576 Ith(y,270)
#define	s_1577 Ith(y,271)
#define	s_1578 Ith(y,272)
#define	s_1579 Ith(y,273)
#define	s_1616 Ith(y,274)
#define	s_1620 Ith(y,275)

	/* *** Definition of the sates derivative *** */

#define	ds_0002 Ith(ydot,0)
#define	ds_0008 Ith(ydot,1)
#define	ds_0009 Ith(ydot,2)
#define	ds_0010 Ith(ydot,3)
#define	ds_0015 Ith(ydot,4)
#define	ds_0016 Ith(ydot,5)
#define	ds_0018 Ith(ydot,6)
#define	ds_0019 Ith(ydot,7)
#define	ds_0025 Ith(ydot,8)
#define	ds_0028 Ith(ydot,9)
#define	ds_0033 Ith(ydot,10)
#define	ds_0037 Ith(ydot,11)
#define	ds_0039 Ith(ydot,12)
#define	ds_0056 Ith(ydot,13)
#define	ds_0061 Ith(ydot,14)
#define	ds_0062 Ith(ydot,15)
#define	ds_0063 Ith(ydot,16)
#define	ds_0066 Ith(ydot,17)
#define	ds_0075 Ith(ydot,18)
#define	ds_0076 Ith(ydot,19)
#define	ds_0077 Ith(ydot,20)
#define	ds_0078 Ith(ydot,21)
#define	ds_0082 Ith(ydot,22)
#define	ds_0086 Ith(ydot,23)
#define	ds_0089 Ith(ydot,24)
#define	ds_0118 Ith(ydot,25)
#define	ds_0120 Ith(ydot,26)
#define	ds_0122 Ith(ydot,27)
#define	ds_0126 Ith(ydot,28)
#define	ds_0141 Ith(ydot,29)
#define	ds_0142 Ith(ydot,30)
#define	ds_0145 Ith(ydot,31)
#define	ds_0146 Ith(ydot,32)
#define	ds_0158 Ith(ydot,33)
#define	ds_0162 Ith(ydot,34)
#define	ds_0165 Ith(ydot,35)
#define	ds_0176 Ith(ydot,36)
#define	ds_0178 Ith(ydot,37)
#define	ds_0180 Ith(ydot,38)
#define	ds_0188 Ith(ydot,39)
#define	ds_0190 Ith(ydot,40)
#define	ds_0201 Ith(ydot,41)
#define	ds_0204 Ith(ydot,42)
#define	ds_0207 Ith(ydot,43)
#define	ds_0209 Ith(ydot,44)
#define	ds_0210 Ith(ydot,45)
#define	ds_0211 Ith(ydot,46)
#define	ds_0218 Ith(ydot,47)
#define	ds_0231 Ith(ydot,48)
#define	ds_0232 Ith(ydot,49)
#define	ds_0258 Ith(ydot,50)
#define	ds_0259 Ith(ydot,51)
#define	ds_0260 Ith(ydot,52)
#define	ds_0261 Ith(ydot,53)
#define	ds_0262 Ith(ydot,54)
#define	ds_0291 Ith(ydot,55)
#define	ds_0295 Ith(ydot,56)
#define	ds_0296 Ith(ydot,57)
#define	ds_0297 Ith(ydot,58)
#define	ds_0298 Ith(ydot,59)
#define	ds_0299 Ith(ydot,60)
#define	ds_0300 Ith(ydot,61)
#define	ds_0301 Ith(ydot,62)
#define	ds_0302 Ith(ydot,63)
#define	ds_0304 Ith(ydot,64)
#define	ds_0306 Ith(ydot,65)
#define	ds_0312 Ith(ydot,66)
#define	ds_0313 Ith(ydot,67)
#define	ds_0314 Ith(ydot,68)
#define	ds_0322 Ith(ydot,69)
#define	ds_0324 Ith(ydot,70)
#define	ds_0325 Ith(ydot,71)
#define	ds_0326 Ith(ydot,72)
#define	ds_0327 Ith(ydot,73)
#define	ds_0328 Ith(ydot,74)
#define	ds_0335 Ith(ydot,75)
#define	ds_0340 Ith(ydot,76)
#define	ds_0349 Ith(ydot,77)
#define	ds_0359 Ith(ydot,78)
#define	ds_0362 Ith(ydot,79)
#define	ds_0367 Ith(ydot,80)
#define	ds_0373 Ith(ydot,81)
#define	ds_0380 Ith(ydot,82)
#define	ds_0386 Ith(ydot,83)
#define	ds_0390 Ith(ydot,84)
#define	ds_0393 Ith(ydot,85)
#define	ds_0394 Ith(ydot,86)
#define	ds_0403 Ith(ydot,87)
#define	ds_0409 Ith(ydot,88)
#define	ds_0419 Ith(ydot,89)
#define	ds_0423 Ith(ydot,90)
#define	ds_0427 Ith(ydot,91)
#define	ds_0434 Ith(ydot,92)
#define	ds_0445 Ith(ydot,93)
#define	ds_0454 Ith(ydot,94)
#define	ds_0455 Ith(ydot,95)
#define	ds_0456 Ith(ydot,96)
#define	ds_0467 Ith(ydot,97)
#define	ds_0471 Ith(ydot,98)
#define	ds_0475 Ith(ydot,99)
#define	ds_0481 Ith(ydot,100)
#define	ds_0493 Ith(ydot,101)
#define	ds_0499 Ith(ydot,102)
#define	ds_0515 Ith(ydot,103)
#define	ds_0516 Ith(ydot,104)
#define	ds_0522 Ith(ydot,105)
#define	ds_0526 Ith(ydot,106)
#define	ds_0529 Ith(ydot,107)
#define	ds_0539 Ith(ydot,108)
#define	ds_0550 Ith(ydot,109)
#define	ds_0551 Ith(ydot,110)
#define	ds_0555 Ith(ydot,111)
#define	ds_0557 Ith(ydot,112)
#define	ds_0563 Ith(ydot,113)
#define	ds_0567 Ith(ydot,114)
#define	ds_0568 Ith(ydot,115)
#define	ds_0573 Ith(ydot,116)
#define	ds_0574 Ith(ydot,117)
#define	ds_0577 Ith(ydot,118)
#define	ds_0581 Ith(ydot,119)
#define	ds_0582 Ith(ydot,120)
#define	ds_0584 Ith(ydot,121)
#define	ds_0586 Ith(ydot,122)
#define	ds_0587 Ith(ydot,123)
#define	ds_0589 Ith(ydot,124)
#define	ds_0595 Ith(ydot,125)
#define	ds_0602 Ith(ydot,126)
#define	ds_0613 Ith(ydot,127)
#define	ds_0615 Ith(ydot,128)
#define	ds_0619 Ith(ydot,129)
#define	ds_0625 Ith(ydot,130)
#define	ds_0629 Ith(ydot,131)
#define	ds_0633 Ith(ydot,132)
#define	ds_0644 Ith(ydot,133)
#define	ds_0645 Ith(ydot,134)
#define	ds_0649 Ith(ydot,135)
#define	ds_0654 Ith(ydot,136)
#define	ds_0656 Ith(ydot,137)
#define	ds_0657 Ith(ydot,138)
#define	ds_0662 Ith(ydot,139)
#define	ds_0664 Ith(ydot,140)
#define	ds_0666 Ith(ydot,141)
#define	ds_0672 Ith(ydot,142)
#define	ds_0680 Ith(ydot,143)
#define	ds_0700 Ith(ydot,144)
#define	ds_0709 Ith(ydot,145)
#define	ds_0710 Ith(ydot,146)
#define	ds_0722 Ith(ydot,147)
#define	ds_0725 Ith(ydot,148)
#define	ds_0739 Ith(ydot,149)
#define	ds_0743 Ith(ydot,150)
#define	ds_0745 Ith(ydot,151)
#define	ds_0750 Ith(ydot,152)
#define	ds_0754 Ith(ydot,153)
#define	ds_0764 Ith(ydot,154)
#define	ds_0765 Ith(ydot,155)
#define	ds_0767 Ith(ydot,156)
#define	ds_0773 Ith(ydot,157)
#define	ds_0782 Ith(ydot,158)
#define	ds_0785 Ith(ydot,159)
#define	ds_0835 Ith(ydot,160)
#define	ds_0836 Ith(ydot,161)
#define	ds_0837 Ith(ydot,162)
#define	ds_0841 Ith(ydot,163)
#define	ds_0849 Ith(ydot,164)
#define	ds_0918 Ith(ydot,165)
#define	ds_0940 Ith(ydot,166)
#define	ds_0943 Ith(ydot,167)
#define	ds_0951 Ith(ydot,168)
#define	ds_0953 Ith(ydot,169)
#define	ds_0955 Ith(ydot,170)
#define	ds_0959 Ith(ydot,171)
#define	ds_0965 Ith(ydot,172)
#define	ds_0969 Ith(ydot,173)
#define	ds_0973 Ith(ydot,174)
#define	ds_0978 Ith(ydot,175)
#define	ds_0979 Ith(ydot,176)
#define	ds_0980 Ith(ydot,177)
#define	ds_0981 Ith(ydot,178)
#define	ds_0991 Ith(ydot,179)
#define	ds_0999 Ith(ydot,180)
#define	ds_1003 Ith(ydot,181)
#define	ds_1006 Ith(ydot,182)
#define	ds_1010 Ith(ydot,183)
#define	ds_1011 Ith(ydot,184)
#define	ds_1012 Ith(ydot,185)
#define	ds_1014 Ith(ydot,186)
#define	ds_1016 Ith(ydot,187)
#define	ds_1020 Ith(ydot,188)
#define	ds_1021 Ith(ydot,189)
#define	ds_1025 Ith(ydot,190)
#define	ds_1029 Ith(ydot,191)
#define	ds_1032 Ith(ydot,192)
#define	ds_1035 Ith(ydot,193)
#define	ds_1038 Ith(ydot,194)
#define	ds_1039 Ith(ydot,195)
#define	ds_1045 Ith(ydot,196)
#define	ds_1048 Ith(ydot,197)
#define	ds_1051 Ith(ydot,198)
#define	ds_1056 Ith(ydot,199)
#define	ds_1059 Ith(ydot,200)
#define	ds_1065 Ith(ydot,201)
#define	ds_1073 Ith(ydot,202)
#define	ds_1084 Ith(ydot,203)
#define	ds_1101 Ith(ydot,204)
#define	ds_1107 Ith(ydot,205)
#define	ds_1151 Ith(ydot,206)
#define	ds_1153 Ith(ydot,207)
#define	ds_1161 Ith(ydot,208)
#define	ds_1176 Ith(ydot,209)
#define	ds_1182 Ith(ydot,210)
#define	ds_1187 Ith(ydot,211)
#define	ds_1191 Ith(ydot,212)
#define	ds_1192 Ith(ydot,213)
#define	ds_1194 Ith(ydot,214)
#define	ds_1195 Ith(ydot,215)
#define	ds_1198 Ith(ydot,216)
#define	ds_1203 Ith(ydot,217)
#define	ds_1207 Ith(ydot,218)
#define	ds_1212 Ith(ydot,219)
#define	ds_1233 Ith(ydot,220)
#define	ds_1234 Ith(ydot,221)
#define	ds_1238 Ith(ydot,222)
#define	ds_1255 Ith(ydot,223)
#define	ds_1266 Ith(ydot,224)
#define	ds_1269 Ith(ydot,225)
#define	ds_1270 Ith(ydot,226)
#define	ds_1271 Ith(ydot,227)
#define	ds_1275 Ith(ydot,228)
#define	ds_1286 Ith(ydot,229)
#define	ds_1302 Ith(ydot,230)
#define	ds_1322 Ith(ydot,231)
#define	ds_1331 Ith(ydot,232)
#define	ds_1337 Ith(ydot,233)
#define	ds_1342 Ith(ydot,234)
#define	ds_1343 Ith(ydot,235)
#define	ds_1346 Ith(ydot,236)
#define	ds_1351 Ith(ydot,237)
#define	ds_1360 Ith(ydot,238)
#define	ds_1364 Ith(ydot,239)
#define	ds_1365 Ith(ydot,240)
#define	ds_1366 Ith(ydot,241)
#define	ds_1376 Ith(ydot,242)
#define	ds_1377 Ith(ydot,243)
#define	ds_1386 Ith(ydot,244)
#define	ds_1399 Ith(ydot,245)
#define	ds_1405 Ith(ydot,246)
#define	ds_1408 Ith(ydot,247)
#define	ds_1413 Ith(ydot,248)
#define	ds_1416 Ith(ydot,249)
#define	ds_1426 Ith(ydot,250)
#define	ds_1427 Ith(ydot,251)
#define	ds_1429 Ith(ydot,252)
#define	ds_1445 Ith(ydot,253)
#define	ds_1447 Ith(ydot,254)
#define	ds_1449 Ith(ydot,255)
#define	ds_1454 Ith(ydot,256)
#define	ds_1467 Ith(ydot,257)
#define	ds_1469 Ith(ydot,258)
#define	ds_1487 Ith(ydot,259)
#define	ds_1520 Ith(ydot,260)
#define	ds_1524 Ith(ydot,261)
#define	ds_1535 Ith(ydot,262)
#define	ds_1537 Ith(ydot,263)
#define	ds_1538 Ith(ydot,264)
#define	ds_1543 Ith(ydot,265)
#define	ds_1545 Ith(ydot,266)
#define	ds_1559 Ith(ydot,267)
#define	ds_1565 Ith(ydot,268)
#define	ds_1569 Ith(ydot,269)
#define	ds_1576 Ith(ydot,270)
#define	ds_1577 Ith(ydot,271)
#define	ds_1578 Ith(ydot,272)
#define	ds_1579 Ith(ydot,273)
#define	ds_1616 Ith(ydot,274)
#define	ds_1620 Ith(ydot,275)

	/* *** Definition of the parameters *** */

#define	Vmax_0001     (*amigo_model).pars[0]
#define	Keq_0001      (*amigo_model).pars[1]
#define	Km0025_0001   (*amigo_model).pars[2]
#define	Km0709_0001   (*amigo_model).pars[3]
#define	Km0710_0001   (*amigo_model).pars[4]
#define	Km1399_0001   (*amigo_model).pars[5]
#define	Vmax_0004     (*amigo_model).pars[6]
#define	Keq_0004      (*amigo_model).pars[7]
#define	Km0063_0004   (*amigo_model).pars[8]
#define	Km0709_0004   (*amigo_model).pars[9]
#define	Km0710_0004   (*amigo_model).pars[10]
#define	Km1399_0004   (*amigo_model).pars[11]
#define	Vmax_0005     (*amigo_model).pars[12]
#define	Keq_0005      (*amigo_model).pars[13]
#define	Km1543_0005   (*amigo_model).pars[14]
#define	Km0002_0005   (*amigo_model).pars[15]
#define	Km1538_0005   (*amigo_model).pars[16]
#define	Vmax_0007     (*amigo_model).pars[17]
#define	Keq_0007      (*amigo_model).pars[18]
#define	Km0077_0007   (*amigo_model).pars[19]
#define	Km0312_0007   (*amigo_model).pars[20]
#define	Vmax_0008     (*amigo_model).pars[21]
#define	Keq_0008      (*amigo_model).pars[22]
#define	Km0082_0008   (*amigo_model).pars[23]
#define	Km0380_0008   (*amigo_model).pars[24]
#define	Km0529_0008   (*amigo_model).pars[25]
#define	Km1331_0008   (*amigo_model).pars[26]
#define	Vmax_0012     (*amigo_model).pars[27]
#define	Keq_0012      (*amigo_model).pars[28]
#define	Km0991_0012   (*amigo_model).pars[29]
#define	Km1203_0012   (*amigo_model).pars[30]
#define	Km0118_0012   (*amigo_model).pars[31]
#define	Km1198_0012   (*amigo_model).pars[32]
#define	Vmax_0014     (*amigo_model).pars[33]
#define	Keq_0014      (*amigo_model).pars[34]
#define	Km0142_0014   (*amigo_model).pars[35]
#define	Km0313_0014   (*amigo_model).pars[36]
#define	Km0419_0014   (*amigo_model).pars[37]
#define	Vmax_0015     (*amigo_model).pars[38]
#define	Keq_0015      (*amigo_model).pars[39]
#define	Km0141_0015   (*amigo_model).pars[40]
#define	Km1212_0015   (*amigo_model).pars[41]
#define	Km0142_0015   (*amigo_model).pars[42]
#define	Km1207_0015   (*amigo_model).pars[43]
#define	Vmax_0016     (*amigo_model).pars[44]
#define	Keq_0016      (*amigo_model).pars[45]
#define	Km0178_0016   (*amigo_model).pars[46]
#define	Km1399_0016   (*amigo_model).pars[47]
#define	Km0039_0016   (*amigo_model).pars[48]
#define	Km0456_0016   (*amigo_model).pars[49]
#define	Vmax_0018     (*amigo_model).pars[50]
#define	Keq_0018      (*amigo_model).pars[51]
#define	Km0176_0018   (*amigo_model).pars[52]
#define	Km0991_0018   (*amigo_model).pars[53]
#define	Km0180_0018   (*amigo_model).pars[54]
#define	Km0953_0018   (*amigo_model).pars[55]
#define	Vmax_0020     (*amigo_model).pars[56]
#define	Keq_0020      (*amigo_model).pars[57]
#define	Km0551_0020   (*amigo_model).pars[58]
#define	Km1360_0020   (*amigo_model).pars[59]
#define	Km0349_0020   (*amigo_model).pars[60]
#define	Km1322_0020   (*amigo_model).pars[61]
#define	Vmax_0023     (*amigo_model).pars[62]
#define	Keq_0023      (*amigo_model).pars[63]
#define	Km0162_0023   (*amigo_model).pars[64]
#define	Km0165_0023   (*amigo_model).pars[65]
#define	Vmax_0024     (*amigo_model).pars[66]
#define	Keq_0024      (*amigo_model).pars[67]
#define	Km0232_0024   (*amigo_model).pars[68]
#define	Km0373_0024   (*amigo_model).pars[69]
#define	Km0162_0024   (*amigo_model).pars[70]
#define	Km0529_0024   (*amigo_model).pars[71]
#define	Vmax_0027     (*amigo_model).pars[72]
#define	Keq_0027      (*amigo_model).pars[73]
#define	Km0835_0027   (*amigo_model).pars[74]
#define	Km0454_0027   (*amigo_model).pars[75]
#define	Vmax_0029     (*amigo_model).pars[76]
#define	Keq_0029      (*amigo_model).pars[77]
#define	Km0010_0029   (*amigo_model).pars[78]
#define	Km0291_0029   (*amigo_model).pars[79]
#define	Km0456_0029   (*amigo_model).pars[80]
#define	Vmax_0032     (*amigo_model).pars[81]
#define	Keq_0032      (*amigo_model).pars[82]
#define	Km0390_0032   (*amigo_model).pars[83]
#define	Km0423_0032   (*amigo_model).pars[84]
#define	Km1322_0032   (*amigo_model).pars[85]
#define	Vmax_0038     (*amigo_model).pars[86]
#define	Keq_0038      (*amigo_model).pars[87]
#define	Km0577_0038   (*amigo_model).pars[88]
#define	Km0158_0038   (*amigo_model).pars[89]
#define	Km0722_0038   (*amigo_model).pars[90]
#define	Vmax_0039     (*amigo_model).pars[91]
#define	Keq_0039      (*amigo_model).pars[92]
#define	Km0210_0039   (*amigo_model).pars[93]
#define	Km0211_0039   (*amigo_model).pars[94]
#define	Vmax_0040     (*amigo_model).pars[95]
#define	Keq_0040      (*amigo_model).pars[96]
#define	Km0349_0040   (*amigo_model).pars[97]
#define	Km0210_0040   (*amigo_model).pars[98]
#define	Km1322_0040   (*amigo_model).pars[99]
#define	Vmax_0041     (*amigo_model).pars[100]
#define	Keq_0041      (*amigo_model).pars[101]
#define	Km0231_0041   (*amigo_model).pars[102]
#define	Km1212_0041   (*amigo_model).pars[103]
#define	Km1207_0041   (*amigo_model).pars[104]
#define	Km1445_0041   (*amigo_model).pars[105]
#define	Vmax_0060     (*amigo_model).pars[106]
#define	Keq_0060      (*amigo_model).pars[107]
#define	Km0165_0060   (*amigo_model).pars[108]
#define	Km0009_0060   (*amigo_model).pars[109]
#define	Vmax_0061     (*amigo_model).pars[110]
#define	Keq_0061      (*amigo_model).pars[111]
#define	Km0009_0061   (*amigo_model).pars[112]
#define	Km1198_0061   (*amigo_model).pars[113]
#define	Km0010_0061   (*amigo_model).pars[114]
#define	Km1203_0061   (*amigo_model).pars[115]
#define	Vmax_0065     (*amigo_model).pars[116]
#define	Keq_0065      (*amigo_model).pars[117]
#define	Km0261_0065   (*amigo_model).pars[118]
#define	Km1360_0065   (*amigo_model).pars[119]
#define	Km0324_0065   (*amigo_model).pars[120]
#define	Km1322_0065   (*amigo_model).pars[121]
#define	Vmax_0079     (*amigo_model).pars[122]
#define	Keq_0079      (*amigo_model).pars[123]
#define	Km0301_0079   (*amigo_model).pars[124]
#define	Km0434_0079   (*amigo_model).pars[125]
#define	Km0999_0079   (*amigo_model).pars[126]
#define	Km0302_0079   (*amigo_model).pars[127]
#define	Km0394_0079   (*amigo_model).pars[128]
#define	Km0991_0079   (*amigo_model).pars[129]
#define	Km1322_0079   (*amigo_model).pars[130]
#define	Vmax_0080     (*amigo_model).pars[131]
#define	Keq_0080      (*amigo_model).pars[132]
#define	Km0306_0080   (*amigo_model).pars[133]
#define	Km1212_0080   (*amigo_model).pars[134]
#define	Km0322_0080   (*amigo_model).pars[135]
#define	Km1207_0080   (*amigo_model).pars[136]
#define	Vmax_0091     (*amigo_model).pars[137]
#define	Keq_0091      (*amigo_model).pars[138]
#define	Km0335_0091   (*amigo_model).pars[139]
#define	Km0340_0091   (*amigo_model).pars[140]
#define	Vmax_0096     (*amigo_model).pars[141]
#define	Keq_0096      (*amigo_model).pars[142]
#define	Km0146_0096   (*amigo_model).pars[143]
#define	Km1212_0096   (*amigo_model).pars[144]
#define	Km0016_0096   (*amigo_model).pars[145]
#define	Km1207_0096   (*amigo_model).pars[146]
#define	Vmax_0097     (*amigo_model).pars[147]
#define	Keq_0097      (*amigo_model).pars[148]
#define	Km1399_0097   (*amigo_model).pars[149]
#define	Km0146_0097   (*amigo_model).pars[150]
#define	Km0456_0097   (*amigo_model).pars[151]
#define	Vmax_0103     (*amigo_model).pars[152]
#define	Keq_0103      (*amigo_model).pars[153]
#define	Km0373_0103   (*amigo_model).pars[154]
#define	Km0367_0103   (*amigo_model).pars[155]
#define	Km0529_0103   (*amigo_model).pars[156]
#define	Vmax_0108     (*amigo_model).pars[157]
#define	Keq_0108      (*amigo_model).pars[158]
#define	Km0373_0108   (*amigo_model).pars[159]
#define	Km0434_0108   (*amigo_model).pars[160]
#define	Km0445_0108   (*amigo_model).pars[161]
#define	Km0394_0108   (*amigo_model).pars[162]
#define	Km1101_0108   (*amigo_model).pars[163]
#define	Km1322_0108   (*amigo_model).pars[164]
#define	Vmax_0111     (*amigo_model).pars[165]
#define	Keq_0111      (*amigo_model).pars[166]
#define	Km0373_0111   (*amigo_model).pars[167]
#define	Km0362_0111   (*amigo_model).pars[168]
#define	Km0529_0111   (*amigo_model).pars[169]
#define	Vmax_0115     (*amigo_model).pars[170]
#define	Keq_0115      (*amigo_model).pars[171]
#define	Km0434_0115   (*amigo_model).pars[172]
#define	Km1192_0115   (*amigo_model).pars[173]
#define	Km0394_0115   (*amigo_model).pars[174]
#define	Km1191_0115   (*amigo_model).pars[175]
#define	Vmax_0118     (*amigo_model).pars[176]
#define	Keq_0118      (*amigo_model).pars[177]
#define	Km0145_0118   (*amigo_model).pars[178]
#define	Km0991_0118   (*amigo_model).pars[179]
#define	Km0180_0118   (*amigo_model).pars[180]
#define	Km1182_0118   (*amigo_model).pars[181]
#define	Vmax_0142     (*amigo_model).pars[182]
#define	Keq_0142      (*amigo_model).pars[183]
#define	Km0386_0142   (*amigo_model).pars[184]
#define	Km0434_0142   (*amigo_model).pars[185]
#define	Km0394_0142   (*amigo_model).pars[186]
#define	Km0423_0142   (*amigo_model).pars[187]
#define	Vmax_0144     (*amigo_model).pars[188]
#define	Keq_0144      (*amigo_model).pars[189]
#define	Km1413_0144   (*amigo_model).pars[190]
#define	Km0386_0144   (*amigo_model).pars[191]
#define	Km1012_0144   (*amigo_model).pars[192]
#define	Vmax_0148     (*amigo_model).pars[193]
#define	Keq_0148      (*amigo_model).pars[194]
#define	Km0423_0148   (*amigo_model).pars[195]
#define	Km0434_0148   (*amigo_model).pars[196]
#define	Km0394_0148   (*amigo_model).pars[197]
#define	Vmax_0151     (*amigo_model).pars[198]
#define	Keq_0151      (*amigo_model).pars[199]
#define	Km0299_0151   (*amigo_model).pars[200]
#define	Km0403_0151   (*amigo_model).pars[201]
#define	Km0725_0151   (*amigo_model).pars[202]
#define	Vmax_0152     (*amigo_model).pars[203]
#define	Keq_0152      (*amigo_model).pars[204]
#define	Km0393_0152   (*amigo_model).pars[205]
#define	Km0423_0152   (*amigo_model).pars[206]
#define	Km0725_0152   (*amigo_model).pars[207]
#define	Vmax_0153     (*amigo_model).pars[208]
#define	Keq_0153      (*amigo_model).pars[209]
#define	Km0785_0153   (*amigo_model).pars[210]
#define	Km0849_0153   (*amigo_model).pars[211]
#define	Km0973_0153   (*amigo_model).pars[212]
#define	Km0393_0153   (*amigo_model).pars[213]
#define	Km0739_0153   (*amigo_model).pars[214]
#define	Km1322_0153   (*amigo_model).pars[215]
#define	Vmax_0154     (*amigo_model).pars[216]
#define	Keq_0154      (*amigo_model).pars[217]
#define	Km0298_0154   (*amigo_model).pars[218]
#define	Km0434_0154   (*amigo_model).pars[219]
#define	Km0201_0154   (*amigo_model).pars[220]
#define	Km0394_0154   (*amigo_model).pars[221]
#define	Vmax_0165     (*amigo_model).pars[222]
#define	Keq_0165      (*amigo_model).pars[223]
#define	Km0359_0165   (*amigo_model).pars[224]
#define	Km1203_0165   (*amigo_model).pars[225]
#define	Km0680_0165   (*amigo_model).pars[226]
#define	Km1198_0165   (*amigo_model).pars[227]
#define	Vmax_0173     (*amigo_model).pars[228]
#define	Keq_0173      (*amigo_model).pars[229]
#define	Km0359_0173   (*amigo_model).pars[230]
#define	Km1207_0173   (*amigo_model).pars[231]
#define	Km0362_0173   (*amigo_model).pars[232]
#define	Km1212_0173   (*amigo_model).pars[233]
#define	Vmax_0195     (*amigo_model).pars[234]
#define	Keq_0195      (*amigo_model).pars[235]
#define	Km0568_0195   (*amigo_model).pars[236]
#define	Km1543_0195   (*amigo_model).pars[237]
#define	Km0409_0195   (*amigo_model).pars[238]
#define	Km1538_0195   (*amigo_model).pars[239]
#define	Vmax_0202     (*amigo_model).pars[240]
#define	Keq_0202      (*amigo_model).pars[241]
#define	Km0427_0202   (*amigo_model).pars[242]
#define	Km1386_0202   (*amigo_model).pars[243]
#define	Km0633_0202   (*amigo_model).pars[244]
#define	Km1187_0202   (*amigo_model).pars[245]
#define	Vmax_0203     (*amigo_model).pars[246]
#define	Keq_0203      (*amigo_model).pars[247]
#define	Km0515_0203   (*amigo_model).pars[248]
#define	Km0999_0203   (*amigo_model).pars[249]
#define	Km0427_0203   (*amigo_model).pars[250]
#define	Km0991_0203   (*amigo_model).pars[251]
#define	Km1399_0203   (*amigo_model).pars[252]
#define	Vmax_0207     (*amigo_model).pars[253]
#define	Keq_0207      (*amigo_model).pars[254]
#define	Km0015_0207   (*amigo_model).pars[255]
#define	Km0725_0207   (*amigo_model).pars[256]
#define	Km0965_0207   (*amigo_model).pars[257]
#define	Vmax_0208     (*amigo_model).pars[258]
#define	Keq_0208      (*amigo_model).pars[259]
#define	Km0434_0208   (*amigo_model).pars[260]
#define	Km0973_0208   (*amigo_model).pars[261]
#define	Km0979_0208   (*amigo_model).pars[262]
#define	Km0015_0208   (*amigo_model).pars[263]
#define	Km0423_0208   (*amigo_model).pars[264]
#define	Km0633_0208   (*amigo_model).pars[265]
#define	Vmax_0211     (*amigo_model).pars[266]
#define	Keq_0211      (*amigo_model).pars[267]
#define	Km0434_0211   (*amigo_model).pars[268]
#define	Km0973_0211   (*amigo_model).pars[269]
#define	Km0999_0211   (*amigo_model).pars[270]
#define	Km0423_0211   (*amigo_model).pars[271]
#define	Km0633_0211   (*amigo_model).pars[272]
#define	Km0969_0211   (*amigo_model).pars[273]
#define	Km0991_0211   (*amigo_model).pars[274]
#define	Vmax_0214     (*amigo_model).pars[275]
#define	Keq_0214      (*amigo_model).pars[276]
#define	Km0455_0214   (*amigo_model).pars[277]
#define	Km0973_0214   (*amigo_model).pars[278]
#define	Km1194_0214   (*amigo_model).pars[279]
#define	Km1322_0214   (*amigo_model).pars[280]
#define	Vmax_0215     (*amigo_model).pars[281]
#define	Keq_0215      (*amigo_model).pars[282]
#define	Km0434_0215   (*amigo_model).pars[283]
#define	Km0973_0215   (*amigo_model).pars[284]
#define	Km0295_0215   (*amigo_model).pars[285]
#define	Km0394_0215   (*amigo_model).pars[286]
#define	Vmax_0216     (*amigo_model).pars[287]
#define	Keq_0216      (*amigo_model).pars[288]
#define	Km0991_0216   (*amigo_model).pars[289]
#define	Km1271_0216   (*amigo_model).pars[290]
#define	Km0180_0216   (*amigo_model).pars[291]
#define	Km0973_0216   (*amigo_model).pars[292]
#define	Vmax_0219     (*amigo_model).pars[293]
#define	Keq_0219      (*amigo_model).pars[294]
#define	Km0295_0219   (*amigo_model).pars[295]
#define	Km1212_0219   (*amigo_model).pars[296]
#define	Km0978_0219   (*amigo_model).pars[297]
#define	Km1207_0219   (*amigo_model).pars[298]
#define	Km1322_0219   (*amigo_model).pars[299]
#define	Vmax_0225     (*amigo_model).pars[300]
#define	Keq_0225      (*amigo_model).pars[301]
#define	Km0434_0225   (*amigo_model).pars[302]
#define	Km1386_0225   (*amigo_model).pars[303]
#define	Km0326_0225   (*amigo_model).pars[304]
#define	Km0633_0225   (*amigo_model).pars[305]
#define	Vmax_0226     (*amigo_model).pars[306]
#define	Keq_0226      (*amigo_model).pars[307]
#define	Km0394_0226   (*amigo_model).pars[308]
#define	Km1322_0226   (*amigo_model).pars[309]
#define	Km0434_0226   (*amigo_model).pars[310]
#define	Vmax_0231     (*amigo_model).pars[311]
#define	Keq_0231      (*amigo_model).pars[312]
#define	Km0262_0231   (*amigo_model).pars[313]
#define	Km1212_0231   (*amigo_model).pars[314]
#define	Km0122_0231   (*amigo_model).pars[315]
#define	Km1207_0231   (*amigo_model).pars[316]
#define	Vmax_0233     (*amigo_model).pars[317]
#define	Keq_0233      (*amigo_model).pars[318]
#define	Km0664_0233   (*amigo_model).pars[319]
#define	Km1212_0233   (*amigo_model).pars[320]
#define	Km1275_0233   (*amigo_model).pars[321]
#define	Km0662_0233   (*amigo_model).pars[322]
#define	Km1207_0233   (*amigo_model).pars[323]
#define	Vmax_0234     (*amigo_model).pars[324]
#define	Keq_0234      (*amigo_model).pars[325]
#define	Km1207_0234   (*amigo_model).pars[326]
#define	Km1578_0234   (*amigo_model).pars[327]
#define	Km0456_0234   (*amigo_model).pars[328]
#define	Km1212_0234   (*amigo_model).pars[329]
#define	Km1579_0234   (*amigo_model).pars[330]
#define	Vmax_0235     (*amigo_model).pars[331]
#define	Keq_0235      (*amigo_model).pars[332]
#define	Km0297_0235   (*amigo_model).pars[333]
#define	Km1198_0235   (*amigo_model).pars[334]
#define	Km0209_0235   (*amigo_model).pars[335]
#define	Km0456_0235   (*amigo_model).pars[336]
#define	Km1203_0235   (*amigo_model).pars[337]
#define	Vmax_0236     (*amigo_model).pars[338]
#define	Keq_0236      (*amigo_model).pars[339]
#define	Km0209_0236   (*amigo_model).pars[340]
#define	Km1212_0236   (*amigo_model).pars[341]
#define	Km0296_0236   (*amigo_model).pars[342]
#define	Km1207_0236   (*amigo_model).pars[343]
#define	Vmax_0237     (*amigo_model).pars[344]
#define	Keq_0237      (*amigo_model).pars[345]
#define	Km1212_0237   (*amigo_model).pars[346]
#define	Km1579_0237   (*amigo_model).pars[347]
#define	Km1207_0237   (*amigo_model).pars[348]
#define	Km1569_0237   (*amigo_model).pars[349]
#define	Vmax_0238     (*amigo_model).pars[350]
#define	Keq_0238      (*amigo_model).pars[351]
#define	Km0296_0238   (*amigo_model).pars[352]
#define	Km1212_0238   (*amigo_model).pars[353]
#define	Km1275_0238   (*amigo_model).pars[354]
#define	Km1207_0238   (*amigo_model).pars[355]
#define	Km1576_0238   (*amigo_model).pars[356]
#define	Vmax_0239     (*amigo_model).pars[357]
#define	Keq_0239      (*amigo_model).pars[358]
#define	Km1212_0239   (*amigo_model).pars[359]
#define	Km1275_0239   (*amigo_model).pars[360]
#define	Km1576_0239   (*amigo_model).pars[361]
#define	Km1207_0239   (*amigo_model).pars[362]
#define	Km1577_0239   (*amigo_model).pars[363]
#define	Vmax_0240     (*amigo_model).pars[364]
#define	Keq_0240      (*amigo_model).pars[365]
#define	Km1212_0240   (*amigo_model).pars[366]
#define	Km1275_0240   (*amigo_model).pars[367]
#define	Km1577_0240   (*amigo_model).pars[368]
#define	Km1207_0240   (*amigo_model).pars[369]
#define	Km1578_0240   (*amigo_model).pars[370]
#define	Vmax_0241     (*amigo_model).pars[371]
#define	Keq_0241      (*amigo_model).pars[372]
#define	Km0122_0241   (*amigo_model).pars[373]
#define	Km1212_0241   (*amigo_model).pars[374]
#define	Km1275_0241   (*amigo_model).pars[375]
#define	Km0297_0241   (*amigo_model).pars[376]
#define	Km1207_0241   (*amigo_model).pars[377]
#define	Vmax_0242     (*amigo_model).pars[378]
#define	Keq_0242      (*amigo_model).pars[379]
#define	Km0657_0242   (*amigo_model).pars[380]
#define	Km1212_0242   (*amigo_model).pars[381]
#define	Km1275_0242   (*amigo_model).pars[382]
#define	Km0664_0242   (*amigo_model).pars[383]
#define	Km1207_0242   (*amigo_model).pars[384]
#define	Vmax_0243     (*amigo_model).pars[385]
#define	Keq_0243      (*amigo_model).pars[386]
#define	Km0700_0243   (*amigo_model).pars[387]
#define	Km0657_0243   (*amigo_model).pars[388]
#define	Vmax_0244     (*amigo_model).pars[389]
#define	Keq_0244      (*amigo_model).pars[390]
#define	Km0662_0244   (*amigo_model).pars[391]
#define	Km1212_0244   (*amigo_model).pars[392]
#define	Km0666_0244   (*amigo_model).pars[393]
#define	Km1207_0244   (*amigo_model).pars[394]
#define	Vmax_0250     (*amigo_model).pars[395]
#define	Keq_0250      (*amigo_model).pars[396]
#define	Km0434_0250   (*amigo_model).pars[397]
#define	Km0445_0250   (*amigo_model).pars[398]
#define	Km0999_0250   (*amigo_model).pars[399]
#define	Km0394_0250   (*amigo_model).pars[400]
#define	Km0455_0250   (*amigo_model).pars[401]
#define	Km0991_0250   (*amigo_model).pars[402]
#define	Km1322_0250   (*amigo_model).pars[403]
#define	Vmax_0257     (*amigo_model).pars[404]
#define	Keq_0257      (*amigo_model).pars[405]
#define	Km0539_0257   (*amigo_model).pars[406]
#define	Km1331_0257   (*amigo_model).pars[407]
#define	Km0471_0257   (*amigo_model).pars[408]
#define	Km0633_0257   (*amigo_model).pars[409]
#define	Vmax_0259     (*amigo_model).pars[410]
#define	Keq_0259      (*amigo_model).pars[411]
#define	Km0475_0259   (*amigo_model).pars[412]
#define	Km1212_0259   (*amigo_model).pars[413]
#define	Km1275_0259   (*amigo_model).pars[414]
#define	Km0481_0259   (*amigo_model).pars[415]
#define	Km1207_0259   (*amigo_model).pars[416]
#define	Vmax_0267     (*amigo_model).pars[417]
#define	Keq_0267      (*amigo_model).pars[418]
#define	Km0481_0267   (*amigo_model).pars[419]
#define	Km1212_0267   (*amigo_model).pars[420]
#define	Km1275_0267   (*amigo_model).pars[421]
#define	Km0493_0267   (*amigo_model).pars[422]
#define	Km1207_0267   (*amigo_model).pars[423]
#define	Vmax_0269     (*amigo_model).pars[424]
#define	Keq_0269      (*amigo_model).pars[425]
#define	Km0493_0269   (*amigo_model).pars[426]
#define	Km1212_0269   (*amigo_model).pars[427]
#define	Km1275_0269   (*amigo_model).pars[428]
#define	Km0499_0269   (*amigo_model).pars[429]
#define	Km1207_0269   (*amigo_model).pars[430]
#define	Vmax_0278     (*amigo_model).pars[431]
#define	Keq_0278      (*amigo_model).pars[432]
#define	Km0515_0278   (*amigo_model).pars[433]
#define	Km1377_0278   (*amigo_model).pars[434]
#define	Vmax_0279     (*amigo_model).pars[435]
#define	Keq_0279      (*amigo_model).pars[436]
#define	Km0324_0279   (*amigo_model).pars[437]
#define	Km0515_0279   (*amigo_model).pars[438]
#define	Km1322_0279   (*amigo_model).pars[439]
#define	Vmax_0280     (*amigo_model).pars[440]
#define	Keq_0280      (*amigo_model).pars[441]
#define	Km0516_0280   (*amigo_model).pars[442]
#define	Km0940_0280   (*amigo_model).pars[443]
#define	Vmax_0300     (*amigo_model).pars[444]
#define	Keq_0300      (*amigo_model).pars[445]
#define	Km0373_0300   (*amigo_model).pars[446]
#define	Km1271_0300   (*amigo_model).pars[447]
#define	Km0522_0300   (*amigo_model).pars[448]
#define	Km0529_0300   (*amigo_model).pars[449]
#define	Vmax_0302     (*amigo_model).pars[450]
#define	Keq_0302      (*amigo_model).pars[451]
#define	Km0522_0302   (*amigo_model).pars[452]
#define	Km0516_0302   (*amigo_model).pars[453]
#define	Vmax_0307     (*amigo_model).pars[454]
#define	Keq_0307      (*amigo_model).pars[455]
#define	Km0419_0307   (*amigo_model).pars[456]
#define	Km0434_0307   (*amigo_model).pars[457]
#define	Km1559_0307   (*amigo_model).pars[458]
#define	Km0394_0307   (*amigo_model).pars[459]
#define	Km0539_0307   (*amigo_model).pars[460]
#define	Km1322_0307   (*amigo_model).pars[461]
#define	Vmax_0309     (*amigo_model).pars[462]
#define	Keq_0309      (*amigo_model).pars[463]
#define	Km1012_0309   (*amigo_model).pars[464]
#define	Km1039_0309   (*amigo_model).pars[465]
#define	Km0980_0309   (*amigo_model).pars[466]
#define	Vmax_0310     (*amigo_model).pars[467]
#define	Keq_0310      (*amigo_model).pars[468]
#define	Km0980_0310   (*amigo_model).pars[469]
#define	Km0178_0310   (*amigo_model).pars[470]
#define	Km0419_0310   (*amigo_model).pars[471]
#define	Km0981_0310   (*amigo_model).pars[472]
#define	Vmax_0311     (*amigo_model).pars[473]
#define	Keq_0311      (*amigo_model).pars[474]
#define	Km0981_0311   (*amigo_model).pars[475]
#define	Km1233_0311   (*amigo_model).pars[476]
#define	Km0362_0311   (*amigo_model).pars[477]
#define	Km0980_0311   (*amigo_model).pars[478]
#define	Vmax_0312     (*amigo_model).pars[479]
#define	Keq_0312      (*amigo_model).pars[480]
#define	Km0841_0312   (*amigo_model).pars[481]
#define	Km1234_0312   (*amigo_model).pars[482]
#define	Km0362_0312   (*amigo_model).pars[483]
#define	Km0981_0312   (*amigo_model).pars[484]
#define	Vmax_0317     (*amigo_model).pars[485]
#define	Keq_0317      (*amigo_model).pars[486]
#define	Km1059_0317   (*amigo_model).pars[487]
#define	Km1212_0317   (*amigo_model).pars[488]
#define	Km1275_0317   (*amigo_model).pars[489]
#define	Km0262_0317   (*amigo_model).pars[490]
#define	Km0722_0317   (*amigo_model).pars[491]
#define	Km1207_0317   (*amigo_model).pars[492]
#define	Vmax_0326     (*amigo_model).pars[493]
#define	Keq_0326      (*amigo_model).pars[494]
#define	Km0419_0326   (*amigo_model).pars[495]
#define	Km0654_0326   (*amigo_model).pars[496]
#define	Km0589_0326   (*amigo_model).pars[497]
#define	Vmax_0330     (*amigo_model).pars[498]
#define	Keq_0330      (*amigo_model).pars[499]
#define	Km0394_0330   (*amigo_model).pars[500]
#define	Km0613_0330   (*amigo_model).pars[501]
#define	Km0434_0330   (*amigo_model).pars[502]
#define	Km0615_0330   (*amigo_model).pars[503]
#define	Vmax_0336     (*amigo_model).pars[504]
#define	Keq_0336      (*amigo_model).pars[505]
#define	Km0529_0336   (*amigo_model).pars[506]
#define	Km1524_0336   (*amigo_model).pars[507]
#define	Km0380_0336   (*amigo_model).pars[508]
#define	Km0619_0336   (*amigo_model).pars[509]
#define	Vmax_0337     (*amigo_model).pars[510]
#define	Keq_0337      (*amigo_model).pars[511]
#define	Km1331_0337   (*amigo_model).pars[512]
#define	Km0619_0337   (*amigo_model).pars[513]
#define	Km1322_0337   (*amigo_model).pars[514]
#define	Vmax_0339     (*amigo_model).pars[515]
#define	Keq_0339      (*amigo_model).pars[516]
#define	Km0061_0339   (*amigo_model).pars[517]
#define	Km1275_0339   (*amigo_model).pars[518]
#define	Km0837_0339   (*amigo_model).pars[519]
#define	Km1269_0339   (*amigo_model).pars[520]
#define	Vmax_0340     (*amigo_model).pars[521]
#define	Keq_0340      (*amigo_model).pars[522]
#define	Km1084_0340   (*amigo_model).pars[523]
#define	Km1445_0340   (*amigo_model).pars[524]
#define	Km0475_0340   (*amigo_model).pars[525]
#define	Vmax_0344     (*amigo_model).pars[526]
#define	Keq_0344      (*amigo_model).pars[527]
#define	Km0625_0344   (*amigo_model).pars[528]
#define	Km1212_0344   (*amigo_model).pars[529]
#define	Km1207_0344   (*amigo_model).pars[530]
#define	Km1487_0344   (*amigo_model).pars[531]
#define	Vmax_0349     (*amigo_model).pars[532]
#define	Keq_0349      (*amigo_model).pars[533]
#define	Km1194_0349   (*amigo_model).pars[534]
#define	Km0061_0349   (*amigo_model).pars[535]
#define	Vmax_0352     (*amigo_model).pars[536]
#define	Keq_0352      (*amigo_model).pars[537]
#define	Km0016_0352   (*amigo_model).pars[538]
#define	Km0232_0352   (*amigo_model).pars[539]
#define	Vmax_0353     (*amigo_model).pars[540]
#define	Keq_0353      (*amigo_model).pars[541]
#define	Km0008_0353   (*amigo_model).pars[542]
#define	Km0056_0353   (*amigo_model).pars[543]
#define	Vmax_0355     (*amigo_model).pars[544]
#define	Keq_0355      (*amigo_model).pars[545]
#define	Km0943_0355   (*amigo_model).pars[546]
#define	Km1376_0355   (*amigo_model).pars[547]
#define	Km0633_0355   (*amigo_model).pars[548]
#define	Km0745_0355   (*amigo_model).pars[549]
#define	Vmax_0361     (*amigo_model).pars[550]
#define	Keq_0361      (*amigo_model).pars[551]
#define	Km0645_0361   (*amigo_model).pars[552]
#define	Km0743_0361   (*amigo_model).pars[553]
#define	Km0644_0361   (*amigo_model).pars[554]
#define	Km0739_0361   (*amigo_model).pars[555]
#define	Vmax_0362     (*amigo_model).pars[556]
#define	Keq_0362      (*amigo_model).pars[557]
#define	Km0644_0362   (*amigo_model).pars[558]
#define	Km0645_0362   (*amigo_model).pars[559]
#define	Km1107_0362   (*amigo_model).pars[560]
#define	Vmax_0364     (*amigo_model).pars[561]
#define	Keq_0364      (*amigo_model).pars[562]
#define	Km0656_0364   (*amigo_model).pars[563]
#define	Km0633_0364   (*amigo_model).pars[564]
#define	Km0654_0364   (*amigo_model).pars[565]
#define	Vmax_0366     (*amigo_model).pars[566]
#define	Keq_0366      (*amigo_model).pars[567]
#define	Km0188_0366   (*amigo_model).pars[568]
#define	Km1360_0366   (*amigo_model).pars[569]
#define	Vmax_0386     (*amigo_model).pars[570]
#define	Keq_0386      (*amigo_model).pars[571]
#define	Km0595_0386   (*amigo_model).pars[572]
#define	Km1101_0386   (*amigo_model).pars[573]
#define	Km1212_0386   (*amigo_model).pars[574]
#define	Km0456_0386   (*amigo_model).pars[575]
#define	Km0529_0386   (*amigo_model).pars[576]
#define	Km1065_0386   (*amigo_model).pars[577]
#define	Km1207_0386   (*amigo_model).pars[578]
#define	Vmax_0387     (*amigo_model).pars[579]
#define	Keq_0387      (*amigo_model).pars[580]
#define	Km1065_0387   (*amigo_model).pars[581]
#define	Km1101_0387   (*amigo_model).pars[582]
#define	Km1212_0387   (*amigo_model).pars[583]
#define	Km0456_0387   (*amigo_model).pars[584]
#define	Km0529_0387   (*amigo_model).pars[585]
#define	Km1161_0387   (*amigo_model).pars[586]
#define	Km1207_0387   (*amigo_model).pars[587]
#define	Vmax_0389     (*amigo_model).pars[588]
#define	Keq_0389      (*amigo_model).pars[589]
#define	Km1101_0389   (*amigo_model).pars[590]
#define	Km1161_0389   (*amigo_model).pars[591]
#define	Km1212_0389   (*amigo_model).pars[592]
#define	Km0456_0389   (*amigo_model).pars[593]
#define	Km0529_0389   (*amigo_model).pars[594]
#define	Km1207_0389   (*amigo_model).pars[595]
#define	Km1286_0389   (*amigo_model).pars[596]
#define	Vmax_0391     (*amigo_model).pars[597]
#define	Keq_0391      (*amigo_model).pars[598]
#define	Km1101_0391   (*amigo_model).pars[599]
#define	Km1212_0391   (*amigo_model).pars[600]
#define	Km1286_0391   (*amigo_model).pars[601]
#define	Km0456_0391   (*amigo_model).pars[602]
#define	Km0529_0391   (*amigo_model).pars[603]
#define	Km1207_0391   (*amigo_model).pars[604]
#define	Km1449_0391   (*amigo_model).pars[605]
#define	Vmax_0393     (*amigo_model).pars[606]
#define	Keq_0393      (*amigo_model).pars[607]
#define	Km1101_0393   (*amigo_model).pars[608]
#define	Km1212_0393   (*amigo_model).pars[609]
#define	Km1449_0393   (*amigo_model).pars[610]
#define	Km0456_0393   (*amigo_model).pars[611]
#define	Km0529_0393   (*amigo_model).pars[612]
#define	Km1084_0393   (*amigo_model).pars[613]
#define	Km1207_0393   (*amigo_model).pars[614]
#define	Vmax_0397     (*amigo_model).pars[615]
#define	Keq_0397      (*amigo_model).pars[616]
#define	Km1101_0397   (*amigo_model).pars[617]
#define	Km1212_0397   (*amigo_model).pars[618]
#define	Km1255_0397   (*amigo_model).pars[619]
#define	Km0456_0397   (*amigo_model).pars[620]
#define	Km0529_0397   (*amigo_model).pars[621]
#define	Km0602_0397   (*amigo_model).pars[622]
#define	Km1207_0397   (*amigo_model).pars[623]
#define	Vmax_0398     (*amigo_model).pars[624]
#define	Keq_0398      (*amigo_model).pars[625]
#define	Km0373_0398   (*amigo_model).pars[626]
#define	Km1101_0398   (*amigo_model).pars[627]
#define	Km1212_0398   (*amigo_model).pars[628]
#define	Km0456_0398   (*amigo_model).pars[629]
#define	Km0529_0398   (*amigo_model).pars[630]
#define	Km1207_0398   (*amigo_model).pars[631]
#define	Km1255_0398   (*amigo_model).pars[632]
#define	Vmax_0399     (*amigo_model).pars[633]
#define	Keq_0399      (*amigo_model).pars[634]
#define	Km0423_0399   (*amigo_model).pars[635]
#define	Km0602_0399   (*amigo_model).pars[636]
#define	Km0633_0399   (*amigo_model).pars[637]
#define	Km0434_0399   (*amigo_model).pars[638]
#define	Km0529_0399   (*amigo_model).pars[639]
#define	Km0595_0399   (*amigo_model).pars[640]
#define	Vmax_0407     (*amigo_model).pars[641]
#define	Keq_0407      (*amigo_model).pars[642]
#define	Km0423_0407   (*amigo_model).pars[643]
#define	Km0633_0407   (*amigo_model).pars[644]
#define	Km1454_0407   (*amigo_model).pars[645]
#define	Km0434_0407   (*amigo_model).pars[646]
#define	Km0529_0407   (*amigo_model).pars[647]
#define	Km1449_0407   (*amigo_model).pars[648]
#define	Vmax_0432     (*amigo_model).pars[649]
#define	Keq_0432      (*amigo_model).pars[650]
#define	Km0602_0432   (*amigo_model).pars[651]
#define	Km1101_0432   (*amigo_model).pars[652]
#define	Km1212_0432   (*amigo_model).pars[653]
#define	Km0456_0432   (*amigo_model).pars[654]
#define	Km0529_0432   (*amigo_model).pars[655]
#define	Km1073_0432   (*amigo_model).pars[656]
#define	Km1207_0432   (*amigo_model).pars[657]
#define	Vmax_0433     (*amigo_model).pars[658]
#define	Keq_0433      (*amigo_model).pars[659]
#define	Km1073_0433   (*amigo_model).pars[660]
#define	Km1101_0433   (*amigo_model).pars[661]
#define	Km1212_0433   (*amigo_model).pars[662]
#define	Km0456_0433   (*amigo_model).pars[663]
#define	Km0529_0433   (*amigo_model).pars[664]
#define	Km1176_0433   (*amigo_model).pars[665]
#define	Km1207_0433   (*amigo_model).pars[666]
#define	Vmax_0434     (*amigo_model).pars[667]
#define	Keq_0434      (*amigo_model).pars[668]
#define	Km1101_0434   (*amigo_model).pars[669]
#define	Km1176_0434   (*amigo_model).pars[670]
#define	Km1212_0434   (*amigo_model).pars[671]
#define	Km0456_0434   (*amigo_model).pars[672]
#define	Km0529_0434   (*amigo_model).pars[673]
#define	Km1207_0434   (*amigo_model).pars[674]
#define	Km1302_0434   (*amigo_model).pars[675]
#define	Vmax_0435     (*amigo_model).pars[676]
#define	Keq_0435      (*amigo_model).pars[677]
#define	Km1101_0435   (*amigo_model).pars[678]
#define	Km1212_0435   (*amigo_model).pars[679]
#define	Km1302_0435   (*amigo_model).pars[680]
#define	Km0456_0435   (*amigo_model).pars[681]
#define	Km0529_0435   (*amigo_model).pars[682]
#define	Km1207_0435   (*amigo_model).pars[683]
#define	Km1454_0435   (*amigo_model).pars[684]
#define	Vmax_0438     (*amigo_model).pars[685]
#define	Keq_0438      (*amigo_model).pars[686]
#define	Km0710_0438   (*amigo_model).pars[687]
#define	Km1275_0438   (*amigo_model).pars[688]
#define	Km0709_0438   (*amigo_model).pars[689]
#define	Vmax_0439     (*amigo_model).pars[690]
#define	Keq_0439      (*amigo_model).pars[691]
#define	Km0709_0439   (*amigo_model).pars[692]
#define	Km1535_0439   (*amigo_model).pars[693]
#define	Km0710_0439   (*amigo_model).pars[694]
#define	Km1537_0439   (*amigo_model).pars[695]
#define	Vmax_0445     (*amigo_model).pars[696]
#define	Keq_0445      (*amigo_model).pars[697]
#define	Km0722_0445   (*amigo_model).pars[698]
#define	Km1198_0445   (*amigo_model).pars[699]
#define	Km0456_0445   (*amigo_model).pars[700]
#define	Km1203_0445   (*amigo_model).pars[701]
#define	Vmax_0446     (*amigo_model).pars[702]
#define	Keq_0446      (*amigo_model).pars[703]
#define	Km0120_0446   (*amigo_model).pars[704]
#define	Km0394_0446   (*amigo_model).pars[705]
#define	Km1322_0446   (*amigo_model).pars[706]
#define	Km0434_0446   (*amigo_model).pars[707]
#define	Km0722_0446   (*amigo_model).pars[708]
#define	Km1487_0446   (*amigo_model).pars[709]
#define	Vmax_0450     (*amigo_model).pars[710]
#define	Keq_0450      (*amigo_model).pars[711]
#define	Km0555_0450   (*amigo_model).pars[712]
#define	Km0629_0450   (*amigo_model).pars[713]
#define	Km0764_0450   (*amigo_model).pars[714]
#define	Vmax_0451     (*amigo_model).pars[715]
#define	Keq_0451      (*amigo_model).pars[716]
#define	Km0725_0451   (*amigo_model).pars[717]
#define	Km0066_0451   (*amigo_model).pars[718]
#define	Vmax_0462     (*amigo_model).pars[719]
#define	Keq_0462      (*amigo_model).pars[720]
#define	Km0745_0462   (*amigo_model).pars[721]
#define	Km0943_0462   (*amigo_model).pars[722]
#define	Km0190_0462   (*amigo_model).pars[723]
#define	Km0633_0462   (*amigo_model).pars[724]
#define	Vmax_0466     (*amigo_model).pars[725]
#define	Keq_0466      (*amigo_model).pars[726]
#define	Km0568_0466   (*amigo_model).pars[727]
#define	Km1207_0466   (*amigo_model).pars[728]
#define	Km0335_0466   (*amigo_model).pars[729]
#define	Km1212_0466   (*amigo_model).pars[730]
#define	Vmax_0467     (*amigo_model).pars[731]
#define	Keq_0467      (*amigo_model).pars[732]
#define	Km0568_0467   (*amigo_model).pars[733]
#define	Km0557_0467   (*amigo_model).pars[734]
#define	Vmax_0470     (*amigo_model).pars[735]
#define	Keq_0470      (*amigo_model).pars[736]
#define	Km0180_0470   (*amigo_model).pars[737]
#define	Km0419_0470   (*amigo_model).pars[738]
#define	Km1203_0470   (*amigo_model).pars[739]
#define	Km0991_0470   (*amigo_model).pars[740]
#define	Km1198_0470   (*amigo_model).pars[741]
#define	Vmax_0471     (*amigo_model).pars[742]
#define	Keq_0471      (*amigo_model).pars[743]
#define	Km0180_0471   (*amigo_model).pars[744]
#define	Km0419_0471   (*amigo_model).pars[745]
#define	Km1212_0471   (*amigo_model).pars[746]
#define	Km0991_0471   (*amigo_model).pars[747]
#define	Km1207_0471   (*amigo_model).pars[748]
#define	Vmax_0476     (*amigo_model).pars[749]
#define	Keq_0476      (*amigo_model).pars[750]
#define	Km0419_0476   (*amigo_model).pars[751]
#define	Km0434_0476   (*amigo_model).pars[752]
#define	Km0991_0476   (*amigo_model).pars[753]
#define	Km0394_0476   (*amigo_model).pars[754]
#define	Km0999_0476   (*amigo_model).pars[755]
#define	Km1322_0476   (*amigo_model).pars[756]
#define	Vmax_0481     (*amigo_model).pars[757]
#define	Keq_0481      (*amigo_model).pars[758]
#define	Km0754_0481   (*amigo_model).pars[759]
#define	Km1212_0481   (*amigo_model).pars[760]
#define	Km0750_0481   (*amigo_model).pars[761]
#define	Km1207_0481   (*amigo_model).pars[762]
#define	Vmax_0483     (*amigo_model).pars[763]
#define	Keq_0483      (*amigo_model).pars[764]
#define	Km0750_0483   (*amigo_model).pars[765]
#define	Km0837_0483   (*amigo_model).pars[766]
#define	Km0754_0483   (*amigo_model).pars[767]
#define	Vmax_0486     (*amigo_model).pars[768]
#define	Keq_0486      (*amigo_model).pars[769]
#define	Km0764_0486   (*amigo_model).pars[770]
#define	Km1198_0486   (*amigo_model).pars[771]
#define	Km1322_0486   (*amigo_model).pars[772]
#define	Km0075_0486   (*amigo_model).pars[773]
#define	Km1203_0486   (*amigo_model).pars[774]
#define	Vmax_0489     (*amigo_model).pars[775]
#define	Keq_0489      (*amigo_model).pars[776]
#define	Km0767_0489   (*amigo_model).pars[777]
#define	Km0765_0489   (*amigo_model).pars[778]
#define	Km1322_0489   (*amigo_model).pars[779]
#define	Vmax_0491     (*amigo_model).pars[780]
#define	Keq_0491      (*amigo_model).pars[781]
#define	Km0629_0491   (*amigo_model).pars[782]
#define	Km1203_0491   (*amigo_model).pars[783]
#define	Km0767_0491   (*amigo_model).pars[784]
#define	Km1198_0491   (*amigo_model).pars[785]
#define	Vmax_0495     (*amigo_model).pars[786]
#define	Keq_0495      (*amigo_model).pars[787]
#define	Km0380_0495   (*amigo_model).pars[788]
#define	Km0767_0495   (*amigo_model).pars[789]
#define	Km0082_0495   (*amigo_model).pars[790]
#define	Km0529_0495   (*amigo_model).pars[791]
#define	Vmax_0499     (*amigo_model).pars[792]
#define	Keq_0499      (*amigo_model).pars[793]
#define	Km0120_0499   (*amigo_model).pars[794]
#define	Km0325_0499   (*amigo_model).pars[795]
#define	Km0301_0499   (*amigo_model).pars[796]
#define	Km1487_0499   (*amigo_model).pars[797]
#define	Vmax_0502     (*amigo_model).pars[798]
#define	Keq_0502      (*amigo_model).pars[799]
#define	Km1039_0502   (*amigo_model).pars[800]
#define	Km1487_0502   (*amigo_model).pars[801]
#define	Km0306_0502   (*amigo_model).pars[802]
#define	Km1003_0502   (*amigo_model).pars[803]
#define	Vmax_0510     (*amigo_model).pars[804]
#define	Keq_0510      (*amigo_model).pars[805]
#define	Km1543_0510   (*amigo_model).pars[806]
#define	Km0773_0510   (*amigo_model).pars[807]
#define	Km1538_0510   (*amigo_model).pars[808]
#define	Vmax_0514     (*amigo_model).pars[809]
#define	Keq_0514      (*amigo_model).pars[810]
#define	Km0434_0514   (*amigo_model).pars[811]
#define	Km0999_0514   (*amigo_model).pars[812]
#define	Km1565_0514   (*amigo_model).pars[813]
#define	Km0423_0514   (*amigo_model).pars[814]
#define	Km0633_0514   (*amigo_model).pars[815]
#define	Km0782_0514   (*amigo_model).pars[816]
#define	Km0991_0514   (*amigo_model).pars[817]
#define	Vmax_0525     (*amigo_model).pars[818]
#define	Keq_0525      (*amigo_model).pars[819]
#define	Km0785_0525   (*amigo_model).pars[820]
#define	Km0141_0525   (*amigo_model).pars[821]
#define	Km0633_0525   (*amigo_model).pars[822]
#define	Km0722_0525   (*amigo_model).pars[823]
#define	Vmax_0528     (*amigo_model).pars[824]
#define	Keq_0528      (*amigo_model).pars[825]
#define	Km0434_0528   (*amigo_model).pars[826]
#define	Km0782_0528   (*amigo_model).pars[827]
#define	Km0394_0528   (*amigo_model).pars[828]
#define	Km0739_0528   (*amigo_model).pars[829]
#define	Vmax_0529     (*amigo_model).pars[830]
#define	Keq_0529      (*amigo_model).pars[831]
#define	Km0586_0529   (*amigo_model).pars[832]
#define	Km0782_0529   (*amigo_model).pars[833]
#define	Km0582_0529   (*amigo_model).pars[834]
#define	Km0739_0529   (*amigo_model).pars[835]
#define	Vmax_0534     (*amigo_model).pars[836]
#define	Keq_0534      (*amigo_model).pars[837]
#define	Km0434_0534   (*amigo_model).pars[838]
#define	Km0563_0534   (*amigo_model).pars[839]
#define	Km0394_0534   (*amigo_model).pars[840]
#define	Km0568_0534   (*amigo_model).pars[841]
#define	Vmax_0536     (*amigo_model).pars[842]
#define	Keq_0536      (*amigo_model).pars[843]
#define	Km1010_0536   (*amigo_model).pars[844]
#define	Km1198_0536   (*amigo_model).pars[845]
#define	Km1006_0536   (*amigo_model).pars[846]
#define	Km1203_0536   (*amigo_model).pars[847]
#define	Vmax_0537     (*amigo_model).pars[848]
#define	Keq_0537      (*amigo_model).pars[849]
#define	Km1011_0537   (*amigo_model).pars[850]
#define	Km1010_0537   (*amigo_model).pars[851]
#define	Km1322_0537   (*amigo_model).pars[852]
#define	Vmax_0538     (*amigo_model).pars[853]
#define	Keq_0538      (*amigo_model).pars[854]
#define	Km0207_0538   (*amigo_model).pars[855]
#define	Km0991_0538   (*amigo_model).pars[856]
#define	Km0180_0538   (*amigo_model).pars[857]
#define	Km1011_0538   (*amigo_model).pars[858]
#define	Vmax_0542     (*amigo_model).pars[859]
#define	Keq_0542      (*amigo_model).pars[860]
#define	Km0454_0542   (*amigo_model).pars[861]
#define	Km0836_0542   (*amigo_model).pars[862]
#define	Vmax_0543     (*amigo_model).pars[863]
#define	Keq_0543      (*amigo_model).pars[864]
#define	Km0180_0543   (*amigo_model).pars[865]
#define	Km0373_0543   (*amigo_model).pars[866]
#define	Km0529_0543   (*amigo_model).pars[867]
#define	Km0835_0543   (*amigo_model).pars[868]
#define	Vmax_0545     (*amigo_model).pars[869]
#define	Keq_0545      (*amigo_model).pars[870]
#define	Km0836_0545   (*amigo_model).pars[871]
#define	Km1198_0545   (*amigo_model).pars[872]
#define	Km0176_0545   (*amigo_model).pars[873]
#define	Km1203_0545   (*amigo_model).pars[874]
#define	Km0456_0545   (*amigo_model).pars[875]
#define	Vmax_0547     (*amigo_model).pars[876]
#define	Keq_0547      (*amigo_model).pars[877]
#define	Km0978_0547   (*amigo_model).pars[878]
#define	Km1212_0547   (*amigo_model).pars[879]
#define	Km1014_0547   (*amigo_model).pars[880]
#define	Km1207_0547   (*amigo_model).pars[881]
#define	Vmax_0548     (*amigo_model).pars[882]
#define	Keq_0548      (*amigo_model).pars[883]
#define	Km0434_0548   (*amigo_model).pars[884]
#define	Km1014_0548   (*amigo_model).pars[885]
#define	Km0394_0548   (*amigo_model).pars[886]
#define	Km1238_0548   (*amigo_model).pars[887]
#define	Vmax_0549     (*amigo_model).pars[888]
#define	Keq_0549      (*amigo_model).pars[889]
#define	Km0373_0549   (*amigo_model).pars[890]
#define	Km1014_0549   (*amigo_model).pars[891]
#define	Km0529_0549   (*amigo_model).pars[892]
#define	Km1233_0549   (*amigo_model).pars[893]
#define	Vmax_0550     (*amigo_model).pars[894]
#define	Keq_0550      (*amigo_model).pars[895]
#define	Km0837_0550   (*amigo_model).pars[896]
#define	Km1616_0550   (*amigo_model).pars[897]
#define	Km1620_0550   (*amigo_model).pars[898]
#define	Vmax_0553     (*amigo_model).pars[899]
#define	Keq_0553      (*amigo_model).pars[900]
#define	Km0033_0553   (*amigo_model).pars[901]
#define	Km0025_0553   (*amigo_model).pars[902]
#define	Km0750_0553   (*amigo_model).pars[903]
#define	Vmax_0558     (*amigo_model).pars[904]
#define	Keq_0558      (*amigo_model).pars[905]
#define	Km0218_0558   (*amigo_model).pars[906]
#define	Km1212_0558   (*amigo_model).pars[907]
#define	Km0028_0558   (*amigo_model).pars[908]
#define	Km0529_0558   (*amigo_model).pars[909]
#define	Km1207_0558   (*amigo_model).pars[910]
#define	Vmax_0559     (*amigo_model).pars[911]
#define	Keq_0559      (*amigo_model).pars[912]
#define	Km0367_0559   (*amigo_model).pars[913]
#define	Km0373_0559   (*amigo_model).pars[914]
#define	Km0218_0559   (*amigo_model).pars[915]
#define	Km0529_0559   (*amigo_model).pars[916]
#define	Vmax_0563     (*amigo_model).pars[917]
#define	Keq_0563      (*amigo_model).pars[918]
#define	Km0312_0563   (*amigo_model).pars[919]
#define	Km0999_0563   (*amigo_model).pars[920]
#define	Km0403_0563   (*amigo_model).pars[921]
#define	Km0550_0563   (*amigo_model).pars[922]
#define	Km0991_0563   (*amigo_model).pars[923]
#define	Vmax_0564     (*amigo_model).pars[924]
#define	Keq_0564      (*amigo_model).pars[925]
#define	Km0550_0564   (*amigo_model).pars[926]
#define	Km0207_0564   (*amigo_model).pars[927]
#define	Vmax_0565     (*amigo_model).pars[928]
#define	Keq_0565      (*amigo_model).pars[929]
#define	Km0849_0565   (*amigo_model).pars[930]
#define	Km1198_0565   (*amigo_model).pars[931]
#define	Km1203_0565   (*amigo_model).pars[932]
#define	Km1565_0565   (*amigo_model).pars[933]
#define	Vmax_0566     (*amigo_model).pars[934]
#define	Keq_0566      (*amigo_model).pars[935]
#define	Km0076_0566   (*amigo_model).pars[936]
#define	Km0086_0566   (*amigo_model).pars[937]
#define	Km0456_0566   (*amigo_model).pars[938]
#define	Vmax_0568     (*amigo_model).pars[939]
#define	Keq_0568      (*amigo_model).pars[940]
#define	Km0633_0568   (*amigo_model).pars[941]
#define	Km1322_0568   (*amigo_model).pars[942]
#define	Vmax_0570     (*amigo_model).pars[943]
#define	Keq_0570      (*amigo_model).pars[944]
#define	Km1365_0570   (*amigo_model).pars[945]
#define	Km0849_0570   (*amigo_model).pars[946]
#define	Vmax_0594     (*amigo_model).pars[947]
#define	Keq_0594      (*amigo_model).pars[948]
#define	Km0089_0594   (*amigo_model).pars[949]
#define	Km0499_0594   (*amigo_model).pars[950]
#define	Km0619_0594   (*amigo_model).pars[951]
#define	Km0918_0594   (*amigo_model).pars[952]
#define	Vmax_0658     (*amigo_model).pars[953]
#define	Keq_0658      (*amigo_model).pars[954]
#define	Km0940_0658   (*amigo_model).pars[955]
#define	Km1198_0658   (*amigo_model).pars[956]
#define	Km0180_0658   (*amigo_model).pars[957]
#define	Km0456_0658   (*amigo_model).pars[958]
#define	Km1203_0658   (*amigo_model).pars[959]
#define	Vmax_0661     (*amigo_model).pars[960]
#define	Keq_0661      (*amigo_model).pars[961]
#define	Km0940_0661   (*amigo_model).pars[962]
#define	Km1207_0661   (*amigo_model).pars[963]
#define	Km0180_0661   (*amigo_model).pars[964]
#define	Km0456_0661   (*amigo_model).pars[965]
#define	Km1212_0661   (*amigo_model).pars[966]
#define	Vmax_0663     (*amigo_model).pars[967]
#define	Keq_0663      (*amigo_model).pars[968]
#define	Km0056_0663   (*amigo_model).pars[969]
#define	Km0991_0663   (*amigo_model).pars[970]
#define	Km0180_0663   (*amigo_model).pars[971]
#define	Km1016_0663   (*amigo_model).pars[972]
#define	Vmax_0667     (*amigo_model).pars[973]
#define	Keq_0667      (*amigo_model).pars[974]
#define	Km0943_0667   (*amigo_model).pars[975]
#define	Km1376_0667   (*amigo_model).pars[976]
#define	Vmax_0669     (*amigo_model).pars[977]
#define	Keq_0669      (*amigo_model).pars[978]
#define	Km0039_0669   (*amigo_model).pars[979]
#define	Km1212_0669   (*amigo_model).pars[980]
#define	Km0008_0669   (*amigo_model).pars[981]
#define	Km1207_0669   (*amigo_model).pars[982]
#define	Vmax_0670     (*amigo_model).pars[983]
#define	Keq_0670      (*amigo_model).pars[984]
#define	Km1020_0670   (*amigo_model).pars[985]
#define	Km0427_0670   (*amigo_model).pars[986]
#define	Km0955_0670   (*amigo_model).pars[987]
#define	Vmax_0674     (*amigo_model).pars[988]
#define	Keq_0674      (*amigo_model).pars[989]
#define	Km0991_0674   (*amigo_model).pars[990]
#define	Km1399_0674   (*amigo_model).pars[991]
#define	Km0180_0674   (*amigo_model).pars[992]
#define	Km0955_0674   (*amigo_model).pars[993]
#define	Vmax_0678     (*amigo_model).pars[994]
#define	Keq_0678      (*amigo_model).pars[995]
#define	Km0953_0678   (*amigo_model).pars[996]
#define	Km1212_0678   (*amigo_model).pars[997]
#define	Km0959_0678   (*amigo_model).pars[998]
#define	Km1207_0678   (*amigo_model).pars[999]
#define	Vmax_0688     (*amigo_model).pars[1000]
#define	Keq_0688      (*amigo_model).pars[1001]
#define	Km1151_0688   (*amigo_model).pars[1002]
#define	Km1212_0688   (*amigo_model).pars[1003]
#define	Km0062_0688   (*amigo_model).pars[1004]
#define	Km1207_0688   (*amigo_model).pars[1005]
#define	Vmax_0694     (*amigo_model).pars[1006]
#define	Keq_0694      (*amigo_model).pars[1007]
#define	Km1048_0694   (*amigo_model).pars[1008]
#define	Km1275_0694   (*amigo_model).pars[1009]
#define	Km1195_0694   (*amigo_model).pars[1010]
#define	Vmax_0696     (*amigo_model).pars[1011]
#define	Keq_0696      (*amigo_model).pars[1012]
#define	Km0062_0696   (*amigo_model).pars[1013]
#define	Km1198_0696   (*amigo_model).pars[1014]
#define	Km0063_0696   (*amigo_model).pars[1015]
#define	Km1203_0696   (*amigo_model).pars[1016]
#define	Vmax_0697     (*amigo_model).pars[1017]
#define	Keq_0697      (*amigo_model).pars[1018]
#define	Km0750_0697   (*amigo_model).pars[1019]
#define	Km1151_0697   (*amigo_model).pars[1020]
#define	Km0033_0697   (*amigo_model).pars[1021]
#define	Vmax_0698     (*amigo_model).pars[1022]
#define	Keq_0698      (*amigo_model).pars[1023]
#define	Km0037_0698   (*amigo_model).pars[1024]
#define	Km1059_0698   (*amigo_model).pars[1025]
#define	Vmax_0699     (*amigo_model).pars[1026]
#define	Keq_0699      (*amigo_model).pars[1027]
#define	Km0291_0699   (*amigo_model).pars[1028]
#define	Km0991_0699   (*amigo_model).pars[1029]
#define	Km0180_0699   (*amigo_model).pars[1030]
#define	Km1021_0699   (*amigo_model).pars[1031]
#define	Vmax_0713     (*amigo_model).pars[1032]
#define	Keq_0713      (*amigo_model).pars[1033]
#define	Km0066_0713   (*amigo_model).pars[1034]
#define	Km1198_0713   (*amigo_model).pars[1035]
#define	Km1203_0713   (*amigo_model).pars[1036]
#define	Km1271_0713   (*amigo_model).pars[1037]
#define	Vmax_0722     (*amigo_model).pars[1038]
#define	Keq_0722      (*amigo_model).pars[1039]
#define	Km0573_0722   (*amigo_model).pars[1040]
#define	Km0785_0722   (*amigo_model).pars[1041]
#define	Km0633_0722   (*amigo_model).pars[1042]
#define	Km0743_0722   (*amigo_model).pars[1043]
#define	Vmax_0723     (*amigo_model).pars[1044]
#define	Keq_0723      (*amigo_model).pars[1045]
#define	Km0557_0723   (*amigo_model).pars[1046]
#define	Km0574_0723   (*amigo_model).pars[1047]
#define	Vmax_0724     (*amigo_model).pars[1048]
#define	Keq_0724      (*amigo_model).pars[1049]
#define	Km0304_0724   (*amigo_model).pars[1050]
#define	Km0120_0724   (*amigo_model).pars[1051]
#define	Vmax_0726     (*amigo_model).pars[1052]
#define	Keq_0726      (*amigo_model).pars[1053]
#define	Km0434_0726   (*amigo_model).pars[1054]
#define	Km1029_0726   (*amigo_model).pars[1055]
#define	Km0633_0726   (*amigo_model).pars[1056]
#define	Km1322_0726   (*amigo_model).pars[1057]
#define	Km1416_0726   (*amigo_model).pars[1058]
#define	Vmax_0727     (*amigo_model).pars[1059]
#define	Keq_0727      (*amigo_model).pars[1060]
#define	Km0322_0727   (*amigo_model).pars[1061]
#define	Km1012_0727   (*amigo_model).pars[1062]
#define	Km1029_0727   (*amigo_model).pars[1063]
#define	Km1487_0727   (*amigo_model).pars[1064]
#define	Vmax_0731     (*amigo_model).pars[1065]
#define	Keq_0731      (*amigo_model).pars[1066]
#define	Km0306_0731   (*amigo_model).pars[1067]
#define	Km1198_0731   (*amigo_model).pars[1068]
#define	Km0304_0731   (*amigo_model).pars[1069]
#define	Km1203_0731   (*amigo_model).pars[1070]
#define	Vmax_0732     (*amigo_model).pars[1071]
#define	Keq_0732      (*amigo_model).pars[1072]
#define	Km0306_0732   (*amigo_model).pars[1073]
#define	Km1207_0732   (*amigo_model).pars[1074]
#define	Km0304_0732   (*amigo_model).pars[1075]
#define	Km1212_0732   (*amigo_model).pars[1076]
#define	Vmax_0736     (*amigo_model).pars[1077]
#define	Keq_0736      (*amigo_model).pars[1078]
#define	Km0028_0736   (*amigo_model).pars[1079]
#define	Km0539_0736   (*amigo_model).pars[1080]
#define	Km0019_0736   (*amigo_model).pars[1081]
#define	Km0467_0736   (*amigo_model).pars[1082]
#define	Vmax_0739     (*amigo_model).pars[1083]
#define	Keq_0739      (*amigo_model).pars[1084]
#define	Km0018_0739   (*amigo_model).pars[1085]
#define	Km0434_0739   (*amigo_model).pars[1086]
#define	Km0394_0739   (*amigo_model).pars[1087]
#define	Km0456_0739   (*amigo_model).pars[1088]
#define	Km0943_0739   (*amigo_model).pars[1089]
#define	Km1322_0739   (*amigo_model).pars[1090]
#define	Vmax_0757     (*amigo_model).pars[1091]
#define	Keq_0757      (*amigo_model).pars[1092]
#define	Km0126_0757   (*amigo_model).pars[1093]
#define	Km1153_0757   (*amigo_model).pars[1094]
#define	Km1322_0757   (*amigo_model).pars[1095]
#define	Vmax_0758     (*amigo_model).pars[1096]
#define	Keq_0758      (*amigo_model).pars[1097]
#define	Km0568_0758   (*amigo_model).pars[1098]
#define	Km0126_0758   (*amigo_model).pars[1099]
#define	Vmax_0759     (*amigo_model).pars[1100]
#define	Keq_0759      (*amigo_model).pars[1101]
#define	Km1191_0759   (*amigo_model).pars[1102]
#define	Km1212_0759   (*amigo_model).pars[1103]
#define	Km0145_0759   (*amigo_model).pars[1104]
#define	Km1207_0759   (*amigo_model).pars[1105]
#define	Km1322_0759   (*amigo_model).pars[1106]
#define	Vmax_0762     (*amigo_model).pars[1107]
#define	Keq_0762      (*amigo_model).pars[1108]
#define	Km1195_0762   (*amigo_model).pars[1109]
#define	Km0722_0762   (*amigo_model).pars[1110]
#define	Km1020_0762   (*amigo_model).pars[1111]
#define	Vmax_0770     (*amigo_model).pars[1112]
#define	Keq_0770      (*amigo_model).pars[1113]
#define	Km1203_0770   (*amigo_model).pars[1114]
#define	Km1537_0770   (*amigo_model).pars[1115]
#define	Km1198_0770   (*amigo_model).pars[1116]
#define	Km1535_0770   (*amigo_model).pars[1117]
#define	Vmax_0792     (*amigo_model).pars[1118]
#define	Keq_0792      (*amigo_model).pars[1119]
#define	Km0467_0792   (*amigo_model).pars[1120]
#define	Km0526_0792   (*amigo_model).pars[1121]
#define	Km1322_0792   (*amigo_model).pars[1122]
#define	Vmax_0800     (*amigo_model).pars[1123]
#define	Keq_0800      (*amigo_model).pars[1124]
#define	Km0434_0800   (*amigo_model).pars[1125]
#define	Km0739_0800   (*amigo_model).pars[1126]
#define	Km0394_0800   (*amigo_model).pars[1127]
#define	Km0785_0800   (*amigo_model).pars[1128]
#define	Vmax_0806     (*amigo_model).pars[1129]
#define	Keq_0806      (*amigo_model).pars[1130]
#define	Km0539_0806   (*amigo_model).pars[1131]
#define	Km0467_0806   (*amigo_model).pars[1132]
#define	Km1322_0806   (*amigo_model).pars[1133]
#define	Vmax_0811     (*amigo_model).pars[1134]
#define	Keq_0811      (*amigo_model).pars[1135]
#define	Km0434_0811   (*amigo_model).pars[1136]
#define	Km1538_0811   (*amigo_model).pars[1137]
#define	Km0394_0811   (*amigo_model).pars[1138]
#define	Km1559_0811   (*amigo_model).pars[1139]
#define	Vmax_0813     (*amigo_model).pars[1140]
#define	Keq_0813      (*amigo_model).pars[1141]
#define	Km0841_0813   (*amigo_model).pars[1142]
#define	Km1233_0813   (*amigo_model).pars[1143]
#define	Km0362_0813   (*amigo_model).pars[1144]
#define	Km1012_0813   (*amigo_model).pars[1145]
#define	Vmax_0816     (*amigo_model).pars[1146]
#define	Keq_0816      (*amigo_model).pars[1147]
#define	Km0455_0816   (*amigo_model).pars[1148]
#define	Km1266_0816   (*amigo_model).pars[1149]
#define	Km0979_0816   (*amigo_model).pars[1150]
#define	Km1322_0816   (*amigo_model).pars[1151]
#define	Vmax_0818     (*amigo_model).pars[1152]
#define	Keq_0818      (*amigo_model).pars[1153]
#define	Km0991_0818   (*amigo_model).pars[1154]
#define	Km1182_0818   (*amigo_model).pars[1155]
#define	Km1192_0818   (*amigo_model).pars[1156]
#define	Km1266_0818   (*amigo_model).pars[1157]
#define	Vmax_0820     (*amigo_model).pars[1158]
#define	Keq_0820      (*amigo_model).pars[1159]
#define	Km1269_0820   (*amigo_model).pars[1160]
#define	Km1386_0820   (*amigo_model).pars[1161]
#define	Km0633_0820   (*amigo_model).pars[1162]
#define	Km1270_0820   (*amigo_model).pars[1163]
#define	Vmax_0821     (*amigo_model).pars[1164]
#define	Keq_0821      (*amigo_model).pars[1165]
#define	Km1270_0821   (*amigo_model).pars[1166]
#define	Km0456_0821   (*amigo_model).pars[1167]
#define	Km1545_0821   (*amigo_model).pars[1168]
#define	Vmax_0851     (*amigo_model).pars[1169]
#define	Keq_0851      (*amigo_model).pars[1170]
#define	Km0951_0851   (*amigo_model).pars[1171]
#define	Km0991_0851   (*amigo_model).pars[1172]
#define	Km0180_0851   (*amigo_model).pars[1173]
#define	Km1032_0851   (*amigo_model).pars[1174]
#define	Vmax_0855     (*amigo_model).pars[1175]
#define	Keq_0855      (*amigo_model).pars[1176]
#define	Km0302_0855   (*amigo_model).pars[1177]
#define	Km0434_0855   (*amigo_model).pars[1178]
#define	Km0300_0855   (*amigo_model).pars[1179]
#define	Km0394_0855   (*amigo_model).pars[1180]
#define	Km1322_0855   (*amigo_model).pars[1181]
#define	Vmax_0858     (*amigo_model).pars[1182]
#define	Keq_0858      (*amigo_model).pars[1183]
#define	Km1351_0858   (*amigo_model).pars[1184]
#define	Km1416_0858   (*amigo_model).pars[1185]
#define	Km1343_0858   (*amigo_model).pars[1186]
#define	Km1413_0858   (*amigo_model).pars[1187]
#define	Vmax_0874     (*amigo_model).pars[1188]
#define	Keq_0874      (*amigo_model).pars[1189]
#define	Km0471_0874   (*amigo_model).pars[1190]
#define	Km1153_0874   (*amigo_model).pars[1191]
#define	Km0089_0874   (*amigo_model).pars[1192]
#define	Km0526_0874   (*amigo_model).pars[1193]
#define	Vmax_0877     (*amigo_model).pars[1194]
#define	Keq_0877      (*amigo_model).pars[1195]
#define	Km1337_0877   (*amigo_model).pars[1196]
#define	Km0456_0877   (*amigo_model).pars[1197]
#define	Km1351_0877   (*amigo_model).pars[1198]
#define	Vmax_0880     (*amigo_model).pars[1199]
#define	Keq_0880      (*amigo_model).pars[1200]
#define	Km0471_0880   (*amigo_model).pars[1201]
#define	Km1039_0880   (*amigo_model).pars[1202]
#define	Km0526_0880   (*amigo_model).pars[1203]
#define	Km1337_0880   (*amigo_model).pars[1204]
#define	Vmax_0883     (*amigo_model).pars[1205]
#define	Keq_0883      (*amigo_model).pars[1206]
#define	Km0201_0883   (*amigo_model).pars[1207]
#define	Km1616_0883   (*amigo_model).pars[1208]
#define	Km0390_0883   (*amigo_model).pars[1209]
#define	Km1469_0883   (*amigo_model).pars[1210]
#define	Km1620_0883   (*amigo_model).pars[1211]
#define	Vmax_0886     (*amigo_model).pars[1212]
#define	Keq_0886      (*amigo_model).pars[1213]
#define	Km0434_0886   (*amigo_model).pars[1214]
#define	Km0557_0886   (*amigo_model).pars[1215]
#define	Km0394_0886   (*amigo_model).pars[1216]
#define	Km0555_0886   (*amigo_model).pars[1217]
#define	Vmax_0887     (*amigo_model).pars[1218]
#define	Keq_0887      (*amigo_model).pars[1219]
#define	Km0434_0887   (*amigo_model).pars[1220]
#define	Km1427_0887   (*amigo_model).pars[1221]
#define	Km0394_0887   (*amigo_model).pars[1222]
#define	Km1426_0887   (*amigo_model).pars[1223]
#define	Vmax_0888     (*amigo_model).pars[1224]
#define	Keq_0888      (*amigo_model).pars[1225]
#define	Km0568_0888   (*amigo_model).pars[1226]
#define	Km0567_0888   (*amigo_model).pars[1227]
#define	Vmax_0889     (*amigo_model).pars[1228]
#define	Keq_0889      (*amigo_model).pars[1229]
#define	Km0340_0889   (*amigo_model).pars[1230]
#define	Km1207_0889   (*amigo_model).pars[1231]
#define	Km0456_0889   (*amigo_model).pars[1232]
#define	Km0577_0889   (*amigo_model).pars[1233]
#define	Km1212_0889   (*amigo_model).pars[1234]
#define	Vmax_0891     (*amigo_model).pars[1235]
#define	Keq_0891      (*amigo_model).pars[1236]
#define	Km0260_0891   (*amigo_model).pars[1237]
#define	Km1198_0891   (*amigo_model).pars[1238]
#define	Km0258_0891   (*amigo_model).pars[1239]
#define	Km1203_0891   (*amigo_model).pars[1240]
#define	Vmax_0892     (*amigo_model).pars[1241]
#define	Keq_0892      (*amigo_model).pars[1242]
#define	Km0075_0892   (*amigo_model).pars[1243]
#define	Km0394_0892   (*amigo_model).pars[1244]
#define	Km0260_0892   (*amigo_model).pars[1245]
#define	Km0434_0892   (*amigo_model).pars[1246]
#define	Vmax_0893     (*amigo_model).pars[1247]
#define	Keq_0893      (*amigo_model).pars[1248]
#define	Km0260_0893   (*amigo_model).pars[1249]
#define	Km0188_0893   (*amigo_model).pars[1250]
#define	Vmax_0900     (*amigo_model).pars[1251]
#define	Keq_0900      (*amigo_model).pars[1252]
#define	Km1342_0900   (*amigo_model).pars[1253]
#define	Km1416_0900   (*amigo_model).pars[1254]
#define	Km1346_0900   (*amigo_model).pars[1255]
#define	Km1413_0900   (*amigo_model).pars[1256]
#define	Vmax_0901     (*amigo_model).pars[1257]
#define	Keq_0901      (*amigo_model).pars[1258]
#define	Km1343_0901   (*amigo_model).pars[1259]
#define	Km1416_0901   (*amigo_model).pars[1260]
#define	Km1342_0901   (*amigo_model).pars[1261]
#define	Km1413_0901   (*amigo_model).pars[1262]
#define	Vmax_0902     (*amigo_model).pars[1263]
#define	Keq_0902      (*amigo_model).pars[1264]
#define	Km0574_0902   (*amigo_model).pars[1265]
#define	Km0573_0902   (*amigo_model).pars[1266]
#define	Vmax_0904     (*amigo_model).pars[1267]
#define	Keq_0904      (*amigo_model).pars[1268]
#define	Km0019_0904   (*amigo_model).pars[1269]
#define	Km0434_0904   (*amigo_model).pars[1270]
#define	Km0018_0904   (*amigo_model).pars[1271]
#define	Km0394_0904   (*amigo_model).pars[1272]
#define	Vmax_0908     (*amigo_model).pars[1273]
#define	Keq_0908      (*amigo_model).pars[1274]
#define	Km0434_0908   (*amigo_model).pars[1275]
#define	Km0973_0908   (*amigo_model).pars[1276]
#define	Km1364_0908   (*amigo_model).pars[1277]
#define	Km0299_0908   (*amigo_model).pars[1278]
#define	Km0394_0908   (*amigo_model).pars[1279]
#define	Km1322_0908   (*amigo_model).pars[1280]
#define	Vmax_0909     (*amigo_model).pars[1281]
#define	Keq_0909      (*amigo_model).pars[1282]
#define	Km0078_0909   (*amigo_model).pars[1283]
#define	Km0077_0909   (*amigo_model).pars[1284]
#define	Vmax_0910     (*amigo_model).pars[1285]
#define	Keq_0910      (*amigo_model).pars[1286]
#define	Km0326_0910   (*amigo_model).pars[1287]
#define	Km0078_0910   (*amigo_model).pars[1288]
#define	Km0633_0910   (*amigo_model).pars[1289]
#define	Vmax_0911     (*amigo_model).pars[1290]
#define	Keq_0911      (*amigo_model).pars[1291]
#define	Km0300_0911   (*amigo_model).pars[1292]
#define	Km0456_0911   (*amigo_model).pars[1293]
#define	Km0434_0911   (*amigo_model).pars[1294]
#define	Km1364_0911   (*amigo_model).pars[1295]
#define	Km0394_0911   (*amigo_model).pars[1296]
#define	Km1322_0911   (*amigo_model).pars[1297]
#define	Vmax_0912     (*amigo_model).pars[1298]
#define	Keq_0912      (*amigo_model).pars[1299]
#define	Km0120_0912   (*amigo_model).pars[1300]
#define	Km0403_0912   (*amigo_model).pars[1301]
#define	Km1365_0912   (*amigo_model).pars[1302]
#define	Km1487_0912   (*amigo_model).pars[1303]
#define	Vmax_0913     (*amigo_model).pars[1304]
#define	Keq_0913      (*amigo_model).pars[1305]
#define	Km1187_0913   (*amigo_model).pars[1306]
#define	Km0076_0913   (*amigo_model).pars[1307]
#define	Vmax_0914     (*amigo_model).pars[1308]
#define	Keq_0914      (*amigo_model).pars[1309]
#define	Km0327_0914   (*amigo_model).pars[1310]
#define	Km0434_0914   (*amigo_model).pars[1311]
#define	Km1003_0914   (*amigo_model).pars[1312]
#define	Km0325_0914   (*amigo_model).pars[1313]
#define	Km0394_0914   (*amigo_model).pars[1314]
#define	Km1322_0914   (*amigo_model).pars[1315]
#define	Vmax_0915     (*amigo_model).pars[1316]
#define	Keq_0915      (*amigo_model).pars[1317]
#define	Km0999_0915   (*amigo_model).pars[1318]
#define	Km1386_0915   (*amigo_model).pars[1319]
#define	Km0327_0915   (*amigo_model).pars[1320]
#define	Km0633_0915   (*amigo_model).pars[1321]
#define	Km0991_0915   (*amigo_model).pars[1322]
#define	Vmax_0916     (*amigo_model).pars[1323]
#define	Keq_0916      (*amigo_model).pars[1324]
#define	Km0434_0916   (*amigo_model).pars[1325]
#define	Km1408_0916   (*amigo_model).pars[1326]
#define	Km0423_0916   (*amigo_model).pars[1327]
#define	Km1386_0916   (*amigo_model).pars[1328]
#define	Vmax_0917     (*amigo_model).pars[1329]
#define	Keq_0917      (*amigo_model).pars[1330]
#define	Km0259_0917   (*amigo_model).pars[1331]
#define	Km1039_0917   (*amigo_model).pars[1332]
#define	Km1322_0917   (*amigo_model).pars[1333]
#define	Vmax_0918     (*amigo_model).pars[1334]
#define	Keq_0918      (*amigo_model).pars[1335]
#define	Km0258_0918   (*amigo_model).pars[1336]
#define	Km0991_0918   (*amigo_model).pars[1337]
#define	Km0180_0918   (*amigo_model).pars[1338]
#define	Km0259_0918   (*amigo_model).pars[1339]
#define	Vmax_0919     (*amigo_model).pars[1340]
#define	Keq_0919      (*amigo_model).pars[1341]
#define	Km1084_0919   (*amigo_model).pars[1342]
#define	Km1366_0919   (*amigo_model).pars[1343]
#define	Km0481_0919   (*amigo_model).pars[1344]
#define	Vmax_0922     (*amigo_model).pars[1345]
#define	Keq_0922      (*amigo_model).pars[1346]
#define	Km1212_0922   (*amigo_model).pars[1347]
#define	Km1275_0922   (*amigo_model).pars[1348]
#define	Km1445_0922   (*amigo_model).pars[1349]
#define	Km1207_0922   (*amigo_model).pars[1350]
#define	Km1366_0922   (*amigo_model).pars[1351]
#define	Vmax_0938     (*amigo_model).pars[1352]
#define	Keq_0938      (*amigo_model).pars[1353]
#define	Km1377_0938   (*amigo_model).pars[1354]
#define	Km0456_0938   (*amigo_model).pars[1355]
#define	Km0951_0938   (*amigo_model).pars[1356]
#define	Vmax_0939     (*amigo_model).pars[1357]
#define	Keq_0939      (*amigo_model).pars[1358]
#define	Km1207_0939   (*amigo_model).pars[1359]
#define	Km1377_0939   (*amigo_model).pars[1360]
#define	Km0204_0939   (*amigo_model).pars[1361]
#define	Km0456_0939   (*amigo_model).pars[1362]
#define	Km1212_0939   (*amigo_model).pars[1363]
#define	Vmax_0957     (*amigo_model).pars[1364]
#define	Keq_0957      (*amigo_model).pars[1365]
#define	Km0118_0957   (*amigo_model).pars[1366]
#define	Km1212_0957   (*amigo_model).pars[1367]
#define	Km1035_0957   (*amigo_model).pars[1368]
#define	Km1207_0957   (*amigo_model).pars[1369]
#define	Vmax_0958     (*amigo_model).pars[1370]
#define	Keq_0958      (*amigo_model).pars[1371]
#define	Km0434_0958   (*amigo_model).pars[1372]
#define	Km0445_0958   (*amigo_model).pars[1373]
#define	Km1399_0958   (*amigo_model).pars[1374]
#define	Km0394_0958   (*amigo_model).pars[1375]
#define	Km1271_0958   (*amigo_model).pars[1376]
#define	Km1322_0958   (*amigo_model).pars[1377]
#define	Vmax_0959     (*amigo_model).pars[1378]
#define	Keq_0959      (*amigo_model).pars[1379]
#define	Km1399_0959   (*amigo_model).pars[1380]
#define	Km0359_0959   (*amigo_model).pars[1381]
#define	Km0456_0959   (*amigo_model).pars[1382]
#define	Vmax_0961     (*amigo_model).pars[1383]
#define	Keq_0961      (*amigo_model).pars[1384]
#define	Km0529_0961   (*amigo_model).pars[1385]
#define	Km1198_0961   (*amigo_model).pars[1386]
#define	Km1399_0961   (*amigo_model).pars[1387]
#define	Km0373_0961   (*amigo_model).pars[1388]
#define	Km0456_0961   (*amigo_model).pars[1389]
#define	Km1203_0961   (*amigo_model).pars[1390]
#define	Vmax_0962     (*amigo_model).pars[1391]
#define	Keq_0962      (*amigo_model).pars[1392]
#define	Km0394_0962   (*amigo_model).pars[1393]
#define	Km1360_0962   (*amigo_model).pars[1394]
#define	Km0434_0962   (*amigo_model).pars[1395]
#define	Km1399_0962   (*amigo_model).pars[1396]
#define	Vmax_0967     (*amigo_model).pars[1397]
#define	Keq_0967      (*amigo_model).pars[1398]
#define	Km0158_0967   (*amigo_model).pars[1399]
#define	Km0314_0967   (*amigo_model).pars[1400]
#define	Km0328_0967   (*amigo_model).pars[1401]
#define	Km1322_0967   (*amigo_model).pars[1402]
#define	Vmax_0968     (*amigo_model).pars[1403]
#define	Keq_0968      (*amigo_model).pars[1404]
#define	Km0328_0968   (*amigo_model).pars[1405]
#define	Km0314_0968   (*amigo_model).pars[1406]
#define	Km1405_0968   (*amigo_model).pars[1407]
#define	Vmax_0970     (*amigo_model).pars[1408]
#define	Keq_0970      (*amigo_model).pars[1409]
#define	Km0434_0970   (*amigo_model).pars[1410]
#define	Km1616_0970   (*amigo_model).pars[1411]
#define	Km0586_0970   (*amigo_model).pars[1412]
#define	Km1620_0970   (*amigo_model).pars[1413]
#define	Vmax_0973     (*amigo_model).pars[1414]
#define	Keq_0973      (*amigo_model).pars[1415]
#define	Km1559_0973   (*amigo_model).pars[1416]
#define	Km1616_0973   (*amigo_model).pars[1417]
#define	Km0656_0973   (*amigo_model).pars[1418]
#define	Km1620_0973   (*amigo_model).pars[1419]
#define	Vmax_0974     (*amigo_model).pars[1420]
#define	Keq_0974      (*amigo_model).pars[1421]
#define	Km0394_0974   (*amigo_model).pars[1422]
#define	Km1616_0974   (*amigo_model).pars[1423]
#define	Km0582_0974   (*amigo_model).pars[1424]
#define	Km1620_0974   (*amigo_model).pars[1425]
#define	Vmax_0976     (*amigo_model).pars[1426]
#define	Keq_0976      (*amigo_model).pars[1427]
#define	Km0467_0976   (*amigo_model).pars[1428]
#define	Km1616_0976   (*amigo_model).pars[1429]
#define	Km0587_0976   (*amigo_model).pars[1430]
#define	Km1620_0976   (*amigo_model).pars[1431]
#define	Vmax_0978     (*amigo_model).pars[1432]
#define	Keq_0978      (*amigo_model).pars[1433]
#define	Km0739_0978   (*amigo_model).pars[1434]
#define	Km1616_0978   (*amigo_model).pars[1435]
#define	Km0613_0978   (*amigo_model).pars[1436]
#define	Km1620_0978   (*amigo_model).pars[1437]
#define	Vmax_0982     (*amigo_model).pars[1438]
#define	Keq_0982      (*amigo_model).pars[1439]
#define	Km0577_0982   (*amigo_model).pars[1440]
#define	Km1408_0982   (*amigo_model).pars[1441]
#define	Vmax_0984     (*amigo_model).pars[1442]
#define	Keq_0984      (*amigo_model).pars[1443]
#define	Km0577_0984   (*amigo_model).pars[1444]
#define	Km0581_0984   (*amigo_model).pars[1445]
#define	Vmax_0986     (*amigo_model).pars[1446]
#define	Keq_0986      (*amigo_model).pars[1447]
#define	Km1416_0986   (*amigo_model).pars[1448]
#define	Km1569_0986   (*amigo_model).pars[1449]
#define	Km0700_0986   (*amigo_model).pars[1450]
#define	Km1413_0986   (*amigo_model).pars[1451]
#define	Vmax_0988     (*amigo_model).pars[1452]
#define	Keq_0988      (*amigo_model).pars[1453]
#define	Km1038_0988   (*amigo_model).pars[1454]
#define	Km1198_0988   (*amigo_model).pars[1455]
#define	Km0180_0988   (*amigo_model).pars[1456]
#define	Km1025_0988   (*amigo_model).pars[1457]
#define	Km1203_0988   (*amigo_model).pars[1458]
#define	Vmax_0989     (*amigo_model).pars[1459]
#define	Keq_0989      (*amigo_model).pars[1460]
#define	Km0959_0989   (*amigo_model).pars[1461]
#define	Km0991_0989   (*amigo_model).pars[1462]
#define	Km1212_0989   (*amigo_model).pars[1463]
#define	Km1038_0989   (*amigo_model).pars[1464]
#define	Km1207_0989   (*amigo_model).pars[1465]
#define	Vmax_0990     (*amigo_model).pars[1466]
#define	Keq_0990      (*amigo_model).pars[1467]
#define	Km1426_0990   (*amigo_model).pars[1468]
#define	Km0551_0990   (*amigo_model).pars[1469]
#define	Km0629_0990   (*amigo_model).pars[1470]
#define	Vmax_0992     (*amigo_model).pars[1471]
#define	Keq_0992      (*amigo_model).pars[1472]
#define	Km0373_0992   (*amigo_model).pars[1473]
#define	Km1039_0992   (*amigo_model).pars[1474]
#define	Km0529_0992   (*amigo_model).pars[1475]
#define	Km1234_0992   (*amigo_model).pars[1476]
#define	Vmax_0993     (*amigo_model).pars[1477]
#define	Keq_0993      (*amigo_model).pars[1478]
#define	Km1039_0993   (*amigo_model).pars[1479]
#define	Km1302_0993   (*amigo_model).pars[1480]
#define	Km0231_0993   (*amigo_model).pars[1481]
#define	Km0456_0993   (*amigo_model).pars[1482]
#define	Km0529_0993   (*amigo_model).pars[1483]
#define	Vmax_0996     (*amigo_model).pars[1484]
#define	Keq_0996      (*amigo_model).pars[1485]
#define	Km0211_0996   (*amigo_model).pars[1486]
#define	Km1212_0996   (*amigo_model).pars[1487]
#define	Km1207_0996   (*amigo_model).pars[1488]
#define	Km1429_0996   (*amigo_model).pars[1489]
#define	Vmax_0997     (*amigo_model).pars[1490]
#define	Keq_0997      (*amigo_model).pars[1491]
#define	Km0434_0997   (*amigo_model).pars[1492]
#define	Km1429_0997   (*amigo_model).pars[1493]
#define	Km0261_0997   (*amigo_model).pars[1494]
#define	Km0394_0997   (*amigo_model).pars[1495]
#define	Vmax_1010     (*amigo_model).pars[1496]
#define	Keq_1010      (*amigo_model).pars[1497]
#define	Km1203_1010   (*amigo_model).pars[1498]
#define	Km1275_1010   (*amigo_model).pars[1499]
#define	Km1447_1010   (*amigo_model).pars[1500]
#define	Km0037_1010   (*amigo_model).pars[1501]
#define	Km1198_1010   (*amigo_model).pars[1502]
#define	Vmax_1011     (*amigo_model).pars[1503]
#define	Keq_1011      (*amigo_model).pars[1504]
#define	Km1212_1011   (*amigo_model).pars[1505]
#define	Km1275_1011   (*amigo_model).pars[1506]
#define	Km1447_1011   (*amigo_model).pars[1507]
#define	Km0037_1011   (*amigo_model).pars[1508]
#define	Km1207_1011   (*amigo_model).pars[1509]
#define	Vmax_1012     (*amigo_model).pars[1510]
#define	Keq_1012      (*amigo_model).pars[1511]
#define	Km0190_1012   (*amigo_model).pars[1512]
#define	Km1212_1012   (*amigo_model).pars[1513]
#define	Km0633_1012   (*amigo_model).pars[1514]
#define	Km1207_1012   (*amigo_model).pars[1515]
#define	Km1447_1012   (*amigo_model).pars[1516]
#define	Vmax_1014     (*amigo_model).pars[1517]
#define	Keq_1014      (*amigo_model).pars[1518]
#define	Km0666_1014   (*amigo_model).pars[1519]
#define	Km0595_1014   (*amigo_model).pars[1520]
#define	Km0672_1014   (*amigo_model).pars[1521]
#define	Vmax_1026     (*amigo_model).pars[1522]
#define	Keq_1026      (*amigo_model).pars[1523]
#define	Km0394_1026   (*amigo_model).pars[1524]
#define	Km1467_1026   (*amigo_model).pars[1525]
#define	Km0298_1026   (*amigo_model).pars[1526]
#define	Km1322_1026   (*amigo_model).pars[1527]
#define	Vmax_1027     (*amigo_model).pars[1528]
#define	Keq_1027      (*amigo_model).pars[1529]
#define	Km1212_1027   (*amigo_model).pars[1530]
#define	Km1469_1027   (*amigo_model).pars[1531]
#define	Km0841_1027   (*amigo_model).pars[1532]
#define	Km1207_1027   (*amigo_model).pars[1533]
#define	Vmax_1038     (*amigo_model).pars[1534]
#define	Keq_1038      (*amigo_model).pars[1535]
#define	Km1212_1038   (*amigo_model).pars[1536]
#define	Km1620_1038   (*amigo_model).pars[1537]
#define	Km1207_1038   (*amigo_model).pars[1538]
#define	Km1616_1038   (*amigo_model).pars[1539]
#define	Vmax_1041     (*amigo_model).pars[1540]
#define	Keq_1041      (*amigo_model).pars[1541]
#define	Km1238_1041   (*amigo_model).pars[1542]
#define	Km1045_1041   (*amigo_model).pars[1543]
#define	Km1322_1041   (*amigo_model).pars[1544]
#define	Vmax_1045     (*amigo_model).pars[1545]
#define	Keq_1045      (*amigo_model).pars[1546]
#define	Km0306_1045   (*amigo_model).pars[1547]
#define	Km0654_1045   (*amigo_model).pars[1548]
#define	Km0625_1045   (*amigo_model).pars[1549]
#define	Km0649_1045   (*amigo_model).pars[1550]
#define	Vmax_1049     (*amigo_model).pars[1551]
#define	Keq_1049      (*amigo_model).pars[1552]
#define	Km0581_1049   (*amigo_model).pars[1553]
#define	Km1408_1049   (*amigo_model).pars[1554]
#define	Km0764_1049   (*amigo_model).pars[1555]
#define	Km1427_1049   (*amigo_model).pars[1556]
#define	Vmax_1050     (*amigo_model).pars[1557]
#define	Keq_1050      (*amigo_model).pars[1558]
#define	Km0551_1050   (*amigo_model).pars[1559]
#define	Km0581_1050   (*amigo_model).pars[1560]
#define	Km0557_1050   (*amigo_model).pars[1561]
#define	Km0764_1050   (*amigo_model).pars[1562]
#define	Vmax_1051     (*amigo_model).pars[1563]
#define	Keq_1051      (*amigo_model).pars[1564]
#define	Km0409_1051   (*amigo_model).pars[1565]
#define	Km1322_1051   (*amigo_model).pars[1566]
#define	Km1520_1051   (*amigo_model).pars[1567]
#define	Vmax_1052     (*amigo_model).pars[1568]
#define	Keq_1052      (*amigo_model).pars[1569]
#define	Km0619_1052   (*amigo_model).pars[1570]
#define	Km0595_1052   (*amigo_model).pars[1571]
#define	Km1524_1052   (*amigo_model).pars[1572]
#define	Vmax_1054     (*amigo_model).pars[1573]
#define	Keq_1054      (*amigo_model).pars[1574]
#define	Km0764_1054   (*amigo_model).pars[1575]
#define	Km0629_1054   (*amigo_model).pars[1576]
#define	Vmax_1055     (*amigo_model).pars[1577]
#define	Keq_1055      (*amigo_model).pars[1578]
#define	Km0086_1055   (*amigo_model).pars[1579]
#define	Km1039_1055   (*amigo_model).pars[1580]
#define	Km0764_1055   (*amigo_model).pars[1581]
#define	Km1048_1055   (*amigo_model).pars[1582]
#define	Vmax_1063     (*amigo_model).pars[1583]
#define	Keq_1063      (*amigo_model).pars[1584]
#define	Km0204_1063   (*amigo_model).pars[1585]
#define	Km0991_1063   (*amigo_model).pars[1586]
#define	Km0180_1063   (*amigo_model).pars[1587]
#define	Km1051_1063   (*amigo_model).pars[1588]
#define	Vmax_1072     (*amigo_model).pars[1589]
#define	Keq_1072      (*amigo_model).pars[1590]
#define	Km0434_1072   (*amigo_model).pars[1591]
#define	Km1545_1072   (*amigo_model).pars[1592]
#define	Km0394_1072   (*amigo_model).pars[1593]
#define	Km1538_1072   (*amigo_model).pars[1594]
#define	Vmax_1084     (*amigo_model).pars[1595]
#define	Keq_1084      (*amigo_model).pars[1596]
#define	Km0567_1084   (*amigo_model).pars[1597]
#define	Km1559_1084   (*amigo_model).pars[1598]
#define	Km0633_1084   (*amigo_model).pars[1599]
#define	Km1543_1084   (*amigo_model).pars[1600]
#define	Vmax_1087     (*amigo_model).pars[1601]
#define	Keq_1087      (*amigo_model).pars[1602]
#define	Km0232_1087   (*amigo_model).pars[1603]
#define	Km0991_1087   (*amigo_model).pars[1604]
#define	Km0180_1087   (*amigo_model).pars[1605]
#define	Km1056_1087   (*amigo_model).pars[1606]
#define	Vmax_1106     (*amigo_model).pars[1607]
#define	Km0362_1106   (*amigo_model).pars[1608]
#define	Vmax_1115     (*amigo_model).pars[1609]
#define	Km0420_1115   (*amigo_model).pars[1610]
#define	Km0419_1115   (*amigo_model).pars[1611]
#define	Vmax_1166     (*amigo_model).pars[1612]
#define	Km0565_1166   (*amigo_model).pars[1613]
#define	Km0563_1166   (*amigo_model).pars[1614]
#define	Vmax_1172     (*amigo_model).pars[1615]
#define	Km0765_1172   (*amigo_model).pars[1616]
#define	Vmax_1244     (*amigo_model).pars[1617]
#define	Km1324_1244   (*amigo_model).pars[1618]
#define	Km1322_1244   (*amigo_model).pars[1619]
#define	Vmax_1266     (*amigo_model).pars[1620]
#define	Km1468_1266   (*amigo_model).pars[1621]
#define	Km1467_1266   (*amigo_model).pars[1622]
#define	Vmax_1664     (*amigo_model).pars[1623]
#define	Keq_1664      (*amigo_model).pars[1624]
#define	Km0456_1664   (*amigo_model).pars[1625]
#define	Km0445_1664   (*amigo_model).pars[1626]
#define	Vmax_1697     (*amigo_model).pars[1627]
#define	Km0456_1697   (*amigo_model).pars[1628]
#define	Vmax_1704     (*amigo_model).pars[1629]
#define	Keq_1704      (*amigo_model).pars[1630]
#define	Km0394_1704   (*amigo_model).pars[1631]
#define	Km0587_1704   (*amigo_model).pars[1632]
#define	Km0434_1704   (*amigo_model).pars[1633]
#define	Km0589_1704   (*amigo_model).pars[1634]
#define	Vmax_1729     (*amigo_model).pars[1635]
#define	Keq_1729      (*amigo_model).pars[1636]
#define	Km0394_1729   (*amigo_model).pars[1637]
#define	Km0582_1729   (*amigo_model).pars[1638]
#define	Km0434_1729   (*amigo_model).pars[1639]
#define	Km0584_1729   (*amigo_model).pars[1640]
#define	Vmax_1762     (*amigo_model).pars[1641]
#define	Km0680_1762   (*amigo_model).pars[1642]
#define	Vmax_1936     (*amigo_model).pars[1643]
#define	Keq_1936      (*amigo_model).pars[1644]
#define	Km0629_1936   (*amigo_model).pars[1645]
#define	Km1151_1936   (*amigo_model).pars[1646]
#define	Km1322_1936   (*amigo_model).pars[1647]
#define	Vmax_1979     (*amigo_model).pars[1648]
#define	Km1277_1979   (*amigo_model).pars[1649]
#define	Km1275_1979   (*amigo_model).pars[1650]
#define	Vmax_2030     (*amigo_model).pars[1651]
#define	Keq_2030      (*amigo_model).pars[1652]
#define	Km0313_2030   (*amigo_model).pars[1653]
#define	Km0314_2030   (*amigo_model).pars[1654]
#define	Km1322_2030   (*amigo_model).pars[1655]
#define	Vmax_2079     (*amigo_model).pars[1656]
#define	Km1520_2079   (*amigo_model).pars[1657]
#define	V0_2111       (*amigo_model).pars[1658]
#define	ic0002_2111   (*amigo_model).pars[1659]
#define	ep0002_2111   (*amigo_model).pars[1660]
#define	ic0423_2111   (*amigo_model).pars[1661]
#define	ep0423_2111   (*amigo_model).pars[1662]
#define	ic0434_2111   (*amigo_model).pars[1663]
#define	ep0434_2111   (*amigo_model).pars[1664]
#define	ic0526_2111   (*amigo_model).pars[1665]
#define	ep0526_2111   (*amigo_model).pars[1666]
#define	ic0584_2111   (*amigo_model).pars[1667]
#define	ep0584_2111   (*amigo_model).pars[1668]
#define	ic0589_2111   (*amigo_model).pars[1669]
#define	ep0589_2111   (*amigo_model).pars[1670]
#define	ic0615_2111   (*amigo_model).pars[1671]
#define	ep0615_2111   (*amigo_model).pars[1672]
#define	ic0649_2111   (*amigo_model).pars[1673]
#define	ep0649_2111   (*amigo_model).pars[1674]
#define	ic0773_2111   (*amigo_model).pars[1675]
#define	ep0773_2111   (*amigo_model).pars[1676]
#define	ic0782_2111   (*amigo_model).pars[1677]
#define	ep0782_2111   (*amigo_model).pars[1678]
#define	ic0955_2111   (*amigo_model).pars[1679]
#define	ep0955_2111   (*amigo_model).pars[1680]
#define	ic0965_2111   (*amigo_model).pars[1681]
#define	ep0965_2111   (*amigo_model).pars[1682]
#define	ic0969_2111   (*amigo_model).pars[1683]
#define	ep0969_2111   (*amigo_model).pars[1684]
#define	ic0973_2111   (*amigo_model).pars[1685]
#define	ep0973_2111   (*amigo_model).pars[1686]
#define	ic0981_2111   (*amigo_model).pars[1687]
#define	ep0981_2111   (*amigo_model).pars[1688]
#define	ic0991_2111   (*amigo_model).pars[1689]
#define	ep0991_2111   (*amigo_model).pars[1690]
#define	ic0999_2111   (*amigo_model).pars[1691]
#define	ep0999_2111   (*amigo_model).pars[1692]
#define	ic1003_2111   (*amigo_model).pars[1693]
#define	ep1003_2111   (*amigo_model).pars[1694]
#define	ic1006_2111   (*amigo_model).pars[1695]
#define	ep1006_2111   (*amigo_model).pars[1696]
#define	ic1016_2111   (*amigo_model).pars[1697]
#define	ep1016_2111   (*amigo_model).pars[1698]
#define	ic1021_2111   (*amigo_model).pars[1699]
#define	ep1021_2111   (*amigo_model).pars[1700]
#define	ic1025_2111   (*amigo_model).pars[1701]
#define	ep1025_2111   (*amigo_model).pars[1702]
#define	ic1029_2111   (*amigo_model).pars[1703]
#define	ep1029_2111   (*amigo_model).pars[1704]
#define	ic1032_2111   (*amigo_model).pars[1705]
#define	ep1032_2111   (*amigo_model).pars[1706]
#define	ic1035_2111   (*amigo_model).pars[1707]
#define	ep1035_2111   (*amigo_model).pars[1708]
#define	ic1039_2111   (*amigo_model).pars[1709]
#define	ep1039_2111   (*amigo_model).pars[1710]
#define	ic1045_2111   (*amigo_model).pars[1711]
#define	ep1045_2111   (*amigo_model).pars[1712]
#define	ic1048_2111   (*amigo_model).pars[1713]
#define	ep1048_2111   (*amigo_model).pars[1714]
#define	ic1051_2111   (*amigo_model).pars[1715]
#define	ep1051_2111   (*amigo_model).pars[1716]
#define	ic1056_2111   (*amigo_model).pars[1717]
#define	ep1056_2111   (*amigo_model).pars[1718]
#define	ic1107_2111   (*amigo_model).pars[1719]
#define	ep1107_2111   (*amigo_model).pars[1720]
#define	ic1405_2111   (*amigo_model).pars[1721]
#define	ep1405_2111   (*amigo_model).pars[1722]
#define	ic1467_2111   (*amigo_model).pars[1723]
#define	ep1467_2111   (*amigo_model).pars[1724]
#define	ic1520_2111   (*amigo_model).pars[1725]
#define	ep1520_2111   (*amigo_model).pars[1726]
#define	ic1545_2111   (*amigo_model).pars[1727]
#define	ep1545_2111   (*amigo_model).pars[1728]
#define	ic0089_2111   (*amigo_model).pars[1729]
#define	ep0089_2111   (*amigo_model).pars[1730]
#define	ic0122_2111   (*amigo_model).pars[1731]
#define	ep0122_2111   (*amigo_model).pars[1732]
#define	ic0918_2111   (*amigo_model).pars[1733]
#define	ep0918_2111   (*amigo_model).pars[1734]
#define	ic0657_2111   (*amigo_model).pars[1735]
#define	ep0657_2111   (*amigo_model).pars[1736]
#define	ic0662_2111   (*amigo_model).pars[1737]
#define	ep0662_2111   (*amigo_model).pars[1738]
#define	ic0666_2111   (*amigo_model).pars[1739]
#define	ep0666_2111   (*amigo_model).pars[1740]
#define	ic0672_2111   (*amigo_model).pars[1741]
#define	ep0672_2111   (*amigo_model).pars[1742]
#define	ic0595_2111   (*amigo_model).pars[1743]
#define	ep0595_2111   (*amigo_model).pars[1744]
#define	ic0700_2111   (*amigo_model).pars[1745]
#define	ep0700_2111   (*amigo_model).pars[1746]
#define	ic1059_2111   (*amigo_model).pars[1747]
#define	ep1059_2111   (*amigo_model).pars[1748]
#define	ic1337_2111   (*amigo_model).pars[1749]
#define	ep1337_2111   (*amigo_model).pars[1750]
#define	ic1346_2111   (*amigo_model).pars[1751]
#define	ep1346_2111   (*amigo_model).pars[1752]
#define	ic1351_2111   (*amigo_model).pars[1753]
#define	ep1351_2111   (*amigo_model).pars[1754]
#define	ic1524_2111   (*amigo_model).pars[1755]
#define	ep1524_2111   (*amigo_model).pars[1756]
#define	ic1569_2111   (*amigo_model).pars[1757]
#define	ep1569_2111   (*amigo_model).pars[1758]
#define	s_0364        (*amigo_model).pars[1759]
#define	s_0420        (*amigo_model).pars[1760]
#define	s_0458        (*amigo_model).pars[1761]
#define	s_0681        (*amigo_model).pars[1762]
#define	s_0766        (*amigo_model).pars[1763]
#define	s_1277        (*amigo_model).pars[1764]
#define	s_1324        (*amigo_model).pars[1765]
#define	s_1468        (*amigo_model).pars[1766]
#define	s_1521        (*amigo_model).pars[1767]
#define	cell          (*amigo_model).pars[1768]
#define	extracellular (*amigo_model).pars[1769]
#define s_0565	((*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim])

/* Right hand side of the system (f(t,x,p))*/
int amigoRHS_B1(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;
	/* *** Definition of the algebraic variables *** */

	double	r_0001;
	double	r_0004;
	double	r_0005;
	double	r_0007;
	double	r_0008;
	double	r_0012;
	double	r_0014;
	double	r_0015;
	double	r_0016;
	double	r_0018;
	double	r_0020;
	double	r_0023;
	double	r_0024;
	double	r_0027;
	double	r_0029;
	double	r_0032;
	double	r_0038;
	double	r_0039;
	double	r_0040;
	double	r_0041;
	double	r_0060;
	double	r_0061;
	double	r_0065;
	double	r_0079;
	double	r_0080;
	double	r_0091;
	double	r_0096;
	double	r_0097;
	double	r_0103;
	double	r_0108;
	double	r_0111;
	double	r_0115;
	double	r_0118;
	double	r_0142;
	double	r_0144;
	double	r_0148;
	double	r_0151;
	double	r_0152;
	double	r_0153;
	double	r_0154;
	double	r_0165;
	double	r_0173;
	double	r_0195;
	double	r_0202;
	double	r_0203;
	double	r_0207;
	double	r_0208;
	double	r_0211;
	double	r_0214;
	double	r_0215;
	double	r_0216;
	double	r_0219;
	double	r_0225;
	double	r_0226;
	double	r_0231;
	double	r_0233;
	double	r_0234;
	double	r_0235;
	double	r_0236;
	double	r_0237;
	double	r_0238;
	double	r_0239;
	double	r_0240;
	double	r_0241;
	double	r_0242;
	double	r_0243;
	double	r_0244;
	double	r_0250;
	double	r_0257;
	double	r_0259;
	double	r_0267;
	double	r_0269;
	double	r_0278;
	double	r_0279;
	double	r_0280;
	double	r_0300;
	double	r_0302;
	double	r_0307;
	double	r_0309;
	double	r_0310;
	double	r_0311;
	double	r_0312;
	double	r_0317;
	double	r_0326;
	double	r_0330;
	double	r_0336;
	double	r_0337;
	double	r_0339;
	double	r_0340;
	double	r_0344;
	double	r_0349;
	double	r_0352;
	double	r_0353;
	double	r_0355;
	double	r_0361;
	double	r_0362;
	double	r_0364;
	double	r_0366;
	double	r_0386;
	double	r_0387;
	double	r_0389;
	double	r_0391;
	double	r_0393;
	double	r_0397;
	double	r_0398;
	double	r_0399;
	double	r_0407;
	double	r_0432;
	double	r_0433;
	double	r_0434;
	double	r_0435;
	double	r_0438;
	double	r_0439;
	double	r_0445;
	double	r_0446;
	double	r_0450;
	double	r_0451;
	double	r_0462;
	double	r_0466;
	double	r_0467;
	double	r_0470;
	double	r_0471;
	double	r_0476;
	double	r_0481;
	double	r_0483;
	double	r_0486;
	double	r_0489;
	double	r_0491;
	double	r_0495;
	double	r_0499;
	double	r_0502;
	double	r_0510;
	double	r_0514;
	double	r_0525;
	double	r_0528;
	double	r_0529;
	double	r_0534;
	double	r_0536;
	double	r_0537;
	double	r_0538;
	double	r_0542;
	double	r_0543;
	double	r_0545;
	double	r_0547;
	double	r_0548;
	double	r_0549;
	double	r_0550;
	double	r_0553;
	double	r_0558;
	double	r_0559;
	double	r_0563;
	double	r_0564;
	double	r_0565;
	double	r_0566;
	double	r_0568;
	double	r_0570;
	double	r_0594;
	double	r_0658;
	double	r_0661;
	double	r_0663;
	double	r_0667;
	double	r_0669;
	double	r_0670;
	double	r_0674;
	double	r_0678;
	double	r_0688;
	double	r_0694;
	double	r_0696;
	double	r_0697;
	double	r_0698;
	double	r_0699;
	double	r_0713;
	double	r_0722;
	double	r_0723;
	double	r_0724;
	double	r_0726;
	double	r_0727;
	double	r_0731;
	double	r_0732;
	double	r_0736;
	double	r_0739;
	double	r_0757;
	double	r_0758;
	double	r_0759;
	double	r_0762;
	double	r_0770;
	double	r_0792;
	double	r_0800;
	double	r_0806;
	double	r_0811;
	double	r_0813;
	double	r_0816;
	double	r_0818;
	double	r_0820;
	double	r_0821;
	double	r_0851;
	double	r_0855;
	double	r_0858;
	double	r_0874;
	double	r_0877;
	double	r_0880;
	double	r_0883;
	double	r_0886;
	double	r_0887;
	double	r_0888;
	double	r_0889;
	double	r_0891;
	double	r_0892;
	double	r_0893;
	double	r_0900;
	double	r_0901;
	double	r_0902;
	double	r_0904;
	double	r_0908;
	double	r_0909;
	double	r_0910;
	double	r_0911;
	double	r_0912;
	double	r_0913;
	double	r_0914;
	double	r_0915;
	double	r_0916;
	double	r_0917;
	double	r_0918;
	double	r_0919;
	double	r_0922;
	double	r_0938;
	double	r_0939;
	double	r_0957;
	double	r_0958;
	double	r_0959;
	double	r_0961;
	double	r_0962;
	double	r_0967;
	double	r_0968;
	double	r_0970;
	double	r_0973;
	double	r_0974;
	double	r_0976;
	double	r_0978;
	double	r_0982;
	double	r_0984;
	double	r_0986;
	double	r_0988;
	double	r_0989;
	double	r_0990;
	double	r_0992;
	double	r_0993;
	double	r_0996;
	double	r_0997;
	double	r_1010;
	double	r_1011;
	double	r_1012;
	double	r_1014;
	double	r_1026;
	double	r_1027;
	double	r_1038;
	double	r_1041;
	double	r_1045;
	double	r_1049;
	double	r_1050;
	double	r_1051;
	double	r_1052;
	double	r_1054;
	double	r_1055;
	double	r_1063;
	double	r_1072;
	double	r_1084;
	double	r_1087;
	double	r_1106;
	double	r_1115;
	double	r_1166;
	double	r_1172;
	double	r_1244;
	double	r_1266;
	double	r_1664;
	double	r_1697;
	double	r_1704;
	double	r_1729;
	double	r_1762;
	double	r_1936;
	double	r_1979;
	double	r_2030;
	double	r_2079;
	double	r_2111;

	/* *** Equations *** */

	r_0001=cell*Vmax_0001*(s_0025*pow(s_0709,2)-pow(s_0710,2)*s_1399/Keq_0001)/(Km0025_0001*pow(Km0709_0001,2))/((1+s_0025/Km0025_0001)*pow(1+s_0709/Km0709_0001,2)+pow(1+s_0710/Km0710_0001,2)*(1+s_1399/Km1399_0001)-1);
	r_0004=cell*Vmax_0004*(s_0063*pow(s_0709,2)-pow(s_0710,2)*s_1399/Keq_0004)/(Km0063_0004*pow(Km0709_0004,2))/((1+s_0063/Km0063_0004)*pow(1+s_0709/Km0709_0004,2)+pow(1+s_0710/Km0710_0004,2)*(1+s_1399/Km1399_0004)-1);
	r_0005=cell*Vmax_0005*(s_1543-s_0002*s_1538/Keq_0005)/Km1543_0005/(1+s_1543/Km1543_0005+(1+s_0002/Km0002_0005)*(1+s_1538/Km1538_0005)-1);
	r_0007=cell*Vmax_0007*(s_0077-s_0312/Keq_0007)/Km0077_0007/(1+s_0077/Km0077_0007+1+s_0312/Km0312_0007-1);
	r_0008=cell*Vmax_0008*(s_0082*s_0380-s_0529*s_1331/Keq_0008)/(Km0082_0008*Km0380_0008)/((1+s_0082/Km0082_0008)*(1+s_0380/Km0380_0008)+(1+s_0529/Km0529_0008)*(1+s_1331/Km1331_0008)-1);
	r_0012=cell*Vmax_0012*(s_0991*s_1203-s_0118*s_1198/Keq_0012)/(Km0991_0012*Km1203_0012)/((1+s_0991/Km0991_0012)*(1+s_1203/Km1203_0012)+(1+s_0118/Km0118_0012)*(1+s_1198/Km1198_0012)-1);
	r_0014=cell*Vmax_0014*(s_0142-s_0313*s_0419/Keq_0014)/Km0142_0014/(1+s_0142/Km0142_0014+(1+s_0313/Km0313_0014)*(1+s_0419/Km0419_0014)-1);
	r_0015=cell*Vmax_0015*(s_0141*s_1212-s_0142*s_1207/Keq_0015)/(Km0141_0015*Km1212_0015)/((1+s_0141/Km0141_0015)*(1+s_1212/Km1212_0015)+(1+s_0142/Km0142_0015)*(1+s_1207/Km1207_0015)-1);
	r_0016=cell*Vmax_0016*(s_0178*s_1399-s_0039*s_0456/Keq_0016)/(Km0178_0016*Km1399_0016)/((1+s_0178/Km0178_0016)*(1+s_1399/Km1399_0016)+(1+s_0039/Km0039_0016)*(1+s_0456/Km0456_0016)-1);
	r_0018=cell*Vmax_0018*(s_0176*s_0991-s_0180*s_0953/Keq_0018)/(Km0176_0018*Km0991_0018)/((1+s_0176/Km0176_0018)*(1+s_0991/Km0991_0018)+(1+s_0180/Km0180_0018)*(1+s_0953/Km0953_0018)-1);
	r_0020=cell*Vmax_0020*(s_0551*s_1360-s_0349*s_1322/Keq_0020)/(Km0551_0020*Km1360_0020)/((1+s_0551/Km0551_0020)*(1+s_1360/Km1360_0020)+(1+s_0349/Km0349_0020)*(1+s_1322/Km1322_0020)-1);
	r_0023=cell*Vmax_0023*(s_0162-s_0165/Keq_0023)/Km0162_0023/(1+s_0162/Km0162_0023+1+s_0165/Km0165_0023-1);
	r_0024=cell*Vmax_0024*(s_0232*s_0373-s_0162*s_0529/Keq_0024)/(Km0232_0024*Km0373_0024)/((1+s_0232/Km0232_0024)*(1+s_0373/Km0373_0024)+(1+s_0162/Km0162_0024)*(1+s_0529/Km0529_0024)-1);
	r_0027=cell*Vmax_0027*(s_0835-s_0454/Keq_0027)/Km0835_0027/(1+s_0835/Km0835_0027+1+s_0454/Km0454_0027-1);
	r_0029=cell*Vmax_0029*(s_0010-s_0291*s_0456/Keq_0029)/Km0010_0029/(1+s_0010/Km0010_0029+(1+s_0291/Km0291_0029)*(1+s_0456/Km0456_0029)-1);
	r_0032=cell*Vmax_0032*(s_0390-s_0423*s_1322/Keq_0032)/Km0390_0032/(1+s_0390/Km0390_0032+(1+s_0423/Km0423_0032)*(1+s_1322/Km1322_0032)-1);
	r_0038=cell*Vmax_0038*(s_0577-s_0158*s_0722/Keq_0038)/Km0577_0038/(1+s_0577/Km0577_0038+(1+s_0158/Km0158_0038)*(1+s_0722/Km0722_0038)-1);
	r_0039=cell*Vmax_0039*(s_0210-s_0211/Keq_0039)/Km0210_0039/(1+s_0210/Km0210_0039+1+s_0211/Km0211_0039-1);
	r_0040=cell*Vmax_0040*(s_0349-s_0210*s_1322/Keq_0040)/Km0349_0040/(1+s_0349/Km0349_0040+(1+s_0210/Km0210_0040)*(1+s_1322/Km1322_0040)-1);
	r_0041=cell*Vmax_0041*(s_0231*s_1212-s_1207*s_1445/Keq_0041)/(Km0231_0041*Km1212_0041)/((1+s_0231/Km0231_0041)*(1+s_1212/Km1212_0041)+(1+s_1207/Km1207_0041)*(1+s_1445/Km1445_0041)-1);
	r_0060=cell*Vmax_0060*(s_0165-s_0009/Keq_0060)/Km0165_0060/(1+s_0165/Km0165_0060+1+s_0009/Km0009_0060-1);
	r_0061=cell*Vmax_0061*(s_0009*s_1198-s_0010*s_1203/Keq_0061)/(Km0009_0061*Km1198_0061)/((1+s_0009/Km0009_0061)*(1+s_1198/Km1198_0061)+(1+s_0010/Km0010_0061)*(1+s_1203/Km1203_0061)-1);
	r_0065=cell*Vmax_0065*(s_0261*s_1360-s_0324*s_1322/Keq_0065)/(Km0261_0065*Km1360_0065)/((1+s_0261/Km0261_0065)*(1+s_1360/Km1360_0065)+(1+s_0324/Km0324_0065)*(1+s_1322/Km1322_0065)-1);
	r_0079=cell*Vmax_0079*(s_0301*s_0434*s_0999-s_0302*s_0394*s_0991*s_1322/Keq_0079)/(Km0301_0079*Km0434_0079*Km0999_0079)/((1+s_0301/Km0301_0079)*(1+s_0434/Km0434_0079)*(1+s_0999/Km0999_0079)+(1+s_0302/Km0302_0079)*(1+s_0394/Km0394_0079)*(1+s_0991/Km0991_0079)*(1+s_1322/Km1322_0079)-1);
	r_0080=cell*Vmax_0080*(s_0306*s_1212-s_0322*s_1207/Keq_0080)/(Km0306_0080*Km1212_0080)/((1+s_0306/Km0306_0080)*(1+s_1212/Km1212_0080)+(1+s_0322/Km0322_0080)*(1+s_1207/Km1207_0080)-1);
	r_0091=cell*Vmax_0091*(s_0335-s_0340/Keq_0091)/Km0335_0091/(1+s_0335/Km0335_0091+1+s_0340/Km0340_0091-1);
	r_0096=cell*Vmax_0096*(s_0146*s_1212-s_0016*s_1207/Keq_0096)/(Km0146_0096*Km1212_0096)/((1+s_0146/Km0146_0096)*(1+s_1212/Km1212_0096)+(1+s_0016/Km0016_0096)*(1+s_1207/Km1207_0096)-1);
	r_0097=cell*Vmax_0097*(pow(s_1399,2)-s_0146*s_0456/Keq_0097)/pow(Km1399_0097,2)/(pow(1+s_1399/Km1399_0097,2)+(1+s_0146/Km0146_0097)*(1+s_0456/Km0456_0097)-1);
	r_0103=cell*Vmax_0103*(pow(s_0373,2)-s_0367*s_0529/Keq_0103)/pow(Km0373_0103,2)/(pow(1+s_0373/Km0373_0103,2)+(1+s_0367/Km0367_0103)*(1+s_0529/Km0529_0103)-1);
	r_0108=cell*Vmax_0108*(s_0373*s_0434*s_0445-s_0394*s_1101*s_1322/Keq_0108)/(Km0373_0108*Km0434_0108*Km0445_0108)/((1+s_0373/Km0373_0108)*(1+s_0434/Km0434_0108)*(1+s_0445/Km0445_0108)+(1+s_0394/Km0394_0108)*(1+s_1101/Km1101_0108)*(1+s_1322/Km1322_0108)-1);
	r_0111=cell*Vmax_0111*(s_0373-s_0362*s_0529/Keq_0111)/Km0373_0111/(1+s_0373/Km0373_0111+(1+s_0362/Km0362_0111)*(1+s_0529/Km0529_0111)-1);
	r_0115=cell*Vmax_0115*(s_0434*s_1192-s_0394*s_1191/Keq_0115)/(Km0434_0115*Km1192_0115)/((1+s_0434/Km0434_0115)*(1+s_1192/Km1192_0115)+(1+s_0394/Km0394_0115)*(1+s_1191/Km1191_0115)-1);
	r_0118=cell*Vmax_0118*(s_0145*s_0991-s_0180*s_1182/Keq_0118)/(Km0145_0118*Km0991_0118)/((1+s_0145/Km0145_0118)*(1+s_0991/Km0991_0118)+(1+s_0180/Km0180_0118)*(1+s_1182/Km1182_0118)-1);
	r_0142=cell*Vmax_0142*(s_0386*s_0434-s_0394*s_0423/Keq_0142)/(Km0386_0142*Km0434_0142)/((1+s_0386/Km0386_0142)*(1+s_0434/Km0434_0142)+(1+s_0394/Km0394_0142)*(1+s_0423/Km0423_0142)-1);
	r_0144=cell*Vmax_0144*(s_1413-s_0386*s_1012/Keq_0144)/Km1413_0144/(1+s_1413/Km1413_0144+(1+s_0386/Km0386_0144)*(1+s_1012/Km1012_0144)-1);
	r_0148=cell*Vmax_0148*(s_0423*s_0434-pow(s_0394,2)/Keq_0148)/(Km0423_0148*Km0434_0148)/((1+s_0423/Km0423_0148)*(1+s_0434/Km0434_0148)+pow(1+s_0394/Km0394_0148,2)-1);
	r_0151=cell*Vmax_0151*(s_0299-s_0403*s_0725/Keq_0151)/Km0299_0151/(1+s_0299/Km0299_0151+(1+s_0403/Km0403_0151)*(1+s_0725/Km0725_0151)-1);
	r_0152=cell*Vmax_0152*(s_0393-s_0423*s_0725/Keq_0152)/Km0393_0152/(1+s_0393/Km0393_0152+(1+s_0423/Km0423_0152)*(1+s_0725/Km0725_0152)-1);
	r_0153=cell*Vmax_0153*(s_0785*s_0849*s_0973-s_0393*s_0739*s_1322/Keq_0153)/(Km0785_0153*Km0849_0153*Km0973_0153)/((1+s_0785/Km0785_0153)*(1+s_0849/Km0849_0153)*(1+s_0973/Km0973_0153)+(1+s_0393/Km0393_0153)*(1+s_0739/Km0739_0153)*(1+s_1322/Km1322_0153)-1);
	r_0154=cell*Vmax_0154*(s_0298*s_0434-s_0201*s_0394/Keq_0154)/(Km0298_0154*Km0434_0154)/((1+s_0298/Km0298_0154)*(1+s_0434/Km0434_0154)+(1+s_0201/Km0201_0154)*(1+s_0394/Km0394_0154)-1);
	r_0165=cell*Vmax_0165*(s_0359*s_1203-s_0680*s_1198/Keq_0165)/(Km0359_0165*Km1203_0165)/((1+s_0359/Km0359_0165)*(1+s_1203/Km1203_0165)+(1+s_0680/Km0680_0165)*(1+s_1198/Km1198_0165)-1);
	r_0173=cell*Vmax_0173*(s_0359*s_1207-s_0362*s_1212/Keq_0173)/(Km0359_0173*Km1207_0173)/((1+s_0359/Km0359_0173)*(1+s_1207/Km1207_0173)+(1+s_0362/Km0362_0173)*(1+s_1212/Km1212_0173)-1);
	r_0195=cell*Vmax_0195*(s_0568*s_1543-s_0409*s_1538/Keq_0195)/(Km0568_0195*Km1543_0195)/((1+s_0568/Km0568_0195)*(1+s_1543/Km1543_0195)+(1+s_0409/Km0409_0195)*(1+s_1538/Km1538_0195)-1);
	r_0202=cell*Vmax_0202*(s_0427*s_1386-s_0633*s_1187/Keq_0202)/(Km0427_0202*Km1386_0202)/((1+s_0427/Km0427_0202)*(1+s_1386/Km1386_0202)+(1+s_0633/Km0633_0202)*(1+s_1187/Km1187_0202)-1);
	r_0203=cell*Vmax_0203*(s_0515*s_0999-s_0427*s_0991*s_1399/Keq_0203)/(Km0515_0203*Km0999_0203)/((1+s_0515/Km0515_0203)*(1+s_0999/Km0999_0203)+(1+s_0427/Km0427_0203)*(1+s_0991/Km0991_0203)*(1+s_1399/Km1399_0203)-1);
	r_0207=cell*Vmax_0207*(s_0015-s_0725*s_0965/Keq_0207)/Km0015_0207/(1+s_0015/Km0015_0207+(1+s_0725/Km0725_0207)*(1+s_0965/Km0965_0207)-1);
	r_0208=cell*Vmax_0208*(s_0434*s_0973*s_0979-s_0015*s_0423*s_0633/Keq_0208)/(Km0434_0208*Km0973_0208*Km0979_0208)/((1+s_0434/Km0434_0208)*(1+s_0973/Km0973_0208)*(1+s_0979/Km0979_0208)+(1+s_0015/Km0015_0208)*(1+s_0423/Km0423_0208)*(1+s_0633/Km0633_0208)-1);
	r_0211=cell*Vmax_0211*(s_0434*s_0973*s_0999-s_0423*s_0633*s_0969*s_0991/Keq_0211)/(Km0434_0211*Km0973_0211*Km0999_0211)/((1+s_0434/Km0434_0211)*(1+s_0973/Km0973_0211)*(1+s_0999/Km0999_0211)+(1+s_0423/Km0423_0211)*(1+s_0633/Km0633_0211)*(1+s_0969/Km0969_0211)*(1+s_0991/Km0991_0211)-1);
	r_0214=cell*Vmax_0214*(s_0455*s_0973-s_1194*s_1322/Keq_0214)/(Km0455_0214*Km0973_0214)/((1+s_0455/Km0455_0214)*(1+s_0973/Km0973_0214)+(1+s_1194/Km1194_0214)*(1+s_1322/Km1322_0214)-1);
	r_0215=cell*Vmax_0215*(s_0434*s_0973-s_0295*s_0394/Keq_0215)/(Km0434_0215*Km0973_0215)/((1+s_0434/Km0434_0215)*(1+s_0973/Km0973_0215)+(1+s_0295/Km0295_0215)*(1+s_0394/Km0394_0215)-1);
	r_0216=cell*Vmax_0216*(s_0991*s_1271-s_0180*s_0973/Keq_0216)/(Km0991_0216*Km1271_0216)/((1+s_0991/Km0991_0216)*(1+s_1271/Km1271_0216)+(1+s_0180/Km0180_0216)*(1+s_0973/Km0973_0216)-1);
	r_0219=cell*Vmax_0219*(s_0295*s_1212-s_0978*s_1207*s_1322/Keq_0219)/(Km0295_0219*Km1212_0219)/((1+s_0295/Km0295_0219)*(1+s_1212/Km1212_0219)+(1+s_0978/Km0978_0219)*(1+s_1207/Km1207_0219)*(1+s_1322/Km1322_0219)-1);
	r_0225=cell*Vmax_0225*(s_0434*s_1386-s_0326*s_0633/Keq_0225)/(Km0434_0225*Km1386_0225)/((1+s_0434/Km0434_0225)*(1+s_1386/Km1386_0225)+(1+s_0326/Km0326_0225)*(1+s_0633/Km0633_0225)-1);
	r_0226=cell*Vmax_0226*(s_0394*s_1322-s_0434/Keq_0226)/(Km0394_0226*Km1322_0226)/((1+s_0394/Km0394_0226)*(1+s_1322/Km1322_0226)+1+s_0434/Km0434_0226-1);
	r_0231=cell*Vmax_0231*(s_0262*s_1212-s_0122*s_1207/Keq_0231)/(Km0262_0231*Km1212_0231)/((1+s_0262/Km0262_0231)*(1+s_1212/Km1212_0231)+(1+s_0122/Km0122_0231)*(1+s_1207/Km1207_0231)-1);
	r_0233=cell*Vmax_0233*(s_0664*s_1212*s_1275-s_0662*s_1207/Keq_0233)/(Km0664_0233*Km1212_0233*Km1275_0233)/((1+s_0664/Km0664_0233)*(1+s_1212/Km1212_0233)*(1+s_1275/Km1275_0233)+(1+s_0662/Km0662_0233)*(1+s_1207/Km1207_0233)-1);
	r_0234=cell*Vmax_0234*(s_1207*s_1578-s_0456*s_1212*s_1579/Keq_0234)/(Km1207_0234*Km1578_0234)/((1+s_1207/Km1207_0234)*(1+s_1578/Km1578_0234)+(1+s_0456/Km0456_0234)*(1+s_1212/Km1212_0234)*(1+s_1579/Km1579_0234)-1);
	r_0235=cell*Vmax_0235*(s_0297*s_1198-s_0209*s_0456*s_1203/Keq_0235)/(Km0297_0235*Km1198_0235)/((1+s_0297/Km0297_0235)*(1+s_1198/Km1198_0235)+(1+s_0209/Km0209_0235)*(1+s_0456/Km0456_0235)*(1+s_1203/Km1203_0235)-1);
	r_0236=cell*Vmax_0236*(s_0209*s_1212-s_0296*s_1207/Keq_0236)/(Km0209_0236*Km1212_0236)/((1+s_0209/Km0209_0236)*(1+s_1212/Km1212_0236)+(1+s_0296/Km0296_0236)*(1+s_1207/Km1207_0236)-1);
	r_0237=cell*Vmax_0237*(s_1212*s_1579-s_1207*s_1569/Keq_0237)/(Km1212_0237*Km1579_0237)/((1+s_1212/Km1212_0237)*(1+s_1579/Km1579_0237)+(1+s_1207/Km1207_0237)*(1+s_1569/Km1569_0237)-1);
	r_0238=cell*Vmax_0238*(s_0296*s_1212*s_1275-s_1207*s_1576/Keq_0238)/(Km0296_0238*Km1212_0238*Km1275_0238)/((1+s_0296/Km0296_0238)*(1+s_1212/Km1212_0238)*(1+s_1275/Km1275_0238)+(1+s_1207/Km1207_0238)*(1+s_1576/Km1576_0238)-1);
	r_0239=cell*Vmax_0239*(s_1212*s_1275*s_1576-s_1207*s_1577/Keq_0239)/(Km1212_0239*Km1275_0239*Km1576_0239)/((1+s_1212/Km1212_0239)*(1+s_1275/Km1275_0239)*(1+s_1576/Km1576_0239)+(1+s_1207/Km1207_0239)*(1+s_1577/Km1577_0239)-1);
	r_0240=cell*Vmax_0240*(s_1212*s_1275*s_1577-s_1207*s_1578/Keq_0240)/(Km1212_0240*Km1275_0240*Km1577_0240)/((1+s_1212/Km1212_0240)*(1+s_1275/Km1275_0240)*(1+s_1577/Km1577_0240)+(1+s_1207/Km1207_0240)*(1+s_1578/Km1578_0240)-1);
	r_0241=cell*Vmax_0241*(s_0122*pow(s_1212,3)*pow(s_1275,3)-s_0297*pow(s_1207,3)/Keq_0241)/(Km0122_0241*pow(Km1212_0241,3)*pow(Km1275_0241,3))/((1+s_0122/Km0122_0241)*pow(1+s_1212/Km1212_0241,3)*pow(1+s_1275/Km1275_0241,3)+(1+s_0297/Km0297_0241)*pow(1+s_1207/Km1207_0241,3)-1);
	r_0242=cell*Vmax_0242*(s_0657*s_1212*s_1275-s_0664*s_1207/Keq_0242)/(Km0657_0242*Km1212_0242*Km1275_0242)/((1+s_0657/Km0657_0242)*(1+s_1212/Km1212_0242)*(1+s_1275/Km1275_0242)+(1+s_0664/Km0664_0242)*(1+s_1207/Km1207_0242)-1);
	r_0243=cell*Vmax_0243*(s_0700-s_0657/Keq_0243)/Km0700_0243/(1+s_0700/Km0700_0243+1+s_0657/Km0657_0243-1);
	r_0244=cell*Vmax_0244*(s_0662*s_1212-s_0666*s_1207/Keq_0244)/(Km0662_0244*Km1212_0244)/((1+s_0662/Km0662_0244)*(1+s_1212/Km1212_0244)+(1+s_0666/Km0666_0244)*(1+s_1207/Km1207_0244)-1);
	r_0250=cell*Vmax_0250*(pow(s_0434,2)*s_0445*s_0999-pow(s_0394,2)*s_0455*s_0991*s_1322/Keq_0250)/(pow(Km0434_0250,2)*Km0445_0250*Km0999_0250)/(pow(1+s_0434/Km0434_0250,2)*(1+s_0445/Km0445_0250)*(1+s_0999/Km0999_0250)+pow(1+s_0394/Km0394_0250,2)*(1+s_0455/Km0455_0250)*(1+s_0991/Km0991_0250)*(1+s_1322/Km1322_0250)-1);
	r_0257=cell*Vmax_0257*(s_0539*s_1331-s_0471*s_0633/Keq_0257)/(Km0539_0257*Km1331_0257)/((1+s_0539/Km0539_0257)*(1+s_1331/Km1331_0257)+(1+s_0471/Km0471_0257)*(1+s_0633/Km0633_0257)-1);
	r_0259=cell*Vmax_0259*(s_0475*s_1212*s_1275-s_0481*s_1207/Keq_0259)/(Km0475_0259*Km1212_0259*Km1275_0259)/((1+s_0475/Km0475_0259)*(1+s_1212/Km1212_0259)*(1+s_1275/Km1275_0259)+(1+s_0481/Km0481_0259)*(1+s_1207/Km1207_0259)-1);
	r_0267=cell*Vmax_0267*(s_0481*s_1212*s_1275-s_0493*s_1207/Keq_0267)/(Km0481_0267*Km1212_0267*Km1275_0267)/((1+s_0481/Km0481_0267)*(1+s_1212/Km1212_0267)*(1+s_1275/Km1275_0267)+(1+s_0493/Km0493_0267)*(1+s_1207/Km1207_0267)-1);
	r_0269=cell*Vmax_0269*(s_0493*s_1212*s_1275-s_0499*s_1207/Keq_0269)/(Km0493_0269*Km1212_0269*Km1275_0269)/((1+s_0493/Km0493_0269)*(1+s_1212/Km1212_0269)*(1+s_1275/Km1275_0269)+(1+s_0499/Km0499_0269)*(1+s_1207/Km1207_0269)-1);
	r_0278=cell*Vmax_0278*(s_0515-s_1377/Keq_0278)/Km0515_0278/(1+s_0515/Km0515_0278+1+s_1377/Km1377_0278-1);
	r_0279=cell*Vmax_0279*(s_0324-s_0515*s_1322/Keq_0279)/Km0324_0279/(1+s_0324/Km0324_0279+(1+s_0515/Km0515_0279)*(1+s_1322/Km1322_0279)-1);
	r_0280=cell*Vmax_0280*(s_0516-s_0940/Keq_0280)/Km0516_0280/(1+s_0516/Km0516_0280+1+s_0940/Km0940_0280-1);
	r_0300=cell*Vmax_0300*(s_0373*s_1271-s_0522*s_0529/Keq_0300)/(Km0373_0300*Km1271_0300)/((1+s_0373/Km0373_0300)*(1+s_1271/Km1271_0300)+(1+s_0522/Km0522_0300)*(1+s_0529/Km0529_0300)-1);
	r_0302=cell*Vmax_0302*(s_0522-s_0516/Keq_0302)/Km0522_0302/(1+s_0522/Km0522_0302+1+s_0516/Km0516_0302-1);
	r_0307=cell*Vmax_0307*(s_0419*s_0434*s_1559-s_0394*s_0539*s_1322/Keq_0307)/(Km0419_0307*Km0434_0307*Km1559_0307)/((1+s_0419/Km0419_0307)*(1+s_0434/Km0434_0307)*(1+s_1559/Km1559_0307)+(1+s_0394/Km0394_0307)*(1+s_0539/Km0539_0307)*(1+s_1322/Km1322_0307)-1);
	r_0309=cell*Vmax_0309*(s_1012*s_1039-s_0980/Keq_0309)/(Km1012_0309*Km1039_0309)/((1+s_1012/Km1012_0309)*(1+s_1039/Km1039_0309)+1+s_0980/Km0980_0309-1);
	r_0310=cell*Vmax_0310*(s_0980-s_0178*s_0419*s_0981/Keq_0310)/Km0980_0310/(1+s_0980/Km0980_0310+(1+s_0178/Km0178_0310)*(1+s_0419/Km0419_0310)*(1+s_0981/Km0981_0310)-1);
	r_0311=cell*Vmax_0311*(s_0981*s_1233-s_0362*s_0980/Keq_0311)/(Km0981_0311*Km1233_0311)/((1+s_0981/Km0981_0311)*(1+s_1233/Km1233_0311)+(1+s_0362/Km0362_0311)*(1+s_0980/Km0980_0311)-1);
	r_0312=cell*Vmax_0312*(s_0841*s_1234-s_0362*s_0981/Keq_0312)/(Km0841_0312*Km1234_0312)/((1+s_0841/Km0841_0312)*(1+s_1234/Km1234_0312)+(1+s_0362/Km0362_0312)*(1+s_0981/Km0981_0312)-1);
	r_0317=cell*Vmax_0317*(s_1059*pow(s_1212,3)*pow(s_1275,3)-s_0262*s_0722*pow(s_1207,3)/Keq_0317)/(Km1059_0317*pow(Km1212_0317,3)*pow(Km1275_0317,3))/((1+s_1059/Km1059_0317)*pow(1+s_1212/Km1212_0317,3)*pow(1+s_1275/Km1275_0317,3)+(1+s_0262/Km0262_0317)*(1+s_0722/Km0722_0317)*pow(1+s_1207/Km1207_0317,3)-1);
	r_0326=cell*Vmax_0326*(s_0419*s_0654-s_0589/Keq_0326)/(Km0419_0326*Km0654_0326)/((1+s_0419/Km0419_0326)*(1+s_0654/Km0654_0326)+1+s_0589/Km0589_0326-1);
	r_0330=cell*Vmax_0330*(s_0394*s_0613-s_0434*s_0615/Keq_0330)/(Km0394_0330*Km0613_0330)/((1+s_0394/Km0394_0330)*(1+s_0613/Km0613_0330)+(1+s_0434/Km0434_0330)*(1+s_0615/Km0615_0330)-1);
	r_0336=cell*Vmax_0336*(s_0529*s_1524-s_0380*s_0619/Keq_0336)/(Km0529_0336*Km1524_0336)/((1+s_0529/Km0529_0336)*(1+s_1524/Km1524_0336)+(1+s_0380/Km0380_0336)*(1+s_0619/Km0619_0336)-1);
	r_0337=cell*Vmax_0337*(s_1331-s_0619*s_1322/Keq_0337)/Km1331_0337/(1+s_1331/Km1331_0337+(1+s_0619/Km0619_0337)*(1+s_1322/Km1322_0337)-1);
	r_0339=cell*Vmax_0339*(s_0061*s_1275-s_0837*s_1269/Keq_0339)/(Km0061_0339*Km1275_0339)/((1+s_0061/Km0061_0339)*(1+s_1275/Km1275_0339)+(1+s_0837/Km0837_0339)*(1+s_1269/Km1269_0339)-1);
	r_0340=cell*Vmax_0340*(s_1084*s_1445-s_0475/Keq_0340)/(Km1084_0340*Km1445_0340)/((1+s_1084/Km1084_0340)*(1+s_1445/Km1445_0340)+1+s_0475/Km0475_0340-1);
	r_0344=cell*Vmax_0344*(s_0625*s_1212-s_1207*s_1487/Keq_0344)/(Km0625_0344*Km1212_0344)/((1+s_0625/Km0625_0344)*(1+s_1212/Km1212_0344)+(1+s_1207/Km1207_0344)*(1+s_1487/Km1487_0344)-1);
	r_0349=cell*Vmax_0349*(s_1194-s_0061/Keq_0349)/Km1194_0349/(1+s_1194/Km1194_0349+1+s_0061/Km0061_0349-1);
	r_0352=cell*Vmax_0352*(s_0016-s_0232/Keq_0352)/Km0016_0352/(1+s_0016/Km0016_0352+1+s_0232/Km0232_0352-1);
	r_0353=cell*Vmax_0353*(s_0008-s_0056/Keq_0353)/Km0008_0353/(1+s_0008/Km0008_0353+1+s_0056/Km0056_0353-1);
	r_0355=cell*Vmax_0355*(s_0943*s_1376-s_0633*s_0745/Keq_0355)/(Km0943_0355*Km1376_0355)/((1+s_0943/Km0943_0355)*(1+s_1376/Km1376_0355)+(1+s_0633/Km0633_0355)*(1+s_0745/Km0745_0355)-1);
	r_0361=cell*Vmax_0361*(s_0645*s_0743-s_0644*s_0739/Keq_0361)/(Km0645_0361*Km0743_0361)/((1+s_0645/Km0645_0361)*(1+s_0743/Km0743_0361)+(1+s_0644/Km0644_0361)*(1+s_0739/Km0739_0361)-1);
	r_0362=cell*Vmax_0362*(s_0644-s_0645*s_1107/Keq_0362)/Km0644_0362/(1+s_0644/Km0644_0362+(1+s_0645/Km0645_0362)*(1+s_1107/Km1107_0362)-1);
	r_0364=cell*Vmax_0364*(s_0656-s_0633*s_0654/Keq_0364)/Km0656_0364/(1+s_0656/Km0656_0364+(1+s_0633/Km0633_0364)*(1+s_0654/Km0654_0364)-1);
	r_0366=cell*Vmax_0366*(s_0188-s_1360/Keq_0366)/Km0188_0366/(1+s_0188/Km0188_0366+1+s_1360/Km1360_0366-1);
	r_0386=cell*Vmax_0386*(s_0595*s_1101*pow(s_1212,2)-s_0456*s_0529*s_1065*pow(s_1207,2)/Keq_0386)/(Km0595_0386*Km1101_0386*pow(Km1212_0386,2))/((1+s_0595/Km0595_0386)*(1+s_1101/Km1101_0386)*pow(1+s_1212/Km1212_0386,2)+(1+s_0456/Km0456_0386)*(1+s_0529/Km0529_0386)*(1+s_1065/Km1065_0386)*pow(1+s_1207/Km1207_0386,2)-1);
	r_0387=cell*Vmax_0387*(s_1065*s_1101*pow(s_1212,2)-s_0456*s_0529*s_1161*pow(s_1207,2)/Keq_0387)/(Km1065_0387*Km1101_0387*pow(Km1212_0387,2))/((1+s_1065/Km1065_0387)*(1+s_1101/Km1101_0387)*pow(1+s_1212/Km1212_0387,2)+(1+s_0456/Km0456_0387)*(1+s_0529/Km0529_0387)*(1+s_1161/Km1161_0387)*pow(1+s_1207/Km1207_0387,2)-1);
	r_0389=cell*Vmax_0389*(s_1101*s_1161*pow(s_1212,2)-s_0456*s_0529*pow(s_1207,2)*s_1286/Keq_0389)/(Km1101_0389*Km1161_0389*pow(Km1212_0389,2))/((1+s_1101/Km1101_0389)*(1+s_1161/Km1161_0389)*pow(1+s_1212/Km1212_0389,2)+(1+s_0456/Km0456_0389)*(1+s_0529/Km0529_0389)*pow(1+s_1207/Km1207_0389,2)*(1+s_1286/Km1286_0389)-1);
	r_0391=cell*Vmax_0391*(s_1101*pow(s_1212,2)*s_1286-s_0456*s_0529*pow(s_1207,2)*s_1449/Keq_0391)/(Km1101_0391*pow(Km1212_0391,2)*Km1286_0391)/((1+s_1101/Km1101_0391)*pow(1+s_1212/Km1212_0391,2)*(1+s_1286/Km1286_0391)+(1+s_0456/Km0456_0391)*(1+s_0529/Km0529_0391)*pow(1+s_1207/Km1207_0391,2)*(1+s_1449/Km1449_0391)-1);
	r_0393=cell*Vmax_0393*(pow(s_1101,3)*pow(s_1212,6)*s_1449-pow(s_0456,3)*pow(s_0529,3)*s_1084*pow(s_1207,6)/Keq_0393)/(pow(Km1101_0393,3)*pow(Km1212_0393,6)*Km1449_0393)/(pow(1+s_1101/Km1101_0393,3)*pow(1+s_1212/Km1212_0393,6)*(1+s_1449/Km1449_0393)+pow(1+s_0456/Km0456_0393,3)*pow(1+s_0529/Km0529_0393,3)*(1+s_1084/Km1084_0393)*pow(1+s_1207/Km1207_0393,6)-1);
	r_0397=cell*Vmax_0397*(s_1101*pow(s_1212,2)*s_1255-s_0456*s_0529*s_0602*pow(s_1207,2)/Keq_0397)/(Km1101_0397*pow(Km1212_0397,2)*Km1255_0397)/((1+s_1101/Km1101_0397)*pow(1+s_1212/Km1212_0397,2)*(1+s_1255/Km1255_0397)+(1+s_0456/Km0456_0397)*(1+s_0529/Km0529_0397)*(1+s_0602/Km0602_0397)*pow(1+s_1207/Km1207_0397,2)-1);
	r_0398=cell*Vmax_0398*(s_0373*pow(s_1101,3)*pow(s_1212,6)-pow(s_0456,3)*pow(s_0529,3)*pow(s_1207,6)*s_1255/Keq_0398)/(Km0373_0398*pow(Km1101_0398,3)*pow(Km1212_0398,6))/((1+s_0373/Km0373_0398)*pow(1+s_1101/Km1101_0398,3)*pow(1+s_1212/Km1212_0398,6)+pow(1+s_0456/Km0456_0398,3)*pow(1+s_0529/Km0529_0398,3)*pow(1+s_1207/Km1207_0398,6)*(1+s_1255/Km1255_0398)-1);
	r_0399=cell*Vmax_0399*(s_0423*s_0602*s_0633-s_0434*s_0529*s_0595/Keq_0399)/(Km0423_0399*Km0602_0399*Km0633_0399)/((1+s_0423/Km0423_0399)*(1+s_0602/Km0602_0399)*(1+s_0633/Km0633_0399)+(1+s_0434/Km0434_0399)*(1+s_0529/Km0529_0399)*(1+s_0595/Km0595_0399)-1);
	r_0407=cell*Vmax_0407*(s_0423*s_0633*s_1454-s_0434*s_0529*s_1449/Keq_0407)/(Km0423_0407*Km0633_0407*Km1454_0407)/((1+s_0423/Km0423_0407)*(1+s_0633/Km0633_0407)*(1+s_1454/Km1454_0407)+(1+s_0434/Km0434_0407)*(1+s_0529/Km0529_0407)*(1+s_1449/Km1449_0407)-1);
	r_0432=cell*Vmax_0432*(s_0602*s_1101*pow(s_1212,2)-s_0456*s_0529*s_1073*pow(s_1207,2)/Keq_0432)/(Km0602_0432*Km1101_0432*pow(Km1212_0432,2))/((1+s_0602/Km0602_0432)*(1+s_1101/Km1101_0432)*pow(1+s_1212/Km1212_0432,2)+(1+s_0456/Km0456_0432)*(1+s_0529/Km0529_0432)*(1+s_1073/Km1073_0432)*pow(1+s_1207/Km1207_0432,2)-1);
	r_0433=cell*Vmax_0433*(s_1073*s_1101*pow(s_1212,2)-s_0456*s_0529*s_1176*pow(s_1207,2)/Keq_0433)/(Km1073_0433*Km1101_0433*pow(Km1212_0433,2))/((1+s_1073/Km1073_0433)*(1+s_1101/Km1101_0433)*pow(1+s_1212/Km1212_0433,2)+(1+s_0456/Km0456_0433)*(1+s_0529/Km0529_0433)*(1+s_1176/Km1176_0433)*pow(1+s_1207/Km1207_0433,2)-1);
	r_0434=cell*Vmax_0434*(s_1101*s_1176*pow(s_1212,2)-s_0456*s_0529*pow(s_1207,2)*s_1302/Keq_0434)/(Km1101_0434*Km1176_0434*pow(Km1212_0434,2))/((1+s_1101/Km1101_0434)*(1+s_1176/Km1176_0434)*pow(1+s_1212/Km1212_0434,2)+(1+s_0456/Km0456_0434)*(1+s_0529/Km0529_0434)*pow(1+s_1207/Km1207_0434,2)*(1+s_1302/Km1302_0434)-1);
	r_0435=cell*Vmax_0435*(s_1101*pow(s_1212,2)*s_1302-s_0456*s_0529*pow(s_1207,2)*s_1454/Keq_0435)/(Km1101_0435*pow(Km1212_0435,2)*Km1302_0435)/((1+s_1101/Km1101_0435)*pow(1+s_1212/Km1212_0435,2)*(1+s_1302/Km1302_0435)+(1+s_0456/Km0456_0435)*(1+s_0529/Km0529_0435)*pow(1+s_1207/Km1207_0435,2)*(1+s_1454/Km1454_0435)-1);
	r_0438=cell*Vmax_0438*(pow(s_0710,4)*s_1275-pow(s_0709,4)/Keq_0438)/(pow(Km0710_0438,4)*Km1275_0438)/(pow(1+s_0710/Km0710_0438,4)*(1+s_1275/Km1275_0438)+pow(1+s_0709/Km0709_0438,4)-1);
	r_0439=cell*Vmax_0439*(pow(s_0709,2)*s_1535-pow(s_0710,2)*s_1537/Keq_0439)/(pow(Km0709_0439,2)*Km1535_0439)/(pow(1+s_0709/Km0709_0439,2)*(1+s_1535/Km1535_0439)+pow(1+s_0710/Km0710_0439,2)*(1+s_1537/Km1537_0439)-1);
	r_0445=cell*Vmax_0445*(s_0722*s_1198-s_0456*s_1203/Keq_0445)/(Km0722_0445*Km1198_0445)/((1+s_0722/Km0722_0445)*(1+s_1198/Km1198_0445)+(1+s_0456/Km0456_0445)*(1+s_1203/Km1203_0445)-1);
	r_0446=cell*Vmax_0446*(s_0120*s_0394*s_1322-s_0434*s_0722*s_1487/Keq_0446)/(Km0120_0446*Km0394_0446*Km1322_0446)/((1+s_0120/Km0120_0446)*(1+s_0394/Km0394_0446)*(1+s_1322/Km1322_0446)+(1+s_0434/Km0434_0446)*(1+s_0722/Km0722_0446)*(1+s_1487/Km1487_0446)-1);
	r_0450=cell*Vmax_0450*(s_0555-s_0629*s_0764/Keq_0450)/Km0555_0450/(1+s_0555/Km0555_0450+(1+s_0629/Km0629_0450)*(1+s_0764/Km0764_0450)-1);
	r_0451=cell*Vmax_0451*(s_0725-s_0066/Keq_0451)/Km0725_0451/(1+s_0725/Km0725_0451+1+s_0066/Km0066_0451-1);
	r_0462=cell*Vmax_0462*(s_0745*s_0943-s_0190*s_0633/Keq_0462)/(Km0745_0462*Km0943_0462)/((1+s_0745/Km0745_0462)*(1+s_0943/Km0943_0462)+(1+s_0190/Km0190_0462)*(1+s_0633/Km0633_0462)-1);
	r_0466=cell*Vmax_0466*(s_0568*s_1207-s_0335*s_1212/Keq_0466)/(Km0568_0466*Km1207_0466)/((1+s_0568/Km0568_0466)*(1+s_1207/Km1207_0466)+(1+s_0335/Km0335_0466)*(1+s_1212/Km1212_0466)-1);
	r_0467=cell*Vmax_0467*(s_0568-s_0557/Keq_0467)/Km0568_0467/(1+s_0568/Km0568_0467+1+s_0557/Km0557_0467-1);
	r_0470=cell*Vmax_0470*(s_0180*s_0419*s_1203-s_0991*s_1198/Keq_0470)/(Km0180_0470*Km0419_0470*Km1203_0470)/((1+s_0180/Km0180_0470)*(1+s_0419/Km0419_0470)*(1+s_1203/Km1203_0470)+(1+s_0991/Km0991_0470)*(1+s_1198/Km1198_0470)-1);
	r_0471=cell*Vmax_0471*(s_0180*s_0419*s_1212-s_0991*s_1207/Keq_0471)/(Km0180_0471*Km0419_0471*Km1212_0471)/((1+s_0180/Km0180_0471)*(1+s_0419/Km0419_0471)*(1+s_1212/Km1212_0471)+(1+s_0991/Km0991_0471)*(1+s_1207/Km1207_0471)-1);
	r_0476=cell*Vmax_0476*(s_0419*s_0434*s_0991-s_0394*s_0999*s_1322/Keq_0476)/(Km0419_0476*Km0434_0476*Km0991_0476)/((1+s_0419/Km0419_0476)*(1+s_0434/Km0434_0476)*(1+s_0991/Km0991_0476)+(1+s_0394/Km0394_0476)*(1+s_0999/Km0999_0476)*(1+s_1322/Km1322_0476)-1);
	r_0481=cell*Vmax_0481*(s_0754*s_1212-pow(s_0750,2)*s_1207/Keq_0481)/(Km0754_0481*Km1212_0481)/((1+s_0754/Km0754_0481)*(1+s_1212/Km1212_0481)+pow(1+s_0750/Km0750_0481,2)*(1+s_1207/Km1207_0481)-1);
	r_0483=cell*Vmax_0483*(pow(s_0750,2)*s_0837-s_0754/Keq_0483)/(pow(Km0750_0483,2)*Km0837_0483)/(pow(1+s_0750/Km0750_0483,2)*(1+s_0837/Km0837_0483)+1+s_0754/Km0754_0483-1);
	r_0486=cell*Vmax_0486*(s_0764*s_1198*s_1322-s_0075*s_1203/Keq_0486)/(Km0764_0486*Km1198_0486*Km1322_0486)/((1+s_0764/Km0764_0486)*(1+s_1198/Km1198_0486)*(1+s_1322/Km1322_0486)+(1+s_0075/Km0075_0486)*(1+s_1203/Km1203_0486)-1);
	r_0489=cell*Vmax_0489*(s_0767-s_0765*s_1322/Keq_0489)/Km0767_0489/(1+s_0767/Km0767_0489+(1+s_0765/Km0765_0489)*(1+s_1322/Km1322_0489)-1);
	r_0491=cell*Vmax_0491*(s_0629*s_1203-s_0767*s_1198/Keq_0491)/(Km0629_0491*Km1203_0491)/((1+s_0629/Km0629_0491)*(1+s_1203/Km1203_0491)+(1+s_0767/Km0767_0491)*(1+s_1198/Km1198_0491)-1);
	r_0495=cell*Vmax_0495*(s_0380*s_0767-s_0082*s_0529/Keq_0495)/(Km0380_0495*Km0767_0495)/((1+s_0380/Km0380_0495)*(1+s_0767/Km0767_0495)+(1+s_0082/Km0082_0495)*(1+s_0529/Km0529_0495)-1);
	r_0499=cell*Vmax_0499*(s_0120*s_0325-s_0301*s_1487/Keq_0499)/(Km0120_0499*Km0325_0499)/((1+s_0120/Km0120_0499)*(1+s_0325/Km0325_0499)+(1+s_0301/Km0301_0499)*(1+s_1487/Km1487_0499)-1);
	r_0502=cell*Vmax_0502*(s_1039*s_1487-s_0306*s_1003/Keq_0502)/(Km1039_0502*Km1487_0502)/((1+s_1039/Km1039_0502)*(1+s_1487/Km1487_0502)+(1+s_0306/Km0306_0502)*(1+s_1003/Km1003_0502)-1);
	r_0510=cell*Vmax_0510*(s_1543-s_0773*s_1538/Keq_0510)/Km1543_0510/(1+s_1543/Km1543_0510+(1+s_0773/Km0773_0510)*(1+s_1538/Km1538_0510)-1);
	r_0514=cell*Vmax_0514*(s_0434*s_0999*s_1565-s_0423*s_0633*s_0782*s_0991/Keq_0514)/(Km0434_0514*Km0999_0514*Km1565_0514)/((1+s_0434/Km0434_0514)*(1+s_0999/Km0999_0514)*(1+s_1565/Km1565_0514)+(1+s_0423/Km0423_0514)*(1+s_0633/Km0633_0514)*(1+s_0782/Km0782_0514)*(1+s_0991/Km0991_0514)-1);
	r_0525=cell*Vmax_0525*(s_0785-s_0141*s_0633*s_0722/Keq_0525)/Km0785_0525/(1+s_0785/Km0785_0525+(1+s_0141/Km0141_0525)*(1+s_0633/Km0633_0525)*(1+s_0722/Km0722_0525)-1);
	r_0528=cell*Vmax_0528*(s_0434*s_0782-s_0394*s_0739/Keq_0528)/(Km0434_0528*Km0782_0528)/((1+s_0434/Km0434_0528)*(1+s_0782/Km0782_0528)+(1+s_0394/Km0394_0528)*(1+s_0739/Km0739_0528)-1);
	r_0529=cell*Vmax_0529*(s_0586*s_0782-s_0582*s_0739/Keq_0529)/(Km0586_0529*Km0782_0529)/((1+s_0586/Km0586_0529)*(1+s_0782/Km0782_0529)+(1+s_0582/Km0582_0529)*(1+s_0739/Km0739_0529)-1);
	r_0534=cell*Vmax_0534*(s_0434*s_0563-s_0394*s_0568/Keq_0534)/(Km0434_0534*Km0563_0534)/((1+s_0434/Km0434_0534)*(1+s_0563/Km0563_0534)+(1+s_0394/Km0394_0534)*(1+s_0568/Km0568_0534)-1);
	r_0536=cell*Vmax_0536*(s_1010*pow(s_1198,2)-s_1006*pow(s_1203,2)/Keq_0536)/(Km1010_0536*pow(Km1198_0536,2))/((1+s_1010/Km1010_0536)*pow(1+s_1198/Km1198_0536,2)+(1+s_1006/Km1006_0536)*pow(1+s_1203/Km1203_0536,2)-1);
	r_0537=cell*Vmax_0537*(s_1011-s_1010*s_1322/Keq_0537)/Km1011_0537/(1+s_1011/Km1011_0537+(1+s_1010/Km1010_0537)*(1+s_1322/Km1322_0537)-1);
	r_0538=cell*Vmax_0538*(s_0207*s_0991-s_0180*s_1011/Keq_0538)/(Km0207_0538*Km0991_0538)/((1+s_0207/Km0207_0538)*(1+s_0991/Km0991_0538)+(1+s_0180/Km0180_0538)*(1+s_1011/Km1011_0538)-1);
	r_0542=cell*Vmax_0542*(s_0454-s_0836/Keq_0542)/Km0454_0542/(1+s_0454/Km0454_0542+1+s_0836/Km0836_0542-1);
	r_0543=cell*Vmax_0543*(s_0180*s_0373-s_0529*s_0835/Keq_0543)/(Km0180_0543*Km0373_0543)/((1+s_0180/Km0180_0543)*(1+s_0373/Km0373_0543)+(1+s_0529/Km0529_0543)*(1+s_0835/Km0835_0543)-1);
	r_0545=cell*Vmax_0545*(s_0836*s_1198-s_0176*s_1203*s_0456/Keq_0545)/(Km0836_0545*Km1198_0545)/((1+s_0836/Km0836_0545)*(1+s_1198/Km1198_0545)+(1+s_0176/Km0176_0545)*(1+s_1203/Km1203_0545)*(1+s_0456/Km0456_0545)-1);
	r_0547=cell*Vmax_0547*(s_0978*s_1212-s_1014*s_1207/Keq_0547)/(Km0978_0547*Km1212_0547)/((1+s_0978/Km0978_0547)*(1+s_1212/Km1212_0547)+(1+s_1014/Km1014_0547)*(1+s_1207/Km1207_0547)-1);
	r_0548=cell*Vmax_0548*(s_0434*s_1014-s_0394*s_1238/Keq_0548)/(Km0434_0548*Km1014_0548)/((1+s_0434/Km0434_0548)*(1+s_1014/Km1014_0548)+(1+s_0394/Km0394_0548)*(1+s_1238/Km1238_0548)-1);
	r_0549=cell*Vmax_0549*(s_0373*s_1014-s_0529*s_1233/Keq_0549)/(Km0373_0549*Km1014_0549)/((1+s_0373/Km0373_0549)*(1+s_1014/Km1014_0549)+(1+s_0529/Km0529_0549)*(1+s_1233/Km1233_0549)-1);
	r_0550=cell*Vmax_0550*(s_0837*s_1616-s_1620/Keq_0550)/(Km0837_0550*Km1616_0550)/((1+s_0837/Km0837_0550)*(1+s_1616/Km1616_0550)+1+s_1620/Km1620_0550-1);
	r_0553=cell*Vmax_0553*(s_0033-s_0025*s_0750/Keq_0553)/Km0033_0553/(1+s_0033/Km0033_0553+(1+s_0025/Km0025_0553)*(1+s_0750/Km0750_0553)-1);
	r_0558=cell*Vmax_0558*(s_0218*pow(s_1212,2)-s_0028*s_0529*pow(s_1207,2)/Keq_0558)/(Km0218_0558*pow(Km1212_0558,2))/((1+s_0218/Km0218_0558)*pow(1+s_1212/Km1212_0558,2)+(1+s_0028/Km0028_0558)*(1+s_0529/Km0529_0558)*pow(1+s_1207/Km1207_0558,2)-1);
	r_0559=cell*Vmax_0559*(s_0367*s_0373-s_0218*s_0529/Keq_0559)/(Km0367_0559*Km0373_0559)/((1+s_0367/Km0367_0559)*(1+s_0373/Km0373_0559)+(1+s_0218/Km0218_0559)*(1+s_0529/Km0529_0559)-1);
	r_0563=cell*Vmax_0563*(s_0312*s_0999-s_0403*s_0550*s_0991/Keq_0563)/(Km0312_0563*Km0999_0563)/((1+s_0312/Km0312_0563)*(1+s_0999/Km0999_0563)+(1+s_0403/Km0403_0563)*(1+s_0550/Km0550_0563)*(1+s_0991/Km0991_0563)-1);
	r_0564=cell*Vmax_0564*(s_0550-s_0207/Keq_0564)/Km0550_0564/(1+s_0550/Km0550_0564+1+s_0207/Km0207_0564-1);
	r_0565=cell*Vmax_0565*(s_0849*s_1198-s_1203*s_1565/Keq_0565)/(Km0849_0565*Km1198_0565)/((1+s_0849/Km0849_0565)*(1+s_1198/Km1198_0565)+(1+s_1203/Km1203_0565)*(1+s_1565/Km1565_0565)-1);
	r_0566=cell*Vmax_0566*(s_0076-s_0086*s_0456/Keq_0566)/Km0076_0566/(1+s_0076/Km0076_0566+(1+s_0086/Km0086_0566)*(1+s_0456/Km0456_0566)-1);
	r_0568=cell*Vmax_0568*(s_0633-pow(s_1322,2)/Keq_0568)/Km0633_0568/(1+s_0633/Km0633_0568+pow(1+s_1322/Km1322_0568,2)-1);
	r_0570=cell*Vmax_0570*(s_1365-s_0849/Keq_0570)/Km1365_0570/(1+s_1365/Km1365_0570+1+s_0849/Km0849_0570-1);
	r_0594=cell*Vmax_0594*(s_0089*s_0499-s_0619*s_0918/Keq_0594)/(Km0089_0594*Km0499_0594)/((1+s_0089/Km0089_0594)*(1+s_0499/Km0499_0594)+(1+s_0619/Km0619_0594)*(1+s_0918/Km0918_0594)-1);
	r_0658=cell*Vmax_0658*(s_0940*s_1198-s_0180*s_0456*s_1203/Keq_0658)/(Km0940_0658*Km1198_0658)/((1+s_0940/Km0940_0658)*(1+s_1198/Km1198_0658)+(1+s_0180/Km0180_0658)*(1+s_0456/Km0456_0658)*(1+s_1203/Km1203_0658)-1);
	r_0661=cell*Vmax_0661*(s_0940*s_1207-s_0180*s_0456*s_1212/Keq_0661)/(Km0940_0661*Km1207_0661)/((1+s_0940/Km0940_0661)*(1+s_1207/Km1207_0661)+(1+s_0180/Km0180_0661)*(1+s_0456/Km0456_0661)*(1+s_1212/Km1212_0661)-1);
	r_0663=cell*Vmax_0663*(s_0056*s_0991-s_0180*s_1016/Keq_0663)/(Km0056_0663*Km0991_0663)/((1+s_0056/Km0056_0663)*(1+s_0991/Km0991_0663)+(1+s_0180/Km0180_0663)*(1+s_1016/Km1016_0663)-1);
	r_0667=cell*Vmax_0667*(s_0943-s_1376/Keq_0667)/Km0943_0667/(1+s_0943/Km0943_0667+1+s_1376/Km1376_0667-1);
	r_0669=cell*Vmax_0669*(s_0039*s_1212-s_0008*s_1207/Keq_0669)/(Km0039_0669*Km1212_0669)/((1+s_0039/Km0039_0669)*(1+s_1212/Km1212_0669)+(1+s_0008/Km0008_0669)*(1+s_1207/Km1207_0669)-1);
	r_0670=cell*Vmax_0670*(s_1020-s_0427*s_0955/Keq_0670)/Km1020_0670/(1+s_1020/Km1020_0670+(1+s_0427/Km0427_0670)*(1+s_0955/Km0955_0670)-1);
	r_0674=cell*Vmax_0674*(s_0991*s_1399-s_0180*s_0955/Keq_0674)/(Km0991_0674*Km1399_0674)/((1+s_0991/Km0991_0674)*(1+s_1399/Km1399_0674)+(1+s_0180/Km0180_0674)*(1+s_0955/Km0955_0674)-1);
	r_0678=cell*Vmax_0678*(s_0953*s_1212-s_0959*s_1207/Keq_0678)/(Km0953_0678*Km1212_0678)/((1+s_0953/Km0953_0678)*(1+s_1212/Km1212_0678)+(1+s_0959/Km0959_0678)*(1+s_1207/Km1207_0678)-1);
	r_0688=cell*Vmax_0688*(s_1151*s_1212-s_0062*s_1207/Keq_0688)/(Km1151_0688*Km1212_0688)/((1+s_1151/Km1151_0688)*(1+s_1212/Km1212_0688)+(1+s_0062/Km0062_0688)*(1+s_1207/Km1207_0688)-1);
	r_0694=cell*Vmax_0694*(s_1048*s_1275-s_1195/Keq_0694)/(Km1048_0694*Km1275_0694)/((1+s_1048/Km1048_0694)*(1+s_1275/Km1275_0694)+1+s_1195/Km1195_0694-1);
	r_0696=cell*Vmax_0696*(s_0062*s_1198-s_0063*s_1203/Keq_0696)/(Km0062_0696*Km1198_0696)/((1+s_0062/Km0062_0696)*(1+s_1198/Km1198_0696)+(1+s_0063/Km0063_0696)*(1+s_1203/Km1203_0696)-1);
	r_0697=cell*Vmax_0697*(s_0750*s_1151-s_0033/Keq_0697)/(Km0750_0697*Km1151_0697)/((1+s_0750/Km0750_0697)*(1+s_1151/Km1151_0697)+1+s_0033/Km0033_0697-1);
	r_0698=cell*Vmax_0698*(s_0037-s_1059/Keq_0698)/Km0037_0698/(1+s_0037/Km0037_0698+1+s_1059/Km1059_0698-1);
	r_0699=cell*Vmax_0699*(s_0291*s_0991-s_0180*s_1021/Keq_0699)/(Km0291_0699*Km0991_0699)/((1+s_0291/Km0291_0699)*(1+s_0991/Km0991_0699)+(1+s_0180/Km0180_0699)*(1+s_1021/Km1021_0699)-1);
	r_0713=cell*Vmax_0713*(s_0066*s_1198-s_1203*s_1271/Keq_0713)/(Km0066_0713*Km1198_0713)/((1+s_0066/Km0066_0713)*(1+s_1198/Km1198_0713)+(1+s_1203/Km1203_0713)*(1+s_1271/Km1271_0713)-1);
	r_0722=cell*Vmax_0722*(s_0573*s_0785-s_0633*s_0743/Keq_0722)/(Km0573_0722*Km0785_0722)/((1+s_0573/Km0573_0722)*(1+s_0785/Km0785_0722)+(1+s_0633/Km0633_0722)*(1+s_0743/Km0743_0722)-1);
	r_0723=cell*Vmax_0723*(s_0557-s_0574/Keq_0723)/Km0557_0723/(1+s_0557/Km0557_0723+1+s_0574/Km0574_0723-1);
	r_0724=cell*Vmax_0724*(s_0304-s_0120/Keq_0724)/Km0304_0724/(1+s_0304/Km0304_0724+1+s_0120/Km0120_0724-1);
	r_0726=cell*Vmax_0726*(s_0434*s_1029-s_0633*s_1322*s_1416/Keq_0726)/(Km0434_0726*Km1029_0726)/((1+s_0434/Km0434_0726)*(1+s_1029/Km1029_0726)+(1+s_0633/Km0633_0726)*(1+s_1322/Km1322_0726)*(1+s_1416/Km1416_0726)-1);
	r_0727=cell*Vmax_0727*(s_0322*s_1012-s_1029*s_1487/Keq_0727)/(Km0322_0727*Km1012_0727)/((1+s_0322/Km0322_0727)*(1+s_1012/Km1012_0727)+(1+s_1029/Km1029_0727)*(1+s_1487/Km1487_0727)-1);
	r_0731=cell*Vmax_0731*(s_0306*s_1198-s_0304*s_1203/Keq_0731)/(Km0306_0731*Km1198_0731)/((1+s_0306/Km0306_0731)*(1+s_1198/Km1198_0731)+(1+s_0304/Km0304_0731)*(1+s_1203/Km1203_0731)-1);
	r_0732=cell*Vmax_0732*(s_0306*s_1207-s_0304*s_1212/Keq_0732)/(Km0306_0732*Km1207_0732)/((1+s_0306/Km0306_0732)*(1+s_1207/Km1207_0732)+(1+s_0304/Km0304_0732)*(1+s_1212/Km1212_0732)-1);
	r_0736=cell*Vmax_0736*(s_0028*s_0539-s_0019*s_0467/Keq_0736)/(Km0028_0736*Km0539_0736)/((1+s_0028/Km0028_0736)*(1+s_0539/Km0539_0736)+(1+s_0019/Km0019_0736)*(1+s_0467/Km0467_0736)-1);
	r_0739=cell*Vmax_0739*(s_0018*s_0434-s_0394*s_0456*s_0943*s_1322/Keq_0739)/(Km0018_0739*Km0434_0739)/((1+s_0018/Km0018_0739)*(1+s_0434/Km0434_0739)+(1+s_0394/Km0394_0739)*(1+s_0456/Km0456_0739)*(1+s_0943/Km0943_0739)*(1+s_1322/Km1322_0739)-1);
	r_0757=cell*Vmax_0757*(s_0126-s_1153*s_1322/Keq_0757)/Km0126_0757/(1+s_0126/Km0126_0757+(1+s_1153/Km1153_0757)*(1+s_1322/Km1322_0757)-1);
	r_0758=cell*Vmax_0758*(s_0568-s_0126/Keq_0758)/Km0568_0758/(1+s_0568/Km0568_0758+1+s_0126/Km0126_0758-1);
	r_0759=cell*Vmax_0759*(s_1191*s_1212-s_0145*s_1207*s_1322/Keq_0759)/(Km1191_0759*Km1212_0759)/((1+s_1191/Km1191_0759)*(1+s_1212/Km1212_0759)+(1+s_0145/Km0145_0759)*(1+s_1207/Km1207_0759)*(1+s_1322/Km1322_0759)-1);
	r_0762=cell*Vmax_0762*(s_1195-s_0722*s_1020/Keq_0762)/Km1195_0762/(1+s_1195/Km1195_0762+(1+s_0722/Km0722_0762)*(1+s_1020/Km1020_0762)-1);
	r_0770=cell*Vmax_0770*(s_1203*s_1537-s_1198*s_1535/Keq_0770)/(Km1203_0770*Km1537_0770)/((1+s_1203/Km1203_0770)*(1+s_1537/Km1537_0770)+(1+s_1198/Km1198_0770)*(1+s_1535/Km1535_0770)-1);
	r_0792=cell*Vmax_0792*(s_0467-s_0526*s_1322/Keq_0792)/Km0467_0792/(1+s_0467/Km0467_0792+(1+s_0526/Km0526_0792)*(1+s_1322/Km1322_0792)-1);
	r_0800=cell*Vmax_0800*(s_0434*s_0739-s_0394*s_0785/Keq_0800)/(Km0434_0800*Km0739_0800)/((1+s_0434/Km0434_0800)*(1+s_0739/Km0739_0800)+(1+s_0394/Km0394_0800)*(1+s_0785/Km0785_0800)-1);
	r_0806=cell*Vmax_0806*(s_0539-s_0467*s_1322/Keq_0806)/Km0539_0806/(1+s_0539/Km0539_0806+(1+s_0467/Km0467_0806)*(1+s_1322/Km1322_0806)-1);
	r_0811=cell*Vmax_0811*(s_0434*s_1538-s_0394*s_1559/Keq_0811)/(Km0434_0811*Km1538_0811)/((1+s_0434/Km0434_0811)*(1+s_1538/Km1538_0811)+(1+s_0394/Km0394_0811)*(1+s_1559/Km1559_0811)-1);
	r_0813=cell*Vmax_0813*(s_0841*s_1233-s_0362*s_1012/Keq_0813)/(Km0841_0813*Km1233_0813)/((1+s_0841/Km0841_0813)*(1+s_1233/Km1233_0813)+(1+s_0362/Km0362_0813)*(1+s_1012/Km1012_0813)-1);
	r_0816=cell*Vmax_0816*(s_0455*s_1266-s_0979*s_1322/Keq_0816)/(Km0455_0816*Km1266_0816)/((1+s_0455/Km0455_0816)*(1+s_1266/Km1266_0816)+(1+s_0979/Km0979_0816)*(1+s_1322/Km1322_0816)-1);
	r_0818=cell*Vmax_0818*(s_0991*s_1182-s_1192*s_1266/Keq_0818)/(Km0991_0818*Km1182_0818)/((1+s_0991/Km0991_0818)*(1+s_1182/Km1182_0818)+(1+s_1192/Km1192_0818)*(1+s_1266/Km1266_0818)-1);
	r_0820=cell*Vmax_0820*(s_1269*s_1386-s_0633*s_1270/Keq_0820)/(Km1269_0820*Km1386_0820)/((1+s_1269/Km1269_0820)*(1+s_1386/Km1386_0820)+(1+s_0633/Km0633_0820)*(1+s_1270/Km1270_0820)-1);
	r_0821=cell*Vmax_0821*(s_1270-s_0456*s_1545/Keq_0821)/Km1270_0821/(1+s_1270/Km1270_0821+(1+s_0456/Km0456_0821)*(1+s_1545/Km1545_0821)-1);
	r_0851=cell*Vmax_0851*(s_0951*s_0991-s_0180*s_1032/Keq_0851)/(Km0951_0851*Km0991_0851)/((1+s_0951/Km0951_0851)*(1+s_0991/Km0991_0851)+(1+s_0180/Km0180_0851)*(1+s_1032/Km1032_0851)-1);
	r_0855=cell*Vmax_0855*(s_0302*s_0434-s_0300*s_0394*s_1322/Keq_0855)/(Km0302_0855*Km0434_0855)/((1+s_0302/Km0302_0855)*(1+s_0434/Km0434_0855)+(1+s_0300/Km0300_0855)*(1+s_0394/Km0394_0855)*(1+s_1322/Km1322_0855)-1);
	r_0858=cell*Vmax_0858*(s_1351*s_1416-s_1343*s_1413/Keq_0858)/(Km1351_0858*Km1416_0858)/((1+s_1351/Km1351_0858)*(1+s_1416/Km1416_0858)+(1+s_1343/Km1343_0858)*(1+s_1413/Km1413_0858)-1);
	r_0874=cell*Vmax_0874*(s_0471*s_1153-s_0089*s_0526/Keq_0874)/(Km0471_0874*Km1153_0874)/((1+s_0471/Km0471_0874)*(1+s_1153/Km1153_0874)+(1+s_0089/Km0089_0874)*(1+s_0526/Km0526_0874)-1);
	r_0877=cell*Vmax_0877*(s_1337-s_0456*s_1351/Keq_0877)/Km1337_0877/(1+s_1337/Km1337_0877+(1+s_0456/Km0456_0877)*(1+s_1351/Km1351_0877)-1);
	r_0880=cell*Vmax_0880*(s_0471*s_1039-s_0526*s_1337/Keq_0880)/(Km0471_0880*Km1039_0880)/((1+s_0471/Km0471_0880)*(1+s_1039/Km1039_0880)+(1+s_0526/Km0526_0880)*(1+s_1337/Km1337_0880)-1);
	r_0883=cell*Vmax_0883*(s_0201*s_1616-s_0390*s_1469*s_1620/Keq_0883)/(Km0201_0883*Km1616_0883)/((1+s_0201/Km0201_0883)*(1+s_1616/Km1616_0883)+(1+s_0390/Km0390_0883)*(1+s_1469/Km1469_0883)*(1+s_1620/Km1620_0883)-1);
	r_0886=cell*Vmax_0886*(s_0434*s_0557-s_0394*s_0555/Keq_0886)/(Km0434_0886*Km0557_0886)/((1+s_0434/Km0434_0886)*(1+s_0557/Km0557_0886)+(1+s_0394/Km0394_0886)*(1+s_0555/Km0555_0886)-1);
	r_0887=cell*Vmax_0887*(s_0434*s_1427-s_0394*s_1426/Keq_0887)/(Km0434_0887*Km1427_0887)/((1+s_0434/Km0434_0887)*(1+s_1427/Km1427_0887)+(1+s_0394/Km0394_0887)*(1+s_1426/Km1426_0887)-1);
	r_0888=cell*Vmax_0888*(s_0568-s_0567/Keq_0888)/Km0568_0888/(1+s_0568/Km0568_0888+1+s_0567/Km0567_0888-1);
	r_0889=cell*Vmax_0889*(s_0340*s_1207-s_0456*s_0577*s_1212/Keq_0889)/(Km0340_0889*Km1207_0889)/((1+s_0340/Km0340_0889)*(1+s_1207/Km1207_0889)+(1+s_0456/Km0456_0889)*(1+s_0577/Km0577_0889)*(1+s_1212/Km1212_0889)-1);
	r_0891=cell*Vmax_0891*(s_0260*s_1198-s_0258*s_1203/Keq_0891)/(Km0260_0891*Km1198_0891)/((1+s_0260/Km0260_0891)*(1+s_1198/Km1198_0891)+(1+s_0258/Km0258_0891)*(1+s_1203/Km1203_0891)-1);
	r_0892=cell*Vmax_0892*(s_0075*s_0394-s_0260*s_0434/Keq_0892)/(Km0075_0892*Km0394_0892)/((1+s_0075/Km0075_0892)*(1+s_0394/Km0394_0892)+(1+s_0260/Km0260_0892)*(1+s_0434/Km0434_0892)-1);
	r_0893=cell*Vmax_0893*(s_0260-s_0188/Keq_0893)/Km0260_0893/(1+s_0260/Km0260_0893+1+s_0188/Km0188_0893-1);
	r_0900=cell*Vmax_0900*(s_1342*s_1416-s_1346*s_1413/Keq_0900)/(Km1342_0900*Km1416_0900)/((1+s_1342/Km1342_0900)*(1+s_1416/Km1416_0900)+(1+s_1346/Km1346_0900)*(1+s_1413/Km1413_0900)-1);
	r_0901=cell*Vmax_0901*(s_1343*s_1416-s_1342*s_1413/Keq_0901)/(Km1343_0901*Km1416_0901)/((1+s_1343/Km1343_0901)*(1+s_1416/Km1416_0901)+(1+s_1342/Km1342_0901)*(1+s_1413/Km1413_0901)-1);
	r_0902=cell*Vmax_0902*(s_0574-s_0573/Keq_0902)/Km0574_0902/(1+s_0574/Km0574_0902+1+s_0573/Km0573_0902-1);
	r_0904=cell*Vmax_0904*(s_0019*s_0434-s_0018*s_0394/Keq_0904)/(Km0019_0904*Km0434_0904)/((1+s_0019/Km0019_0904)*(1+s_0434/Km0434_0904)+(1+s_0018/Km0018_0904)*(1+s_0394/Km0394_0904)-1);
	r_0908=cell*Vmax_0908*(s_0434*s_0973*s_1364-s_0299*s_0394*s_1322/Keq_0908)/(Km0434_0908*Km0973_0908*Km1364_0908)/((1+s_0434/Km0434_0908)*(1+s_0973/Km0973_0908)*(1+s_1364/Km1364_0908)+(1+s_0299/Km0299_0908)*(1+s_0394/Km0394_0908)*(1+s_1322/Km1322_0908)-1);
	r_0909=cell*Vmax_0909*(s_0078-s_0077/Keq_0909)/Km0078_0909/(1+s_0078/Km0078_0909+1+s_0077/Km0077_0909-1);
	r_0910=cell*Vmax_0910*(s_0326-s_0078*s_0633/Keq_0910)/Km0326_0910/(1+s_0326/Km0326_0910+(1+s_0078/Km0078_0910)*(1+s_0633/Km0633_0910)-1);
	r_0911=cell*Vmax_0911*(s_0300*s_0456*s_0434-s_1364*s_0394*s_1322/Keq_0911)/(Km0300_0911*Km0456_0911*Km0434_0911)/((1+s_0300/Km0300_0911)*(1+s_0456/Km0456_0911)*(1+s_0434/Km0434_0911)+(1+s_1364/Km1364_0911)*(1+s_0394/Km0394_0911)*(1+s_1322/Km1322_0911)-1);
	r_0912=cell*Vmax_0912*(s_0120*s_0403-s_1365*s_1487/Keq_0912)/(Km0120_0912*Km0403_0912)/((1+s_0120/Km0120_0912)*(1+s_0403/Km0403_0912)+(1+s_1365/Km1365_0912)*(1+s_1487/Km1487_0912)-1);
	r_0913=cell*Vmax_0913*(s_1187-s_0076/Keq_0913)/Km1187_0913/(1+s_1187/Km1187_0913+1+s_0076/Km0076_0913-1);
	r_0914=cell*Vmax_0914*(s_0327*s_0434*s_1003-s_0325*s_0394*s_1322/Keq_0914)/(Km0327_0914*Km0434_0914*Km1003_0914)/((1+s_0327/Km0327_0914)*(1+s_0434/Km0434_0914)*(1+s_1003/Km1003_0914)+(1+s_0325/Km0325_0914)*(1+s_0394/Km0394_0914)*(1+s_1322/Km1322_0914)-1);
	r_0915=cell*Vmax_0915*(s_0999*s_1386-s_0327*s_0633*s_0991/Keq_0915)/(Km0999_0915*Km1386_0915)/((1+s_0999/Km0999_0915)*(1+s_1386/Km1386_0915)+(1+s_0327/Km0327_0915)*(1+s_0633/Km0633_0915)*(1+s_0991/Km0991_0915)-1);
	r_0916=cell*Vmax_0916*(s_0434*s_1408-s_0423*s_1386/Keq_0916)/(Km0434_0916*Km1408_0916)/((1+s_0434/Km0434_0916)*(1+s_1408/Km1408_0916)+(1+s_0423/Km0423_0916)*(1+s_1386/Km1386_0916)-1);
	r_0917=cell*Vmax_0917*(s_0259-s_1039*s_1322/Keq_0917)/Km0259_0917/(1+s_0259/Km0259_0917+(1+s_1039/Km1039_0917)*(1+s_1322/Km1322_0917)-1);
	r_0918=cell*Vmax_0918*(s_0258*s_0991-s_0180*s_0259/Keq_0918)/(Km0258_0918*Km0991_0918)/((1+s_0258/Km0258_0918)*(1+s_0991/Km0991_0918)+(1+s_0180/Km0180_0918)*(1+s_0259/Km0259_0918)-1);
	r_0919=cell*Vmax_0919*(s_1084*s_1366-s_0481/Keq_0919)/(Km1084_0919*Km1366_0919)/((1+s_1084/Km1084_0919)*(1+s_1366/Km1366_0919)+1+s_0481/Km0481_0919-1);
	r_0922=cell*Vmax_0922*(s_1212*s_1275*s_1445-s_1207*s_1366/Keq_0922)/(Km1212_0922*Km1275_0922*Km1445_0922)/((1+s_1212/Km1212_0922)*(1+s_1275/Km1275_0922)*(1+s_1445/Km1445_0922)+(1+s_1207/Km1207_0922)*(1+s_1366/Km1366_0922)-1);
	r_0938=cell*Vmax_0938*(s_1377-s_0456*s_0951/Keq_0938)/Km1377_0938/(1+s_1377/Km1377_0938+(1+s_0456/Km0456_0938)*(1+s_0951/Km0951_0938)-1);
	r_0939=cell*Vmax_0939*(s_1207*s_1377-s_0204*s_0456*s_1212/Keq_0939)/(Km1207_0939*Km1377_0939)/((1+s_1207/Km1207_0939)*(1+s_1377/Km1377_0939)+(1+s_0204/Km0204_0939)*(1+s_0456/Km0456_0939)*(1+s_1212/Km1212_0939)-1);
	r_0957=cell*Vmax_0957*(s_0118*s_1212-s_1035*s_1207/Keq_0957)/(Km0118_0957*Km1212_0957)/((1+s_0118/Km0118_0957)*(1+s_1212/Km1212_0957)+(1+s_1035/Km1035_0957)*(1+s_1207/Km1207_0957)-1);
	r_0958=cell*Vmax_0958*(s_0434*s_0445*s_1399-s_0394*s_1271*s_1322/Keq_0958)/(Km0434_0958*Km0445_0958*Km1399_0958)/((1+s_0434/Km0434_0958)*(1+s_0445/Km0445_0958)*(1+s_1399/Km1399_0958)+(1+s_0394/Km0394_0958)*(1+s_1271/Km1271_0958)*(1+s_1322/Km1322_0958)-1);
	r_0959=cell*Vmax_0959*(s_1399-s_0359*s_0456/Keq_0959)/Km1399_0959/(1+s_1399/Km1399_0959+(1+s_0359/Km0359_0959)*(1+s_0456/Km0456_0959)-1);
	r_0961=cell*Vmax_0961*(s_0529*s_1198*s_1399-s_0373*s_0456*s_1203/Keq_0961)/(Km0529_0961*Km1198_0961*Km1399_0961)/((1+s_0529/Km0529_0961)*(1+s_1198/Km1198_0961)*(1+s_1399/Km1399_0961)+(1+s_0373/Km0373_0961)*(1+s_0456/Km0456_0961)*(1+s_1203/Km1203_0961)-1);
	r_0962=cell*Vmax_0962*(s_0394*s_1360-s_0434*s_1399/Keq_0962)/(Km0394_0962*Km1360_0962)/((1+s_0394/Km0394_0962)*(1+s_1360/Km1360_0962)+(1+s_0434/Km0434_0962)*(1+s_1399/Km1399_0962)-1);
	r_0967=cell*Vmax_0967*(s_0158*s_0314-s_0328*s_1322/Keq_0967)/(Km0158_0967*Km0314_0967)/((1+s_0158/Km0158_0967)*(1+s_0314/Km0314_0967)+(1+s_0328/Km0328_0967)*(1+s_1322/Km1322_0967)-1);
	r_0968=cell*Vmax_0968*(pow(s_0328,2)-s_0314*s_1405/Keq_0968)/pow(Km0328_0968,2)/(pow(1+s_0328/Km0328_0968,2)+(1+s_0314/Km0314_0968)*(1+s_1405/Km1405_0968)-1);
	r_0970=cell*Vmax_0970*(s_0434*s_1616-s_0586*s_1620/Keq_0970)/(Km0434_0970*Km1616_0970)/((1+s_0434/Km0434_0970)*(1+s_1616/Km1616_0970)+(1+s_0586/Km0586_0970)*(1+s_1620/Km1620_0970)-1);
	r_0973=cell*Vmax_0973*(s_1559*s_1616-s_0656*s_1620/Keq_0973)/(Km1559_0973*Km1616_0973)/((1+s_1559/Km1559_0973)*(1+s_1616/Km1616_0973)+(1+s_0656/Km0656_0973)*(1+s_1620/Km1620_0973)-1);
	r_0974=cell*Vmax_0974*(s_0394*s_1616-s_0582*s_1620/Keq_0974)/(Km0394_0974*Km1616_0974)/((1+s_0394/Km0394_0974)*(1+s_1616/Km1616_0974)+(1+s_0582/Km0582_0974)*(1+s_1620/Km1620_0974)-1);
	r_0976=cell*Vmax_0976*(s_0467*s_1616-s_0587*s_1620/Keq_0976)/(Km0467_0976*Km1616_0976)/((1+s_0467/Km0467_0976)*(1+s_1616/Km1616_0976)+(1+s_0587/Km0587_0976)*(1+s_1620/Km1620_0976)-1);
	r_0978=cell*Vmax_0978*(s_0739*s_1616-s_0613*s_1620/Keq_0978)/(Km0739_0978*Km1616_0978)/((1+s_0739/Km0739_0978)*(1+s_1616/Km1616_0978)+(1+s_0613/Km0613_0978)*(1+s_1620/Km1620_0978)-1);
	r_0982=cell*Vmax_0982*(s_0577-s_1408/Keq_0982)/Km0577_0982/(1+s_0577/Km0577_0982+1+s_1408/Km1408_0982-1);
	r_0984=cell*Vmax_0984*(s_0577-s_0581/Keq_0984)/Km0577_0984/(1+s_0577/Km0577_0984+1+s_0581/Km0581_0984-1);
	r_0986=cell*Vmax_0986*(s_1416*s_1569-s_0700*s_1413/Keq_0986)/(Km1416_0986*Km1569_0986)/((1+s_1416/Km1416_0986)*(1+s_1569/Km1569_0986)+(1+s_0700/Km0700_0986)*(1+s_1413/Km1413_0986)-1);
	r_0988=cell*Vmax_0988*(s_1038*s_1198-s_0180*s_1025*s_1203/Keq_0988)/(Km1038_0988*Km1198_0988)/((1+s_1038/Km1038_0988)*(1+s_1198/Km1198_0988)+(1+s_0180/Km0180_0988)*(1+s_1025/Km1025_0988)*(1+s_1203/Km1203_0988)-1);
	r_0989=cell*Vmax_0989*(s_0959*s_0991*s_1212-s_1038*s_1207/Keq_0989)/(Km0959_0989*Km0991_0989*Km1212_0989)/((1+s_0959/Km0959_0989)*(1+s_0991/Km0991_0989)*(1+s_1212/Km1212_0989)+(1+s_1038/Km1038_0989)*(1+s_1207/Km1207_0989)-1);
	r_0990=cell*Vmax_0990*(s_1426-s_0551*s_0629/Keq_0990)/Km1426_0990/(1+s_1426/Km1426_0990+(1+s_0551/Km0551_0990)*(1+s_0629/Km0629_0990)-1);
	r_0992=cell*Vmax_0992*(s_0373*s_1039-s_0529*s_1234/Keq_0992)/(Km0373_0992*Km1039_0992)/((1+s_0373/Km0373_0992)*(1+s_1039/Km1039_0992)+(1+s_0529/Km0529_0992)*(1+s_1234/Km1234_0992)-1);
	r_0993=cell*Vmax_0993*(s_1039*s_1302-s_0231*s_0456*s_0529/Keq_0993)/(Km1039_0993*Km1302_0993)/((1+s_1039/Km1039_0993)*(1+s_1302/Km1302_0993)+(1+s_0231/Km0231_0993)*(1+s_0456/Km0456_0993)*(1+s_0529/Km0529_0993)-1);
	r_0996=cell*Vmax_0996*(s_0211*s_1212-s_1207*s_1429/Keq_0996)/(Km0211_0996*Km1212_0996)/((1+s_0211/Km0211_0996)*(1+s_1212/Km1212_0996)+(1+s_1207/Km1207_0996)*(1+s_1429/Km1429_0996)-1);
	r_0997=cell*Vmax_0997*(s_0434*s_1429-s_0261*s_0394/Keq_0997)/(Km0434_0997*Km1429_0997)/((1+s_0434/Km0434_0997)*(1+s_1429/Km1429_0997)+(1+s_0261/Km0261_0997)*(1+s_0394/Km0394_0997)-1);
	r_1010=cell*Vmax_1010*(s_1203*s_1275*s_1447-s_0037*s_1198/Keq_1010)/(Km1203_1010*Km1275_1010*Km1447_1010)/((1+s_1203/Km1203_1010)*(1+s_1275/Km1275_1010)*(1+s_1447/Km1447_1010)+(1+s_0037/Km0037_1010)*(1+s_1198/Km1198_1010)-1);
	r_1011=cell*Vmax_1011*(s_1212*s_1275*s_1447-s_0037*s_1207/Keq_1011)/(Km1212_1011*Km1275_1011*Km1447_1011)/((1+s_1212/Km1212_1011)*(1+s_1275/Km1275_1011)*(1+s_1447/Km1447_1011)+(1+s_0037/Km0037_1011)*(1+s_1207/Km1207_1011)-1);
	r_1012=cell*Vmax_1012*(pow(s_0190,2)*s_1212-pow(s_0633,2)*s_1207*s_1447/Keq_1012)/(pow(Km0190_1012,2)*Km1212_1012)/(pow(1+s_0190/Km0190_1012,2)*(1+s_1212/Km1212_1012)+pow(1+s_0633/Km0633_1012,2)*(1+s_1207/Km1207_1012)*(1+s_1447/Km1447_1012)-1);
	r_1014=cell*Vmax_1014*(s_0666*s_0595-s_0672/Keq_1014)/(Km0666_1014*Km0595_1014)/((1+s_0666/Km0666_1014)*(1+s_0595/Km0595_1014)+1+s_0672/Km0672_1014-1);
	r_1026=cell*Vmax_1026*(s_0394*s_1467-s_0298*s_1322/Keq_1026)/(Km0394_1026*Km1467_1026)/((1+s_0394/Km0394_1026)*(1+s_1467/Km1467_1026)+(1+s_0298/Km0298_1026)*(1+s_1322/Km1322_1026)-1);
	r_1027=cell*Vmax_1027*(pow(s_1212,3)*s_1469-s_0841*pow(s_1207,3)/Keq_1027)/(pow(Km1212_1027,3)*Km1469_1027)/(pow(1+s_1212/Km1212_1027,3)*(1+s_1469/Km1469_1027)+(1+s_0841/Km0841_1027)*pow(1+s_1207/Km1207_1027,3)-1);
	r_1038=cell*Vmax_1038*(s_1212*s_1620-s_1207*s_1616/Keq_1038)/(Km1212_1038*Km1620_1038)/((1+s_1212/Km1212_1038)*(1+s_1620/Km1620_1038)+(1+s_1207/Km1207_1038)*(1+s_1616/Km1616_1038)-1);
	r_1041=cell*Vmax_1041*(s_1238-s_1045*s_1322/Keq_1041)/Km1238_1041/(1+s_1238/Km1238_1041+(1+s_1045/Km1045_1041)*(1+s_1322/Km1322_1041)-1);
	r_1045=cell*Vmax_1045*(s_0306*s_0654-s_0625*s_0649/Keq_1045)/(Km0306_1045*Km0654_1045)/((1+s_0306/Km0306_1045)*(1+s_0654/Km0654_1045)+(1+s_0625/Km0625_1045)*(1+s_0649/Km0649_1045)-1);
	r_1049=cell*Vmax_1049*(s_0581*s_1408-s_0764*s_1427/Keq_1049)/(Km0581_1049*Km1408_1049)/((1+s_0581/Km0581_1049)*(1+s_1408/Km1408_1049)+(1+s_0764/Km0764_1049)*(1+s_1427/Km1427_1049)-1);
	r_1050=cell*Vmax_1050*(s_0551*s_0581-s_0557*s_0764/Keq_1050)/(Km0551_1050*Km0581_1050)/((1+s_0551/Km0551_1050)*(1+s_0581/Km0581_1050)+(1+s_0557/Km0557_1050)*(1+s_0764/Km0764_1050)-1);
	r_1051=cell*Vmax_1051*(s_0409-s_1322*s_1520/Keq_1051)/Km0409_1051/(1+s_0409/Km0409_1051+(1+s_1322/Km1322_1051)*(1+s_1520/Km1520_1051)-1);
	r_1052=cell*Vmax_1052*(s_0619*s_0595-s_1524/Keq_1052)/(Km0619_1052*Km0595_1052)/((1+s_0619/Km0619_1052)*(1+s_0595/Km0595_1052)+1+s_1524/Km1524_1052-1);
	r_1054=cell*Vmax_1054*(s_0764-s_0629/Keq_1054)/Km0764_1054/(1+s_0764/Km0764_1054+1+s_0629/Km0629_1054-1);
	r_1055=cell*Vmax_1055*(s_0086*s_1039-s_0764*s_1048/Keq_1055)/(Km0086_1055*Km1039_1055)/((1+s_0086/Km0086_1055)*(1+s_1039/Km1039_1055)+(1+s_0764/Km0764_1055)*(1+s_1048/Km1048_1055)-1);
	r_1063=cell*Vmax_1063*(s_0204*s_0991-s_0180*s_1051/Keq_1063)/(Km0204_1063*Km0991_1063)/((1+s_0204/Km0204_1063)*(1+s_0991/Km0991_1063)+(1+s_0180/Km0180_1063)*(1+s_1051/Km1051_1063)-1);
	r_1072=cell*Vmax_1072*(s_0434*s_1545-s_0394*s_1538/Keq_1072)/(Km0434_1072*Km1545_1072)/((1+s_0434/Km0434_1072)*(1+s_1545/Km1545_1072)+(1+s_0394/Km0394_1072)*(1+s_1538/Km1538_1072)-1);
	r_1084=cell*Vmax_1084*(s_0567*s_1559-s_0633*s_1543/Keq_1084)/(Km0567_1084*Km1559_1084)/((1+s_0567/Km0567_1084)*(1+s_1559/Km1559_1084)+(1+s_0633/Km0633_1084)*(1+s_1543/Km1543_1084)-1);
	r_1087=cell*Vmax_1087*(s_0232*s_0991-s_0180*s_1056/Keq_1087)/(Km0232_1087*Km0991_1087)/((1+s_0232/Km0232_1087)*(1+s_0991/Km0991_1087)+(1+s_0180/Km0180_1087)*(1+s_1056/Km1056_1087)-1);
	r_1106=cell*Vmax_1106*s_0362/Km0362_1106/(1+s_0362/Km0362_1106);
	r_1115=cell*Vmax_1115*(s_0420-s_0419)/Km0420_1115/(1+s_0420/Km0420_1115+1+s_0419/Km0419_1115-1);
	r_1166=cell*Vmax_1166*(s_0565-s_0563)/Km0565_1166/(1+s_0565/Km0565_1166+1+s_0563/Km0563_1166-1);
	r_1172=cell*Vmax_1172*s_0765/Km0765_1172/(1+s_0765/Km0765_1172);
	r_1244=cell*Vmax_1244*(s_1324-s_1322)/Km1324_1244/(1+s_1324/Km1324_1244+1+s_1322/Km1322_1244-1);
	r_1266=cell*Vmax_1266*(s_1468-s_1467)/Km1468_1266/(1+s_1468/Km1468_1266+1+s_1467/Km1467_1266-1);
	r_1664=cell*Vmax_1664*(s_0456-s_0445/Keq_1664)/Km0456_1664/(1+s_0456/Km0456_1664+1+s_0445/Km0445_1664-1);
	r_1697=cell*Vmax_1697*s_0456/Km0456_1697/(1+s_0456/Km0456_1697);
	r_1704=cell*Vmax_1704*(s_0394*s_0587-s_0434*s_0589/Keq_1704)/(Km0394_1704*Km0587_1704)/((1+s_0394/Km0394_1704)*(1+s_0587/Km0587_1704)+(1+s_0434/Km0434_1704)*(1+s_0589/Km0589_1704)-1);
	r_1729=cell*Vmax_1729*(s_0394*s_0582-s_0434*s_0584/Keq_1729)/(Km0394_1729*Km0582_1729)/((1+s_0394/Km0394_1729)*(1+s_0582/Km0582_1729)+(1+s_0434/Km0434_1729)*(1+s_0584/Km0584_1729)-1);
	r_1762=cell*Vmax_1762*s_0680/Km0680_1762/(1+s_0680/Km0680_1762);
	r_1936=cell*Vmax_1936*(s_0629-s_1151*s_1322/Keq_1936)/Km0629_1936/(1+s_0629/Km0629_1936+(1+s_1151/Km1151_1936)*(1+s_1322/Km1322_1936)-1);
	r_1979=cell*Vmax_1979*(s_1277-s_1275)/Km1277_1979/(1+s_1277/Km1277_1979+1+s_1275/Km1275_1979-1);
	r_2030=cell*Vmax_2030*(s_0313-s_0314*s_1322/Keq_2030)/Km0313_2030/(1+s_0313/Km0313_2030+(1+s_0314/Km0314_2030)*(1+s_1322/Km1322_2030)-1);
	r_2079=cell*Vmax_2079*s_1520/Km1520_2079/(1+s_1520/Km1520_2079);
	r_2111=cell*fmax(V0_2111*(1+ep0002_2111*log(s_0002/ic0002_2111)+ep0423_2111*log(s_0423/ic0423_2111)+ep0434_2111*log(s_0434/ic0434_2111)+ep0526_2111*log(s_0526/ic0526_2111)+ep0584_2111*log(s_0584/ic0584_2111)+ep0589_2111*log(s_0589/ic0589_2111)+ep0615_2111*log(s_0615/ic0615_2111)+ep0649_2111*log(s_0649/ic0649_2111)+ep0773_2111*log(s_0773/ic0773_2111)+ep0782_2111*log(s_0782/ic0782_2111)+ep0955_2111*log(s_0955/ic0955_2111)+ep0965_2111*log(s_0965/ic0965_2111)+ep0969_2111*log(s_0969/ic0969_2111)+ep0973_2111*log(s_0973/ic0973_2111)+ep0981_2111*log(s_0981/ic0981_2111)+ep0991_2111*log(s_0991/ic0991_2111)+ep0999_2111*log(s_0999/ic0999_2111)+ep1003_2111*log(s_1003/ic1003_2111)+ep1006_2111*log(s_1006/ic1006_2111)+ep1016_2111*log(s_1016/ic1016_2111)+ep1021_2111*log(s_1021/ic1021_2111)+ep1025_2111*log(s_1025/ic1025_2111)+ep1029_2111*log(s_1029/ic1029_2111)+ep1032_2111*log(s_1032/ic1032_2111)+ep1035_2111*log(s_1035/ic1035_2111)+ep1039_2111*log(s_1039/ic1039_2111)+ep1045_2111*log(s_1045/ic1045_2111)+ep1048_2111*log(s_1048/ic1048_2111)+ep1051_2111*log(s_1051/ic1051_2111)+ep1056_2111*log(s_1056/ic1056_2111)+ep1107_2111*log(s_1107/ic1107_2111)+ep1405_2111*log(s_1405/ic1405_2111)+ep1467_2111*log(s_1467/ic1467_2111)+ep1520_2111*log(s_1520/ic1520_2111)+ep1545_2111*log(s_1545/ic1545_2111)+ep0089_2111*log(s_0089/ic0089_2111)+ep0122_2111*log(s_0122/ic0122_2111)+ep0918_2111*log(s_0918/ic0918_2111)+ep0657_2111*log(s_0657/ic0657_2111)+ep0662_2111*log(s_0662/ic0662_2111)+ep0666_2111*log(s_0666/ic0666_2111)+ep0672_2111*log(s_0672/ic0672_2111)+ep0595_2111*log(s_0595/ic0595_2111)+ep0700_2111*log(s_0700/ic0700_2111)+ep1059_2111*log(s_1059/ic1059_2111)+ep1337_2111*log(s_1337/ic1337_2111)+ep1346_2111*log(s_1346/ic1346_2111)+ep1351_2111*log(s_1351/ic1351_2111)+ep1524_2111*log(s_1524/ic1524_2111)+ep1569_2111*log(s_1569/ic1569_2111)),0);
	ds_0002=(r_0005-1.14*r_2111)/cell;
	ds_0008=(-r_0353+r_0669)/cell;
	ds_0009=(r_0060-r_0061)/cell;
	ds_0010=(-r_0029+r_0061)/cell;
	ds_0015=(-r_0207+r_0208)/cell;
	ds_0016=(r_0096-r_0352)/cell;
	ds_0018=(-r_0739+r_0904)/cell;
	ds_0019=(r_0736-r_0904)/cell;
	ds_0025=(-r_0001+r_0553)/cell;
	ds_0028=(r_0558-r_0736)/cell;
	ds_0033=(-r_0553+r_0697)/cell;
	ds_0037=(-r_0698+r_1010+r_1011)/cell;
	ds_0039=(r_0016-r_0669)/cell;
	ds_0056=(r_0353-r_0663)/cell;
	ds_0061=(-r_0339+r_0349)/cell;
	ds_0062=(r_0688-r_0696)/cell;
	ds_0063=(-r_0004+r_0696)/cell;
	ds_0066=(r_0451-r_0713)/cell;
	ds_0075=(r_0486-r_0892)/cell;
	ds_0076=(-r_0566+r_0913)/cell;
	ds_0077=(-r_0007+r_0909)/cell;
	ds_0078=(-r_0909+r_0910)/cell;
	ds_0082=(-r_0008+r_0495)/cell;
	ds_0086=(r_0566-r_1055)/cell;
	ds_0089=(-r_0594+r_0874-0.00153*r_2111)/cell;
	ds_0118=(r_0012-r_0957)/cell;
	ds_0120=(-r_0446-r_0499+r_0724-r_0912)/cell;
	ds_0122=(r_0231-r_0241-5.6e-05*r_2111)/cell;
	ds_0126=(-r_0757+r_0758)/cell;
	ds_0141=(-r_0015+r_0525)/cell;
	ds_0142=(-r_0014+r_0015)/cell;
	ds_0145=(-r_0118+r_0759)/cell;
	ds_0146=(-r_0096+r_0097)/cell;
	ds_0158=(r_0038-r_0967)/cell;
	ds_0162=(-r_0023+r_0024)/cell;
	ds_0165=(r_0023-r_0060)/cell;
	ds_0176=(-r_0018+r_0545)/cell;
	ds_0178=(-r_0016+r_0310)/cell;
	ds_0180=(r_0018+r_0118+r_0216-r_0470-r_0471+r_0538-r_0543+r_0658+r_0661+r_0663+r_0674+r_0699+r_0851+r_0918+r_0988+r_1063+r_1087)/cell;
	ds_0188=(-r_0366+r_0893)/cell;
	ds_0190=(r_0462-2.0*r_1012)/cell;
	ds_0201=(r_0154-r_0883)/cell;
	ds_0204=(r_0939-r_1063)/cell;
	ds_0207=(-r_0538+r_0564)/cell;
	ds_0209=(r_0235-r_0236)/cell;
	ds_0210=(-r_0039+r_0040)/cell;
	ds_0211=(r_0039-r_0996)/cell;
	ds_0218=(-r_0558+r_0559)/cell;
	ds_0231=(-r_0041+r_0993)/cell;
	ds_0232=(-r_0024+r_0352-r_1087)/cell;
	ds_0258=(r_0891-r_0918)/cell;
	ds_0259=(-r_0917+r_0918)/cell;
	ds_0260=(-r_0891+r_0892-r_0893)/cell;
	ds_0261=(-r_0065+r_0997)/cell;
	ds_0262=(-r_0231+r_0317)/cell;
	ds_0291=(r_0029-r_0699)/cell;
	ds_0295=(r_0215-r_0219)/cell;
	ds_0296=(r_0236-r_0238)/cell;
	ds_0297=(-r_0235+r_0241)/cell;
	ds_0298=(-r_0154+r_1026)/cell;
	ds_0299=(-r_0151+r_0908)/cell;
	ds_0300=(r_0855-r_0911)/cell;
	ds_0301=(-r_0079+r_0499)/cell;
	ds_0302=(r_0079-r_0855)/cell;
	ds_0304=(-r_0724+r_0731+r_0732)/cell;
	ds_0306=(-r_0080+r_0502-r_0731-r_0732-r_1045)/cell;
	ds_0312=(r_0007-r_0563)/cell;
	ds_0313=(r_0014-r_2030)/cell;
	ds_0314=(-r_0967+r_0968+r_2030)/cell;
	ds_0322=(r_0080-r_0727)/cell;
	ds_0324=(r_0065-r_0279)/cell;
	ds_0325=(-r_0499+r_0914)/cell;
	ds_0326=(r_0225-r_0910)/cell;
	ds_0327=(-r_0914+r_0915)/cell;
	ds_0328=(r_0967-2.0*r_0968)/cell;
	ds_0335=(-r_0091+r_0466)/cell;
	ds_0340=(r_0091-r_0889)/cell;
	ds_0349=(r_0020-r_0040)/cell;
	ds_0359=(-r_0165-r_0173+r_0959)/cell;
	ds_0362=(r_0111+r_0173+r_0311+r_0312+r_0813-r_1106)/cell;
	ds_0367=(r_0103-r_0559)/cell;
	ds_0373=(-r_0024-2.0*r_0103-r_0108-r_0111-r_0300-r_0398-r_0543-r_0549-r_0559+r_0961-r_0992)/cell;
	ds_0380=(-r_0008+r_0336-r_0495)/cell;
	ds_0386=(-r_0142+r_0144)/cell;
	ds_0390=(-r_0032+r_0883)/cell;
	ds_0393=(-r_0152+r_0153)/cell;
	ds_0394=(r_0079+r_0108+r_0115+r_0142+2.0*r_0148+r_0154+r_0215-r_0226+2.0*r_0250+r_0307-r_0330-r_0446+r_0476+r_0528+r_0534+r_0548+r_0739+r_0800+r_0811+r_0855+r_0886+r_0887-r_0892+r_0904+r_0908+r_0911+r_0914+r_0958-r_0962-r_0974+r_0997-r_1026+r_1072-r_1704-r_1729+59.3*r_2111)/cell;
	ds_0403=(r_0151+r_0563-r_0912)/cell;
	ds_0409=(r_0195-r_1051)/cell;
	ds_0419=(r_0014-r_0307+r_0310-r_0326-r_0470-r_0471-r_0476+r_1115)/cell;
	ds_0423=(r_0032+r_0142-r_0148+r_0152+r_0208+r_0211-r_0399-r_0407+r_0514+r_0916-0.051*r_2111)/cell;
	ds_0427=(-r_0202+r_0203+r_0670)/cell;
	ds_0434=(-r_0079-r_0108-r_0115-r_0142-r_0148-r_0154-r_0208-r_0211-r_0215-r_0225+r_0226-2.0*r_0250-r_0307+r_0330+r_0399+r_0407+r_0446-r_0476-r_0514-r_0528-r_0534-r_0548-r_0726-r_0739-r_0800-r_0811-r_0855-r_0886-r_0887+r_0892-r_0904-r_0908-r_0911-r_0914-r_0916-r_0958+r_0962-r_0970-r_0997-r_1072+r_1704+r_1729-59.3*r_2111)/cell;
	ds_0445=(-r_0108-r_0250-r_0958+r_1664)/cell;
	ds_0454=(r_0027-r_0542)/cell;
	ds_0455=(-r_0214+r_0250-r_0816)/cell;
	ds_0456=(r_0016+r_0029+r_0097+r_0234+r_0235+r_0386+r_0387+r_0389+r_0391+3.0*r_0393+r_0397+3.0*r_0398+r_0432+r_0433+r_0434+r_0435+r_0445+r_0545+r_0566+r_0658+r_0661+r_0739+r_0821+r_0877+r_0889-r_0911+r_0938+r_0939+r_0959+r_0961+r_0993-r_1664-r_1697)/cell;
	ds_0467=(r_0736-r_0792+r_0806-r_0976)/cell;
	ds_0471=(r_0257-r_0874-r_0880)/cell;
	ds_0475=(-r_0259+r_0340)/cell;
	ds_0481=(r_0259-r_0267+r_0919)/cell;
	ds_0493=(r_0267-r_0269)/cell;
	ds_0499=(r_0269-r_0594)/cell;
	ds_0515=(-r_0203-r_0278+r_0279)/cell;
	ds_0516=(-r_0280+r_0302)/cell;
	ds_0522=(r_0300-r_0302)/cell;
	ds_0526=(r_0792+r_0874+r_0880-0.05*r_2111)/cell;
	ds_0529=(r_0008+r_0024+r_0103+r_0111+r_0300-r_0336+r_0386+r_0387+r_0389+r_0391+3.0*r_0393+r_0397+3.0*r_0398+r_0399+r_0407+r_0432+r_0433+r_0434+r_0435+r_0495+r_0543+r_0549+r_0558+r_0559-r_0961+r_0992+r_0993)/cell;
	ds_0539=(-r_0257+r_0307-r_0736-r_0806)/cell;
	ds_0550=(r_0563-r_0564)/cell;
	ds_0551=(-r_0020+r_0990-r_1050)/cell;
	ds_0555=(-r_0450+r_0886)/cell;
	ds_0557=(r_0467-r_0723-r_0886+r_1050)/cell;
	ds_0563=(-r_0534+r_1166)/cell;
	ds_0567=(r_0888-r_1084)/cell;
	ds_0568=(-r_0195-r_0466-r_0467+r_0534-r_0758-r_0888)/cell;
	ds_0573=(-r_0722+r_0902)/cell;
	ds_0574=(r_0723-r_0902)/cell;
	ds_0577=(-r_0038+r_0889-r_0982-r_0984)/cell;
	ds_0581=(r_0984-r_1049-r_1050)/cell;
	ds_0582=(r_0529+r_0974-r_1729)/cell;
	ds_0584=(r_1729-0.00359*r_2111)/cell;
	ds_0586=(-r_0529+r_0970)/cell;
	ds_0587=(r_0976-r_1704)/cell;
	ds_0589=(r_0326+r_1704-0.00243*r_2111)/cell;
	ds_0595=(-r_0386+r_0399-1.8*r_1014-2.6*r_1052-0.0005356*r_2111)/cell;
	ds_0602=(r_0397-r_0399-r_0432)/cell;
	ds_0613=(-r_0330+r_0978)/cell;
	ds_0615=(r_0330-0.00243*r_2111)/cell;
	ds_0619=(r_0336+r_0337+r_0594-r_1052)/cell;
	ds_0625=(-r_0344+r_1045)/cell;
	ds_0629=(r_0450-r_0491+r_0990+r_1054-r_1936)/cell;
	ds_0633=(r_0202+r_0208+r_0211+r_0225+r_0257+r_0355+r_0364-r_0399-r_0407+r_0462+r_0514+r_0525-r_0568+r_0722+r_0726+r_0820+r_0910+r_0915+2.0*r_1012+r_1084)/cell;
	ds_0644=(r_0361-r_0362)/cell;
	ds_0645=(-r_0361+r_0362)/cell;
	ds_0649=(r_1045-0.00359*r_2111)/cell;
	ds_0654=(-r_0326+r_0364-r_1045)/cell;
	ds_0656=(-r_0364+r_0973)/cell;
	ds_0657=(-r_0242+r_0243-9.6e-05*r_2111)/cell;
	ds_0662=(r_0233-r_0244-0.000125*r_2111)/cell;
	ds_0664=(-r_0233+r_0242)/cell;
	ds_0666=(r_0244-r_1014-0.0056*r_2111)/cell;
	ds_0672=(r_1014-0.000812*r_2111)/cell;
	ds_0680=(r_0165-r_1762)/cell;
	ds_0700=(-r_0243+r_0986-0.000114*r_2111)/cell;
	ds_0709=(-2.0*r_0001-2.0*r_0004+4.0*r_0438-2.0*r_0439)/cell;
	ds_0710=(2.0*r_0001+2.0*r_0004-4.0*r_0438+2.0*r_0439)/cell;
	ds_0722=(r_0038+r_0317-r_0445+r_0446+r_0525+r_0762)/cell;
	ds_0725=(r_0151+r_0152+r_0207-r_0451)/cell;
	ds_0739=(r_0153+r_0361+r_0528+r_0529-r_0800-r_0978)/cell;
	ds_0743=(-r_0361+r_0722)/cell;
	ds_0745=(r_0355-r_0462)/cell;
	ds_0750=(2.0*r_0481-2.0*r_0483+r_0553-r_0697)/cell;
	ds_0754=(-r_0481+r_0483)/cell;
	ds_0764=(r_0450-r_0486+r_1049+r_1050-r_1054+r_1055)/cell;
	ds_0765=(r_0489-r_1172)/cell;
	ds_0767=(-r_0489+r_0491-r_0495)/cell;
	ds_0773=(r_0510-0.519*r_2111)/cell;
	ds_0782=(r_0514-r_0528-r_0529-0.051*r_2111)/cell;
	ds_0785=(-r_0153-r_0525-r_0722+r_0800)/cell;
	ds_0835=(-r_0027+r_0543)/cell;
	ds_0836=(r_0542-r_0545)/cell;
	ds_0837=(r_0339-r_0483-r_0550)/cell;
	ds_0841=(-r_0312-r_0813+r_1027)/cell;
	ds_0849=(-r_0153-r_0565+r_0570)/cell;
	ds_0918=(r_0594-0.000538625*r_2111)/cell;
	ds_0940=(r_0280-r_0658-r_0661)/cell;
	ds_0943=(-r_0355-r_0462-r_0667+r_0739)/cell;
	ds_0951=(-r_0851+r_0938)/cell;
	ds_0953=(r_0018-r_0678)/cell;
	ds_0955=(r_0670+r_0674-0.357*r_2111)/cell;
	ds_0959=(r_0678-r_0989)/cell;
	ds_0965=(r_0207-0.136*r_2111)/cell;
	ds_0969=(r_0211-0.172*r_2111)/cell;
	ds_0973=(-r_0153-r_0208-r_0211-r_0214-r_0215+r_0216-r_0908-0.172*r_2111)/cell;
	ds_0978=(r_0219-r_0547)/cell;
	ds_0979=(-r_0208+r_0816)/cell;
	ds_0980=(r_0309-r_0310+r_0311)/cell;
	ds_0981=(r_0310-r_0311+r_0312-0.0429*r_2111)/cell;
	ds_0991=(-r_0012-r_0018+r_0079-r_0118+r_0203+r_0211-r_0216+r_0250+r_0470+r_0471-r_0476+r_0514-r_0538+r_0563-r_0663-r_0674-r_0699-r_0818-r_0851+r_0915-r_0918-r_0989-r_1063-r_1087-0.268*r_2111)/cell;
	ds_0999=(-r_0079-r_0203-r_0211-r_0250+r_0476-r_0514-r_0563-r_0915-0.268*r_2111)/cell;
	ds_1003=(r_0502-r_0914-0.325*r_2111)/cell;
	ds_1006=(r_0536-0.075*r_2111)/cell;
	ds_1010=(-r_0536+r_0537)/cell;
	ds_1011=(-r_0537+r_0538)/cell;
	ds_1012=(r_0144-r_0309-r_0727+r_0813)/cell;
	ds_1014=(r_0547-r_0548-r_0549)/cell;
	ds_1016=(r_0663-0.172*r_2111)/cell;
	ds_1020=(-r_0670+r_0762)/cell;
	ds_1021=(r_0699-0.25*r_2111)/cell;
	ds_1025=(r_0988-0.239*r_2111)/cell;
	ds_1029=(-r_0726+r_0727-0.05*r_2111)/cell;
	ds_1032=(r_0851-0.114*r_2111)/cell;
	ds_1035=(r_0957-0.129*r_2111)/cell;
	ds_1038=(-r_0988+r_0989)/cell;
	ds_1039=(-r_0309-r_0502-r_0880+r_0917-r_0992-r_0993-r_1055-0.254*r_2111)/cell;
	ds_1045=(r_1041-0.197*r_2111)/cell;
	ds_1048=(-r_0694+r_1055-0.028*r_2111)/cell;
	ds_1051=(r_1063-0.0965*r_2111)/cell;
	ds_1056=(r_1087-0.257*r_2111)/cell;
	ds_1059=(-r_0317+r_0698-3.2e-05*r_2111)/cell;
	ds_1065=(r_0386-r_0387)/cell;
	ds_1073=(r_0432-r_0433)/cell;
	ds_1084=(-r_0340+r_0393-r_0919)/cell;
	ds_1101=(r_0108-r_0386-r_0387-r_0389-r_0391-3.0*r_0393-r_0397-3.0*r_0398-r_0432-r_0433-r_0434-r_0435)/cell;
	ds_1107=(r_0362-0.821*r_2111)/cell;
	ds_1151=(-r_0688-r_0697+r_1936)/cell;
	ds_1153=(r_0757-r_0874)/cell;
	ds_1161=(r_0387-r_0389)/cell;
	ds_1176=(r_0433-r_0434)/cell;
	ds_1182=(r_0118-r_0818)/cell;
	ds_1187=(r_0202-r_0913)/cell;
	ds_1191=(r_0115-r_0759)/cell;
	ds_1192=(-r_0115+r_0818)/cell;
	ds_1194=(r_0214-r_0349)/cell;
	ds_1195=(r_0694-r_0762)/cell;
	ds_1198=(r_0012-r_0061+r_0165-r_0235-r_0445+r_0470-r_0486+r_0491-2.0*r_0536-r_0545-r_0565-r_0658-r_0696-r_0713-r_0731+r_0770-r_0891-r_0961-r_0988+r_1010)/cell;
	ds_1203=(-r_0012+r_0061-r_0165+r_0235+r_0445-r_0470+r_0486-r_0491+2.0*r_0536+r_0545+r_0565+r_0658+r_0696+r_0713+r_0731-r_0770+r_0891+r_0961+r_0988-r_1010)/cell;
	ds_1207=(r_0015+r_0041+r_0080+r_0096-r_0173+r_0219+r_0231+r_0233-r_0234+r_0236+r_0237+r_0238+r_0239+r_0240+3.0*r_0241+r_0242+r_0244+r_0259+r_0267+r_0269+3.0*r_0317+r_0344+2.0*r_0386+2.0*r_0387+2.0*r_0389+2.0*r_0391+6.0*r_0393+2.0*r_0397+6.0*r_0398+2.0*r_0432+2.0*r_0433+2.0*r_0434+2.0*r_0435-r_0466+r_0471+r_0481+r_0547+2.0*r_0558-r_0661+r_0669+r_0678+r_0688-r_0732+r_0759-r_0889+r_0922-r_0939+r_0957+r_0989+r_0996+r_1011+r_1012+3.0*r_1027+r_1038)/cell;
	ds_1212=(-r_0015-r_0041-r_0080-r_0096+r_0173-r_0219-r_0231-r_0233+r_0234-r_0236-r_0237-r_0238-r_0239-r_0240-3.0*r_0241-r_0242-r_0244-r_0259-r_0267-r_0269-3.0*r_0317-r_0344-2.0*r_0386-2.0*r_0387-2.0*r_0389-2.0*r_0391-6.0*r_0393-2.0*r_0397-6.0*r_0398-2.0*r_0432-2.0*r_0433-2.0*r_0434-2.0*r_0435+r_0466-r_0471-r_0481-r_0547-2.0*r_0558+r_0661-r_0669-r_0678-r_0688+r_0732-r_0759+r_0889-r_0922+r_0939-r_0957-r_0989-r_0996-r_1011-r_1012-3.0*r_1027-r_1038)/cell;
	ds_1233=(-r_0311+r_0549-r_0813)/cell;
	ds_1234=(-r_0312+r_0992)/cell;
	ds_1238=(r_0548-r_1041)/cell;
	ds_1255=(-r_0397+r_0398)/cell;
	ds_1266=(-r_0816+r_0818)/cell;
	ds_1269=(r_0339-r_0820)/cell;
	ds_1270=(r_0820-r_0821)/cell;
	ds_1271=(-r_0216-r_0300+r_0713+r_0958)/cell;
	ds_1275=(-r_0233-r_0238-r_0239-r_0240-3.0*r_0241-r_0242-r_0259-r_0267-r_0269-3.0*r_0317-r_0339-r_0438-r_0694-r_0922-r_1010-r_1011+r_1979)/cell;
	ds_1286=(r_0389-r_0391)/cell;
	ds_1302=(r_0434-r_0435-r_0993)/cell;
	ds_1322=(r_0020+r_0032+r_0040+r_0065+r_0079+r_0108+r_0153+r_0214+r_0219-r_0226+r_0250+r_0279+r_0307+r_0337-r_0446+r_0476-r_0486+r_0489+r_0537+2.0*r_0568+r_0726+r_0739+r_0757+r_0759+r_0792+r_0806+r_0816+r_0855+r_0908+r_0911+r_0914+r_0917+r_0958+r_0967+r_1026+r_1041+r_1051+r_1244+r_1936+r_2030+59.3*r_2111)/cell;
	ds_1331=(r_0008-r_0257-r_0337)/cell;
	ds_1337=(-r_0877+r_0880-0.000373*r_2111)/cell;
	ds_1342=(-r_0900+r_0901)/cell;
	ds_1343=(r_0858-r_0901)/cell;
	ds_1346=(r_0900-0.00288*r_2111)/cell;
	ds_1351=(-r_0858+r_0877-0.000697*r_2111)/cell;
	ds_1360=(-r_0020-r_0065+r_0366-r_0962)/cell;
	ds_1364=(-r_0908+r_0911)/cell;
	ds_1365=(-r_0570+r_0912)/cell;
	ds_1366=(-r_0919+r_0922)/cell;
	ds_1376=(-r_0355+r_0667)/cell;
	ds_1377=(r_0278-r_0938-r_0939)/cell;
	ds_1386=(-r_0202-r_0225-r_0820-r_0915+r_0916)/cell;
	ds_1399=(r_0001+r_0004-r_0016-2.0*r_0097+r_0203-r_0674-r_0958-r_0959-r_0961+r_0962)/cell;
	ds_1405=(r_0968-0.0009*r_2111)/cell;
	ds_1408=(-r_0916+r_0982-r_1049)/cell;
	ds_1413=(-r_0144+r_0858+r_0900+r_0901+r_0986)/cell;
	ds_1416=(r_0726-r_0858-r_0900-r_0901-r_0986)/cell;
	ds_1426=(r_0887-r_0990)/cell;
	ds_1427=(-r_0887+r_1049)/cell;
	ds_1429=(r_0996-r_0997)/cell;
	ds_1445=(r_0041-r_0340-r_0922)/cell;
	ds_1447=(-r_1010-r_1011+r_1012)/cell;
	ds_1449=(r_0391-r_0393+r_0407)/cell;
	ds_1454=(-r_0407+r_0435)/cell;
	ds_1467=(-r_1026+r_1266-0.02*r_2111)/cell;
	ds_1469=(r_0883-r_1027)/cell;
	ds_1487=(r_0344+r_0446+r_0499-r_0502+r_0727+r_0912)/cell;
	ds_1520=(r_1051-r_2079-0.0234*r_2111)/cell;
	ds_1524=(-r_0336+r_1052-0.000781*r_2111)/cell;
	ds_1535=(-r_0439+r_0770)/cell;
	ds_1537=(r_0439-r_0770)/cell;
	ds_1538=(r_0005+r_0195+r_0510-r_0811+r_1072)/cell;
	ds_1543=(-r_0005-r_0195-r_0510+r_1084)/cell;
	ds_1545=(r_0821-r_1072-0.067*r_2111)/cell;
	ds_1559=(-r_0307+r_0811-r_0973-r_1084)/cell;
	ds_1565=(-r_0514+r_0565)/cell;
	ds_1569=(r_0237-r_0986-1.5e-05*r_2111)/cell;
	ds_1576=(r_0238-r_0239)/cell;
	ds_1577=(r_0239-r_0240)/cell;
	ds_1578=(-r_0234+r_0240)/cell;
	ds_1579=(r_0234-r_0237)/cell;
	ds_1616=(-r_0550-r_0883-r_0970-r_0973-r_0974-r_0976-r_0978+r_1038)/cell;
	ds_1620=(r_0550+r_0883+r_0970+r_0973+r_0974+r_0976+r_0978-r_1038)/cell;

	return(0);

}


/* Jacobian of the system (dfdx)*/
int amigoJAC_B1(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS_B1(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	return(0);

}

#define	 s_0002 (amigo_model->sim_results[0][j]) 
#define	 s_0008 (amigo_model->sim_results[1][j]) 
#define	 s_0009 (amigo_model->sim_results[2][j]) 
#define	 s_0010 (amigo_model->sim_results[3][j]) 
#define	 s_0015 (amigo_model->sim_results[4][j]) 
#define	 s_0016 (amigo_model->sim_results[5][j]) 
#define	 s_0018 (amigo_model->sim_results[6][j]) 
#define	 s_0019 (amigo_model->sim_results[7][j]) 
#define	 s_0025 (amigo_model->sim_results[8][j]) 
#define	 s_0028 (amigo_model->sim_results[9][j]) 
#define	 s_0033 (amigo_model->sim_results[10][j]) 
#define	 s_0037 (amigo_model->sim_results[11][j]) 
#define	 s_0039 (amigo_model->sim_results[12][j]) 
#define	 s_0056 (amigo_model->sim_results[13][j]) 
#define	 s_0061 (amigo_model->sim_results[14][j]) 
#define	 s_0062 (amigo_model->sim_results[15][j]) 
#define	 s_0063 (amigo_model->sim_results[16][j]) 
#define	 s_0066 (amigo_model->sim_results[17][j]) 
#define	 s_0075 (amigo_model->sim_results[18][j]) 
#define	 s_0076 (amigo_model->sim_results[19][j]) 
#define	 s_0077 (amigo_model->sim_results[20][j]) 
#define	 s_0078 (amigo_model->sim_results[21][j]) 
#define	 s_0082 (amigo_model->sim_results[22][j]) 
#define	 s_0086 (amigo_model->sim_results[23][j]) 
#define	 s_0089 (amigo_model->sim_results[24][j]) 
#define	 s_0118 (amigo_model->sim_results[25][j]) 
#define	 s_0120 (amigo_model->sim_results[26][j]) 
#define	 s_0122 (amigo_model->sim_results[27][j]) 
#define	 s_0126 (amigo_model->sim_results[28][j]) 
#define	 s_0141 (amigo_model->sim_results[29][j]) 
#define	 s_0142 (amigo_model->sim_results[30][j]) 
#define	 s_0145 (amigo_model->sim_results[31][j]) 
#define	 s_0146 (amigo_model->sim_results[32][j]) 
#define	 s_0158 (amigo_model->sim_results[33][j]) 
#define	 s_0162 (amigo_model->sim_results[34][j]) 
#define	 s_0165 (amigo_model->sim_results[35][j]) 
#define	 s_0176 (amigo_model->sim_results[36][j]) 
#define	 s_0178 (amigo_model->sim_results[37][j]) 
#define	 s_0180 (amigo_model->sim_results[38][j]) 
#define	 s_0188 (amigo_model->sim_results[39][j]) 
#define	 s_0190 (amigo_model->sim_results[40][j]) 
#define	 s_0201 (amigo_model->sim_results[41][j]) 
#define	 s_0204 (amigo_model->sim_results[42][j]) 
#define	 s_0207 (amigo_model->sim_results[43][j]) 
#define	 s_0209 (amigo_model->sim_results[44][j]) 
#define	 s_0210 (amigo_model->sim_results[45][j]) 
#define	 s_0211 (amigo_model->sim_results[46][j]) 
#define	 s_0218 (amigo_model->sim_results[47][j]) 
#define	 s_0231 (amigo_model->sim_results[48][j]) 
#define	 s_0232 (amigo_model->sim_results[49][j]) 
#define	 s_0258 (amigo_model->sim_results[50][j]) 
#define	 s_0259 (amigo_model->sim_results[51][j]) 
#define	 s_0260 (amigo_model->sim_results[52][j]) 
#define	 s_0261 (amigo_model->sim_results[53][j]) 
#define	 s_0262 (amigo_model->sim_results[54][j]) 
#define	 s_0291 (amigo_model->sim_results[55][j]) 
#define	 s_0295 (amigo_model->sim_results[56][j]) 
#define	 s_0296 (amigo_model->sim_results[57][j]) 
#define	 s_0297 (amigo_model->sim_results[58][j]) 
#define	 s_0298 (amigo_model->sim_results[59][j]) 
#define	 s_0299 (amigo_model->sim_results[60][j]) 
#define	 s_0300 (amigo_model->sim_results[61][j]) 
#define	 s_0301 (amigo_model->sim_results[62][j]) 
#define	 s_0302 (amigo_model->sim_results[63][j]) 
#define	 s_0304 (amigo_model->sim_results[64][j]) 
#define	 s_0306 (amigo_model->sim_results[65][j]) 
#define	 s_0312 (amigo_model->sim_results[66][j]) 
#define	 s_0313 (amigo_model->sim_results[67][j]) 
#define	 s_0314 (amigo_model->sim_results[68][j]) 
#define	 s_0322 (amigo_model->sim_results[69][j]) 
#define	 s_0324 (amigo_model->sim_results[70][j]) 
#define	 s_0325 (amigo_model->sim_results[71][j]) 
#define	 s_0326 (amigo_model->sim_results[72][j]) 
#define	 s_0327 (amigo_model->sim_results[73][j]) 
#define	 s_0328 (amigo_model->sim_results[74][j]) 
#define	 s_0335 (amigo_model->sim_results[75][j]) 
#define	 s_0340 (amigo_model->sim_results[76][j]) 
#define	 s_0349 (amigo_model->sim_results[77][j]) 
#define	 s_0359 (amigo_model->sim_results[78][j]) 
#define	 s_0362 (amigo_model->sim_results[79][j]) 
#define	 s_0367 (amigo_model->sim_results[80][j]) 
#define	 s_0373 (amigo_model->sim_results[81][j]) 
#define	 s_0380 (amigo_model->sim_results[82][j]) 
#define	 s_0386 (amigo_model->sim_results[83][j]) 
#define	 s_0390 (amigo_model->sim_results[84][j]) 
#define	 s_0393 (amigo_model->sim_results[85][j]) 
#define	 s_0394 (amigo_model->sim_results[86][j]) 
#define	 s_0403 (amigo_model->sim_results[87][j]) 
#define	 s_0409 (amigo_model->sim_results[88][j]) 
#define	 s_0419 (amigo_model->sim_results[89][j]) 
#define	 s_0423 (amigo_model->sim_results[90][j]) 
#define	 s_0427 (amigo_model->sim_results[91][j]) 
#define	 s_0434 (amigo_model->sim_results[92][j]) 
#define	 s_0445 (amigo_model->sim_results[93][j]) 
#define	 s_0454 (amigo_model->sim_results[94][j]) 
#define	 s_0455 (amigo_model->sim_results[95][j]) 
#define	 s_0456 (amigo_model->sim_results[96][j]) 
#define	 s_0467 (amigo_model->sim_results[97][j]) 
#define	 s_0471 (amigo_model->sim_results[98][j]) 
#define	 s_0475 (amigo_model->sim_results[99][j]) 
#define	 s_0481 (amigo_model->sim_results[100][j]) 
#define	 s_0493 (amigo_model->sim_results[101][j]) 
#define	 s_0499 (amigo_model->sim_results[102][j]) 
#define	 s_0515 (amigo_model->sim_results[103][j]) 
#define	 s_0516 (amigo_model->sim_results[104][j]) 
#define	 s_0522 (amigo_model->sim_results[105][j]) 
#define	 s_0526 (amigo_model->sim_results[106][j]) 
#define	 s_0529 (amigo_model->sim_results[107][j]) 
#define	 s_0539 (amigo_model->sim_results[108][j]) 
#define	 s_0550 (amigo_model->sim_results[109][j]) 
#define	 s_0551 (amigo_model->sim_results[110][j]) 
#define	 s_0555 (amigo_model->sim_results[111][j]) 
#define	 s_0557 (amigo_model->sim_results[112][j]) 
#define	 s_0563 (amigo_model->sim_results[113][j]) 
#define	 s_0567 (amigo_model->sim_results[114][j]) 
#define	 s_0568 (amigo_model->sim_results[115][j]) 
#define	 s_0573 (amigo_model->sim_results[116][j]) 
#define	 s_0574 (amigo_model->sim_results[117][j]) 
#define	 s_0577 (amigo_model->sim_results[118][j]) 
#define	 s_0581 (amigo_model->sim_results[119][j]) 
#define	 s_0582 (amigo_model->sim_results[120][j]) 
#define	 s_0584 (amigo_model->sim_results[121][j]) 
#define	 s_0586 (amigo_model->sim_results[122][j]) 
#define	 s_0587 (amigo_model->sim_results[123][j]) 
#define	 s_0589 (amigo_model->sim_results[124][j]) 
#define	 s_0595 (amigo_model->sim_results[125][j]) 
#define	 s_0602 (amigo_model->sim_results[126][j]) 
#define	 s_0613 (amigo_model->sim_results[127][j]) 
#define	 s_0615 (amigo_model->sim_results[128][j]) 
#define	 s_0619 (amigo_model->sim_results[129][j]) 
#define	 s_0625 (amigo_model->sim_results[130][j]) 
#define	 s_0629 (amigo_model->sim_results[131][j]) 
#define	 s_0633 (amigo_model->sim_results[132][j]) 
#define	 s_0644 (amigo_model->sim_results[133][j]) 
#define	 s_0645 (amigo_model->sim_results[134][j]) 
#define	 s_0649 (amigo_model->sim_results[135][j]) 
#define	 s_0654 (amigo_model->sim_results[136][j]) 
#define	 s_0656 (amigo_model->sim_results[137][j]) 
#define	 s_0657 (amigo_model->sim_results[138][j]) 
#define	 s_0662 (amigo_model->sim_results[139][j]) 
#define	 s_0664 (amigo_model->sim_results[140][j]) 
#define	 s_0666 (amigo_model->sim_results[141][j]) 
#define	 s_0672 (amigo_model->sim_results[142][j]) 
#define	 s_0680 (amigo_model->sim_results[143][j]) 
#define	 s_0700 (amigo_model->sim_results[144][j]) 
#define	 s_0709 (amigo_model->sim_results[145][j]) 
#define	 s_0710 (amigo_model->sim_results[146][j]) 
#define	 s_0722 (amigo_model->sim_results[147][j]) 
#define	 s_0725 (amigo_model->sim_results[148][j]) 
#define	 s_0739 (amigo_model->sim_results[149][j]) 
#define	 s_0743 (amigo_model->sim_results[150][j]) 
#define	 s_0745 (amigo_model->sim_results[151][j]) 
#define	 s_0750 (amigo_model->sim_results[152][j]) 
#define	 s_0754 (amigo_model->sim_results[153][j]) 
#define	 s_0764 (amigo_model->sim_results[154][j]) 
#define	 s_0765 (amigo_model->sim_results[155][j]) 
#define	 s_0767 (amigo_model->sim_results[156][j]) 
#define	 s_0773 (amigo_model->sim_results[157][j]) 
#define	 s_0782 (amigo_model->sim_results[158][j]) 
#define	 s_0785 (amigo_model->sim_results[159][j]) 
#define	 s_0835 (amigo_model->sim_results[160][j]) 
#define	 s_0836 (amigo_model->sim_results[161][j]) 
#define	 s_0837 (amigo_model->sim_results[162][j]) 
#define	 s_0841 (amigo_model->sim_results[163][j]) 
#define	 s_0849 (amigo_model->sim_results[164][j]) 
#define	 s_0918 (amigo_model->sim_results[165][j]) 
#define	 s_0940 (amigo_model->sim_results[166][j]) 
#define	 s_0943 (amigo_model->sim_results[167][j]) 
#define	 s_0951 (amigo_model->sim_results[168][j]) 
#define	 s_0953 (amigo_model->sim_results[169][j]) 
#define	 s_0955 (amigo_model->sim_results[170][j]) 
#define	 s_0959 (amigo_model->sim_results[171][j]) 
#define	 s_0965 (amigo_model->sim_results[172][j]) 
#define	 s_0969 (amigo_model->sim_results[173][j]) 
#define	 s_0973 (amigo_model->sim_results[174][j]) 
#define	 s_0978 (amigo_model->sim_results[175][j]) 
#define	 s_0979 (amigo_model->sim_results[176][j]) 
#define	 s_0980 (amigo_model->sim_results[177][j]) 
#define	 s_0981 (amigo_model->sim_results[178][j]) 
#define	 s_0991 (amigo_model->sim_results[179][j]) 
#define	 s_0999 (amigo_model->sim_results[180][j]) 
#define	 s_1003 (amigo_model->sim_results[181][j]) 
#define	 s_1006 (amigo_model->sim_results[182][j]) 
#define	 s_1010 (amigo_model->sim_results[183][j]) 
#define	 s_1011 (amigo_model->sim_results[184][j]) 
#define	 s_1012 (amigo_model->sim_results[185][j]) 
#define	 s_1014 (amigo_model->sim_results[186][j]) 
#define	 s_1016 (amigo_model->sim_results[187][j]) 
#define	 s_1020 (amigo_model->sim_results[188][j]) 
#define	 s_1021 (amigo_model->sim_results[189][j]) 
#define	 s_1025 (amigo_model->sim_results[190][j]) 
#define	 s_1029 (amigo_model->sim_results[191][j]) 
#define	 s_1032 (amigo_model->sim_results[192][j]) 
#define	 s_1035 (amigo_model->sim_results[193][j]) 
#define	 s_1038 (amigo_model->sim_results[194][j]) 
#define	 s_1039 (amigo_model->sim_results[195][j]) 
#define	 s_1045 (amigo_model->sim_results[196][j]) 
#define	 s_1048 (amigo_model->sim_results[197][j]) 
#define	 s_1051 (amigo_model->sim_results[198][j]) 
#define	 s_1056 (amigo_model->sim_results[199][j]) 
#define	 s_1059 (amigo_model->sim_results[200][j]) 
#define	 s_1065 (amigo_model->sim_results[201][j]) 
#define	 s_1073 (amigo_model->sim_results[202][j]) 
#define	 s_1084 (amigo_model->sim_results[203][j]) 
#define	 s_1101 (amigo_model->sim_results[204][j]) 
#define	 s_1107 (amigo_model->sim_results[205][j]) 
#define	 s_1151 (amigo_model->sim_results[206][j]) 
#define	 s_1153 (amigo_model->sim_results[207][j]) 
#define	 s_1161 (amigo_model->sim_results[208][j]) 
#define	 s_1176 (amigo_model->sim_results[209][j]) 
#define	 s_1182 (amigo_model->sim_results[210][j]) 
#define	 s_1187 (amigo_model->sim_results[211][j]) 
#define	 s_1191 (amigo_model->sim_results[212][j]) 
#define	 s_1192 (amigo_model->sim_results[213][j]) 
#define	 s_1194 (amigo_model->sim_results[214][j]) 
#define	 s_1195 (amigo_model->sim_results[215][j]) 
#define	 s_1198 (amigo_model->sim_results[216][j]) 
#define	 s_1203 (amigo_model->sim_results[217][j]) 
#define	 s_1207 (amigo_model->sim_results[218][j]) 
#define	 s_1212 (amigo_model->sim_results[219][j]) 
#define	 s_1233 (amigo_model->sim_results[220][j]) 
#define	 s_1234 (amigo_model->sim_results[221][j]) 
#define	 s_1238 (amigo_model->sim_results[222][j]) 
#define	 s_1255 (amigo_model->sim_results[223][j]) 
#define	 s_1266 (amigo_model->sim_results[224][j]) 
#define	 s_1269 (amigo_model->sim_results[225][j]) 
#define	 s_1270 (amigo_model->sim_results[226][j]) 
#define	 s_1271 (amigo_model->sim_results[227][j]) 
#define	 s_1275 (amigo_model->sim_results[228][j]) 
#define	 s_1286 (amigo_model->sim_results[229][j]) 
#define	 s_1302 (amigo_model->sim_results[230][j]) 
#define	 s_1322 (amigo_model->sim_results[231][j]) 
#define	 s_1331 (amigo_model->sim_results[232][j]) 
#define	 s_1337 (amigo_model->sim_results[233][j]) 
#define	 s_1342 (amigo_model->sim_results[234][j]) 
#define	 s_1343 (amigo_model->sim_results[235][j]) 
#define	 s_1346 (amigo_model->sim_results[236][j]) 
#define	 s_1351 (amigo_model->sim_results[237][j]) 
#define	 s_1360 (amigo_model->sim_results[238][j]) 
#define	 s_1364 (amigo_model->sim_results[239][j]) 
#define	 s_1365 (amigo_model->sim_results[240][j]) 
#define	 s_1366 (amigo_model->sim_results[241][j]) 
#define	 s_1376 (amigo_model->sim_results[242][j]) 
#define	 s_1377 (amigo_model->sim_results[243][j]) 
#define	 s_1386 (amigo_model->sim_results[244][j]) 
#define	 s_1399 (amigo_model->sim_results[245][j]) 
#define	 s_1405 (amigo_model->sim_results[246][j]) 
#define	 s_1408 (amigo_model->sim_results[247][j]) 
#define	 s_1413 (amigo_model->sim_results[248][j]) 
#define	 s_1416 (amigo_model->sim_results[249][j]) 
#define	 s_1426 (amigo_model->sim_results[250][j]) 
#define	 s_1427 (amigo_model->sim_results[251][j]) 
#define	 s_1429 (amigo_model->sim_results[252][j]) 
#define	 s_1445 (amigo_model->sim_results[253][j]) 
#define	 s_1447 (amigo_model->sim_results[254][j]) 
#define	 s_1449 (amigo_model->sim_results[255][j]) 
#define	 s_1454 (amigo_model->sim_results[256][j]) 
#define	 s_1467 (amigo_model->sim_results[257][j]) 
#define	 s_1469 (amigo_model->sim_results[258][j]) 
#define	 s_1487 (amigo_model->sim_results[259][j]) 
#define	 s_1520 (amigo_model->sim_results[260][j]) 
#define	 s_1524 (amigo_model->sim_results[261][j]) 
#define	 s_1535 (amigo_model->sim_results[262][j]) 
#define	 s_1537 (amigo_model->sim_results[263][j]) 
#define	 s_1538 (amigo_model->sim_results[264][j]) 
#define	 s_1543 (amigo_model->sim_results[265][j]) 
#define	 s_1545 (amigo_model->sim_results[266][j]) 
#define	 s_1559 (amigo_model->sim_results[267][j]) 
#define	 s_1565 (amigo_model->sim_results[268][j]) 
#define	 s_1569 (amigo_model->sim_results[269][j]) 
#define	 s_1576 (amigo_model->sim_results[270][j]) 
#define	 s_1577 (amigo_model->sim_results[271][j]) 
#define	 s_1578 (amigo_model->sim_results[272][j]) 
#define	 s_1579 (amigo_model->sim_results[273][j]) 
#define	 s_1616 (amigo_model->sim_results[274][j]) 
#define	 s_1620 (amigo_model->sim_results[275][j]) 



void amigoRHS_get_OBS_B1(void* data){

	int j;
	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;
    double s_0565_substitute;

	 switch (amigo_model->exp_num){

		#define	 s_0075_obs amigo_model->obs_results[0][j] 
		#define	 s_0180_obs amigo_model->obs_results[1][j] 
		#define	 s_0188_obs amigo_model->obs_results[2][j] 
		#define	 s_0259_obs amigo_model->obs_results[3][j] 
		#define	 s_0260_obs amigo_model->obs_results[4][j] 
		#define	 s_0359_obs amigo_model->obs_results[5][j] 
		#define	 s_0362_obs amigo_model->obs_results[6][j] 
		#define	 s_0394_obs amigo_model->obs_results[7][j] 
		#define	 s_0409_obs amigo_model->obs_results[8][j] 
		#define	 s_0423_obs amigo_model->obs_results[9][j] 
		#define	 s_0434_obs amigo_model->obs_results[10][j] 
		#define	 s_0555_obs amigo_model->obs_results[11][j] 
		#define	 s_0557_obs amigo_model->obs_results[12][j] 
		#define	 s_0563_obs amigo_model->obs_results[13][j] 
		#define	 s_0567_obs amigo_model->obs_results[14][j] 
		#define	 s_0568_obs amigo_model->obs_results[15][j] 
		#define	 s_0586_obs amigo_model->obs_results[16][j] 
		#define	 s_0629_obs amigo_model->obs_results[17][j] 
		#define	 s_0680_obs amigo_model->obs_results[18][j] 
		#define	 s_0765_obs amigo_model->obs_results[19][j] 
		#define	 s_0764_obs amigo_model->obs_results[20][j] 
		#define	 s_0767_obs amigo_model->obs_results[21][j] 
		#define	 s_0785_obs amigo_model->obs_results[22][j] 
		#define	 s_0849_obs amigo_model->obs_results[23][j] 
		#define	 s_0955_obs amigo_model->obs_results[24][j] 
		#define	 s_0991_obs amigo_model->obs_results[25][j] 
		#define	 s_1003_obs amigo_model->obs_results[26][j] 
		#define	 s_1039_obs amigo_model->obs_results[27][j] 
		#define	 s_1045_obs amigo_model->obs_results[28][j] 
		#define	 s_1198_obs amigo_model->obs_results[29][j] 
		#define	 s_1203_obs amigo_model->obs_results[30][j] 
		#define	 s_1360_obs amigo_model->obs_results[31][j] 
		#define	 s_1399_obs amigo_model->obs_results[32][j] 
		#define	 s_1520_obs amigo_model->obs_results[33][j] 
		#define	 s_1538_obs amigo_model->obs_results[34][j] 
		#define	 s_1543_obs amigo_model->obs_results[35][j] 
		#define	 s_1559_obs amigo_model->obs_results[36][j] 
		#define	 s_1565_obs amigo_model->obs_results[37][j] 
		#define	 r_1166_obs amigo_model->obs_results[38][j] 
		#define	 r_1697_obs amigo_model->obs_results[39][j] 
		#define	 r_1762_obs amigo_model->obs_results[40][j] 
		#define	 r_1106_obs amigo_model->obs_results[41][j] 
		#define	 r_1172_obs amigo_model->obs_results[42][j] 
		#define	 r_2079_obs amigo_model->obs_results[43][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){
				
                  if(j<10){
                s_0565_substitute=74;
            }else{
                s_0565_substitute=1;
            }   
				s_0075_obs=s_0075;
				s_0180_obs=s_0180;
				s_0188_obs=s_0188;
				s_0259_obs=s_0259;
				s_0260_obs=s_0260;
				s_0359_obs=s_0359;
				s_0362_obs=s_0362;
				s_0394_obs=s_0394;
				s_0409_obs=s_0409;
				s_0423_obs=s_0423;
				s_0434_obs=s_0434;
				s_0555_obs=s_0555;
				s_0557_obs=s_0557;
				s_0563_obs=s_0563;
				s_0567_obs=s_0567;
				s_0568_obs=s_0568;
				s_0586_obs=s_0586;
				s_0629_obs=s_0629;
				s_0680_obs=s_0680;
				s_0765_obs=s_0765;
				s_0764_obs=s_0764;
				s_0767_obs=s_0767;
				s_0785_obs=s_0785;
				s_0849_obs=s_0849;
				s_0955_obs=s_0955;
				s_0991_obs=s_0991;
				s_1003_obs=s_1003;
				s_1039_obs=s_1039;
				s_1045_obs=s_1045;
				s_1198_obs=s_1198;
				s_1203_obs=s_1203;
				s_1360_obs=s_1360;
				s_1399_obs=s_1399;
				s_1520_obs=s_1520;
				s_1538_obs=s_1538;
				s_1543_obs=s_1543;
				s_1559_obs=s_1559;
				s_1565_obs=s_1565;
				r_1166_obs=cell*Vmax_1166*( s_0565_substitute-s_0563)/Km0565_1166/(1+ s_0565_substitute/Km0565_1166+1+s_0563/Km0563_1166-1);
				r_1697_obs=cell*Vmax_1697*s_0456/Km0456_1697/(1+s_0456/Km0456_1697);
				r_1762_obs=cell*Vmax_1762*s_0680/Km0680_1762/(1+s_0680/Km0680_1762);
				r_1106_obs=cell*Vmax_1106*s_0362/Km0362_1106/(1+s_0362/Km0362_1106);
				r_1172_obs=cell*Vmax_1172*s_0765/Km0765_1172/(1+s_0765/Km0765_1172);
				r_2079_obs=cell*Vmax_2079*s_1520/Km1520_2079/(1+s_1520/Km1520_2079);

			}

		 break;

	}

	return(amigo_model);

}

#define	 s_0002 (amigo_model->sens_results[0][j][k]) 
#define	 s_0008 (amigo_model->sens_results[1][j][k]) 
#define	 s_0009 (amigo_model->sens_results[2][j][k]) 
#define	 s_0010 (amigo_model->sens_results[3][j][k]) 
#define	 s_0015 (amigo_model->sens_results[4][j][k]) 
#define	 s_0016 (amigo_model->sens_results[5][j][k]) 
#define	 s_0018 (amigo_model->sens_results[6][j][k]) 
#define	 s_0019 (amigo_model->sens_results[7][j][k]) 
#define	 s_0025 (amigo_model->sens_results[8][j][k]) 
#define	 s_0028 (amigo_model->sens_results[9][j][k]) 
#define	 s_0033 (amigo_model->sens_results[10][j][k]) 
#define	 s_0037 (amigo_model->sens_results[11][j][k]) 
#define	 s_0039 (amigo_model->sens_results[12][j][k]) 
#define	 s_0056 (amigo_model->sens_results[13][j][k]) 
#define	 s_0061 (amigo_model->sens_results[14][j][k]) 
#define	 s_0062 (amigo_model->sens_results[15][j][k]) 
#define	 s_0063 (amigo_model->sens_results[16][j][k]) 
#define	 s_0066 (amigo_model->sens_results[17][j][k]) 
#define	 s_0075 (amigo_model->sens_results[18][j][k]) 
#define	 s_0076 (amigo_model->sens_results[19][j][k]) 
#define	 s_0077 (amigo_model->sens_results[20][j][k]) 
#define	 s_0078 (amigo_model->sens_results[21][j][k]) 
#define	 s_0082 (amigo_model->sens_results[22][j][k]) 
#define	 s_0086 (amigo_model->sens_results[23][j][k]) 
#define	 s_0089 (amigo_model->sens_results[24][j][k]) 
#define	 s_0118 (amigo_model->sens_results[25][j][k]) 
#define	 s_0120 (amigo_model->sens_results[26][j][k]) 
#define	 s_0122 (amigo_model->sens_results[27][j][k]) 
#define	 s_0126 (amigo_model->sens_results[28][j][k]) 
#define	 s_0141 (amigo_model->sens_results[29][j][k]) 
#define	 s_0142 (amigo_model->sens_results[30][j][k]) 
#define	 s_0145 (amigo_model->sens_results[31][j][k]) 
#define	 s_0146 (amigo_model->sens_results[32][j][k]) 
#define	 s_0158 (amigo_model->sens_results[33][j][k]) 
#define	 s_0162 (amigo_model->sens_results[34][j][k]) 
#define	 s_0165 (amigo_model->sens_results[35][j][k]) 
#define	 s_0176 (amigo_model->sens_results[36][j][k]) 
#define	 s_0178 (amigo_model->sens_results[37][j][k]) 
#define	 s_0180 (amigo_model->sens_results[38][j][k]) 
#define	 s_0188 (amigo_model->sens_results[39][j][k]) 
#define	 s_0190 (amigo_model->sens_results[40][j][k]) 
#define	 s_0201 (amigo_model->sens_results[41][j][k]) 
#define	 s_0204 (amigo_model->sens_results[42][j][k]) 
#define	 s_0207 (amigo_model->sens_results[43][j][k]) 
#define	 s_0209 (amigo_model->sens_results[44][j][k]) 
#define	 s_0210 (amigo_model->sens_results[45][j][k]) 
#define	 s_0211 (amigo_model->sens_results[46][j][k]) 
#define	 s_0218 (amigo_model->sens_results[47][j][k]) 
#define	 s_0231 (amigo_model->sens_results[48][j][k]) 
#define	 s_0232 (amigo_model->sens_results[49][j][k]) 
#define	 s_0258 (amigo_model->sens_results[50][j][k]) 
#define	 s_0259 (amigo_model->sens_results[51][j][k]) 
#define	 s_0260 (amigo_model->sens_results[52][j][k]) 
#define	 s_0261 (amigo_model->sens_results[53][j][k]) 
#define	 s_0262 (amigo_model->sens_results[54][j][k]) 
#define	 s_0291 (amigo_model->sens_results[55][j][k]) 
#define	 s_0295 (amigo_model->sens_results[56][j][k]) 
#define	 s_0296 (amigo_model->sens_results[57][j][k]) 
#define	 s_0297 (amigo_model->sens_results[58][j][k]) 
#define	 s_0298 (amigo_model->sens_results[59][j][k]) 
#define	 s_0299 (amigo_model->sens_results[60][j][k]) 
#define	 s_0300 (amigo_model->sens_results[61][j][k]) 
#define	 s_0301 (amigo_model->sens_results[62][j][k]) 
#define	 s_0302 (amigo_model->sens_results[63][j][k]) 
#define	 s_0304 (amigo_model->sens_results[64][j][k]) 
#define	 s_0306 (amigo_model->sens_results[65][j][k]) 
#define	 s_0312 (amigo_model->sens_results[66][j][k]) 
#define	 s_0313 (amigo_model->sens_results[67][j][k]) 
#define	 s_0314 (amigo_model->sens_results[68][j][k]) 
#define	 s_0322 (amigo_model->sens_results[69][j][k]) 
#define	 s_0324 (amigo_model->sens_results[70][j][k]) 
#define	 s_0325 (amigo_model->sens_results[71][j][k]) 
#define	 s_0326 (amigo_model->sens_results[72][j][k]) 
#define	 s_0327 (amigo_model->sens_results[73][j][k]) 
#define	 s_0328 (amigo_model->sens_results[74][j][k]) 
#define	 s_0335 (amigo_model->sens_results[75][j][k]) 
#define	 s_0340 (amigo_model->sens_results[76][j][k]) 
#define	 s_0349 (amigo_model->sens_results[77][j][k]) 
#define	 s_0359 (amigo_model->sens_results[78][j][k]) 
#define	 s_0362 (amigo_model->sens_results[79][j][k]) 
#define	 s_0367 (amigo_model->sens_results[80][j][k]) 
#define	 s_0373 (amigo_model->sens_results[81][j][k]) 
#define	 s_0380 (amigo_model->sens_results[82][j][k]) 
#define	 s_0386 (amigo_model->sens_results[83][j][k]) 
#define	 s_0390 (amigo_model->sens_results[84][j][k]) 
#define	 s_0393 (amigo_model->sens_results[85][j][k]) 
#define	 s_0394 (amigo_model->sens_results[86][j][k]) 
#define	 s_0403 (amigo_model->sens_results[87][j][k]) 
#define	 s_0409 (amigo_model->sens_results[88][j][k]) 
#define	 s_0419 (amigo_model->sens_results[89][j][k]) 
#define	 s_0423 (amigo_model->sens_results[90][j][k]) 
#define	 s_0427 (amigo_model->sens_results[91][j][k]) 
#define	 s_0434 (amigo_model->sens_results[92][j][k]) 
#define	 s_0445 (amigo_model->sens_results[93][j][k]) 
#define	 s_0454 (amigo_model->sens_results[94][j][k]) 
#define	 s_0455 (amigo_model->sens_results[95][j][k]) 
#define	 s_0456 (amigo_model->sens_results[96][j][k]) 
#define	 s_0467 (amigo_model->sens_results[97][j][k]) 
#define	 s_0471 (amigo_model->sens_results[98][j][k]) 
#define	 s_0475 (amigo_model->sens_results[99][j][k]) 
#define	 s_0481 (amigo_model->sens_results[100][j][k]) 
#define	 s_0493 (amigo_model->sens_results[101][j][k]) 
#define	 s_0499 (amigo_model->sens_results[102][j][k]) 
#define	 s_0515 (amigo_model->sens_results[103][j][k]) 
#define	 s_0516 (amigo_model->sens_results[104][j][k]) 
#define	 s_0522 (amigo_model->sens_results[105][j][k]) 
#define	 s_0526 (amigo_model->sens_results[106][j][k]) 
#define	 s_0529 (amigo_model->sens_results[107][j][k]) 
#define	 s_0539 (amigo_model->sens_results[108][j][k]) 
#define	 s_0550 (amigo_model->sens_results[109][j][k]) 
#define	 s_0551 (amigo_model->sens_results[110][j][k]) 
#define	 s_0555 (amigo_model->sens_results[111][j][k]) 
#define	 s_0557 (amigo_model->sens_results[112][j][k]) 
#define	 s_0563 (amigo_model->sens_results[113][j][k]) 
#define	 s_0567 (amigo_model->sens_results[114][j][k]) 
#define	 s_0568 (amigo_model->sens_results[115][j][k]) 
#define	 s_0573 (amigo_model->sens_results[116][j][k]) 
#define	 s_0574 (amigo_model->sens_results[117][j][k]) 
#define	 s_0577 (amigo_model->sens_results[118][j][k]) 
#define	 s_0581 (amigo_model->sens_results[119][j][k]) 
#define	 s_0582 (amigo_model->sens_results[120][j][k]) 
#define	 s_0584 (amigo_model->sens_results[121][j][k]) 
#define	 s_0586 (amigo_model->sens_results[122][j][k]) 
#define	 s_0587 (amigo_model->sens_results[123][j][k]) 
#define	 s_0589 (amigo_model->sens_results[124][j][k]) 
#define	 s_0595 (amigo_model->sens_results[125][j][k]) 
#define	 s_0602 (amigo_model->sens_results[126][j][k]) 
#define	 s_0613 (amigo_model->sens_results[127][j][k]) 
#define	 s_0615 (amigo_model->sens_results[128][j][k]) 
#define	 s_0619 (amigo_model->sens_results[129][j][k]) 
#define	 s_0625 (amigo_model->sens_results[130][j][k]) 
#define	 s_0629 (amigo_model->sens_results[131][j][k]) 
#define	 s_0633 (amigo_model->sens_results[132][j][k]) 
#define	 s_0644 (amigo_model->sens_results[133][j][k]) 
#define	 s_0645 (amigo_model->sens_results[134][j][k]) 
#define	 s_0649 (amigo_model->sens_results[135][j][k]) 
#define	 s_0654 (amigo_model->sens_results[136][j][k]) 
#define	 s_0656 (amigo_model->sens_results[137][j][k]) 
#define	 s_0657 (amigo_model->sens_results[138][j][k]) 
#define	 s_0662 (amigo_model->sens_results[139][j][k]) 
#define	 s_0664 (amigo_model->sens_results[140][j][k]) 
#define	 s_0666 (amigo_model->sens_results[141][j][k]) 
#define	 s_0672 (amigo_model->sens_results[142][j][k]) 
#define	 s_0680 (amigo_model->sens_results[143][j][k]) 
#define	 s_0700 (amigo_model->sens_results[144][j][k]) 
#define	 s_0709 (amigo_model->sens_results[145][j][k]) 
#define	 s_0710 (amigo_model->sens_results[146][j][k]) 
#define	 s_0722 (amigo_model->sens_results[147][j][k]) 
#define	 s_0725 (amigo_model->sens_results[148][j][k]) 
#define	 s_0739 (amigo_model->sens_results[149][j][k]) 
#define	 s_0743 (amigo_model->sens_results[150][j][k]) 
#define	 s_0745 (amigo_model->sens_results[151][j][k]) 
#define	 s_0750 (amigo_model->sens_results[152][j][k]) 
#define	 s_0754 (amigo_model->sens_results[153][j][k]) 
#define	 s_0764 (amigo_model->sens_results[154][j][k]) 
#define	 s_0765 (amigo_model->sens_results[155][j][k]) 
#define	 s_0767 (amigo_model->sens_results[156][j][k]) 
#define	 s_0773 (amigo_model->sens_results[157][j][k]) 
#define	 s_0782 (amigo_model->sens_results[158][j][k]) 
#define	 s_0785 (amigo_model->sens_results[159][j][k]) 
#define	 s_0835 (amigo_model->sens_results[160][j][k]) 
#define	 s_0836 (amigo_model->sens_results[161][j][k]) 
#define	 s_0837 (amigo_model->sens_results[162][j][k]) 
#define	 s_0841 (amigo_model->sens_results[163][j][k]) 
#define	 s_0849 (amigo_model->sens_results[164][j][k]) 
#define	 s_0918 (amigo_model->sens_results[165][j][k]) 
#define	 s_0940 (amigo_model->sens_results[166][j][k]) 
#define	 s_0943 (amigo_model->sens_results[167][j][k]) 
#define	 s_0951 (amigo_model->sens_results[168][j][k]) 
#define	 s_0953 (amigo_model->sens_results[169][j][k]) 
#define	 s_0955 (amigo_model->sens_results[170][j][k]) 
#define	 s_0959 (amigo_model->sens_results[171][j][k]) 
#define	 s_0965 (amigo_model->sens_results[172][j][k]) 
#define	 s_0969 (amigo_model->sens_results[173][j][k]) 
#define	 s_0973 (amigo_model->sens_results[174][j][k]) 
#define	 s_0978 (amigo_model->sens_results[175][j][k]) 
#define	 s_0979 (amigo_model->sens_results[176][j][k]) 
#define	 s_0980 (amigo_model->sens_results[177][j][k]) 
#define	 s_0981 (amigo_model->sens_results[178][j][k]) 
#define	 s_0991 (amigo_model->sens_results[179][j][k]) 
#define	 s_0999 (amigo_model->sens_results[180][j][k]) 
#define	 s_1003 (amigo_model->sens_results[181][j][k]) 
#define	 s_1006 (amigo_model->sens_results[182][j][k]) 
#define	 s_1010 (amigo_model->sens_results[183][j][k]) 
#define	 s_1011 (amigo_model->sens_results[184][j][k]) 
#define	 s_1012 (amigo_model->sens_results[185][j][k]) 
#define	 s_1014 (amigo_model->sens_results[186][j][k]) 
#define	 s_1016 (amigo_model->sens_results[187][j][k]) 
#define	 s_1020 (amigo_model->sens_results[188][j][k]) 
#define	 s_1021 (amigo_model->sens_results[189][j][k]) 
#define	 s_1025 (amigo_model->sens_results[190][j][k]) 
#define	 s_1029 (amigo_model->sens_results[191][j][k]) 
#define	 s_1032 (amigo_model->sens_results[192][j][k]) 
#define	 s_1035 (amigo_model->sens_results[193][j][k]) 
#define	 s_1038 (amigo_model->sens_results[194][j][k]) 
#define	 s_1039 (amigo_model->sens_results[195][j][k]) 
#define	 s_1045 (amigo_model->sens_results[196][j][k]) 
#define	 s_1048 (amigo_model->sens_results[197][j][k]) 
#define	 s_1051 (amigo_model->sens_results[198][j][k]) 
#define	 s_1056 (amigo_model->sens_results[199][j][k]) 
#define	 s_1059 (amigo_model->sens_results[200][j][k]) 
#define	 s_1065 (amigo_model->sens_results[201][j][k]) 
#define	 s_1073 (amigo_model->sens_results[202][j][k]) 
#define	 s_1084 (amigo_model->sens_results[203][j][k]) 
#define	 s_1101 (amigo_model->sens_results[204][j][k]) 
#define	 s_1107 (amigo_model->sens_results[205][j][k]) 
#define	 s_1151 (amigo_model->sens_results[206][j][k]) 
#define	 s_1153 (amigo_model->sens_results[207][j][k]) 
#define	 s_1161 (amigo_model->sens_results[208][j][k]) 
#define	 s_1176 (amigo_model->sens_results[209][j][k]) 
#define	 s_1182 (amigo_model->sens_results[210][j][k]) 
#define	 s_1187 (amigo_model->sens_results[211][j][k]) 
#define	 s_1191 (amigo_model->sens_results[212][j][k]) 
#define	 s_1192 (amigo_model->sens_results[213][j][k]) 
#define	 s_1194 (amigo_model->sens_results[214][j][k]) 
#define	 s_1195 (amigo_model->sens_results[215][j][k]) 
#define	 s_1198 (amigo_model->sens_results[216][j][k]) 
#define	 s_1203 (amigo_model->sens_results[217][j][k]) 
#define	 s_1207 (amigo_model->sens_results[218][j][k]) 
#define	 s_1212 (amigo_model->sens_results[219][j][k]) 
#define	 s_1233 (amigo_model->sens_results[220][j][k]) 
#define	 s_1234 (amigo_model->sens_results[221][j][k]) 
#define	 s_1238 (amigo_model->sens_results[222][j][k]) 
#define	 s_1255 (amigo_model->sens_results[223][j][k]) 
#define	 s_1266 (amigo_model->sens_results[224][j][k]) 
#define	 s_1269 (amigo_model->sens_results[225][j][k]) 
#define	 s_1270 (amigo_model->sens_results[226][j][k]) 
#define	 s_1271 (amigo_model->sens_results[227][j][k]) 
#define	 s_1275 (amigo_model->sens_results[228][j][k]) 
#define	 s_1286 (amigo_model->sens_results[229][j][k]) 
#define	 s_1302 (amigo_model->sens_results[230][j][k]) 
#define	 s_1322 (amigo_model->sens_results[231][j][k]) 
#define	 s_1331 (amigo_model->sens_results[232][j][k]) 
#define	 s_1337 (amigo_model->sens_results[233][j][k]) 
#define	 s_1342 (amigo_model->sens_results[234][j][k]) 
#define	 s_1343 (amigo_model->sens_results[235][j][k]) 
#define	 s_1346 (amigo_model->sens_results[236][j][k]) 
#define	 s_1351 (amigo_model->sens_results[237][j][k]) 
#define	 s_1360 (amigo_model->sens_results[238][j][k]) 
#define	 s_1364 (amigo_model->sens_results[239][j][k]) 
#define	 s_1365 (amigo_model->sens_results[240][j][k]) 
#define	 s_1366 (amigo_model->sens_results[241][j][k]) 
#define	 s_1376 (amigo_model->sens_results[242][j][k]) 
#define	 s_1377 (amigo_model->sens_results[243][j][k]) 
#define	 s_1386 (amigo_model->sens_results[244][j][k]) 
#define	 s_1399 (amigo_model->sens_results[245][j][k]) 
#define	 s_1405 (amigo_model->sens_results[246][j][k]) 
#define	 s_1408 (amigo_model->sens_results[247][j][k]) 
#define	 s_1413 (amigo_model->sens_results[248][j][k]) 
#define	 s_1416 (amigo_model->sens_results[249][j][k]) 
#define	 s_1426 (amigo_model->sens_results[250][j][k]) 
#define	 s_1427 (amigo_model->sens_results[251][j][k]) 
#define	 s_1429 (amigo_model->sens_results[252][j][k]) 
#define	 s_1445 (amigo_model->sens_results[253][j][k]) 
#define	 s_1447 (amigo_model->sens_results[254][j][k]) 
#define	 s_1449 (amigo_model->sens_results[255][j][k]) 
#define	 s_1454 (amigo_model->sens_results[256][j][k]) 
#define	 s_1467 (amigo_model->sens_results[257][j][k]) 
#define	 s_1469 (amigo_model->sens_results[258][j][k]) 
#define	 s_1487 (amigo_model->sens_results[259][j][k]) 
#define	 s_1520 (amigo_model->sens_results[260][j][k]) 
#define	 s_1524 (amigo_model->sens_results[261][j][k]) 
#define	 s_1535 (amigo_model->sens_results[262][j][k]) 
#define	 s_1537 (amigo_model->sens_results[263][j][k]) 
#define	 s_1538 (amigo_model->sens_results[264][j][k]) 
#define	 s_1543 (amigo_model->sens_results[265][j][k]) 
#define	 s_1545 (amigo_model->sens_results[266][j][k]) 
#define	 s_1559 (amigo_model->sens_results[267][j][k]) 
#define	 s_1565 (amigo_model->sens_results[268][j][k]) 
#define	 s_1569 (amigo_model->sens_results[269][j][k]) 
#define	 s_1576 (amigo_model->sens_results[270][j][k]) 
#define	 s_1577 (amigo_model->sens_results[271][j][k]) 
#define	 s_1578 (amigo_model->sens_results[272][j][k]) 
#define	 s_1579 (amigo_model->sens_results[273][j][k]) 
#define	 s_1616 (amigo_model->sens_results[274][j][k]) 
#define	 s_1620 (amigo_model->sens_results[275][j][k]) 



void amigoRHS_get_sens_OBS_B1(void* data){
	int j,k;

	AMIGO_model* amigo_model=(AMIGO_model*)data;
    double s_0565_substitute=0;

	 switch (amigo_model->exp_num){


		 case 0:
 
             
		#define	 s_0075_obs amigo_model->sens_results[0][j][k] 
		#define	 s_0180_obs amigo_model->sens_results[1][j][k] 
		#define	 s_0188_obs amigo_model->sens_results[2][j][k] 
		#define	 s_0259_obs amigo_model->sens_results[3][j][k] 
		#define	 s_0260_obs amigo_model->sens_results[4][j][k] 
		#define	 s_0359_obs amigo_model->sens_results[5][j][k] 
		#define	 s_0362_obs amigo_model->sens_results[6][j][k] 
		#define	 s_0394_obs amigo_model->sens_results[7][j][k] 
		#define	 s_0409_obs amigo_model->sens_results[8][j][k] 
		#define	 s_0423_obs amigo_model->sens_results[9][j][k] 
		#define	 s_0434_obs amigo_model->sens_results[10][j][k] 
		#define	 s_0555_obs amigo_model->sens_results[11][j][k] 
		#define	 s_0557_obs amigo_model->sens_results[12][j][k] 
		#define	 s_0563_obs amigo_model->sens_results[13][j][k] 
		#define	 s_0567_obs amigo_model->sens_results[14][j][k] 
		#define	 s_0568_obs amigo_model->sens_results[15][j][k] 
		#define	 s_0586_obs amigo_model->sens_results[16][j][k] 
		#define	 s_0629_obs amigo_model->sens_results[17][j][k] 
		#define	 s_0680_obs amigo_model->sens_results[18][j][k] 
		#define	 s_0765_obs amigo_model->sens_results[19][j][k] 
		#define	 s_0764_obs amigo_model->sens_results[20][j][k] 
		#define	 s_0767_obs amigo_model->sens_results[21][j][k] 
		#define	 s_0785_obs amigo_model->sens_results[22][j][k] 
		#define	 s_0849_obs amigo_model->sens_results[23][j][k] 
		#define	 s_0955_obs amigo_model->sens_results[24][j][k] 
		#define	 s_0991_obs amigo_model->sens_results[25][j][k] 
		#define	 s_1003_obs amigo_model->sens_results[26][j][k] 
		#define	 s_1039_obs amigo_model->sens_results[27][j][k] 
		#define	 s_1045_obs amigo_model->sens_results[28][j][k] 
		#define	 s_1198_obs amigo_model->sens_results[29][j][k] 
		#define	 s_1203_obs amigo_model->sens_results[30][j][k] 
		#define	 s_1360_obs amigo_model->sens_results[31][j][k] 
		#define	 s_1399_obs amigo_model->sens_results[32][j][k] 
		#define	 s_1520_obs amigo_model->sens_results[33][j][k] 
		#define	 s_1538_obs amigo_model->sens_results[34][j][k] 
		#define	 s_1543_obs amigo_model->sens_results[35][j][k] 
		#define	 s_1559_obs amigo_model->sens_results[36][j][k] 
		#define	 s_1565_obs amigo_model->sens_results[37][j][k] 
		#define	 r_1166_obs amigo_model->sens_results[38][j][k] 
		#define	 r_1697_obs amigo_model->sens_results[39][j][k] 
		#define	 r_1762_obs amigo_model->sens_results[40][j][k] 
		#define	 r_1106_obs amigo_model->sens_results[41][j][k] 
		#define	 r_1172_obs amigo_model->sens_results[42][j][k] 
		#define	 r_2079_obs amigo_model->sens_results[43][j][k] 

			 for (j = 0; j < amigo_model->n_times; ++j){
				 for (k = 0; k < amigo_model->n_total_x; ++k){
					
					s_0565_substitute=0;
					s_0075_obs=s_0075;
					s_0180_obs=s_0180;
					s_0188_obs=s_0188;
					s_0259_obs=s_0259;
					s_0260_obs=s_0260;
					s_0359_obs=s_0359;
					s_0362_obs=s_0362;
					s_0394_obs=s_0394;
					s_0409_obs=s_0409;
					s_0423_obs=s_0423;
					s_0434_obs=s_0434;
					s_0555_obs=s_0555;
					s_0557_obs=s_0557;
					s_0563_obs=s_0563;
					s_0567_obs=s_0567;
					s_0568_obs=s_0568;
					s_0586_obs=s_0586;
					s_0629_obs=s_0629;
					s_0680_obs=s_0680;
					s_0765_obs=s_0765;
					s_0764_obs=s_0764;
					s_0767_obs=s_0767;
					s_0785_obs=s_0785;
					s_0849_obs=s_0849;
					s_0955_obs=s_0955;
					s_0991_obs=s_0991;
					s_1003_obs=s_1003;
					s_1039_obs=s_1039;
					s_1045_obs=s_1045;
					s_1198_obs=s_1198;
					s_1203_obs=s_1203;
					s_1360_obs=s_1360;
					s_1399_obs=s_1399;
					s_1520_obs=s_1520;
					s_1538_obs=s_1538;
					s_1543_obs=s_1543;
					s_1559_obs=s_1559;
					s_1565_obs=s_1565;
                    
					r_1166_obs=cell*Vmax_1166*(s_0565_substitute-s_0563)/Km0565_1166/(1+s_0565_substitute/Km0565_1166+1+s_0563/Km0563_1166-1);
					r_1697_obs=cell*Vmax_1697*s_0456/Km0456_1697/(1+s_0456/Km0456_1697);
					r_1762_obs=cell*Vmax_1762*s_0680/Km0680_1762/(1+s_0680/Km0680_1762);
					r_1106_obs=cell*Vmax_1106*s_0362/Km0362_1106/(1+s_0362/Km0362_1106);
					r_1172_obs=cell*Vmax_1172*s_0765/Km0765_1172/(1+s_0765/Km0765_1172);
					r_2079_obs=cell*Vmax_2079*s_1520/Km1520_2079/(1+s_1520/Km1520_2079);
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon_B1(void* data, realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
}
