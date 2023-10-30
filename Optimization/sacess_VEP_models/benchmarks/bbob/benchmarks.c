
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "benchmarksdeclare.h"
#include "benchmarkshelper.h"
#include "bbobStructures.h"

#define NHIGHPEAKS21 101
#define NHIGHPEAKS22 21

static double * tmpvect;
static double * tmx;
static double ** rotation;
static double ** rot2;
static double ** linearTF;
static double * peaks21;
static double * peaks22;
static int * rperm;
static int * rperm21;
static int * rperm22;
static double ** Xlocal;
static double ** Xlocal21;
static double ** Xlocal22;
static double ** arrScales;
static double ** arrScales21;
static double ** arrScales22;

/*
 * Noiseless functions testbed. All functions are ranged in [-5, 5]^DIM.
 */
/*isInitDone status changes when either DIM or trialid change.*/
/*it also changes when a new initialisation has been done (since it rewrites the values of Xopt, Fopt...)*/

TwoDoubles f1BBOB(double* x) {
    /*Sphere function*/

    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 1;
    double Fadd, r, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        isInitDone = 1;
    }

    Fadd = Fopt;
    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        r = x[i] - Xopt[i];
        Ftrue += r * r;
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    res.Ftrue = Ftrue;
    res.Fval = Fval;
    return res;
}

TwoDoubles f2BBOB(double* x) {
    /* separable ellipsoid with monotone transformation, condition 1e6*/

    int i, rseed; /*Loop over dim*/
    static double condition = 1e6;
    static unsigned int funcId = 2;
    double Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        
        computeXopt(rseed, DIM);
        isInitDone = 1;
    }

    Fadd = Fopt;

    for (i = 0; i < DIM; i++)
    {
        tmx[i] = x[i] - Xopt[i];
    }

    monotoneTFosc(tmx);

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += powl(condition, ((double)i)/((double)(DIM-1))) * tmx[i] * tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f3BBOB(double* x) {
    /* Rastrigin with monotone transformation separable "condition" 10*/
    int i, rseed; /*Loop over dim*/

    static unsigned int funcId = 3;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp, tmp2, Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        isInitDone = 1;
    }

    Fadd = Fopt;
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = x[i] - Xopt[i];
    }

    monotoneTFosc(tmx);
    for (i = 0; i < DIM; i++)
    {
        tmp = ((double)i)/((double)(DIM-1));
        if (tmx[i] > 0)
            tmx[i] = powl(tmx[i], 1 + beta * tmp * sqrtl(tmx[i]));
        tmx[i] = powl(sqrtl(condition), tmp) * tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < DIM; i++)
    {
        tmp += cosl(2*M_PI*tmx[i]);
        tmp2 += tmx[i]*tmx[i];
    }
    Ftrue = 10 * (DIM - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f4BBOB(double* x) {
    /* skew Rastrigin-Bueche, condition 10, skew-"condition" 100*/

    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 4;
    static double condition = 10.;
    static double alpha = 100.;
    double tmp, tmp2, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = 3 + 10000 * trialid; /* Not the same as before.*/
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        for (i = 0; i < DIM; i += 2)
            Xopt[i] = fabsl(Xopt[i]); /*Skew*/
        isInitDone = 1;
    }
    Fadd = Fopt;

    for (i = 0; i < DIM; i++) {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
            Fpen += tmp * tmp;
    }
    Fpen *= 1e2;
    Fadd += Fpen;

    for (i = 0; i < DIM; i++)
    {
        tmx[i] = x[i] - Xopt[i];
    }

    monotoneTFosc(tmx);
    for (i = 0; i < DIM; i++)
    {
        if (i % 2 == 0 && tmx[i] > 0)
            tmx[i] = sqrtl(alpha) * tmx[i];
        tmx[i] = powl(sqrtl(condition), ((double)i)/((double)(DIM-1))) * tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < DIM; i++)
    {
        tmp += cosl(2*M_PI*tmx[i]);
        tmp2 += tmx[i]*tmx[i];
    }
    Ftrue = 10 * (DIM - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f5BBOB(double* x) {
    /* linear slope*/
    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 5;
    static double alpha = 100.;
    static double Fadd; /*Treatment is different from other functions.*/
    double tmp, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        Fadd = Fopt;
        computeXopt(rseed, DIM);
        for (i = 0; i < DIM; i ++)
        {
            tmp = powl(sqrtl(alpha), ((double)i)/((double)(DIM-1)));
            if (Xopt[i] > 0)
            {
                Xopt[i] = 5.;
            }
            else if (Xopt[i] < 0)
            {
                Xopt[i] = -5.;
            }
            Fadd += 5. * tmp;
        }
        isInitDone = 1;
    }

    /* BOUNDARY HANDLING*/
    /* move "too" good coordinates back into domain*/
    for (i = 0; i < DIM; i++) {
        if ((Xopt[i] == 5.) && (x[i] > 5))
            tmx[i] = 5.;
        else if ((Xopt[i] == -5.) && (x[i] < -5))
            tmx[i] = -5.;
        else
            tmx[i] = x[i];
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        if (Xopt[i] > 0) {
            Ftrue -= powl(sqrtl(alpha), ((double)i)/((double)(DIM-1))) * tmx[i];
        } else {
            Ftrue += powl(sqrtl(alpha), ((double)i)/((double)(DIM-1))) * tmx[i];
        }
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f6BBOB(double* x) {
    /* attractive sector function*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 6;
    static double alpha = 100.;
    double Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        static double condition = 10.;
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        /* decouple scaling from function definition*/
        for (i = 0; i < DIM; i ++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++) {
                    linearTF[i][j] += rotation[i][k] * powl(sqrtl(condition), ((double)k)/((double)(DIM-1))) * rot2[k][j];
                }
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {

        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += linearTF[i][j] * (x[j] - Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        if (tmx[i] * Xopt[i] > 0)
            tmx[i] *= alpha;
        Ftrue += tmx[i] * tmx[i];
    }

    /*MonotoneTFosc...*/
    if (Ftrue > 0)
    {
        Ftrue = powl(expl(logl(Ftrue)/0.1 + 0.49*(sinl(logl(Ftrue)/0.1) + sinl(0.79*logl(Ftrue)/0.1))), 0.1);
    }
    else if (Ftrue < 0)
    {
        Ftrue = -powl(expl(logl(-Ftrue)/0.1 + 0.49*(sinl(0.55 * logl(-Ftrue)/0.1) + sinl(0.31*logl(-Ftrue)/0.1))), 0.1);
    }
    Ftrue = powl(Ftrue, 0.9);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f7BBOB(double* x) {
    /* step-ellipsoid, condition 100*/

    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {

        tmpvect[i] = 0.;
        tmp = sqrtl(powl(condition/10., ((double)i)/((double)(DIM-1))));
        for (j = 0; j < DIM; j++) {
            tmpvect[i] += tmp * rot2[i][j] * (x[j] - Xopt[j]);
        }

    }
    x1 = tmpvect[0];

    for (i = 0; i < DIM; i++) {
        if (fabsl(tmpvect[i]) > 0.5)
            tmpvect[i] = round1(tmpvect[i]);
        else
            tmpvect[i] = round1(alpha * tmpvect[i])/alpha;
    }

    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += rotation[i][j] * tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += powl(condition, ((double)i)/((double)(DIM-1))) * tmx[i] * tmx[i];
    }
    Ftrue = 0.1 * fmax1(1e-4 * fabsl(x1), Ftrue);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f8BBOB(double* x) {
    /* Rosenbrock, non-rotated*/
    static unsigned int funcId = 8;
    int i, rseed; /*Loop over dim*/
    static double scales;
    double tmp, Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;

        scales = fmax1(1., sqrtl((double)DIM) / 8.);
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        for (i = 0; i < DIM; i ++)
            Xopt[i] *= 0.75;
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp = (tmx[i] * tmx[i] - tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < DIM - 1; i ++)
    {
        tmp = (tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f9BBOB(double* x) {
    /* Rosenbrock, rotated*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 9;
    double scales, tmp, Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        /* computeXopt(rseed, DIM);*/
        computeRotation(rotation, rseed, DIM);
        scales = fmax1(1., sqrtl((double)DIM) / 8.);
        for (i = 0; i < DIM; i ++)
        {
            for (j = 0; j < DIM; j++)
                linearTF[i][j] = scales * rotation[i][j];
        }
/*         for (i = 0; i < DIM; i++)
           {
               Xopt[i] = 0.;
               for (j = 0; j < DIM; j++)
               {
                   Xopt[i] += linearTF[j][i] * 0.5/scales/scales;
                   //computed only if Xopt is returned which is not the case at this point.
               }
            }*/
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            tmx[i] += linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp = (tmx[i] * tmx[i] - tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < DIM - 1; i ++)
    {
       tmp = (tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f10BBOB(double* x) {
    /* ellipsoid with monotone transformation, condition 1e6*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 10;
    static double condition = 1e6;
    double Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
    }

    monotoneTFosc(tmx);
    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += powl(condition, ((double)i)/((double)(DIM-1))) * tmx[i] * tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f11BBOB(double* x) {
    /* discus (tablet) with monotone transformation, condition 1e6*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 11;
    static double condition = 1e6;
    double Fadd, Fval, Ftrue;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
    }

    monotoneTFosc(tmx);

    /* COMPUTATION core*/
    Ftrue = condition * tmx[0] * tmx[0];
    for (i = 1; i < DIM; i++)
    {
        Ftrue += tmx[i] * tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f12BBOB(double* x) {
    /* bent cigar with asymmetric space distortion, condition 1e6*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 12;
    static double condition = 1e6;
    static double beta = 0.5;
    double Fadd, Fval, Ftrue;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed + 1000000, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmpvect[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
        if (tmpvect[i] > 0)
        {
            tmpvect[i] = powl(tmpvect[i], 1 + beta * ((double)i)/((double)(DIM-1)) * sqrtl(tmpvect[i]));
        }
    }

    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += rotation[i][j] * tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    Ftrue = tmx[0] * tmx[0];
    for (i = 1; i < DIM; i++)
    {
        Ftrue += condition * tmx[i] * tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f13BBOB(double* x) {
    /* sharp ridge*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 13;
    static double condition = 10.;
    static double alpha = 100.;
    double Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++)
                {
                    linearTF[i][j] += rotation[i][k] * powl(sqrtl(condition), ((double)k)/((double)(DIM-1))) * rot2[k][j];
                }
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += linearTF[i][j] * (x[j] - Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 1; i < DIM; i++)
    {
        Ftrue += tmx[i] * tmx[i];
    }
    Ftrue = alpha * sqrtl(Ftrue);
    Ftrue += tmx[0] * tmx[0];

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f14BBOB(double* x) {
    /* sum of different powlers, between x^2 and x^6*/
    int i, j, rseed;
    static unsigned int funcId = 14;
    static double alpha = 4.;
    double Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++) {
            tmx[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += powl(fabsl(tmx[i]), 2. + alpha * ((double)i)/((double)(DIM-1)));
    }
    Ftrue = sqrtl(Ftrue);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f15BBOB(double* x) {
    /* Rastrigin with asymmetric non-linear distortion, "condition" 10*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 15;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp = 0., tmp2 = 0., Fadd, Fval, Ftrue;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++) {
                    linearTF[i][j] += rotation[i][k] * powl(sqrtl(condition), ((double)k)/((double)(DIM-1))) * rot2[k][j];
                }
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmpvect[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
    }

    monotoneTFosc(tmpvect);
    for (i = 0; i < DIM; i++)
    {
        if (tmpvect[i] > 0)
            tmpvect[i] = powl(tmpvect[i], 1 + beta * ((double)i)/((double)(DIM-1)) * sqrtl(tmpvect[i]));
    }
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += linearTF[i][j] * tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp += cosl(2. * M_PI * tmx[i]);
        tmp2 += tmx[i] * tmx[i];
    }
    Ftrue = 10. * ((double)DIM - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f16BBOB(double* x) {
    /* Weierstrass, condition 100*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 16;
    static double condition = 100.;
    static double aK[12];
    static double bK[12];
    static double F0;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++) {
                    linearTF[i][j] += rotation[i][k] * powl(1./sqrtl(condition), ((double)k)/((double)(DIM-1))) * rot2[k][j];
                }
            }
        }

        F0 = 0.;
        for (i = 0; i < 12; i ++) /* number of summands, 20 in CEC2005, 10/12 saves 30% of time*/
        {
            aK[i] = powl(0.5, (double)i);
            bK[i] = powl(3., (double)i);
            F0 += aK[i] * cosl(2 * M_PI * bK[i] * 0.5);
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10./(double)DIM * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmpvect[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
    }

    monotoneTFosc(tmpvect);
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += linearTF[i][j] * tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = 0.;
        for (j = 0; j < 12; j++)
        {
            tmp += cosl(2 * M_PI * (tmx[i] + 0.5) * bK[j]) * aK[j];
        }
        Ftrue += tmp;
    }
    Ftrue = 10. * powl(Ftrue/(double)DIM - F0, 3.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f17BBOB(double* x) {
    /* Schaffers F7 with asymmetric non-linear transformation, condition 10*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmpvect[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
        if (tmpvect[i] > 0)
            tmpvect[i] = powl(tmpvect[i], 1 + beta * ((double)i)/((double)(DIM-1)) * sqrtl(tmpvect[i]));
    }

    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        tmp = powl(sqrtl(condition), ((double)i)/((double)(DIM-1)));
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += tmp * rot2[i][j] * tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp = tmx[i] * tmx[i] + tmx[i+1] * tmx[i+1];
        Ftrue += powl(tmp, 0.25) * (powl(sinl(50 * powl(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = powl(Ftrue/(double)(DIM - 1), 2.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f18BBOB(double* x) {
    /* Schaffers F7 with asymmetric non-linear transformation, condition 1000*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 18;
    static double condition = 1e3;
    static double beta = 0.5;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = 17 + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmpvect[i] += rotation[i][j] * (x[j] - Xopt[j]);
        }
        if (tmpvect[i] > 0)
            tmpvect[i] = powl(tmpvect[i], 1. + beta * ((double)i)/((double)(DIM-1)) * sqrtl(tmpvect[i]));
    }

    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        tmp = powl(sqrtl(condition), ((double)i)/((double)(DIM-1)));
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += tmp * rot2[i][j] * tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp = tmx[i] * tmx[i] + tmx[i+1] * tmx[i+1];
        Ftrue += powl(tmp, 0.25) * (powl(sinl(50. * powl(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = powl(Ftrue/(double)(DIM - 1), 2.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f19BBOB(double* x) {
    /* F8F2 sum of Griewank-Rosenbrock 2-D blocks*/
    int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 19;
    double scales, F2, tmp = 0., tmp2, Fadd, Fval, Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        /* computeXopt(rseed, DIM); Xopt is not used.*/
        scales = fmax1(1., sqrtl((double)DIM) / 8.);
        computeRotation(rotation, rseed, DIM);
        for (i = 0; i < DIM; i ++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = scales * rotation[i][j];
            }
        }
        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.;
            for (j = 0; j < DIM; j++)
            {
                Xopt[i] += linearTF[j][i] * 0.5/scales/scales;
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            tmx[i] += linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp2 = tmx[i] * tmx[i] -tmx[i+1];
        F2 = 100. * tmp2 * tmp2;
        tmp2 = 1 - tmx[i];
        F2 += tmp2 * tmp2;
        tmp += F2 / 4000. - cosl(F2);
    }
    Ftrue = 10. + 10. * tmp / (double)(DIM - 1);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f20BBOB(double* x) {
    /* Schwefel with tridiagonal variable transformation*/
    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 20;
    static double condition = 10.;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        unif(tmpvect, DIM, rseed);
        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.5 * 4.2096874633;
            if (tmpvect[i] - 0.5 < 0)
                Xopt[i] *= -1.;
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmpvect[i] = 2. * x[i];
        if (Xopt[i] < 0.)
            tmpvect[i] *= -1.;
    }

    tmx[0] = tmpvect[0];
    for (i = 1; i < DIM; i++)
    {
        tmx[i] = tmpvect[i] + 0.25 * (tmpvect[i-1] - 2. * fabsl(Xopt[i-1]));
    }

    for (i = 0; i < DIM; i++)
    {
        tmx[i] -= 2 * fabsl(Xopt[i]);
        tmx[i] *= powl(sqrtl(condition), ((double)i)/((double)(DIM-1)));
        tmx[i] = 100. * (tmx[i] + 2 * fabsl(Xopt[i]));
    }

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(tmx[i]) - 500.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 0.01 * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += tmx[i] * sinl(sqrtl(fabsl(tmx[i])));
    }
    Ftrue = 0.01 * ((418.9828872724339) - Ftrue / (double)DIM);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f21BBOB(double* x) {
    /* Gallagher with 101 Gaussian peaks, condition up to 1000, one global rotation*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 21;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
    static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)DIM;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeRotation(rotation, rseed, DIM);
        peaks = peaks21;
        unif(peaks, NHIGHPEAKS21 - 1, rseed);
        rperm = rperm21;
        for (i = 0; i < NHIGHPEAKS21 - 1; i++)
            rperm[i] = i;
        qsort(rperm, NHIGHPEAKS21 - 1, sizeof(int), compare_doubles);

        /* Random permutation*/
        arrCondition[0] = sqrtl(maxcondition);
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS21; i++)
        {
            arrCondition[i] = powl(maxcondition, (double)(rperm[i-1])/((double)(NHIGHPEAKS21-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        arrScales = arrScales21;
        for (i = 0; i < NHIGHPEAKS21; i++)
        {
            unif(peaks, DIM, rseed + 1000 * i);
            for (j = 0; j < DIM; j++)
                rperm[j] = j;
            qsort(rperm, DIM, sizeof(int), compare_doubles);
            for (j = 0; j < DIM; j++)
            {
                arrScales[i][j] = powl(arrCondition[i], ((double)rperm[j])/((double)(DIM-1)) - 0.5);
            }
        }

        unif(peaks, DIM * NHIGHPEAKS21, rseed);
        Xlocal = Xlocal21;
        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.8 * (10. * peaks[i] -5.);
            for (j = 0; j < NHIGHPEAKS21; j++)
            {
                Xlocal[i][j] = 0.;
                for (k = 0; k < DIM; k++)
                {
                    Xlocal[i][j] += rotation[i][k] * (10. * peaks[j * DIM + k] -5.);
                }
                if (j == 0)
                    Xlocal[i][j] *= 0.8;
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp = (tmx[j] - Xlocal[j][i]);
            tmp2 += arrScales[i][j] * tmp * tmp;
        }
        tmp2 = peakvalues[i] * expl(fac * tmp2);
        f = fmax1(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = logl(f)/a;
        Ftrue = powl(expl(Ftrue + 0.49*(sinl(Ftrue) + sinl(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = logl(-f)/a;
        Ftrue = -powl(expl(Ftrue + 0.49*(sinl(0.55 * Ftrue) + sinl(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f22BBOB(double* x) {
    /* Gallagher with 21 Gaussian peaks, condition up to 1000, one global rotation*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 22;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS22];
    static double peakvalues[NHIGHPEAKS22];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)DIM;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeRotation(rotation, rseed, DIM);
        peaks = peaks22;
        unif(peaks, NHIGHPEAKS22 - 1, rseed);
        rperm = rperm22;
        for (i = 0; i < NHIGHPEAKS22 - 1; i++)
            rperm[i] = i;
        qsort(rperm, NHIGHPEAKS22 - 1, sizeof(int), compare_doubles);
        /* Random permutation*/
        arrCondition[0] = maxcondition;
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS22; i++)
        {
            arrCondition[i] = powl(maxcondition, (double)(rperm[i-1])/((double)(NHIGHPEAKS22-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS22-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        arrScales = arrScales22;
        for (i = 0; i < NHIGHPEAKS22; i++)
        {
            unif(peaks, DIM, rseed + 1000 * i);
            for (j = 0; j < DIM; j++)
                rperm[j] = j;
            qsort(rperm, DIM, sizeof(int), compare_doubles);
            for (j = 0; j < DIM; j++)
            {
                arrScales[i][j] = powl(arrCondition[i], ((double)rperm[j])/((double)(DIM-1)) - 0.5);
            }
        }

        unif(peaks, DIM * NHIGHPEAKS22, rseed);
        Xlocal = Xlocal22;
        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.8 * (9.8 * peaks[i] -4.9);
            for (j = 0; j < NHIGHPEAKS22; j++)
            {
                Xlocal[i][j] = 0.;
                for (k = 0; k < DIM; k++)
                {
                    Xlocal[i][j] += rotation[i][k] * (9.8 * peaks[j * DIM + k] -4.9);
                }
                if (j == 0)
                    Xlocal[i][j] *= 0.8;
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmx[i] += rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS22; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp = (tmx[j] - Xlocal[j][i]);
            tmp2 += arrScales[i][j] * tmp * tmp;
        }
        tmp2 = peakvalues[i] * expl(fac * tmp2);
        f = fmax1(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = logl(f)/a;
        Ftrue = powl(expl(Ftrue + 0.49*(sinl(Ftrue) + sinl(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = logl(-f)/a;
        Ftrue = -powl(expl(Ftrue + 0.49*(sinl(0.55 * Ftrue) + sinl(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    /* free(Xopt); //Not used!*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f23BBOB(double* x) {
    /* Katsuura function*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 23;
    static double condition = 100.;
    double Fadd = 0., Fpen = 0., tmp, Ftrue = 0., arr, prod = 1., tmp2, Fval;
    double *ptmx, *plinTF, *ptmp;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++)
                {
                    linearTF[i][j] += rotation[i][k] * powl(sqrtl(condition), ((double)k)/(double)(DIM - 1)) * rot2[k][j];
                }
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    /* write rotated difference vector into tmx*/
    for (j = 0; j < DIM; j++)  /* store difference vector*/
        tmpvect[j] = x[j] - Xopt[j];
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.;
        ptmx = &tmx[i];
        plinTF = linearTF[i];
        ptmp = tmpvect;
        for (j = 0; j < DIM; j++) {
            *ptmx += *plinTF++ * *ptmp++;
        }
    }

/*     for (i = 0; i < DIM; i++) {
           tmx[i] = 0.;
           for (j = 0; j < DIM; j++) {
               tmx[i] += linearTF[i][j] * (x[j] - Xopt[j]);
           }
       }*/

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = 0.;
        for (j = 1; j < 33; j++)
        {
            tmp2 = powl(2., (double)j);
            arr = tmx[i] * tmp2;
            tmp += fabsl(arr - round1(arr)) / tmp2;
        }
        tmp = 1. + tmp * (double)(i + 1);
        prod *= tmp;
    }
    Ftrue = 10./(double)DIM/(double)DIM * (-1. + powl(prod, 10./powl((double)DIM, 1.2)));
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f24BBOB(double* x) {
    /* Lunacek bi-Rastrigin, condition 100*/
    /* in PPSN 2008, Rastrigin part rotated and scaled*/
    int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 24;
    static double condition = 100.;
    static double mu1 = 2.5;
    double Fadd, Fpen = 0., tmp, Ftrue = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0., Fval;
    double s = 1. - 0.5 / (sqrtl((double)(DIM + 20)) - 4.1);
    static double d = 1.;
    double mu2 = -sqrtl((mu1 * mu1 - d) / s);
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = funcId + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeRotation(rotation, rseed + 1000000, DIM);
        computeRotation(rot2, rseed, DIM);
        gauss(tmpvect, DIM, rseed);
        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.5 * mu1;
            if (tmpvect[i] < 0.)
                Xopt[i] *= -1.;
        }

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                linearTF[i][j] = 0.;
                for (k = 0; k < DIM; k++) {
                    linearTF[i][j] += rotation[i][k] * powl(sqrtl(condition), ((double)k)/((double)(DIM-1))) * rot2[k][j];
                }
            }
        }
        isInitDone = 1;
    }
    Fadd = Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabsl(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 1e4 * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        tmx[i] = 2. * x[i];
        if (Xopt[i] < 0.)
            tmx[i] *= -1.;
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < DIM; i++)
    {
        tmp2 += (tmx[i] - mu1) * (tmx[i] - mu1);
        tmp3 += (tmx[i] - mu2) * (tmx[i] - mu2);
        tmp4 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp4 += linearTF[i][j] * (tmx[j] - mu1);
        }
        tmp += cosl(2 * M_PI * tmp4);
    }
    Ftrue = fmin1(tmp2, d * (double)DIM + s * tmp3) + 10. * ((double)DIM - tmp);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}


void initbenchmarks(void) /*ParamStruct params)*/
{

    int i;
    tmpvect = (double*)malloc(sizeof(double) * DIM);
    tmx = (double*)malloc(sizeof(double) * DIM);
    rotation = (double**)malloc(sizeof(double*) * DIM);
    rot2 = (double**)malloc(sizeof(double*) * DIM);
    linearTF = (double**)malloc(sizeof(double*) * DIM);
    peaks21 = (double*)malloc(sizeof(double) * DIM * NHIGHPEAKS21);
    rperm21 = (int*)malloc(sizeof(int) * (int)fmax1(DIM, NHIGHPEAKS21 - 1));
    arrScales21 = (double**)malloc(sizeof(double*) * NHIGHPEAKS21);
    Xlocal21 = (double**)malloc(sizeof(double*) * DIM);
    peaks22 = (double*)malloc(sizeof(double) * DIM * NHIGHPEAKS22);
    rperm22 = (int*)malloc(sizeof(int) * (int)fmax1(DIM,  NHIGHPEAKS22 - 1));
    arrScales22 = (double**)malloc(sizeof(double*) * NHIGHPEAKS22);
    Xlocal22 = (double**)malloc(sizeof(double*) * DIM);

    for (i = 0; i < DIM; i++)
    {
        rotation[i] = (double*)malloc(sizeof(double) * DIM);
        rot2[i] = (double*)malloc(sizeof(double) * DIM);
        linearTF[i] = (double*)malloc(sizeof(double) * DIM);
        Xlocal21[i] = (double*)malloc(sizeof(double) * NHIGHPEAKS21);
        Xlocal22[i] = (double*)malloc(sizeof(double) * NHIGHPEAKS22);
    }
    for (i = 0; i < NHIGHPEAKS21; i++)
        arrScales21[i] = (double*)malloc(sizeof(double) * DIM);
    for (i = 0; i < NHIGHPEAKS22; i++)
        arrScales22[i] = (double*)malloc(sizeof(double) * DIM);

    return;
}

void finibenchmarks(void)
{
    free(tmpvect);
    free(tmx);
    freeStarStar(rotation, DIM);
    freeStarStar(rot2, DIM);
    freeStarStar(linearTF, DIM);
    free(peaks21);
    free(rperm21);
    freeStarStar(arrScales21, NHIGHPEAKS21);
    freeStarStar(Xlocal21, DIM);
    free(peaks22);
    free(rperm22);
    freeStarStar(arrScales22, NHIGHPEAKS22);
    freeStarStar(Xlocal22, DIM);
    
    tmpvect = NULL;
    tmx =  NULL;
    rotation =  NULL;
    rot2 =  NULL;
    linearTF =  NULL;
    peaks21 =  NULL;
    rperm21 =  NULL;
    arrScales21 =  NULL;
    Xlocal21 =  NULL;
    peaks22 =  NULL;
    rperm22 =  NULL;
    arrScales22 =  NULL;
    Xlocal22 =  NULL;
    return;
}

bbobFunction handles[24] = { &f1BBOB, &f2BBOB, &f3BBOB, &f4BBOB, &f5BBOB, &f6BBOB, &f7BBOB, &f8BBOB, &f9BBOB, &f10BBOB, 
&f11BBOB, &f12BBOB, &f13BBOB, &f14BBOB, &f15BBOB, &f16BBOB, &f17BBOB, &f18BBOB, &f19BBOB, &f20BBOB, &f21BBOB, &f22BBOB, &f23BBOB, &f24BBOB};
unsigned int handlesLength = 24;

