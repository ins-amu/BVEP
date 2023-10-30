
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "benchmarksdeclare.h"
#include "benchmarkshelper.h"
#include "bbobStructures.h"

#define NHIGHPEAKS21 101

static double * tmx;
static double * tmpvect;
static double ** rotation;
static double ** rot2;
static double ** linearTF;
static double * peaks21;
static int * rperm;
static int * rperm21;
static double ** Xlocal;
static double ** Xlocal21;
static double ** arrScales;
static double ** arrScales21;

/*
 * Noisy functions testbed. All functions are ranged in [-5, 5]^DIM.
 */

TwoDoubles f101(double* x) {
    /*sphere with moderate Gauss noise*/
    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 101;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }

    Fval = FGauss(Ftrue, 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f102(double* x) {
    /* sphere with moderate uniform noise*/
    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 102;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.01 * (0.49 + 1./DIM), 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f103(double* x) {
    /* sphere with moderate Cauchy noise*/
    int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 103;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 0.01, 0.05);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f104(double* x) {
    /* Rosenbrock non-rotated with moderate Gauss noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 104;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    Ftrue = 0.;
    for (i = 0; i < DIM - 1; i++)
    {
        tmp = (tmx[i] * tmx[i] - tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < DIM - 1; i ++)
    {
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FGauss(Ftrue, 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f105(double* x) {
    /* Rosenbrock non-rotated with moderate uniform noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 105;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
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
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FUniform(Ftrue, 0.01 * (0.49 + 1./DIM), 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f106(double* x) {
    /* Rosenbrock non-rotated with moderate Cauchy noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 106;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
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
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FCauchy(Ftrue, 0.01, 0.05);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f107(double* x) {
    /*sphere with Gauss noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 107;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FGauss(Ftrue, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f108(double* x) {
    /*sphere with uniform noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 108;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f109(double* x) {
    /*sphere with Cauchy noise*/
    int i, rseed; /*Loop over dim*/
    static int funcId = 109;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
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
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = (x[i] - Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 1., 0.2);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f110(double* x) {
    /* Rosenbrock non-rotated with Gauss noise*/
    int i, rseed; /*Loop over dim */
    static int funcId = 110;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
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
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FGauss(Ftrue, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f111(double* x) {
    /* Rosenbrock non-rotated with moderate uniform noise*/
    int i, rseed; /*Loop over dim */
    static int funcId = 111;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
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
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f112(double* x) {
    /* Rosenbrock non-rotated with moderate Cauchy noise*/
    int i, rseed; /*Loop over dim */
    static int funcId = 112;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        scales = fmax(1., sqrtl((double)DIM) / 8.);
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = scales * (x[i] - 0.75 * Xopt[i]) + 1;
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
        tmp = (tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 1., 0.2);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f113(double* x) {
    /* step-ellipsoid with gauss noise, condition 100*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 113;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        {
            tmpvect[i] = round1(tmpvect[i]);
        }
        else
        {
            tmpvect[i] = round1(alpha * tmpvect[i])/alpha;
        }
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
    Ftrue = 0.1 * fmax(1e-4 * fabsl(x1), Ftrue);
    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f114(double* x) {
    /* step-ellipsoid with uniform noise, condition 100*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 114;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        {
            tmpvect[i] = round1(tmpvect[i]);
        }
        else
        {
            tmpvect[i] = round1(alpha * tmpvect[i])/alpha;
        }
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
    Ftrue = 0.1 * fmax(1e-4 * fabsl(x1), Ftrue);
    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f115(double* x) {
    /* step-ellipsoid with Cauchy noise, condition 100*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 115;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        {
            tmpvect[i] = round1(tmpvect[i]);
        }
        else
        {
            tmpvect[i] = round1(alpha * tmpvect[i])/alpha;
        }
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
    Ftrue = 0.1 * fmax(1e-4 * fabsl(x1), Ftrue);
    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f116(double* x) {
    /* ellipsoid with Gauss noise, monotone x-transformation, condition 1e4*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 116;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f117(double* x) {
    /* ellipsoid with uniform noise, monotone x-transformation, condition 1e4*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 117;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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

    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f118(double* x) {
    /* ellipsoid with Cauchy noise, monotone x-transformation, condition 1e4*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 118;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f119(double* x) {
    /* sum of different powlers with Gauss Noise, between x^2 and x^6*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 119;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(fabsl(tmx[i]), 2 + alpha * ((double)i)/((double)(DIM-1)));
    }
    Ftrue = sqrtl(Ftrue);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f120(double* x) {
    /* sum of different powlers with uniform noise, between x^2 and x^6*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 120;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(fabsl(tmx[i]), 2 + alpha * ((double)i)/((double)(DIM-1)));
    }
    Ftrue = sqrtl(Ftrue);

    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f121(double* x) {
    /* sum of different powlers with seldom Cauchy Noise, between x^2 and x^6*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 121;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        computeXopt(rseed, DIM);
        computeRotation(rotation, rseed + 1000000, DIM);
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(fabsl(tmx[i]), 2 + alpha * ((double)i)/((double)(DIM-1)));
    }
    Ftrue = sqrtl(Ftrue);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f122(double* x) {
    /* Schaffers F7 with Gauss noise, with asymmetric non-linear transformation, condition 10*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 122;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(tmp, 0.25) * (powl(sinl(50. * powl(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = powl(Ftrue/(double)(DIM - 1), 2.);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f123(double* x) {
    /* Schaffers F7 with uniform noise, with asymmetric non-linear transformation, condition 10*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 123;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(tmp, 0.25) * (powl(sinl(50. * powl(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = powl(Ftrue/(double)(DIM - 1), 2.);

    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f124(double* x) {
    /* Schaffers F7 with seldom Cauchy noise, with asymmetric non-linear transformation, condition 10*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 124;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
        Ftrue += powl(tmp, 0.25) * (powl(sinl(50. * powl(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = powl(Ftrue/(double)(DIM - 1), 2.);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f125(double* x) {
    /* F8F2 sum of Griewank-Rosenbrock 2-D blocks with Gauss noise*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 125;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        /* computeXopt(rseed, DIM);*/
        scales = fmax(1., sqrtl((double)DIM) / 8.);
        computeRotation(rotation, rseed, DIM);
/*        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.;
            for (j = 0; j < DIM; j++)
            {
                Xopt[i] += rotation[j][i] * 0.5/scales;
                //computed only if Xopt is returned which is not the case at this point.
            }
        }*/
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            tmx[i] += scales * rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < DIM - 1; i++)
    {
        F2 = 100. * (tmx[i] * tmx[i] - tmx[i+1]) * (tmx[i] * tmx[i] - tmx[i+1]) + (1 - tmx[i]) * (1 - tmx[i]);
        tmp += F2 / 4000. - cosl(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(DIM - 1);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f126(double* x) {
    /* F8F2 sum of Griewank-Rosenbrock 2-D blocks with uniform noise*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 126;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        /* computeXopt(rseed, DIM);*/
        scales = fmax(1., sqrtl((double)DIM) / 8.);
        computeRotation(rotation, rseed, DIM);
/*        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.;
            for (j = 0; j < DIM; j++)
            {
                Xopt[i] += rotation[j][i] * 0.5/scales;
                //computed only if Xopt is returned which is not the case at this point.
            }
        }*/
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            tmx[i] += scales * rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < DIM - 1; i++)
    {
        F2 = 100. * (tmx[i] * tmx[i] - tmx[i+1]) * (tmx[i] * tmx[i] - tmx[i+1]) + (1 - tmx[i]) * (1 - tmx[i]);
        tmp += F2 / 4000. - cosl(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(DIM - 1);

    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f127(double* x) {
    /* F8F2 sum of Griewank-Rosenbrock 2-D blocks with seldom Cauchy noise*/
    int i, j, rseed; /*Loop over dim*/
    static int funcId = 127;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    TwoDoubles res;

    if (!isInitDone)
    {
        rseed = rrseed + 10000 * trialid;
        /*INITIALIZATION*/
        Fopt = computeFopt(funcId, trialid);
        /* computeXopt(rseed, DIM);*/
        scales = fmax(1., sqrtl((double)DIM) / 8.);
        computeRotation(rotation, rseed, DIM);
/*        for (i = 0; i < DIM; i++)
        {
            Xopt[i] = 0.;
            for (j = 0; j < DIM; j++)
            {
                Xopt[i] += rotation[j][i] * 0.5/scales;
                //computed only if Xopt is returned which is not the case at this point.
            }
        }*/
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
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            tmx[i] += scales * rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < DIM - 1; i++)
    {
        F2 = 100. * (tmx[i] * tmx[i] - tmx[i+1]) * (tmx[i] * tmx[i] - tmx[i+1]) + (1 - tmx[i]) * (1 - tmx[i]);
        tmp += F2 / 4000. - cosl(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(DIM - 1);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f128(double* x) {
    /* Gallagher with 101 Gaussian peaks with Gauss noise, condition up to 1000, one global rotation*/
    int i, j, k, rseed; /*Loop over dim*/
    static int funcId = 128;
    static int rrseed = 21;
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
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
            tmp2 += arrScales[i][j] * (tmx[j] - Xlocal[j][i]) * (tmx[j] - Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * expl(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
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

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;
    /* free(Xopt); //Not used!*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f129(double* x) {
    /* Gallagher with 101 Gaussian peaks with uniform noise, condition up to 1000, one global rotation*/
    int i, j, k, rseed; /*Loop over dim*/
    static int funcId = 129;
    static int rrseed = 21;
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
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
            tmp2 += arrScales[i][j] * (tmx[j] - Xlocal[j][i]) * (tmx[j] - Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * expl(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
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

    Fval = FUniform(Ftrue, 0.49 + 1./DIM, 1.);

    Ftrue += Fadd;
    Fval += Fadd;
    /* free(Xopt); //Not used!*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

TwoDoubles f130(double* x) {
    /* Gallagher with 101 Gaussian peaks with Cauchy noise, condition up to 1000, one global rotation*/
    int i, j, k, rseed; /*Loop over dim*/
    static int funcId = 130;
    static int rrseed = 21;
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
        rseed = rrseed + 10000 * trialid;
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
    Fadd += 100. * Fpen;

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
            tmp2 += arrScales[i][j] * (tmx[j] - Xlocal[j][i]) * (tmx[j] - Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * expl(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
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

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;
    /* free(Xopt); //Not used!*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

void initbenchmarksnoisy(void)
{
    int i;

    tmpvect = (double*)malloc(sizeof(double) * DIM);
    tmx = (double*)malloc(sizeof(double) * DIM);
    rotation = (double**)malloc(sizeof(double*) * DIM);
    rot2 = (double**)malloc(sizeof(double*) * DIM);
    linearTF = (double**)malloc(sizeof(double*) * DIM);
    peaks21 = (double*)malloc(sizeof(double) * DIM * NHIGHPEAKS21);
    rperm21 = (int*)malloc(sizeof(int) * fmax(DIM, NHIGHPEAKS21 - 1));
    arrScales21 = (double**)malloc(sizeof(double*) * NHIGHPEAKS21);
    Xlocal21 = (double**)malloc(sizeof(double*) * DIM);

    for (i = 0; i < DIM; i++)
    {
        rotation[i] = (double*)malloc(sizeof(double) * DIM);
        rot2[i] = (double*)malloc(sizeof(double) * DIM);
        linearTF[i] = (double*)malloc(sizeof(double) * DIM);
        Xlocal21[i] = (double*)malloc(sizeof(double) * NHIGHPEAKS21);
    }
    for (i = 0; i < NHIGHPEAKS21; i++)
        arrScales21[i] = (double*)malloc(sizeof(double) * DIM);

    return;
}

void finibenchmarksnoisy(void)
{
    free(tmpvect);
    free(tmx);
    freeStarStar(rotation, DIM);
    freeStarStar(rot2, DIM);
    freeStarStar(linearTF, DIM);
    free(peaks21);
    //free(rperm);    
    free(rperm21);
    //freeStarStar(arrScales, NHIGHPEAKS21);    
    freeStarStar(arrScales21, NHIGHPEAKS21);
    //freeStarStar(Xlocal, DIM);    
    freeStarStar(Xlocal21, DIM);
    return;
}


bbobFunction handlesNoisy[30] = { &f101, &f102, &f103, &f104, &f105, &f106, &f107, &f108, &f109, &f110, &f111, &f112, &f113, &f114, &f115, &f116, &f117, &f118, &f119, &f120, &f121, &f122, &f123, &f124, &f125, &f126, &f127, &f128, &f129, &f130};
unsigned int handlesNoisyLength = 30;

