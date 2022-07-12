functions {
    /* form of epileptor x1 dynamics */
    real f1(real x, real z) {
        if (x<0) {
            return x*x*x + 2*x*x;
        } else {
            real z4;
            z4 <- z - 4.0;
            return -0.6*z4*z4*x + 5*x*x;
        }
    }

    real gamma_a_from_u_v(real u, real v) {
        return (u*sqrt(u * u + 4 * v) + u * u + 2 * v) / (2 * v);
    }

    real gamma_b_from_a_u(real a, real u) {
        return (a - 1.0) / u;
    }

}

data {
    int nn;
    int nt;
    int ns;
    real I1;
    real tau0;
    real dt;
    matrix[ns, nn] gain;
    matrix[nt, ns] seeg_log_power;
    vector[nn] Ic;
    matrix<lower=0.0>[nn, nn] SC;
    real SC_var; // prior on variance of connectivity strengths
    real K_lo;
    real K_u;
    real K_v;
    real x0_lo;
    real x0_hi;
    real eps_hi;
    real sig_hi;
    real zlim[2];
}

transformed data {
    real xlim[2];
    // real zlim[2];
    real K_hi;
    xlim[1] <- -2.5;
    xlim[2] <- 1.5;
    /*
    zlim[1] <- 2.0;
    zlim[2] <- 5.0;
    */
    K_hi <- 10.0;
}

parameters {
    matrix<lower=xlim[1], upper=xlim[2]>[nn, nt] x;
    matrix<lower=zlim[1], upper=zlim[2]>[nn, nt] z;
    vector<lower=x0_lo, upper=x0_hi>[nn] x0;
    real<lower=K_lo, upper=K_hi> K;
    real<lower=0.0, upper=eps_hi> eps;
    real<lower=2.0, upper=10.0> amp;
    real<lower=1.25, upper=2.0> offset;
    real<lower=0.0, upper=sig_hi> sig;
    matrix<lower=0.0>[nn, nn] FC;
    real<lower=0.0, upper=1.0> tt;
}

model {
    real dx;
    real dz;
    real gx;
    real fc_a;
    real fc_b;
    real fc_u;
    real fc_v;
    real K_a;
    real K_b;

    K_a <- gamma_a_from_u_v(K_u, K_v);
    K_b <- gamma_b_from_a_u(K_a, K_u);
    K ~ gamma(K_a, K_b) T[K_lo, K_hi];

    /* functional connectivity */
    for (i in 1:nn) {
        for (j in 1:nn) {
            if (i>=j) {
                if ((i==j) || (SC[i, j]==0.0)) {
                    fc_u <- 1e-6;
                    fc_v <- 1e-3;
                } else {
                    fc_u <- SC[i, j];
                    fc_v <- fc_u * SC_var;
                }
                fc_a <- gamma_a_from_u_v(fc_u, fc_v);
                fc_b <- gamma_b_from_a_u(fc_a, fc_u);
                FC[i, j] ~ gamma(fc_a, fc_b);
                FC[j, i] ~ gamma(fc_a, fc_b);
            }
        }
    }

    /* integrate & predict */
    for (t in 1:(nt-1)) {
        for (i in 1:nn) {
            gx <- Ic[i] * (-1.8 - x[i, t]);
            for (j in 1:nn)
                if (i!=j)
                    gx <- gx + FC[i, j]*(x[j, t] - x[i, t]);
            // dx <- 1.0 - f1(x[i, t], z[i, t]) - z[i, t] + I1;
            dx <- 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
            dz <- (1/tau0)*(4*(x[i, t] - x0[i]) - z[i, t] - K*gx);
            x[i, t+1] ~ normal(x[i, t] + dt*tt*dx, sig); // T[xlim[1], xlim[2]];
            z[i, t+1] ~ normal(z[i, t] + dt*tt*dz, sig); // T[zlim[1], zlim[2]];
        }
        seeg_log_power[t] ~ normal(amp*(log(gain*exp(col(x, t))) + offset), eps);
    }
}
