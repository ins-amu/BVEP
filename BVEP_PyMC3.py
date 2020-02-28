"""
Written in INS Marseille.

"""
import pymc3 as pm
import numpy as np
import theano
import theano.tensor as tt

####################################################################################################

def dx(x, z, I1):
    dx_eval = 1 - x**3 - 2 * x**2 - z + 3.1
    return dx_eval

####################################################################################################

def dz(x, z, SC, K, eta, tau0):
    nn = x.size
    x_diff = np.repeat(x[:, np.newaxis], nn, axis=1) - x
    gx = K * SC * x_diff.T
    dz_eval = (1 / tau0) * (4 * (x - eta) - z - gx.sum(axis=1))
    return dz_eval

####################################################################################################

def step_ode(x_prev, z_prev, dt, SC, K, x0, I1, tau0):
    x_next = x_prev + dt*dx(x_prev, z_prev, I1)
    z_next = z_prev + dt*dz(x_prev, z_prev, SC, K, x0, tau0)
    return x_next, z_next

####################################################################################################

def step_sde(x_eta, z_eta, x_prev, z_prev, dt, SC, K, eta, I1, tau0, sig):
    x_next = x_prev + dt*dx(x_prev, z_prev, I1) + tt.sqrt(dt)*x_eta*sig
    z_next = z_prev + dt*dz(x_prev, z_prev, SC, K, eta, tau0)+ tt.sqrt(dt)*z_eta*sig
    return x_next, z_nex

####################################################################################################


class BVEP_cen:
    def __init__(self, consts, prior_mu, obs):
        self.consts = consts
        self.prior_mu = prior_mu
        self.obs = obs
        self.model = pm.Model()
        with self.model:
            BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
            eta = pm.Normal('eta', mu=self.prior_mu['eta'] , sd=1.0, shape=self.consts['nn'])
            amplitude = pm.TruncatedNormal('amplitude', mu=self.prior_mu['amplitude'] , sd=1.0, lower=0)
            offset = pm.TruncatedNormal('offset', mu=self.prior_mu['offset'], sd=1.0)
            K = pm.TruncatedNormal('K', mu=self.prior_mu['K'] , sd= 1., lower=0.0)
            eps = pm.TruncatedNormal('eps', mu=.0, sd=1.0, lower=0.0)
            sig = pm.TruncatedNormal('sig', mu=.0, sd=1.0, lower=0.0)
            x_init = pm.Normal('x_init', mu=self.prior_mu['x_init'],  sd=1., shape=self.consts['nn'])
            z_init = pm.Normal('z_init', mu=self.prior_mu['z_init'],  sd=1., shape=self.consts['nn'])

            # Cast constants in the model as tensors using theano shared variables
            time_step = theano.shared(self.consts['time_step'], 'time_step')
            SC = theano.shared(self.consts['SC'], 'SC')
            I1 = theano.shared(self.consts['I1'], 'I1')
            tau0 = theano.shared(self.consts['tau0'], 'tau0')

            output, updates = theano.scan(fn=step_ode, outputs_info=[x_init, z_init], non_sequences=[time_step, SC, K, eta, I1, tau0],  n_steps=self.consts['nt'])
            
            x_sym = output[0]
            z_sym = output[1]
            
            x = pm.Normal('x', mu=x_sym, sd=tt.sqrt(time_step)*sig, shape=(self.consts['nt'], self.consts['nn']))
            z = pm.Normal('z', mu=z_sym, sd=tt.sqrt(time_step)*sig, shape=(self.consts['nt'], self.consts['nn']))

            xhat = pm.Deterministic('xhat', amplitude * x + offset)

            xs = pm.Normal('xs', mu=xhat, sd=eps, shape=(self.consts['nt'], self.consts['nn']), observed=self.obs['xs'])

####################################################################################################

class BVEP_noncen:
    def __init__(self, consts, prior_mu, obs):
        self.consts = consts
        self.prior_mu = prior_mu
        self.obs = obs
        self.model = pm.Model()
        with self.model:
            BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
            eta_star = pm.Normal('eta_star', mu=0.0, sd=1.0,  shape=self.consts['nn'])
            eta = pm.Deterministic('eta', self.prior_mu['eta'] + eta_star)
            amplitude_star = pm.Normal('amplitude_star', mu=.0, sd=1.0)
            amplitude = pm.Deterministic('amplitude', self.prior_mu['amplitude'] + amplitude_star)
            offset_star = pm.Normal('offset_star', mu=0.0, sd=1.0)
            offset = pm.Deterministic('offset', self.prior_mu['offset']+ offset_star)
            K_star = pm.Normal('K_star', mu=0.0, sd=1.0)
            K = pm.Deterministic('K', self.prior_mu['K'] + K_star)
            eps = BoundedNormal('eps', mu=.0, sd=1.0)
            sig = BoundedNormal('sig', mu=.0, sd=1.0)
            
            x_init_star = pm.Normal('x_init_star',  mu=0, sd=1., shape=self.consts['nn'])
            x_init = pm.Deterministic('x_init', self.prior_mu['x_init'] + x_init_star)
            z_init_star = pm.Normal('z_init_star',  mu=0, sd=1.,  shape=self.consts['nn'])
            z_init = pm.Deterministic('z_init', self.prior_mu['z_init'] + z_init_star)
            x_eta = pm.Normal('x_eta', mu=0, sd=1.0, shape=(self.consts['nt'], self.consts['nn']))
            z_eta = pm.Normal('z_eta', mu=0, sd=1.0, shape=(self.consts['nt'], self.consts['nn']))

            # Cast constants in the model as tensors using theano shared variables
            time_step = theano.shared(self.consts['time_step'], 'time_step')
            SC = theano.shared(self.consts['SC'], 'SC')
            I1 = theano.shared(self.consts['I1'], 'I1')
            tau0 = theano.shared(self.consts['tau0'], 'tau0')

            output, updates = theano.scan(fn=step_sde, sequences=[x_eta, z_eta], outputs_info=[x_init, z_init], non_sequences=[time_step, SC, K, eta, I1, tau0, sig],  n_steps=self.consts['nt'])
            
            x_sym = output[0]
            z_sym = output[1]

            x = pm.Deterministic('x', x_sym)
            z = pm.Deterministic('z', z_sym)

            xhat = pm.Deterministic('xhat', amplitude * x_sym + offset)

            xs = pm.Normal('xs', mu=xhat, sd=eps, shape=(self.consts['nt'], self.consts['nn']), observed=self.obs['xs'])


####################################################################################################
####################################################################################################

