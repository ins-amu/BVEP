import os
import sys
import time
import errno
import time
import timeit
import numpy as np
from numba import jit
import warnings
warnings.simplefilter("ignore")


@jit(nopython=False) 
def VEP2Dmodel(params, constants, init_conditions, SC, dt, ts):
    
    nt=ts.shape[0]
    nn=SC.shape[0]
    
    #parameters
    eta=params[0:nn]
    K=params[-1]
    
    # fixed parameters
    tau0, I1, sigma=constants[0], constants[1], constants[2]


    # simulation from initial point
    x = np.zeros((nn, nt))  # fast voltage
    z = np.zeros((nn, nt))  # slow voltage

    # initial conditions
    x_init, z_init=init_conditions[0], init_conditions[1]
    for i in range(0, nn):
          x[i, 0] = x_init
          z[i, 0] = z_init
    

    # integrate ODE
    for t in range(0, nt-1):
        for i in range(0, nn):
            gx = 0;
            for j in range(0, nn):
                    gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
            dx = 1.0 -  np.power(x[i, t], 3) - 2.0*np.power(x[i, t], 2) - z[i, t] + I1;
            dz = (1./tau0)*(4.*(x[i, t] - eta[i]) - z[i, t] - K*gx);
            x[i, t+1] = x[i, t] + dt*dx + np.sqrt(dt) * sigma * np.random.randn() 
            z[i, t+1] = z[i, t] + dt*dz + np.sqrt(dt) * sigma * np.random.randn()  
  
    
    return np.concatenate((x.reshape(-1) , z.reshape(-1) )) 
    
cwd = os.getcwd() + '/benchmarks/python'
Res_dir='data_output_files'
weights = np.loadtxt(os.path.join(cwd+"/ExperimentalData/connectivity", "weights.txt"))


# normalize Connectome
weights = weights/np.max(weights)
num_regions = len(weights)

hz_val=-3.5 
pz_val=-2.4
ez_val=-1.6

ez_idx = np.array([6, 34],  dtype=np.int32)
pz_idx = np.array([5, 11, 27],  dtype=np.int32)

weights[np.ix_(np.array([5,11]), ez_idx)] = 4.0
weights[np.ix_(np.array([27]), ez_idx)] = 1.0

SC=weights
SC.shape

T = 14.0
dt=0.1
ts = np.arange(0, T + dt, dt)

nt=ts.shape[0]
nn=SC.shape[0]

tau0=10.
I1=3.1    
sigma=0. # we assume ODE so no noise in the system
constants = np.array([tau0, I1, sigma])

x_init=-2.5
z_init=3.5
init_conditions = np.array([x_init, z_init])

eta_true = np.ones(nn)*hz_val
eta_true[ez_idx] = ez_val
eta_true[pz_idx] = pz_val


K_true=1. # global copuling parameters

params_true = np.append(eta_true, K_true)

Sim_true = VEP2Dmodel(params_true, constants, init_conditions, SC, dt, ts)

X_true=Sim_true[0:nn*nt].reshape(nn, nt)
Z_true=Sim_true[nn*nt:2*nn*nt].reshape(nn, nt)

print("PYTHON P3");


def cost_function(params):

    Sim =VEP2Dmodel(params, constants, init_conditions, SC, dt, ts)
    X_model=Sim[0:nn*nt].reshape(nn, nt)

    RMSE=np.sqrt(np.sum((X_model-X_true)**2))

    return RMSE
    
    
#from random import  uniform

#eta_r = [2 for i in range(84)]
#K_r = [1 for i in range(1)]

#params_guess=np.append(eta_r, K_r)
#print(cost_function(params_guess))
#from random import  uniform
#eta_r = [0.3 for i in range(84)]
#K_r = [1 for i in range(1)]
#params_guess=np.append(eta_r, K_r)
#print(Cost_function(params_guess))
