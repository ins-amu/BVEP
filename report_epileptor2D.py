#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 21:00:17 2017

@author: meysamhashemi
"""

import os
import sys
import numpy as np
import scipy as scp
import scipy as scp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import fsolve
from scipy.optimize import root
import warnings
warnings.filterwarnings("ignore")
##########################################################################################################
from parse_csv import parse_csv
##########################################################################################################
print ("Start of report!")
##########################################################################################################
cwd = os.getcwd()
print('Report directory:', cwd )
Mainpath=cwd+'/'
##########################################################################################################
data_input ='data_input_epileptor'
##########################################################################################################
data_input_npz = np.load(Mainpath+'/data_input_files/'+str(data_input)+'.R.npz')
##########################################################################################################
##########################################################################################################
opt_nutsplot=1
opt_loglik_plot=0
##########################################################################################################
##########################################################################################################
xs = data_input_npz['xs']
zs = data_input_npz['zs']
nt = data_input_npz['nt']
dt = data_input_npz['dt']
eta_true= data_input_npz['eta_true']
tau0= data_input_npz['tau0']
I1= data_input_npz['I1']
##########################################################################################################
print ("-"*60)
if '_ode' in sys.argv[1]:
        dynamic_type='ode'
        print('dynamic_type: ', dynamic_type)
elif '_deter' in sys.argv[1]:
        dynamic_type='deterministic'
        print('dynamic_type: ', dynamic_type)
elif '_noncen' in sys.argv[1]:
        dynamic_type='stochastic_noncen'
        print('dynamic_type: ', dynamic_type)  
elif '_cen' in sys.argv[1]:
        dynamic_type='stochastic_cen'
        print('dynamic_type: ', dynamic_type)      
else:
        dynamic_type='unknown'
        print('dynamic_type: ', dynamic_type)
#####################################################
#####################################################        
print ("-"*60)
print ('Report starts!')
#####################################################
#####################################################
print('csv file directory:', sys.argv[1])
fit = parse_csv(sys.argv[1])
#fit = parse_csv('output.csv')
#####################################################
#####################################################
my_file, file_extension = os.path.splitext(sys.argv[1])
script_dir = os.path.dirname(__file__)
print('script_dir:', script_dir)
repor_filename=sys.argv[1].rsplit('.',1)[0].split("/")[-1]
print('repor_filename:', repor_filename)

report_dir = os.path.join(cwd, 'report_'+repor_filename+'/')
if not os.path.isdir(report_dir):
    os.makedirs(report_dir)
#####################################################
#####################################################
numsamples=1000
Burnin=0
nsample=numsamples-Burnin
#####################################################
print ("-"*60)
print ("parameters and their shapes: ")
for key, val in fit.items():
    print(key, val.shape)
print ("-"*60)
##########################################################################################################
amplitude, offset = [fit[_].reshape((-1, 1, 1)) for _ in 'amplitude offset'.split()]
##########################################################################################################
##########################################################################################################
eta_est=fit['eta'][-nsample:]
amplitude_est=fit['amplitude'][-nsample:]
offset_est=fit['offset'][-nsample:]
eps_est=fit['eps'][-nsample:]
##########################################################################################################
if dynamic_type=='ode':
     sig_est=0.0
elif dynamic_type=='deterministic':
     sig_est=0.0
else:
     sig_est=fit['sig'][-nsample:]
##########################################################################################################
eta_mean=np.mean(eta_est)
amplitude_mean= np.mean(amplitude_est)
offset_mean= np.mean(offset_est)
eps_mean=np.mean(eps_est)
sig_mean= np.mean(sig_est)    
##########################################################################################################
##########################################################################################################
if dynamic_type=='ode':
    fit_x=fit['y_hat'][:,0,:]
    fit_z=fit['y_hat'][:,1,:]
else:
    fit_x=fit['x']
    fit_z=fit['z']
##########################################################################################################
##########################################################################################################
plt.figure(figsize=(12, 6))
xlo, xhi = np.percentile(fit_x, [5, 95], axis=0)
t = np.r_[:xlo.shape[0]]

plt.subplot(221)
plt.plot(xs, '--b.',linewidth=2, alpha=0.4)
plt.xlabel('Time (ms)'); plt.ylabel('observed x(t)');
zlo, zhi = np.percentile(fit_z, [5, 95], axis=0)
plt.subplot(222)
plt.fill_between(t, xlo.T, xhi.T, alpha=1., facecolor='b', edgecolor='b')
plt.xlabel('Time (ms)'); plt.ylabel('fitted x(t)');
plt.subplot(223)
plt.plot(zs, '--r.',linewidth=2, alpha=0.4)
plt.xlabel('Time (ms)'); plt.ylabel('non-observed z(t)');
plt.subplot(224)
plt.fill_between(t, zlo.T, zhi.T, alpha=1., facecolor='r', edgecolor='r')
plt.xlabel('Time (ms)'); plt.ylabel('inferred z(t)');
plt.savefig(os.path.join(report_dir, 'figure-hiddenstates.png'),dpi=300) 
##########################################################################################################
if dynamic_type=='ode':
    xs_est = fit['xhat_q'][:,0,:]
    xs_ppc=fit['x_ppc'][:,0,:]
else:
    xs_est = fit['xhat_q']
    xs_ppc=fit['x_ppc']
##########################################################################################################
##########################################################################################################
xest_lo, xest_hi = np.percentile(xs_est, [5, 95], axis=0)
xppc_lo, xppc_hi = np.percentile(xs_ppc, [5, 95], axis=0)
t = np.r_[:xest_lo.shape[0]]
##########################################################################################################
plt.figure(figsize=(12, 4))
plt.plot(xs, '--.',linewidth=2, color='cyan', alpha=0.4, label='$x_{obs}$', zorder=6)
plt.fill_between(t, xppc_lo, xppc_hi, linewidth=4, alpha=0.4, facecolor='r', edgecolor='r', label='$x_{ppc}$', zorder=3)
plt.fill_between(t, xest_lo, xest_hi, linewidth=3, alpha=0.4, facecolor='b', edgecolor='b',  label='$x_{fit}$', zorder=4)
plt.xlabel('Time (ms)'); plt.ylabel('x(t)'); 
plt.legend()
plt.savefig(os.path.join(report_dir, 'figure-fit.png'),dpi=300) 
##########################################################################################################
plt.figure(figsize=(10.2, 4))
for i, variables in enumerate([eta_est]):
    plt.subplot(1, 2, 1)
    plt.plot(variables)
    plt.axhline(y=variables.mean(), linewidth=2, color='cyan')
    plt.axhline(y=eta_true, linestyle='dashed', linewidth=2, color='r')
    plt.xlabel('iterations'); plt.ylabel("$\eta$") 
    
plt.subplot(1, 2, 2)
plt.violinplot(eta_est, widths=0.7, showmeans=True, showextrema=True);
plt.axhline(y=eta_true, linestyle='dashed', linewidth=2, color='r')
#axhline(y=eta_est.mean(), linewidth=2, color='g')
plt.ylabel("$\eta$")       
plt.tight_layout(pad=0.4, w_pad=0.6, h_pad=1.0)
plt.savefig(os.path.join(report_dir, 'figure-etafit.png'),dpi=300) 
##########################################################################################################
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]
##########################################################################################################
def LSE(x1, x2):
    return np.sum((x1 - x2)**2)
##########################################################################################################
def RMSE(x1, x2):
     return np.sqrt(((x1 - x2) ** 2))
##########################################################################################################
def nullclines(X,Y, eta, tau0):   
    dX =  1.0 - X**3 - 2.0*X**2 - Z + I1 
    dY =  (1.0/tau0)*(4*(X-eta) - Z)
    plt.contour(X, Y, dX, 0, colors='r')
    plt.contour(X, Y, dY, 0, colors='b') 
##########################################################################################################
def root_fuc(roots, eta, tau0):   
    x, z = np.empty(1), np.empty(1)
    F = np.empty(2)
    x = roots[0]
    z = roots[1]
    F[0] =1.0 - x**3 - 2.0*x**2 - z + I1 
    F[1] =(1.0/tau0)*(4*(x-eta) - z)
    return F
##########################################################################################################
rGuess=np.array([[-1.,3.]])
true_roots = fsolve(root_fuc,rGuess, args=(eta_true, tau0))
estimated_roots = fsolve(root_fuc,rGuess, args=(eta_mean, tau0))
##########################################################################################################
def true_phase_space(Obs, X, Y, eta, tau0, roots):
     xs, zs = Obs['xs'], Obs['zs']
     plt.plot(xs, zs, '--.g')
     nullclines(X,Z, eta, tau0)
     plt.plot(roots[0],roots[1] , marker="o", markersize=8, color='k')
     plt.xlabel('x'); plt.ylabel('z'); 
     plt.title('True phase-plane')
##########################################################################################################
if dynamic_type=='ode':
    def estimated_phase_space(csv, X, Y, eta, tau0, roots):
         x_lo, x_hi = np.percentile(csv['y_hat'][:,0,:], [5, 95], axis=0)
         z_lo, z_hi = np.percentile(csv['y_hat'][:,1,:], [5, 95], axis=0)
         plt.plot(x_lo, z_lo, '-y')
         plt.plot(x_hi, z_hi, '-y')
         x, z = np.mean(csv['y_hat'][:,0,:], axis=0), mean(csv['y_hat'][:,1,:], axis=0)
         plt.plot(x, z, '--.m')
         nullclines(X,Z, eta, tau0)
         plt.plot(roots[0],roots[1] , marker="o", markersize=8, color='k')
         plt.xlabel('x'); plt.ylabel('z'); 
         plt.title('Estimated phase-plane')
else:
    def estimated_phase_space(csv, X, Y, eta, tau0, roots):
     x_lo, x_hi = np.percentile(csv['x'], [5, 95], axis=0)
     z_lo, z_hi = np.percentile(csv['z'], [5, 95], axis=0)
     plt.plot(x_lo, z_lo, '-y')
     plt.plot(x_hi, z_hi, '-y')
     x, z = np.mean(csv['x'], axis=0), np.mean(csv['z'], axis=0)
     plt.plot(x, z, '--.m')
     nullclines(X,Z, eta, tau0)
     plt.plot(roots[0],roots[1] , marker="o", markersize=8, color='k')
     plt.xlabel('x'); plt.ylabel('z'); 
     plt.title('Estimated phase-plane')
##########################################################################################################
X, Z = np.mgrid[-5.0:5.0:100j, 0.0:10.0:100j]

plt.figure(figsize=(12.2, 5))
plt.subplot(121)
true_phase_space(data_input_npz, X, Z, eta_true, tau0, true_roots)
plt.xlim((-3,2)), plt.ylim((1,6))
plt.subplot(122)
estimated_phase_space(fit, X, Z, eta_est.mean(), tau0, estimated_roots)
plt.xlim((-3,2)), plt.ylim((1,6))
plt.savefig(os.path.join(report_dir, 'figure-phasespace.png'),dpi=300) 
##########################################################################################################
LSE_val=np.average([LSE(xs, xest_lo), LSE(xs, xest_hi)])
RMS_val=RMSE(eta_true, eta_mean)
##########################################################################################################
##########################################################################################################
def trace_nuts(csv, skip=0):
        i=1
        for key, val in csv.items():
             if key[-2:] == '__':
                #print(key, val.shape)
                plt.subplot(4, 2, i)
                plt.plot(csv[key][skip:], alpha=0.5)
                if key in ('stepsize__', ):
                    plt.gca().set_yscale('log')
                plt.title(key)
                plt.grid(1)
             i += 1
##########################################################################################################
if opt_nutsplot==1:
    plt.figure(figsize=(10.3, 6))  
    trace_nuts(fit)
    plt.tight_layout()
    plt.savefig(os.path.join(report_dir, 'figure-Nuts.png'),dpi=300) 
##########################################################################################################
if dynamic_type=='stochastic_noncen':
    opt_eta_plot=1
else:
    opt_eta_plot=0
##########################################################################################################
if opt_eta_plot==1:
    plt.figure(figsize=(8, 6))  
    plt.subplot(211)
    plt.plot(np.mean(fit['x_eta'],axis=0), label='$x_eta(t)$')
    plt.xlabel('t'); plt.ylabel('x_eta'); 
    plt.subplot(212)
    plt.plot(np.mean(fit['z_eta'],axis=0), label='$z_eta(t)$')
    plt.xlabel('t'); plt.ylabel('z_eta'); 
    plt.savefig(os.path.join(report_dir, 'figure-estimatednoise.png'),dpi=300)  
##########################################################################################################
if opt_loglik_plot==1:
    plt.figure(figsize=[12,4])
    for i in range(0,200):
        sns.kdeplot((-1*fit['log_lik'][i,:]))
    #plt.axvline(x=np.max(fit['log_lik']), ymin=0.0, ymax=0.7, linewidth=1.5, linestyle='-', color='r')
    plo=sns.kdeplot(-1*np.mean(fit['log_lik'], axis=0), color="k")
    plt.ylabel('log_lik'); 
    # plo.set(ylim=(0, 4))
    # plo.set(xlim=(-2, 2))
print('max_log_lik:', np.max(fit['log_lik']))
##########################################################################################################
with open(report_dir+'/'+ 'report_fit.txt', 'w') as fd:
    fd.write("eta_mean, eta_std, amplitude, offset, sig, eps, LSE, RMS\n")
    fmt = "%f,%f,%f,%f,%f,%f,%f,%s\n"
    args =  np.mean(eta_est), np.std(eta_est), np.mean(amplitude_est), np.mean(offset_est), np.mean(sig_est), np.mean(eps_est),  LSE_val, RMS_val
    fd.write(fmt % args)
##########################################################################################################
print ("-"*60)
print ('estimated eta is:', np.mean(eta_est))
print ('estimated amplitude is:', np.mean(amplitude_est))
print ('estimated offset is:', np.mean(offset_est))
print ('estimated eps is:', np.mean(eps_est))
print ('estimated sig is:', np.mean(sig_est))
print ("-"*60)
print ('LSE=',LSE_val)                               
print ('RMS=',RMS_val)
print ("-"*60)
print ("-"*60)
##########################################################################################################
plt.show()
##########################################################################################################
print ("-"*60)
print ("End of report!")
print ("-"*60)
print ("-"*60)
print ("-"*60)
##########################################################################################################