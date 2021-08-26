"""
The code written by Stan team (from Michael Betancourt)
Adapted for BVEP model by meysamhashemi INS Marseille  

"""

import os
import sys
import numpy as np
import pandas as pd
import json
import pickle
import subprocess
import matplotlib.pyplot as plt

from psis import psisloo
##########################################################################################################
##########################################################################################################
np.seterr(divide='ignore', invalid='ignore')
plt.rcParams.update({'figure.max_open_warning': 0})
##########################################################################################################
def check_div(divergent):
    n = sum(divergent)
    N = len(divergent)
    divergence=100 * n / N
    return divergence
##########################################################################################################
def check_treedepth(depths, max_depth):
    n = sum(1 for x in depths if x == max_depth)
    N = len(depths)
    treedepth=100 * n / N
    return treedepth
##########################################################################################################
def check_energy(energies):
    """Checks the energy Bayesian fraction of missing information (E-BFMI)"""
    numer = sum((energies[i] - energies[i - 1])**2 for i in range(1, len(energies))) / len(energies)
    denom = np.var(energies)
    energy=numer / denom
    return energy
##########################################################################################################
def check_n_eff(lp__, fit_summary):
    """Checks the effective sample size per iteration"""
    N_Eff = fit_summary['N_Eff']
    n_effs=pd.to_numeric(N_Eff, errors='coerce')
    names = fit_summary['name']
    n_iter = len(lp__)

    i=0
    j=0
    ratios=[]
    for n_eff, name in zip(n_effs, names):
        j=j+1 
        ratio = n_eff / n_iter
        ratios.append(ratio)        
        if (ratio <0.001):
            i=i+1
    ratios_sum=sum(np.asarray(ratios))        
    n_eff_low_percent=i/j     
    return n_eff_low_percent
##########################################################################################################
def check_rhat(fit_summary):
    """Checks the potential scale reduction factors"""
    from math import isnan
    from math import isinf

    Rhat_values = fit_summary['R_hat']
    rhats=pd.to_numeric(Rhat_values, errors='coerce')
    rhats = rhats[~np.isnan(rhats)]   
    names = fit_summary['name']
    prefixes = ["log_lik", "x_ppc", "xhat_q", "num_data", "num_params", "lp__", "divergent__",  "treedepth__", "energy__"]
    names= names.loc[~names.str.startswith(tuple(prefixes))] 
    i=0
    j=0
    for rhat, name in zip(rhats, names):
        j=j+1 
        if (rhat > 1.1):
            i=i+1
    Rhat_large_percent=i/j        
    return Rhat_large_percent
##########################################################################################################
def run_summary(filepath, inputfile):
    head, tail = os.path.split(os.path.split(inputfile)[0])
    fit_filefoder= inputfile.split("/")[-2]
    fit_filename = inputfile.split("/")[-1]
    print('fit_filefoder:', fit_filefoder )
    print('fit_filename:', fit_filename )
    mycsvfile=filepath+str(fit_filefoder)+'/'+str(fit_filename)
    print('mycsvfile:', mycsvfile)
    os.chdir(filepath+str(fit_filefoder)+'/')
    run_stan_summary=True
    for fname in os.listdir(filepath+str(fit_filefoder)+'/'):
        if  os.path.exists("summary"+"_"+str(fit_filename)):
            print('summary already exisited!')
            run_stan_summary=False
    if run_stan_summary==True:
        print('running stan summary!')
        Input = subprocess.getoutput("/home/meysam/cmdstan/bin/stansummary --csv_file=summary"+"_"+str(fit_filename)+" " +str(mycsvfile))
        print('finished stan summary!')
    return pd.read_csv(str("summary"+"_"+str(fit_filename)),  comment='#')
##########################################################################################################
##########################################################################################################
def maxlike(log_lik):
    MLE=np.max(log_lik)
    return MLE
##########################################################################################################
def waic_values(log_lik):
    lppd=np.sum(np.log(np.mean(np.exp(log_lik), axis=0)))
    p_waic = np.sum(np.var(log_lik, axis=0, dtype=np.float64))
    elpd_waic=lppd - p_waic
    waic = -2 * elpd_waic
    return p_waic, lppd, elpd_waic, waic
##########################################################################################################
def loo_values(log_lik):
    LOO_PSIS= psisloo(log_lik)
    loo, loos, ks= LOO_PSIS[0], LOO_PSIS[1], LOO_PSIS[2]
    ks_thr=1.0
    ks_large_percent=((ks > ks_thr).sum()/ks.size)*100
    return loo, ks, ks_large_percent
##########################################################################################################
def Nuts_plot(dict_samples_diagnostics, fit_summary):

    Rhat_values = fit_summary['R_hat']
    rhats=pd.to_numeric(Rhat_values, errors='coerce')
    rhats = rhats[~np.isnan(rhats)]   
    names = fit_summary['name']
    prefixes = ["log_lik", "log_lik_sum",  "Seeg_qqc", "Seeg_ppc", "Seeg_qqc_vect", "x_ppc", "xhat_q",  "num_data", "num_params", "lp__", "accept_stat__",  "stepsize__",  "treedepth__", "n_leapfrog__",  "divergent__",  "energy__"]
    names= names.loc[~names.str.startswith(tuple(prefixes))]  
    names.index = pd.to_numeric(names.index)

    i=0
    j=0
    for rhat, name in zip(rhats, names):
        j=j+1 
        if (rhat > 1.1):
            i=i+1

    Rhat_large_percent=i/j
    Rhat_large_sum=((rhats > 1.1).sum())
    Rhat_large_sum_percent=((rhats > 1.1).sum()/len(rhats.reindex(names.index)))*100
    
    grid = plt.GridSpec(4, 2)  
    i=0
    for key, val in dict_samples_diagnostics.items():
        plt.subplot(grid[int(i/2), int(i % 2)])  
        plt.plot(dict_samples_diagnostics[key], alpha=0.5)
        if key in ('stepsize__', ):
           plt.gca().set_yscale('log')
        plt.title(key)
        plt.grid(1)
        i += 1                     
    plt.subplot(grid[3:, 0])
    plt.plot(rhats.reindex(names.index))
    plt.hlines(y=1.1, xmin=0., xmax=rhats.reindex(names.index).shape[0], linewidth=1.5, color = 'red', linestyle='--', zorder=5)
    #plt.text(-1, 1.4,  r'$\hat R > 1.1:\ \ $'+str(Rhat_large_percent)+r'$/$' + str( len(rhats[names.index])) , size=8, color = 'red', zorder=4)
    plt.text(-1, 1.11,  r'$\hat R > 1.1:\ \ $'+"{:.2f}".format(Rhat_large_percent)+r'$\%$' , size=8, color = 'red', zorder=4)
    plt.grid()
    plt.xlabel('Parameters')
    plt.ylabel(r'$\hat R$')
    plt.subplot(grid[3:, 1])
    plt.hist(rhats.dropna().values, bins=100)
    #plt.hist(rhats.dropna().values, bins=100, range=(1.,1.5))
    plt.xlabel(r'$\hat R$ values')
    plt.ylabel('Count')   
##########################################################################################################
def Elbo_plot(Iter, Elbo):
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    plt.plot(Elbo)
    my_yticks = ax.get_yticks()
    plt.yticks([])
    plt.xlabel('100*iter')
    plt.ylabel('ELBO')
#########################################################################################################
def ks_plot(loglik, len_seizure, num_seizure):
    loo, ks, k_large_percent=loo_values(loglik)
    y=ks
    y_thr=1.0
    plot = plt.scatter(np.r_[:y.shape[0]] , y, c = y, cmap = 'seismic', lw=.01)
    yall=y.flatten()
    plotall = plt.scatter(np.r_[:yall.shape[0]] , yall, c = yall, cmap = 'seismic', lw=.01)
    plt.colorbar(plotall)       
    plt.hlines(y=y_thr, xmin=0., xmax=y.shape[0], linewidth=1.5, color = 'white', linestyle='--' , zorder=2)
    plt.vlines(x= np.linspace(0, len_seizure*num_seizure, num=num_seizure), ymin=0., ymax=np.max(yall), linewidth=.5, color = 'y', linestyle='--', zorder=3)
    plt.text(-1, 1.1,  r'$ks > 1.0:\ \ $'+str(k_large_percent)+r'$\%$' , size=20, color = 'red', zorder=5)
    plt.margins(0.11, 0.1)
    #norm_colorbar(y.flatten())
    plt.ylabel('ks')
    plt.xlabel('data points')
##########################################################################################################

