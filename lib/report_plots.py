#!/usr/bin/env python3
"""
@author: meysamhashemi INS Marseille

"""
        

import os
import sys
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import json
import pickle

import re
import glob
import itertools
from itertools import chain
from operator import itemgetter
from pandas.plotting import scatter_matrix

from report_metrics import LSE, Err, RMSE, LSE_obs, z_score, shrinkage


##########################################################################################################
##########################################################################################################
def plot_showpicks(ts, Obs, Obs_est, showpicks):
        Obs_lo, Obs_hi = np.percentile(Obs_est, [5, 95], axis=0)
        t_Obs = np.r_[:Obs_lo.shape[0]]
        for i in range(len(showpicks)):
              plt.subplot(len(showpicks), 1, i + 1)   
              plt.plot(ts, Obs[:, showpicks[i]],     'r',   alpha=0.8, linewidth=2.0, label=f'activity of node {showpicks[i]+1}')
              plt.plot(ts, Obs_lo[:, showpicks[i]],  'b',   alpha=0.8, linewidth=2.0, label='fitted percentiles')
              plt.plot(ts, Obs_hi[:, showpicks[i]],  'b',   alpha=0.8, linewidth=2.0)
              plt.ylabel('Amp'); 
              plt.ylim([-2.5, 2.5])
              if i==0:
                  plt.title('Fitted Source level data')
              if i==len(showpicks)-1:
                  plt.xlabel('Time []'); 
              plt.legend()
##########################################################################################################  
def plot_prior(nodes, Prior):
    plt.violinplot(Prior, widths=0.7, showmeans=True)
    plt.xticks(np.r_[1:len(nodes)+1], np.r_[1:len(nodes)+1], rotation=30)
    plt.ylabel("$Prior~(x_0)$")
    plt.tight_layout()   
##########################################################################################################      
def plot_posterior(nodes, eta_c, delta_eta, eta_hz, eta_ez, eta_pz, Hz_idx, Ez_idx, Pz_idx, eta_est): 
        parts= plt.violinplot(eta_est, widths=0.7, showmeans=True, showextrema=True);
        plt.plot(Ez_idx+1, eta_ez*np.ones((len(Ez_idx))) ,'o', color='k', alpha=0.9, markersize=3)
        plt.plot(Pz_idx+1, eta_pz*np.ones((len(Pz_idx))) ,'o', color='k', alpha=0.9, markersize=3)
        plt.plot(Hz_idx+1, eta_hz*np.ones((len(Hz_idx))),'o' , color='k', alpha=0.9, markersize=3)
        plt.axhline(y=eta_c, linewidth=.8, color = 'r', linestyle='--')
        plt.axhline(y=eta_c-delta_eta, linewidth=.8, color = 'y', linestyle='--')
        text(-2.85, -1.7, "$EZ$", size=12, color = 'r')
        text(-2.85, -2.6, "$PZ$", size=12, color = 'y')
        text(-2.85, -3.7, "$HZ$", size=12, color = 'g')
        #plt.xticks(np.r_[1:len(nodes)+1], np.r_[1:len(nodes)+1], rotation=90, fontsize=8)  
        plt.xticks(np.arange(1,len(nodes)+2, step=2),np.arange(1, len(nodes)+1, step=2), fontsize=10, rotation=0)
        plt.ylabel(' Posterior ' +r'${(\eta_i)}$', fontsize=14);  
        plt.xlabel('Brain nodes', fontsize=14); 
        plt.tight_layout()
        #plt.text(-9.2, -.9, "B" ,fontsize=24, fontweight='bold')

        for pc in parts['bodies'][0:nodes.shape[0]]:
            pc.set_facecolor('g')
            pc.set_edgecolor('g')
            pc.set_alpha(0.5)
        i = 0
        while i < len(Ez_idx):
            for pc in parts['bodies'][Ez_idx[i]:Ez_idx[i]+1]:
                pc.set_facecolor('r')
                pc.set_edgecolor('r')
                pc.set_alpha(0.8)
            i += 1

        j = 0
        while j < len(Pz_idx):
            for pc in parts['bodies'][Pz_idx[j]:Pz_idx[j]+1]:
                pc.set_facecolor('y')
                pc.set_edgecolor('y')
                pc.set_alpha(0.8)
            j += 1
        
##########################################################################################################
def plot_features(ts, Obs, Obs_est, showpicks):

       Obs_lo, Obs_hi = np.percentile(Obs_est, [5, 95], axis=0)
       t_Obs = np.r_[:Obs_lo.shape[0]]

       plt.fill_between(ts, Obs_lo[:, showpicks[0]], Obs_hi[:, showpicks[0]], color='green', alpha=.6, zorder=8) 
       plt.fill_between(ts, Obs_lo[:, showpicks[1]], Obs_hi[:, showpicks[1]], color='darkgoldenrod', alpha=.9, zorder=8) 
       plt.fill_between(ts, Obs_lo[:, showpicks[2]], Obs_hi[:, showpicks[2]], color='darkred', alpha=.5, zorder=8) 

       plt.plot(ts, Obs[:, showpicks[0]], color='g',alpha=.4, linestyle='--', marker= '.' , markersize=14, linewidth=1, label=f'node {showpicks[0]+1} (HZ)', zorder=4)
       plt.plot(ts, Obs[:, showpicks[1]], color='y',alpha=.7,  linestyle='--', marker= '.' , markersize=14, linewidth=1, label=f'node {showpicks[1]+1} (PZ)', zorder=4)
       plt.plot(ts, Obs[:, showpicks[2]], color='red', alpha=.8, linestyle='--', marker= '.' , markersize=14, linewidth=1, label=f'node {showpicks[2]+1} (EZ)', zorder=4)


       #for i, (pl_lo, pl_hi) in enumerate(zip(Obs_lo.T, Obs_hi.T)):
       #     if i in showpicks:
            # plt.fill_between(ts, pl_lo, pl_hi, alpha=0.95, linewidth=12.0)               
       # leg = []     
       # for i in  showpicks :
       #      leg.append(f'node {i+1}')
       # plt.legend(leg)

       plt.legend(fontsize=9)
       plt.ylabel('Obs. vs fitted '+r'$ x_{1,i}(t)$', fontsize=14);        
       plt.xlabel('Time (s)', fontsize=14); 
       plt.xticks(fontsize=14)
       plt.yticks(fontsize=14)
       #plt.text(-23.0, 1.5, "A" ,fontsize=24, fontweight='bold')
##########################################################################################################
def ppplot(eta_est,  eta_hz, eta_ez, eta_pz, showpicks):
       prior_sd =1.0
       prior_mu=array([-2.5])
       Prior = [np.random.normal(mu, prior_sd, size=10000) for mu in prior_mu ]
       prior_plot=sns.distplot(Prior[0],bins=20, hist=True, kde=True, color="blue", label=f' Prior')
       sns.distplot(eta_est[:, showpicks[0]] ,hist=True, kde=True, bins=20, color = 'g', hist_kws={'edgecolor':'g'}, kde_kws={'linewidth': 2},  label=f' Posterior node {showpicks[0]+1} (HZ)')
       sns.distplot(eta_est[:, showpicks[1]] ,hist=True, kde=True, bins=20, color = 'y', hist_kws={'edgecolor':'y'}, kde_kws={'linewidth': 2},  label=f' Posterior node {showpicks[1]+1} (PZ)')
       sns.distplot(eta_est[:, showpicks[2]] ,hist=True, kde=True, bins=20, color = 'r', hist_kws={'edgecolor':'r'}, kde_kws={'linewidth': 2},  label=f' Posterior node {showpicks[2]+1} (EZ)')
       plt.axvline(x=eta_hz, ymin=0.0, ymax=.65, linewidth=3, linestyle='--', color='gray')
       plt.axvline(x=eta_pz, ymin=0.0, ymax=.7, linewidth=3, linestyle='--', color='gray')
       plt.axvline(x=eta_ez, ymin=0.0, ymax=.5, linewidth=3, linestyle='--', color='gray')
       prior_plot.set(xlim=(-6, 1))
       prior_plot.set(ylim=(0, 2.5))
       plt.legend(fontsize=8)
       plt.xlabel(r'$ (x_{0,i})$', fontsize=10);
       plt.xticks(fontsize=10)
       plt.yticks(fontsize=10)
       #plt.text(-6.75, 2.5, "B" ,fontsize=16, fontweight='bold')
##########################################################################################################
def plot_obs_x(obs_x):
    ts = np.r_[:obs_x.shape[0]]
    plt.plot(ts, obs_x, alpha=0.5, marker= '.' , markersize=4)
    plt.xticks(np.arange(0,121,20), np.arange(0,121,20),fontsize=10)
    plt.yticks(np.arange(3,4,0.5), np.arange(3,4,0.5),fontsize=10)
    plt.ylabel(' Obs.' +r'$ x_{1,i}(t)$', fontsize=10);  
    plt.xlabel('Time (s)', fontsize=10); 
##########################################################################################################
def plot_obs_z(obs_z):
    ts = np.r_[:obs_z.shape[0]]
    plt.plot(ts, obs_z, alpha=0.5, marker= '.' , markersize=4)
    plt.xticks(np.arange(0,121,20), np.arange(0,121,20),fontsize=10)
    plt.yticks(np.arange(3,4,0.5), np.arange(3,4,0.5),fontsize=10)
    plt.ylabel(' Unobserved ' +r'$z_{i}(t)$', fontsize=10);  
    plt.xlabel('Time (s)', fontsize=10); 
##########################################################################################################
def plot_hiddenstates_x(fit_x):
    x_lw, x_hi = percentile(fit_x, [5, 95], axis=0)
    t = np.r_[:x_lw.shape[0]]
    for x_l, x_h in zip(x_lw.T, x_hi.T):
        plt.fill_between(t, x_l, x_h, alpha=0.8,  edgecolor='none',  label='x')
    plt.ylabel(' fitted ' +r'$ x_{1,i}(t)$', fontsize=10);  
    plt.xlabel('Time (s)', fontsize=10); 
##########################################################################################################
def plot_hiddenstates_z(fit_z):
    z_lw, z_hi = percentile(fit_z, [5, 95], axis=0)
    t = np.r_[:z_lo.shape[0]]
    for z_l, z_h in zip(z_lo.T, z_hi.T):
        plt.fill_between(t, z_l, z_h, alpha=0.8, edgecolor='none',  label='z')
    plt.ylabel(' Inferred ' +r'$ z_{i}(t)$', fontsize=10);  
    plt.xlabel('Time (s)', fontsize=10); 
##########################################################################################################
def plot_zscore_shrinkage(nodes, eta_true_mu, eta_est_mu, eta_est_std, prior_std):
    z_score_eta=z_score(eta_true_mu, eta_est_mu, eta_est_std)
    colors= np.random.rand(z_score_eta.shape[0])
    plt.scatter(shrinkage([prior_std]*nodes.shape[0], eta_est_std), z_score_eta ,s=120, c='blue')
    plt.xlabel("Posterior shrinkages", fontsize=14)
    plt.ylabel("Posterior z-scores", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.axis((0,1.1,0,10))
    #plt.text(-.4, 10.4, "C" ,fontsize=24, fontweight='bold')
##########################################################################################################
def pair_plots_params(csv, keys, skip=0):
    n = len(keys)
    if isinstance(csv, dict):
        csv = [csv]  # following assumes list of chains' results
    for i, key_i in enumerate(keys):
        for j, key_j in enumerate(keys):
            plt.subplot(n, n, i*n+j+1)
            for csvi in csv:
                if i==j:
                    plt.hist(csvi[key_i][skip:], 20, log=False)
                    plt.xticks(fontsize = 10)   
                    plt.yticks(fontsize = 10)  
                else:
                    plt.plot(csvi[key_j][skip:], csvi[key_i][skip:], '.')
                    plt.xticks(fontsize = 10)   
                    plt.yticks(fontsize = 10)  
            if i==0:
                plt.title(key_j, fontsize = 14)
            if j==0:
                plt.ylabel(key_i, fontsize = 14)

##########################################################################################################
def pair_plots(samples, params, figname='', sampler=None):
    import numpy as np
    import matplotlib.pyplot as plt
    div_iters = np.where(samples['divergent__'] == 1)[0] if sampler == 'HMC' else []
    plt.figure(figsize=(23, 13))
    nParams = len(params)
    for i in range(nParams):
        for j in range(nParams):
            plt.subplot(nParams, nParams, i * nParams + (j + 1))
            if (i == j):
                plt.hist(samples[params[i]].flatten(), bins=20, color='black')
            else:
                xvals = np.mean(
                    samples[params[j]], axis=1) if len(
                        samples[params[j]].shape) > 1 else samples[params[j]]
                yvals = np.mean(
                    samples[params[i]], axis=1) if len(
                        samples[params[i]].shape) > 1 else samples[params[i]]
                for k in range(xvals.shape[0]):
                    if (k in div_iters):
                        plt.plot(xvals[k], yvals[k], 'ro', alpha=0.8)
                    else:
                        plt.plot(xvals[k], yvals[k], 'ko', alpha=0.1)
            if (i == 0):
                plt.title(params[j], fontsize=13)
            if (j == 0):
                plt.ylabel(params[i], fontsize=13, rotation=90)
    plt.tight_layout()
    if (figname):
        plt.savefig(figname)                
##########################################################################################################
def func_2DepileptorEqs(y, eta, K, I, SC):  
    num_nodes= int(SC.shape[0])
    x = np.empty(num_nodes)
    z = np.empty(num_nodes)
    F = np.empty(2*num_nodes)
    for i in np.r_[0:num_nodes]:
        x[i] = y[2*i]
        z[i] = y[2*i+1]
    for i in np.r_[0:num_nodes]: 
        gx=0
        for j in np.r_[0:num_nodes]:
               gx=gx+SC[i,j]*(x[j]-x[i])  
        F[2*i] = -x[i]**3-2*x[i]**2-z[i]+I
        F[2*i+1] = 4*(x[i]-eta[i])-z[i]-K*gx   
    return F                
##########################################################################################################
def plot_nullclines(i, npz, nodes, SC, eta, K):
    tau0, Ks, I1 = npz['tau0'], npz['Ks'], npz['I1']
    X, Z = np.mgrid[-10.0:10.0:1000j, -10.0:10.0:1000j]
    for ii, idx in enumerate(nodes):
            gx=0
            for jj, idx in enumerate(nodes):
                  gx=gx+SC[ii,jj]*(X[:,jj]-X[:,ii])
            dX = (I1 + 1.0) - X**3.0 - 2.0*X**2.0 - Z
            dZ = (1.0/tau0) * ((4.0 * (X - eta[i])) - Z -K*gx )        
    CS = plt.contour(X, Z, dX, 0, colors='b', linewidths=(1.5,))
    #CS = plt.contour(X, Z, dZ, 0, colors='b', linewidths=(1.,)) 

    x = np.r_[-5:5:0.01]
    z = np.r_[0:10:0.01]
    xx, zz = np.meshgrid(x,z)
    dx = np.zeros([len(z),len(x)])
    dz = np.zeros([len(z),len(x)])

    for ii in range(len(z)):
       for jj in range(len(x)):
            dx[ii,jj] = 1 - x[jj]**3 - 2*x[jj]**2 - z[ii] + I1
            dz[ii,jj] = 1/tau0*(4*(x[jj] - eta[i]) - z[ii])
    plt.streamplot(x,z,dx,dz,density=2.0, color='deepskyblue', linewidth=.5 , zorder=1)
############################################################################################################
def plot_znullcline(i, nodes, roots):
    Nt=10000
    xz_range=np.linspace(-4,5.0, Nt)     
    xs = np.empty((len(nodes),Nt))
    ys = np.empty((len(nodes),Nt))    
    xs[:,]=xz_range
    for ii in (nodes): 
            ys[ii] =4.0*(xs[ii]-roots[2*i])+roots[2*i+1]
    plt.plot(xz_range,ys[i,:], 'blue', linewidth=1.5)
############################################################################################################
def plot_phasespace(csv, npz, nodes, showpicks, SC, eta, K, true_roots, true_roots_K0, estimated_roots, estimated_roots_K0):
    #fig=plt.figure(figsize=(12, 5))

    if 'Obs_seeg' in npz.files :
            obs_x, obs_z = npz['x_source'], npz['z_source'] 
    else:
            obs_x, obs_z = npz['Obs'], npz['Obs2'] 


    plt.subplot(2, len(showpicks), 1)   
    plt.plot(obs_x[:, showpicks[0]]-.28, obs_z[:, showpicks[0]]+1.44, 'g', alpha=0.4, linewidth=6., zorder=4, label=f'node {showpicks[0]+1} (HZ)')    
    plot_znullcline(showpicks[0], nodes, estimated_roots)    
    plot_nullclines(showpicks[0], npz, nodes, SC, eta, K)
    plt.scatter(estimated_roots[2*showpicks[0]],estimated_roots[2*showpicks[0]+1], s=18, facecolors='k', edgecolors='k', linewidth=1., zorder=3)
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.ylabel(r'$ z_{i}$', fontsize=18);   
    plt.title('Simulated HZ', fontsize=12)

    plt.subplot(2, len(showpicks), 2)   
    plt.plot(obs_x[:, showpicks[1]], obs_z[:, showpicks[1]], 'y', alpha=.7, linewidth=2., zorder=3, label=f'node {showpicks[1]+1} (PZ)')
    plot_znullcline(showpicks[1], nodes, true_roots)    
    plot_nullclines(showpicks[1], npz, nodes, SC, eta, K)
    plt.scatter(true_roots[2*showpicks[1]],true_roots[2*showpicks[1]+1] , s=28, facecolors='none', edgecolors='k', linewidth=1., zorder=4)   
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.title('Simulated PZ', fontsize=12)

    plt.subplot(2, len(showpicks),  3)   
    plt.plot(obs_x[:, showpicks[2]], obs_z[:, showpicks[2]], 'r', alpha=0.7, linewidth=2., zorder=3, label=f'node {showpicks[2]+1} (EZ)')   
    plot_znullcline(showpicks[2], nodes, true_roots)    
    plot_nullclines(showpicks[2], npz, nodes, SC, eta, K)
    plt.scatter(true_roots[2*showpicks[2]],true_roots[2*showpicks[2]+1] , s=28, facecolors='none', edgecolors='k', linewidth=1., zorder=4) 
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.title('Simulated EZ', fontsize=12)

    x, z = csv['x'], csv['z'] 
    Nsamples=4
    plt.subplot(2, len(showpicks), 4)   
    for j in range(Nsamples):
           plt.plot(x[-j, :, showpicks[0]], z[-j, :, showpicks[0]], 'darkgreen', alpha=1., linewidth=1., zorder=3)
    plt.plot(x[-j, :, showpicks[0]], z[-j, :, showpicks[0]], 'darkgreen', alpha=1., linewidth=2., zorder=3, label=f'node {showpicks[0]+1} (HZ)')       
    plot_znullcline(showpicks[0], nodes, estimated_roots)    
    plot_nullclines(showpicks[0], npz, nodes, SC, eta, K)
    plt.scatter(estimated_roots[2*showpicks[0]],estimated_roots[2*showpicks[0]+1] , s=24, facecolors='k', linewidth=1.,  edgecolors='k', zorder=4)
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.ylabel(r'$ z_{i}$', fontsize=18);   
    plt.xlabel(r'$ x_{1,i}$', fontsize=18);     
    plt.title('Inferred HZ', fontsize=12)

    plt.subplot(2, len(showpicks), 5)   
    for j in range(Nsamples):
           plt.plot(x[-j, :, showpicks[1]], z[-j, :, showpicks[1]], 'darkgoldenrod', alpha=1., linewidth=1.5, zorder=3) 
    plt.plot(x[-j, :, showpicks[1]], z[-j, :, showpicks[1]], 'darkgoldenrod', alpha=1., linewidth=1.5, zorder=3, label=f'node {showpicks[1]+1} (PZ)') 
    plot_znullcline(showpicks[1], nodes, estimated_roots)    
    plot_nullclines(showpicks[1], npz, nodes, SC, eta, K)
    plt.scatter(estimated_roots[2*showpicks[1]] ,estimated_roots[2*showpicks[1]+1] , s=28, facecolors='none',  linewidth=1., edgecolors='k', zorder=4)   
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.xlabel(r'$ x_{1,i}$', fontsize=18);   
    plt.title('Inferred PZ', fontsize=12)

    plt.subplot(2, len(showpicks),  6)   
    for j in range(Nsamples):
           plt.plot(x[-j, :, showpicks[2]], z[-j, :, showpicks[2]], 'darkred', alpha=1., linewidth=1.5, zorder=3) 
    plt.plot(x[-j, :, showpicks[2]], z[-j, :, showpicks[2]], 'darkred', alpha=1., linewidth=1.5, zorder=3, label=f'node {showpicks[2]+1} (EZ)')        
    plot_znullcline(showpicks[2], nodes, estimated_roots)    
    plot_nullclines(showpicks[2], npz, nodes, SC, eta, K)
    plt.scatter(estimated_roots[2*showpicks[2]],estimated_roots[2*showpicks[2]+1] , s=28, facecolors='none', linewidth=1., edgecolors='k', zorder=4) 
    plt.axis((-4,3,0,7))
    plt.legend(loc=1, fontsize=8)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.xlabel(r'$ x_{1,i}$', fontsize=18);   
    plt.title('Inferred EZ', fontsize=12)

    plt.tight_layout()

    # for i, label in enumerate(('A', 'B', 'C', 'D', 'E', 'F' )):
    #     ax = fig.add_subplot(2,3,i+1)
    #     ax.text(-0.05, 1.13, label, transform=ax.transAxes,
    #       fontsize=24, fontweight='bold', va='top', ha='right')



##########################################################################################################
# fig.tight_layout(pad=2.5, w_pad=2, h_pad=3.)
# for i, idx in enumerate(showpicks):
#         plt.subplot(len(showpicks), 3, len(showpicks)+i + 1)    
#         for j in range(Nsamples):
#             plt.plot(x[-j, :, idx], z[-j, :, idx], 'g', alpha=0.2, linewidth=1.)
#         plt.plot(x[-j, :, idx], z[-j, :, idx], 'g', alpha=0.2, linewidth=1., label=f'node {idx+1}')
#         plot_znullcline(idx, nodes, roots)    
#         plot_nullclines(idx, npz, nodes, SC, eta, K)
#         plot(roots_K0[2*idx], roots_K0[2*idx+1] , marker="o", markersize=4, color='k')
#         plt.plot(roots[2*idx], roots[2*idx+1] , marker="o", markersize=3, color='y')
#         plt.axis((-4,4,-2,8))
#         plt.grid(True)
#         if i>1: plt.xlabel('x(t)')
#         if i % 2 == 0: plt.ylabel('z(t)')
#         plt.legend()
#         plt.title(f'node {idx+1}')   
##########################################################################################################
##########################################################################################################        
#From scikit-learn
def plot_confusion_matrix(cm, target_names, title='Confusion matrix', cmap=None, normalize=False):
    """
    given a sklearn confusion matrix (cm), make a nice plot

    Arguments
    ---------
    cm:           confusion matrix from sklearn.metrics.confusion_matrix

    target_names: given classification classes such as [0, 1, 2]
                  the class names, for example: ['high', 'medium', 'low']

    title:        the text to display at the top of the matrix

    cmap:         the gradient of the values displayed from matplotlib.pyplot.cm
                  see http://matplotlib.org/examples/color/colormaps_reference.html
                  plt.get_cmap('jet') or plt.cm.Blues

    normalize:    If False, plot the raw numbers
                  If True, plot the proportions

    Usage
    -----
    plot_confusion_matrix(cm           = cm,                  # confusion matrix created by
                                                              # sklearn.metrics.confusion_matrix
                          normalize    = True,                # show proportions
                          target_names = y_labels_vals,       # list of names of the classes
                          title        = best_estimator_name) # title of graph

    Citiation
    ---------
    http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html

    """


    accuracy = np.trace(cm) / float(np.sum(cm))
    misclass = 1 - accuracy

    if cmap is None:
        cmap = plt.get_cmap('Blues')

    #plt.figure(figsize=(8, 6))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()

    if target_names is not None:
        tick_marks = np.arange(len(target_names))
        plt.xticks(tick_marks, target_names, rotation=45)
        plt.yticks(tick_marks, target_names)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]


    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        if normalize:
            plt.text(j, i, "{:0.4f}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
        else:
            plt.text(j, i, "{:,}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

    #plt.tight_layout()
    plt.title('Confusion matrix', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('True nodes', fontsize=12)
    plt.xlabel('Predicted nodes\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass), fontsize=12)
    #plt.text(-1.8, -.6, "D" ,fontsize=24, fontweight='bold')
##########################################################################################################        
##########################################################################################################        
