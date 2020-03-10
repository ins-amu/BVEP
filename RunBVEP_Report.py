#!/usr/bin/env python3
"""
@author: meysamhashemi INS Marseille

"""

import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import json
import pickle

from scipy.optimize import fsolve
from scipy.optimize import root

import scikitplot as skplt
from sklearn.metrics import confusion_matrix
import sklearn.preprocessing
from pandas.plotting import scatter_matrix


import re
import glob
from itertools import chain
from operator import itemgetter



from read_csvs import read_samples
from read_csvs import parse_csv
##########################################################################################################
np.seterr(divide='ignore', invalid='ignore')
plt.rcParams.update({'figure.max_open_warning': 0})
##########################################################################################################
import time
starttime = time.time()
##########################################################################################################
##########################################################################################################
print ("-"*60)
print ("Report starts!")
print ("-"*60)
##########################################################################################################
cwd = os.getcwd()
print("The results directory is:", cwd)
##########################################################################################################
Mainpath =cwd+"/"
##########################################################################################################
opt_print_node_abbrev=0
##########################################################################################################
opt_plot_figure_fit=1
opt_plot_figure_diagnostics=1
opt_pairplot_fit=1
##########################################################################################################
#####################################################
#####################################################
Reportmetrics_dir = os.path.join(cwd, "Report_Figs/")

if not os.path.isdir(Reportmetrics_dir):
    os.makedirs(Reportmetrics_dir) 
#####################################################  
if opt_print_node_abbrev==1:
    num_nodes=84
    nodes = np.r_[0:num_nodes]   

    reg_names = []
    with open(Mainpath+'ExperimentalData/centers.txt', 'r') as fd:
        for i, line in enumerate(fd.readlines()):
            reg_names.append(( line.strip().split()[0]))

    print ("The brain nodes are: ") 
    node_name=[]
    node_region=[]
    node_abbrev = []
    for i in nodes:
        rn = reg_names[i]
        node_region.append(rn)
        ra = ''.join([p[0] for p in rn.split('-')[1:] if p]).upper()
        node_abbrev.append(ra)
        print (ra, ' '.join(rn.split('-')))
        node_name.append(ra)
#####################################################
def f_int(x):
    return np.int(x)

f_vector_int = np.vectorize(f_int)
#####################################################
numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
##########################################################################################################
##########################################################################################################
all_folders = sorted(os.listdir(cwd), key=numericalSort)
##########################################################################################################
all_csvfile_list=[]
all_csvfile_convergence_list=[]
all_csvfile_converged_list=[]


all_LSE_list=[]
all_err_eta_list=[]
all_z_score_eta_list=[]
all_shrinkage_eta_list=[]
all_score_std_list=[]
all_score_extrm_list=[]
all_score_ez_list=[]
all_score_pz_list=[]
all_accuracy_list=[]


all_divergence_list=[]
all_treedepth_list=[]
all_energy_list=[]
all_n_eff_low_percent_list=[]
all_Rhat_large_percent_list=[]
all_Rhat_large_percent_list=[]
all_Elbo_list=[]

all_maxlike_list=[]
all_p_waic_list=[]
all_elpd_waic_list=[]
all_waic_list=[]
all_loo_list=[]
all_ks_large_percent_list=[]
##########################################################################################################
I1=3.1
eta_c=-2.05        
delta_eta=1.0
eta_hi=1.0
eta_lo=-5.0                  
prior_std =1.0
sd_thr=0.8
##########################################################################################################

counter_converged_csv = 0 
counter_Notconverged_csv = 0 


for folder in all_folders:
            if folder.startswith("data_output"):
                
                            cwd_=os.path.join(cwd, folder)
                            os.chdir(cwd_)


                            print ("-"*60)
                            print ("csv files found at folder :\t", folder)
                            patient_name, datafile_name, alg_name, measure_level =folder.split("_")[-1], folder.split("_")[-2], folder.split("_")[2], folder.split("_")[3]

                            print ("algorithm:", alg_name)
                            print ("measure level:", measure_level)
                            print ("patient:", patient_name)
                            print ("data file name:", datafile_name)

                            if measure_level=='SL':
                                level_name='Sourcelevel'
                            else:
                            	level_name='Seeglevel'

                            print ("-"*60)
                            print ("loading ground-truth...")

                            eta_true_inx, eta_true_mu, eta_true_sd = np.loadtxt(Mainpath+"ExperimentalData"+"/"+"eta_"+str(measure_level)+'_'+str(patient_name)+".txt", delimiter=',', usecols=(0,1,2), skiprows=1, ).T 

                            with open(Mainpath+"ExperimentalData"+"/"+"true_labels_nodes_"+str(measure_level)+'_'+str(patient_name)+".txt", 'r') as filehandle:  
                                    true_labels_nodes = json.load(filehandle)

                            print ("loading data-input...")
                            npz = np.load(Mainpath+"data_input_files/"+str(patient_name)+"/"+str(level_name)+"/"+str(datafile_name)+"_"+str(patient_name)+".npz")

                            print ("-"*60)
                           
                            data_inut_list_items=[]
                            for key, val in npz.items():
                                data_inut_list_items.append(key)

                  
                            Obs = npz['Obs']
                            Obs2 = npz['Obs2']
                            SC = npz['SC']
                            tau0 = npz['tau0']
                            Ks = npz['Ks']
                            ts = npz['ts']
                            num_nodes=npz['nn']

                            if 'nodes' in  data_inut_list_items:
                                    nodes=npz['nodes']
                            else:
                                    nodes = np.r_[0:num_nodes]


                            if patient_name=='patient1':
                                eta_ez= -1.5
                                eta_pz= -2.6
                                eta_hz= -3.4

                            if patient_name=='patient2':
                                eta_ez= -1.6
                                eta_pz= -2.4
                                eta_hz= -3.65

                            print('eta_ez:', eta_ez)
                            print('eta_pz:', eta_pz)
                            print('eta_hz:', eta_hz)

                            from report_nodetype import node_EzPzHz_idex
                            Ez_indices, Pz_indices, Hz_indices, Ez_arr, Pz_arr, Hz_arr= node_EzPzHz_idex(true_labels_nodes)

                            print ("-"*60)
                            file_csvs = sorted(glob.glob(cwd_ + "/*.csv"), key=numericalSort)
                            file_outs = sorted(glob.glob(cwd_ + "/*.out"), key=numericalSort)
                            print ("-"*60)



                            for filecsv, fileout in zip(file_csvs, file_outs):
                                                    
                                                    print ('For', os.path.basename(os.path.normpath(filecsv)), ':')
                                                    csvfile_name=str(os.path.basename(os.path.normpath(filecsv)).rsplit('.',1)[0])
                                                    all_csvfile_list.append(csvfile_name)


                                                    from report_algconfigs import convergence_opt, convergence_advi, convergence_hmc
                                                    if alg_name=='opt':
                                                            print ("...checking convergence of opt...")
                                                            if convergence_opt(fileout) == 1:
                                                                print("...opt converged.")
                                                            else:
                                                                print("...opt did not converged!")
                                                    elif alg_name=='advi':
                                                            print ("...checking convergence of advi...")
                                                            if convergence_advi(fileout) == 1:
                                                                print("...advi converged.")
                                                            else:
                                                                print("...advi did not converged!")
                                                    elif alg_name=="hmc":
                                                            print ("...checking convergence of hmc...")
                                                            if convergence_hmc(fileout) == 1:
                                                                print("...hmc converged.")
                                                            else:
                                                                print("...hmc did not converged!")            
                                                    else:
                                                                print("unknown input")            





                                                    if  convergence_advi(fileout)==1 or convergence_hmc(fileout)==1: 

                                                            all_csvfile_convergence_list.append('converged')
                                                            all_csvfile_converged_list.append(csvfile_name)

                                                            print ("-"*60)
                                                            print ("parsing values from csv file...")

                                                            if alg_name=="hmc":
                                                                    from report_algconfigs import hmc_config
                                                                    num_samples, num_warmup, save_warmup,  adapt_delta, max_depth, total_time= hmc_config(fileout)
                                                                    num_samples=int(num_samples)
                                                                    dict_fit=read_samples(filecsv, nwarmup=0, nsampling=num_samples,    variables_of_interest=['x', 'z', 'amplitude', 'offset'])
                                                                    dict_eta=read_samples(filecsv,  nwarmup=0, nsampling=num_samples,    variables_of_interest=['eta']) 
                                                                    dict_K=read_samples(filecsv,  nwarmup=0, nsampling=num_samples,    variables_of_interest=['K']) 
                                                                    dict_samples_loglik=read_samples(filecsv, nwarmup=0, nsampling=num_samples, variables_of_interest=['log_lik'])

                                                            if alg_name=="advi":
                                                                    from report_algconfigs import advi_config, advi_elbo
                                                                    grad_samples, output_samples,  output_samples, tol_rel_obj = advi_config(fileout)
                                                                    Iter, Elbo= advi_elbo(fileout)
                                                                    num_iter=int(np.max(Iter))
                                                                    dict_fit=read_samples(filecsv, nsampling=num_iter,  variables_of_interest=['x', 'z', 'amplitude', 'offset'])
                                                                    dict_eta=read_samples(filecsv,  nsampling=num_iter,  variables_of_interest=['eta']) 
                                                                    dict_K=read_samples(filecsv,  nsampling=num_iter,  variables_of_interest=['K']) 
                                                                    dict_samples_loglik=read_samples(filecsv,  nsampling=num_iter,  variables_of_interest=['log_lik'])



                                                            dict_fit_x, dict_fit_z, dict_fit_amp, dict_fit_offset=dict_fit['x'], dict_fit['z'], dict_fit['amplitude'], dict_fit['offset']
                                                            x, z, amplitude, offset = np.transpose(dict_fit_x, (0, 2, 1)),  np.transpose(dict_fit_z, (0, 2, 1)), dict_fit_amp.reshape((-1, 1, 1)),  dict_fit_offset.reshape((-1, 1, 1))
                                                            
                                                            len_seizure, num_seizure =x.shape[1], x.shape[2]

                                                            K_est=np.mean(np.array(list(chain.from_iterable(dict_K.values()))), axis=0)


                                                            if measure_level=='SL':                                                           
                                                                    Obs_est =  amplitude*x + offset
                                                            else:                            
                                                                    gain = npz['gain']
                                                                    Obs_est =  amplitude*(x.dot(gain.T)) + offset

                                                            print ("-"*60)
                                                            print ("computing error metrics...")

                                                            from report_metrics import LSE, Err, RMSE, LSE_obs, z_score, shrinkage

                                                            Obs_low, Obs_hi = np.percentile(Obs_est, [5, 95], axis=0)
                                                            LSE_val=LSE_obs(Obs, Obs_low, Obs_hi)

                                                            chain_eta = chain.from_iterable(dict_eta.values())
                                                            eta_est=np.array(list(chain_eta))
                                                            eta_est_mu=np.mean(eta_est, axis=0)
                                                            eta_est_std=np.std(eta_est, axis=0)
                               
                                                            err_eta=Err(eta_true_mu, eta_est_mu)
                                                            z_score_eta=np.sum(z_score(eta_true_mu, eta_est_mu, eta_est_std))
                                                            shrinkage_eta=np.mean(shrinkage([prior_std]*nodes.shape[0], eta_est_std))



                                                            print ("-"*60)
                                                            print ("computing accuracy...")

                                                            from report_nodetype import node_score_std_extrm, node_score_ezpz
                                                            score_std, score_extrm= node_score_std_extrm(nodes, eta_est, eta_true_mu, sd_thr)
                                                            score_ez, score_pz, ce, cp, nodes_est, est_labels_nodes= node_score_ezpz(nodes, eta_est, true_labels_nodes, eta_c, delta_eta, eta_hi)
                                                            print ('Prediction about Ez and Pz nodes is=', int(score_ez+score_pz),' out of', len(Ez_indices)+len(Pz_indices), ' EzPznodes')                                              

                                                            y_true=true_labels_nodes
                                                            y_pred=est_labels_nodes                                           
                                                            confusion_matrix_values= confusion_matrix(y_true, y_pred, labels=["EZ", "PZ", "HZ"])
                                                            accuracy=100*np.trace(confusion_matrix_values)/np.sum(confusion_matrix_values)
                                                            print ('Prediction about all the node types is=', int(np.trace(confusion_matrix_values)),' out of', len(nodes), ' nodes')


                                                            all_LSE_list.append(LSE_val)
                                                            all_err_eta_list.append(err_eta)
                                                            all_z_score_eta_list.append(z_score_eta)
                                                            all_shrinkage_eta_list.append(shrinkage_eta)
                                                            all_score_std_list.append((score_std))
                                                            all_score_extrm_list.append((score_extrm))
                                                            all_score_ez_list.append((score_ez))
                                                            all_score_pz_list.append((score_pz))
                                                            all_accuracy_list.append(accuracy)


                                                            print ("-"*60)
                                                            print ("computing waic and loo...")


                                                            chain_loglik= chain.from_iterable(dict_samples_loglik.values())
                                                            loglik=np.array(list(chain_loglik))

                                                            from  report_diagnostics_stan import  maxlike, waic_values, loo_values
                                                            maxlike=maxlike(loglik)
                                                            p_waic, lppd, elpd_waic, waic=waic_values(-1*loglik)
                                                            loo, ks, ks_large_percent=loo_values(-1*loglik)


                                                            all_maxlike_list.append(maxlike)
                                                            all_p_waic_list.append(p_waic)
                                                            all_elpd_waic_list.append(elpd_waic)
                                                            all_waic_list.append(waic)
                                                            all_loo_list.append(loo)
                                                            all_ks_large_percent_list.append(ks_large_percent)


                                                            print ("-"*60)
                                                            print ("computing diagnostics...")

                                                            if alg_name=="hmc":
                                                                    from  report_diagnostics_stan import  run_summary, check_rhat, check_div, check_treedepth, check_energy, check_n_eff
                                                                    fit_summary=run_summary(Mainpath, filecsv)

                                                                    dict_samples_diagnostics=read_samples(filecsv, nwarmup=0, nsampling=num_samples,  variables_of_interest=['lp__', 'energy__', 'accept_stat__', 'divergent__', 'treedepth__', 'n_leapfrog__'  ])
                                                                    dict_samples_lp__, dict_samples_divergent__, dict_samples_treedepth__, dict_samples_energy__,  = dict_samples_diagnostics['lp__'], dict_samples_diagnostics['divergent__'], dict_samples_diagnostics['treedepth__'], dict_samples_diagnostics['energy__']

                                                                    divergence__=check_div(dict_samples_divergent__)
                                                                    treedepth__=check_treedepth(dict_samples_treedepth__, int(max_depth))
                                                                    energy__=check_energy(dict_samples_energy__)
                                                                    n_eff_low=check_n_eff(dict_samples_lp__, fit_summary)
                                                                    rhat_large=check_rhat(fit_summary)
                                                                    elbo_val=np.nan 
                                                                    
                                                                    if opt_plot_figure_diagnostics==1:
                                                                            from  report_diagnostics_stan import  Nuts_plot
                                                                            plt.figure(figsize=(12, 16))                                                                    
                                                                            Nuts_plot(dict_samples_diagnostics, fit_summary)
                                                                            plt.savefig(os.path.join(Reportmetrics_dir, 'figure-nuts_'+csvfile_name+'.png'), dpi=300) 
     
                                                                            from  report_diagnostics_stan import  ks_plot
                                                                            plt.figure(figsize=(12, 8))  
                                                                            ks_plot(loglik, len_seizure, num_seizure)
                                                                            plt.savefig(os.path.join(Reportmetrics_dir, 'figure-ks_'+csvfile_name+'.png'), dpi=300) 
         

                                                            if alg_name=="advi":
                                                                    from  report_diagnostics_stan import Elbo_plot
                                                                    divergence__=np.nan
                                                                    treedepth__=np.nan
                                                                    energy__=np.nan
                                                                    n_eff_low=np.nan
                                                                    rhat_large=np.nan
                                                                    elbo_val=np.max(Elbo)    

                                                                    if opt_plot_figure_diagnostics==1:
                                                                            plt.figure(figsize=(12, 4)) 
                                                                            Elbo_plot(Iter, Elbo)                                                                  
                                                                            plt.savefig(os.path.join(Reportmetrics_dir, 'figure-elbo_'+csvfile_name+'.png'), dpi=300) 
                                                                            
                                                                            from  report_diagnostics_stan import  ks_plot
                                                                            plt.figure(figsize=(12, 8))  
                                                                            ks_plot(loglik, len_seizure, num_seizure)
                                                                            plt.savefig(os.path.join(Reportmetrics_dir, 'figure-ks_'+csvfile_name+'.png'), dpi=300) 
         

                                                            all_divergence_list.append(divergence__)
                                                            all_treedepth_list.append(treedepth__)
                                                            all_energy_list.append(energy__)
                                                            all_n_eff_low_percent_list.append(n_eff_low)
                                                            all_Rhat_large_percent_list.append(rhat_large)
                                                            all_Elbo_list.append(elbo_val)
     

                                                            if opt_plot_figure_fit==1:
                                                                    from report_plots import plot_features, plot_zscore_shrinkage,  plot_posterior, ppplot, plot_obs_x, plot_obs_z, plot_hiddenstates_x, plot_hiddenstates_z, plot_confusion_matrix

                                                                    if patient_name=='patient1':
                                                                        showpicks = f_vector_int(np.r_[int(Hz_indices[0]), int(Pz_indices[3]), int(Ez_indices[0])])
                                                                        showpicks=showpicks.tolist()
                                                                        print("Shown nodes are:", [(x+1) for x in showpicks])
                                                                    if patient_name=='patient2':
                                                                        showpicks = f_vector_int(np.r_[int(Hz_indices[0]), int(Pz_indices[0]), int(Ez_indices[0])])
                                                                        showpicks=showpicks.tolist()
                                                                        print("Shown nodes are:", [(x+1) for x in showpicks])

                                

                                                                    plt.figure(figsize=(12, 5))  
                                                                    grid = plt.GridSpec(2, 4)  
                                                                    plt.subplot(grid[0, :2]) 
                                                                    plot_features(0.1*ts, Obs, Obs_est, showpicks)
                                                                    plt.subplot(grid[0, 2]) 
                                                                    plot_zscore_shrinkage(nodes, eta_true_mu, eta_est_mu, eta_est_std, prior_std)
                                                                    plt.subplot(grid[0, 3]) 
                                                                    plot_confusion_matrix(confusion_matrix_values, ["EZ", "PZ", "HZ"])
                                                                    plt.subplot(grid[1, 0:])            
                                                                    plot_posterior(nodes, eta_c, delta_eta, eta_hz, eta_ez, eta_pz, np.asarray(Hz_indices), np.asarray(Ez_indices), np.asarray(Pz_indices), eta_est)
                                                                    plt.savefig(os.path.join(Reportmetrics_dir, 'figure-fit_'+csvfile_name+'.png'), dpi=300) 
                                                                    
                                                                    # ppplot(eta_est,  eta_hz, eta_ez, eta_pz, showpicks)

                                                                    from report_plots import func_2DepileptorEqs,  plot_phasespace
                                                                    from read_csvs import parse_csv
                                                                    fit = parse_csv(filecsv)

                                                                    I=I1+1;
                                                                    InitialGuess=np.array([[-1.0, 3.0]])
                                                                    rGuess = np.repeat(InitialGuess, int(SC.shape[0]), axis=0)

                                                                    true_roots_K0 = fsolve(func_2DepileptorEqs,rGuess, args=(eta_true_mu, 0.0, I, SC))
                                                                    true_roots = fsolve(func_2DepileptorEqs,rGuess, args=(eta_true_mu, Ks, I, SC))

                                                                    estimated_roots_K0 = fsolve(func_2DepileptorEqs,rGuess, args=(eta_est_mu, 0.0, I, SC))
                                                                    estimated_roots = fsolve(func_2DepileptorEqs,rGuess, args=(eta_est_mu, K_est, I, SC))

                                                                    plt.figure()           
                                                                    plot_phasespace(fit, npz, nodes, showpicks, SC,  eta_true_mu, Ks,  true_roots, true_roots_K0, estimated_roots, estimated_roots_K0)
                                                                    plt.savefig(os.path.join(Reportmetrics_dir, 'figure-phaseplane_'+csvfile_name+'.png'), dpi=300) 

                                                        
                                                            if opt_pairplot_fit==1:
                                                                    from read_csvs import parse_csv
                                                                    from report_plots import pair_plots_params
                                                                    fit = parse_csv(filecsv)
                                                                    extras = ' amplitude offset K eps sig'.split()              
                                                                    plt.figure(figsize=(8, 6)) 
                                                                    pair_plots_params(fit, extras)  
                                                                    plt.tight_layout()
                                                                    plt.savefig(os.path.join(Reportmetrics_dir, 'figure-pairplot_'+csvfile_name+'.png'), dpi=300) 

                                                            counter_converged_csv = counter_converged_csv+1

                                                    else:

                                                            print ('The algorithm has not been converged!')

                                                            all_csvfile_convergence_list.append('Notconverged')

                                                            err_eta, z_score_eta, LSE_val=np.inf, np.inf, np.inf
                                                            shrinkage_eta, score_ez, score_pz , score_std, score_extrm, accuracy =np.nan, 0., 0., 0., 0., 0.

                                                            all_LSE_list.append(LSE_val)
                                                            all_err_eta_list.append(err_eta)
                                                            all_z_score_eta_list.append(z_score_eta)
                                                            all_shrinkage_eta_list.append(shrinkage_eta)
                                                            all_score_std_list.append(int(score_std))
                                                            all_score_extrm_list.append((score_extrm))
                                                            all_score_ez_list.append(int(score_ez))
                                                            all_score_pz_list.append(int(score_pz))
                                                            all_accuracy_list.append(accuracy)


                                                            counter_Notconverged_csv = counter_Notconverged_csv+1




counter_all_csv=counter_converged_csv+counter_Notconverged_csv
                                                                       
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
print ("-"*80)
print ('The total csv files that did Not converged are:', counter_Notconverged_csv)
print ('The total csv files that converged are:', counter_converged_csv)
print ('The total csv files are:', counter_all_csv)
print ("-"*80)
####################################################################################################################################################################################################################
number_allcsv_converged=range(counter_converged_csv)
number_allcsv=range(counter_all_csv)
####################################################################################################################################################################################################################with open(Reportmetrics_dir+'Reportmetrics_all.txt', 'w') as fd:
with open(Reportmetrics_dir+'ReportMetrics_allcsvfiles.txt', 'w') as fd:
    fd.write("Index, csvname, convergence, LSE, Err, zscore, shrinkage, score_std, score_extrm, scoreEZ, scorePZ, Accuracy\n")
    fmt = "%d,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
    for idx in range(counter_all_csv):
            args=idx+1, all_csvfile_list[idx], all_csvfile_convergence_list[idx], all_LSE_list[idx], all_err_eta_list[idx], all_z_score_eta_list[idx], all_shrinkage_eta_list[idx], all_score_std_list[idx], all_score_extrm_list[idx], all_score_ez_list[idx], all_score_ez_list[idx], all_accuracy_list[idx]
            fd.write(fmt % args)
####################################################################################################################################################################################################################
with open(Reportmetrics_dir+'Reportdiagnostics_allcsvfiles.txt', 'w') as fd:
    fd.write("Index, csvname, divergence, treedepth, energy, n_eff, Rhat_large_percent, ELBO\n")
    fmt = "%d,%s,%f,%f,%f,%f,%f,%f\n"
    for idx in range(counter_converged_csv):
            args=idx+1,  all_csvfile_converged_list[idx], all_divergence_list[idx], all_treedepth_list[idx], all_energy_list[idx], all_n_eff_low_percent_list[idx], all_Rhat_large_percent_list[idx], all_Elbo_list[idx]
            fd.write(fmt % args)
####################################################################################################################################################################################################################
with open(Reportmetrics_dir+'ReportLikelihood_allcsvfiles.txt', 'w') as fd:
    fd.write("Index, csvname, maxlike, p_waic, elpd_waic, waic, loo, ks_large_percent\n")
    fmt = "%d,%s,%f,%f,%f,%f,%f,%f\n"
    for idx, idx_nonconv in zip(number_allcsv_converged, number_allcsv):
            args=idx+1, all_csvfile_converged_list[idx], all_maxlike_list[idx], all_p_waic_list[idx], all_elpd_waic_list[idx], all_waic_list[idx], all_loo_list[idx], all_ks_large_percent_list[idx]
            fd.write(fmt % args)
####################################################################################################################################################################################################################
endtime = time.time()
print('Report running time:', endtime - starttime)
####################################################################################################################################################################################################################

