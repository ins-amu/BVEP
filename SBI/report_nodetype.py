#!/usr/bin/env python3
"""
@author: meysamhashemi  INS Marseille

"""
import os
import sys
import numpy as np
import re

##########################################################################################################
def common_elements(list1, list2):
    return list(set(list1) & set(list2))
##########################################################################################################
def noncommon_elements(list1, list2):
    return list(set(list1) ^ set(list2))
##########################################################################################################
def f_int(x):
    return np.int(x)

f_vector_int = np.vectorize(f_int)  
##########################################################################################################
def node_EzPzHz_idex(true_labels_nodes):
   
    Ez_indices=[index for index, value in enumerate(true_labels_nodes) if value == 'EZ']
    Pz_indices=[index for index, value in enumerate(true_labels_nodes) if value == 'PZ']
    Hz_indices=[index for index, value in enumerate(true_labels_nodes) if value == 'HZ']

    Ez_arr = np.asarray(Ez_indices,dtype=float)  
    Pz_arr = np.asarray(Pz_indices,dtype=float) 
    Hz_arr = np.asarray(Hz_indices,dtype=float) 

    return  Ez_indices, Pz_indices, Hz_indices, Ez_arr, Pz_arr, Hz_arr
##########################################################################################################
def node_EzPzHz_list(npz, nodes):

    Ez_list_idx = npz['Ez_idx'] 
    Pz_list_idx = npz['Pz_idx'] 
    Hz_list_idx=noncommon_elements((noncommon_elements(nodes, Ez_list_idx)), Pz_list_idx)
    Hz_list_idx=f_vector_int(np.asarray(Hz_list_idx,dtype=float))

    return Ez_list_idx, Pz_list_idx, Hz_list_idx    
##########################################################################################################
def node_score_std_extrm(nodes, x0est, x0true_mu, sd_thr):

    InferredNodes_std=np.array([0.0])  
    InferredNodes_extrm=np.array([0.0])

    nodes_est=[]
    x0_est=[]
                           
    for i, idx in enumerate(nodes):
        x0_i = x0est[:, i]
        nodes_est.append(float(idx))
        x0_est.append(float(np.mean(x0_i)))

        if (np.std(x0_i)<sd_thr) and (np.mean(x0_i)-1*np.sqrt(np.std(x0_i))<=x0true_mu[i]<= np.mean(x0_i)+1*np.sqrt(np.std(x0_i))): 
                InferredNodes_std+=1.0

        if  np.min(x0_i)<=x0true_mu[i]<= np.max(x0_i): 
                InferredNodes_extrm+=1.0
     
   
    score_std=InferredNodes_std/len(nodes)
    score_extrm=InferredNodes_extrm/len(nodes)

    return score_std, score_extrm

##########################################################################################################
def node_score_ezpz(nodes, x0est, true_labels_nodes, x0c, deltax0, x0_hi):

     Ez_indices, Pz_indices, Hz_indices, Ez_arr, Pz_arr, Hz_arr=node_EzPzHz_idex(true_labels_nodes)

     nodes_est=[]
     x0_est=[]
                           
     for i, idx in enumerate(nodes):
        x0_i = x0est[:, i]
        nodes_est.append(float(idx))
        x0_est.append(float(np.mean(x0_i)))

     nodes_est=np.asarray(nodes_est)
     x0_est=np.asarray(x0_est)

     Xr_est=np.array([nodes_est, x0_est]).transpose()

     ch=[]
     ce=[]
     cp=[]
     ci=[]

     est_labels_nodes=[]
     score_ez=np.array([0.0])  
     score_pz=np.array([0.0])     


     for row in Xr_est:
        indexEz = np.where(Ez_arr== row[0])[0]
        indexPz = np.where(Pz_arr== row[0])[0]

        if np.size(indexEz) != 0 and x0c<=row[1]<= x0_hi:
            ch.append([row[0], row[1], 'EZ'])
            ce.append(int(row[0]))
            est_labels_nodes.append('EZ')
            score_ez+=1.0
                                                                         
        elif np.size(indexPz) != 0 and x0c-deltax0<=row[1]<x0c:
            ch.append([row[0], row[1], 'PZ'])
            cp.append(int(row[0]))
            est_labels_nodes.append('PZ')
            score_pz+=1.0    
        else: 
            ch.append([row[0], row[1], 'HZ'])
            ci.append(int(row[0]))
            est_labels_nodes.append('HZ')

     return score_ez, score_pz, ce, cp, nodes_est, est_labels_nodes

##########################################################################################################
##########################################################################################################
  