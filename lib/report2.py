#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 16:14:32 2017

@author: meysamhashemi
"""

import sys
from pylab import *
from parse_csv import parse_csv
import common
import matplotlib.pyplot as plt

npz = np.load('data.R.npz')
Ic = npz['Ic']
ps = npz['seeg_log_power']
W_ = npz['SC']
nodes = npz['nodes']
gain = npz['gain']



reg_names = []
with open('centers.txt', 'r') as fd:
    for i, line in enumerate(fd.readlines()):
        reg_names.append(( line.strip().split()[0]))
        
Ezlist=[ "ctx-lh-posteriorcingulate", "ctx-lh-precuneus"]
Pzlist=[ "ctx-lh-postcentral", "ctx-lh-superiorparietal" ]

Ez_indices = [reg_names.index(c) for c in Ezlist]
Pz_indices = [reg_names.index(c) for c in Pzlist]


Ez = np.asarray(Ez_indices,dtype=float)   
Pz = np.asarray(Pz_indices,dtype=float)   

#Ez=np.array([4, 5, 31, 39, 40, 53], dtype=float)
#Pz=np.array([6,13, 14, 27, 33] , dtype=float) 


x0_lo= -4.0
x0_hi= -1.0
x0c= -2.05
deltax0= 0.4



reportx0 = loadtxt("reportx0.txt", comments="#", delimiter=",",usecols=(0, 1, 2), unpack=False)

Nodes=reportx0[:,0]
x0=reportx0[:,1]
Xr=reportx0[:,0:2]    

def common_elements(list1, list2):
    return list(set(list1) & set(list2))
    
common_nodes=[]
common_nodes=common_elements(Nodes, Ez)

common_nodes.extend(common_elements(Nodes, Pz))

ch = []
Predectivescore=np.array([0.0])  
 
for row in Xr:
    indexEz = np.where(Ez== row[0])[0]
    if np.size(indexEz) != 0:
        if x0c+deltax0<=row[1]<= x0_hi: #EZ: [-1.65, -1.0]
           ch.append([row[0], row[1], 'Ez'])
           Predectivescore+=1.0
           
    indexPz = np.where(Pz== row[0])[0]
    if np.size(indexPz) != 0:
        if x0c<=row[1]<=x0c+deltax0:  #PZ: [-2.05, -1.65]
           ch.append([row[0], row[1], 'Pz'])
           Predectivescore+=1.0     
     
        
        
print 'The Predective Score is', int(Predectivescore),' out of', len(common_nodes), ' common nodes'

pehi=load('pehi.npy')
pelo=load('pelo.npy')
t=load('t.npy')

#showpicks = r_[0, 2]

#figure(figsize=(6, 4))
#for i, (lo, hi) in enumerate(zip(pelo.T, pehi.T)):
#    if i in showpicks:
#        fill_between(t, lo, hi, alpha=0.5, facecolor='b', edgecolor='none')
#plot(ps[:, showpicks], 'k', alpha=0.5)
#xlabel('Time (s)'); ylabel('sEEG(t)'); 
#
#leg = []


y_data=ps[:, 0]

y_model_hi=pehi[:,0]
y_model_lo=pelo[:,0]


def LSE(y1, y2):
    return sum((y1 - y2)**2)


LSE_hi=LSE(y_data, y_model_hi)
LSE_lo=LSE(y_data, y_model_lo)

LSE_val=average([LSE_lo, LSE_hi])

print 'LSE is', LSE_val

show()
