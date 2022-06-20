#!/usr/bin/env python3
"""
@author: meysamhashemi  INS Marseille

"""
import os
import sys
import numpy as np
import scipy as scp

from scipy import signal
from scipy.signal import hilbert

from scipy.signal import savgol_filter
from scipy.signal import find_peaks

from scipy import stats as spstats
from scipy.stats import moment
from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import mode


from numba import jit
######################################################
@jit(nopython=True) 
def normalized(x_array):
    x_array /= np.max(np.abs(x_array))
    return x_array
######################################################
def softmax(x):
    return np.exp(x) / np.sum(np.exp(x), axis=0)    
######################################################

############################################################################################################
def calculate_summary_statistics_features(x, nn, nt, dt, ts, features):
    """Calculate summary statistics

    Parameters
    ----------
    x : output of the simulator

    Returns
    -------
    np.array, summary statistics
    """
    

    X=x.reshape(nn, nt)

    n_summary = 100*nn

    sum_stats_vec = np.concatenate((np.mean(X, axis=1), 
                                    np.median(X, axis=1),
                                    np.std(X, axis=1),
                                    skew(X, axis=1), 
                                    kurtosis(X, axis=1),
                                    ))


    for item in features:


            if item is 'higher_moments':

                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                                moment(X, moment=2, axis=1),
                                                moment(X, moment=3, axis=1),
                                                moment(X, moment=4, axis=1),
                                                moment(X, moment=5, axis=1),
                                                moment(X, moment=6, axis=1),
                                                moment(X, moment=7, axis=1),
                                                moment(X, moment=8, axis=1),
                                                moment(X, moment=9, axis=1), 
                                                moment(X, moment=10, axis=1),        
                                                                   ))


            if item is 'power_envelope':

                        X_area = np.trapz(X, dx=0.0001)
                        X_pwr = np.sum((X*X), axis=1)
                        X_pwr_n = (X_pwr/ X_pwr.max())

                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                                        X_area, X_pwr, X_pwr_n,            
                                                                    ))



            if item is 'spectral_power':

                        fs = 10e3 

                        f, Pxx_den =  signal.periodogram(X, fs)

                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                                        np.max(Pxx_den, axis=1), 
                                                        np.mean(Pxx_den, axis=1),
                                                        np.median(Pxx_den, axis=1),
                                                        np.std(Pxx_den, axis=1),
                                                        skew(Pxx_den, axis=1), 
                                                        kurtosis(Pxx_den, axis=1), 
                                                        moment(Pxx_den, moment=2, axis=1),
                                                        moment(Pxx_den, moment=3, axis=1),
                                                        moment(Pxx_den, moment=4, axis=1),
                                                        moment(Pxx_den, moment=5, axis=1),
                                                        moment(Pxx_den, moment=6, axis=1),
                                                        moment(Pxx_den, moment=2, axis=1),
                                                        moment(Pxx_den, moment=3, axis=1),
                                                        moment(Pxx_den, moment=4, axis=1),
                                                        moment(Pxx_den, moment=5, axis=1),
                                                        moment(Pxx_den, moment=6, axis=1),
                                                        moment(Pxx_den, moment=7, axis=1),
                                                        moment(Pxx_den, moment=8, axis=1),
                                                        moment(Pxx_den, moment=9, axis=1),
                                                        moment(Pxx_den, moment=10, axis=1),
                                                        np.diag(np.dot(Pxx_den, Pxx_den.transpose())),
                                                                           ))

            if item is 'amplitude_phase':

                        analytic_signal = hilbert(X)
                        amplitude_envelope = np.abs(analytic_signal)
                        instantaneous_phase = np.unwrap(np.angle(analytic_signal))

                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                                        np.mean(amplitude_envelope, axis=1),
                                                        np.median(amplitude_envelope, axis=1),
                                                        np.std(amplitude_envelope, axis=1),
                                                        skew(amplitude_envelope, axis=1), 
                                                        kurtosis(amplitude_envelope, axis=1),                 
                                                        np.mean(instantaneous_phase, axis=1),
                                                        np.median(instantaneous_phase, axis=1),
                                                        np.std(instantaneous_phase, axis=1),
                                                        skew(instantaneous_phase, axis=1), 
                                                        kurtosis(instantaneous_phase, axis=1), 
                                                                           ))



            if item is 'seizures_onset':

                        seizures_num=[]
                        seizures_on=[]
                        
                        for i in np.r_[0:nn]:
                                    v=np.zeros(nt)
                                    v= np.array(X[i,:])

                                    v_th=-1
                                    ind = np.where(v < v_th)
                                    v[ind] = v_th

                                    ind = np.where(np.diff(v) < 0)

                                    seizure_times = np.arange(0, ts.shape[0], dt)[ind]
                                    #seizure_times = np.array(ts)[ind]
                                    seizure_times_stim = seizure_times

                                    if seizure_times_stim.shape[0] > 0:
                                                seizure_times_stim = seizure_times_stim[np.append(1, np.diff(seizure_times_stim)) > .5]
                                                seizures_on.append(seizure_times_stim[0])
                                    else:            
                                                seizures_on.append(100000.)

                                    seizures_num.append(seizure_times_stim.shape[0])  


                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                        np.array(seizures_num),
                                        np.array(seizures_on),
                                                       ))


            if item is 'seizures_offset':

                        seizures_num=[]
                        seizures_off=[]
                        
                        for i in np.r_[0:nn]:
                                    v=np.zeros(nt)
                                    v= np.array(X[i,:])

                                    v_th=-1
                                    ind = np.where(v < v_th)
                                    v[ind] = v_th

                                    ind = np.where(np.diff(v) < 0)

                                    seizure_times = np.arange(0, ts.shape[0], dt)[ind]
                                    #seizure_times = np.array(ts)[ind]
                                    seizure_times_stim = seizure_times

                                    if seizure_times_stim.shape[0] > 0:
                                                seizures_off.append(seizure_times_stim[-1])
                                                seizure_times_stim = seizure_times_stim[np.append(1, np.diff(seizure_times_stim)) > .5]
                                    else:            
                                                seizures_off.append(100000.)

                                    seizures_num.append(seizure_times_stim.shape[0])  


                        sum_stats_vec = np.concatenate((sum_stats_vec,
                                        np.array(seizures_num),
                                        np.array(seizures_off),
                                                       ))




        
    sum_stats_vec = sum_stats_vec[0:n_summary]        


    return sum_stats_vec


############################################################################################################



def calculate_summary_statistics(x, nn, ns, nt): 
    """Calculate summary statistics

    Parameters
    ----------
    x : output of the simulator

    Returns
    -------
    np.array, summary statistics
    """
    

    X=x.reshape(ns, nt)

    n_summary = 100*nn
    
           
    X_area = np.trapz(X, dx=0.0001)
    X_pwr = np.sum((X*X), axis=1)
    X_pwr_n = (X_pwr/ X_pwr.max())
        
    sum_stats_vec = np.concatenate((np.mean(X, axis=1), 
                                    np.median(X, axis=1),
                                    np.std(X, axis=1),
                                    skew(X, axis=1), 
                                    kurtosis(X, axis=1),
                                    moment(X, moment=2, axis=1),
                                    moment(X, moment=3, axis=1),
                                    moment(X, moment=4, axis=1),
                                    moment(X, moment=5, axis=1),
                                    moment(X, moment=6, axis=1),
                                    moment(X, moment=7, axis=1),
                                    moment(X, moment=8, axis=1),
                                    moment(X, moment=9, axis=1),
                                    moment(X, moment=10, axis=1),
                                    X_area, X_pwr, X_pwr_n, 
                                    ))
   
    seizures_num=[]
    seizures_on=[]


    for i in np.r_[0:ns]:
                    v=np.zeros(nt)
                    Xhat=np.zeros(nt)
                    Xhat = savgol_filter(X[i,:], 11, 3)

                    v= np.array(Xhat)

                    peaks, _ = find_peaks(v, height=0, rel_height=.3, width=5)


                    if peaks.shape[0] > 0:
                                seizures_on.append(peaks[0])
                    else:            
                                seizures_on.append(100000.)

                    seizures_num.append(peaks.shape[0])  

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    np.array(seizures_num),
                                    np.array(seizures_on),
                                     ))


    fs = 10e3 

    f, Pxx_den =  signal.periodogram(X, fs)

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    np.max(Pxx_den, axis=1), 
                                    np.mean(Pxx_den, axis=1),
                                    np.median(Pxx_den, axis=1),
                                    np.std(Pxx_den, axis=1),
                                    skew(Pxx_den, axis=1), 
                                    kurtosis(Pxx_den, axis=1), 
                                    moment(Pxx_den, moment=2, axis=1),
                                    moment(Pxx_den, moment=3, axis=1),
                                    moment(Pxx_den, moment=4, axis=1),
                                    moment(Pxx_den, moment=5, axis=1),
                                    moment(Pxx_den, moment=6, axis=1),
                                    moment(Pxx_den, moment=2, axis=1),
                                    moment(Pxx_den, moment=3, axis=1),
                                    moment(Pxx_den, moment=4, axis=1),
                                    moment(Pxx_den, moment=5, axis=1),
                                    moment(Pxx_den, moment=6, axis=1),
                                    moment(Pxx_den, moment=7, axis=1),
                                    moment(Pxx_den, moment=8, axis=1),
                                    moment(Pxx_den, moment=9, axis=1),
                                    moment(Pxx_den, moment=10, axis=1),
                                    np.diag(np.dot(Pxx_den, Pxx_den.transpose())),
                                                       ))

            
    analytic_signal = hilbert(X)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    np.mean(amplitude_envelope, axis=1),
                                    np.median(amplitude_envelope, axis=1),
                                    np.std(amplitude_envelope, axis=1),
                                    skew(amplitude_envelope, axis=1), 
                                    kurtosis(amplitude_envelope, axis=1),                 
                                    np.mean(instantaneous_phase, axis=1),
                                    np.median(instantaneous_phase, axis=1),
                                    np.std(instantaneous_phase, axis=1),
                                    skew(instantaneous_phase, axis=1), 
                                    kurtosis(instantaneous_phase, axis=1), 
                                                       ))
                            
            
    sum_stats_vec = sum_stats_vec[0:n_summary]        


    return sum_stats_vec
############################################################################################################
############################################################################################################
def calculate_summary_statistics_normalized(x, nn, ns, nt): 
    """Calculate summary statistics

    Parameters
    ----------
    x : output of the simulator

    Returns
    -------
    np.array, summary statistics
    """
    
    from scipy import signal
    from scipy.signal import hilbert

    from scipy.signal import savgol_filter
    from scipy.signal import find_peaks

    from scipy import stats as spstats
    from scipy.stats import moment
    from scipy.stats import kurtosis
    from scipy.stats import skew
    from scipy.stats import mode
    

    X=x.reshape(ns, nt)

    n_summary = 100*nn
    
           
    X_area = np.trapz(X, dx=0.0001)
    X_pwr = np.sum((X*X), axis=1)
    X_pwr_n = (X_pwr/ X_pwr.max())
        
       
    sum_stats_vec = np.concatenate((normalized(np.mean(X, axis=1)), 
                                    normalized(np.median(X, axis=1)), 
                                    normalized(np.std(X, axis=1)), 
                                    normalized(skew(X, axis=1)), 
                                    normalized(kurtosis(X, axis=1)), 
                                    normalized(moment(X, moment=2, axis=1)), 
                                    normalized(moment(X, moment=3, axis=1)), 
                                    normalized(moment(X, moment=4, axis=1)), 
                                    normalized(moment(X, moment=5, axis=1)), 
                                    normalized(moment(X, moment=6, axis=1)), 
                                    normalized(moment(X, moment=7, axis=1)), 
                                    normalized(moment(X, moment=8, axis=1)), 
                                    normalized(moment(X, moment=9, axis=1)), 
                                    normalized(moment(X, moment=10, axis=1)), 
                                    normalized(X_area), 
                                    normalized(X_pwr), 
                                    normalized(X_pwr_n), 
                                    ))
   
    seizures_num=[]
    seizures_on=[]


    for i in np.r_[0:ns]:
                    v=np.zeros(nt)
                    Xhat=np.zeros(nt)
                    Xhat = savgol_filter(X[i,:], 11, 3)

                    v= np.array(Xhat)

                    peaks, _ = find_peaks(v, height=0, rel_height=.3, width=5)


                    if peaks.shape[0] > 0:
                                seizures_on.append(peaks[0])
                    else:            
                                seizures_on.append(100000.)

                    seizures_num.append(peaks.shape[0])  

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    np.array(seizures_num),
                                    np.array(seizures_on),
                                     ))


    fs = 10e3 

    f, Pxx_den =  signal.periodogram(X, fs)

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    normalized(np.max(Pxx_den, axis=1)), 
                                    normalized(np.mean(Pxx_den, axis=1)),
                                    normalized(np.median(Pxx_den, axis=1)),
                                    normalized(np.std(Pxx_den, axis=1)),
                                    normalized(skew(Pxx_den, axis=1)), 
                                    normalized(kurtosis(Pxx_den, axis=1)), 
                                    normalized(moment(Pxx_den, moment=2, axis=1)),
                                    normalized(moment(Pxx_den, moment=3, axis=1)),
                                    normalized(moment(Pxx_den, moment=4, axis=1)),
                                    normalized(moment(Pxx_den, moment=5, axis=1)),
                                    normalized(moment(Pxx_den, moment=6, axis=1)),
                                    normalized(moment(Pxx_den, moment=2, axis=1)),
                                    normalized(moment(Pxx_den, moment=3, axis=1)),
                                    normalized(moment(Pxx_den, moment=4, axis=1)),
                                    normalized(moment(Pxx_den, moment=5, axis=1)),
                                    normalized(moment(Pxx_den, moment=6, axis=1)),
                                    normalized(moment(Pxx_den, moment=7, axis=1)),
                                    normalized(moment(Pxx_den, moment=8, axis=1)),
                                    normalized(moment(Pxx_den, moment=9, axis=1)),
                                    normalized(moment(Pxx_den, moment=10, axis=1)),
                                    np.diag(np.dot(Pxx_den, Pxx_den.transpose())),
                                                       ))

            
    analytic_signal = hilbert(X)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    normalized(np.mean(amplitude_envelope, axis=1)),
                                    normalized(np.median(amplitude_envelope, axis=1)),
                                    normalized(np.std(amplitude_envelope, axis=1)),
                                    normalized(skew(amplitude_envelope, axis=1)), 
                                    normalized(kurtosis(amplitude_envelope, axis=1)),                 
                                    normalized(np.mean(instantaneous_phase, axis=1)),
                                    normalized(np.median(instantaneous_phase, axis=1)),
                                    normalized(np.std(instantaneous_phase, axis=1)),
                                    normalized(skew(instantaneous_phase, axis=1)), 
                                    normalized(kurtosis(instantaneous_phase, axis=1)), 
                                                       ))
                            
            
    sum_stats_vec = sum_stats_vec[0:n_summary]        


    return sum_stats_vec

############################################################################################################
def calculate_summary_statistics_softmax(x, nn, ns, nt): 
    """Calculate summary statistics

    Parameters
    ----------
    x : output of the simulator

    Returns
    -------
    np.array, summary statistics
    """
    
    from scipy import signal
    from scipy.signal import hilbert

    from scipy.signal import savgol_filter
    from scipy.signal import find_peaks

    from scipy import stats as spstats
    from scipy.stats import moment
    from scipy.stats import kurtosis
    from scipy.stats import skew
    from scipy.stats import mode
    

    X=x.reshape(ns, nt)

    n_summary = 100*nn
    
           
    X_area = np.trapz(X, dx=0.0001)
    X_pwr = np.sum((X*X), axis=1)
    X_pwr_n = (X_pwr/ X_pwr.max())
        
       
    sum_stats_vec = np.concatenate((softmax(np.mean(X, axis=1)), 
                                    softmax(np.median(X, axis=1)), 
                                    softmax(np.std(X, axis=1)), 
                                    softmax(skew(X, axis=1)), 
                                    softmax(kurtosis(X, axis=1)), 
                                    (moment(X, moment=2, axis=1)), 
                                    (moment(X, moment=3, axis=1)), 
                                    (moment(X, moment=4, axis=1)), 
                                    (moment(X, moment=5, axis=1)), 
                                    (moment(X, moment=6, axis=1)), 
                                    (moment(X, moment=7, axis=1)), 
                                    (moment(X, moment=8, axis=1)), 
                                    (moment(X, moment=9, axis=1)), 
                                    (moment(X, moment=10, axis=1)), 
                                    softmax(X_area), 
                                    (X_pwr), 
                                    softmax(X_pwr_n), 
                                    ))
   
    seizures_num=[]
    seizures_on=[]


    for i in np.r_[0:ns]:
                    v=np.zeros(nt)
                    Xhat=np.zeros(nt)
                    Xhat = savgol_filter(X[i,:], 11, 3)

                    v= np.array(Xhat)

                    peaks, _ = find_peaks(v, height=0, rel_height=.3, width=5)


                    if peaks.shape[0] > 0:
                                seizures_on.append(peaks[0])
                    else:            
                                seizures_on.append(100000.)

                    seizures_num.append(peaks.shape[0])  

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    softmax(np.array(seizures_num)),
                                    softmax(np.array(seizures_on)),
                                     ))


    fs = 10e3 

    f, Pxx_den =  signal.periodogram(X, fs)

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    softmax(np.max(Pxx_den, axis=1)), 
                                    softmax(np.mean(Pxx_den, axis=1)),
                                    softmax(np.median(Pxx_den, axis=1)),
                                    softmax(np.std(Pxx_den, axis=1)),
                                    softmax(skew(Pxx_den, axis=1)), 
                                    softmax(kurtosis(Pxx_den, axis=1)), 
                                    softmax(moment(Pxx_den, moment=2, axis=1)),
                                    softmax(moment(Pxx_den, moment=3, axis=1)),
                                    softmax(moment(Pxx_den, moment=4, axis=1)),
                                    softmax(moment(Pxx_den, moment=5, axis=1)),
                                    softmax(moment(Pxx_den, moment=6, axis=1)),
                                    softmax(moment(Pxx_den, moment=2, axis=1)),
                                    softmax(moment(Pxx_den, moment=3, axis=1)),
                                    softmax(moment(Pxx_den, moment=4, axis=1)),
                                    softmax(moment(Pxx_den, moment=5, axis=1)),
                                    softmax(moment(Pxx_den, moment=6, axis=1)),
                                    softmax(moment(Pxx_den, moment=7, axis=1)),
                                    softmax(moment(Pxx_den, moment=8, axis=1)),
                                    softmax(moment(Pxx_den, moment=9, axis=1)),
                                    softmax(moment(Pxx_den, moment=10, axis=1)),
                                    np.diag(np.dot(Pxx_den, Pxx_den.transpose())),
                                                       ))

            
    analytic_signal = hilbert(X)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))

    sum_stats_vec = np.concatenate((sum_stats_vec,
                                    softmax(np.mean(amplitude_envelope, axis=1)),
                                    softmax(np.median(amplitude_envelope, axis=1)),
                                    softmax(np.std(amplitude_envelope, axis=1)),
                                    softmax(skew(amplitude_envelope, axis=1)), 
                                    softmax(kurtosis(amplitude_envelope, axis=1)),                 
                                    softmax(np.mean(instantaneous_phase, axis=1)),
                                    softmax(np.median(instantaneous_phase, axis=1)),
                                    softmax(np.std(instantaneous_phase, axis=1)),
                                    softmax(skew(instantaneous_phase, axis=1)), 
                                    softmax(kurtosis(instantaneous_phase, axis=1)), 
                                                       ))
                            
            
    sum_stats_vec = sum_stats_vec[0:n_summary]        


    return sum_stats_vec    
