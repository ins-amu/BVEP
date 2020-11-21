#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:27:30 2017

@author: meysamhashemi
"""

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from statsmodels.graphics import tsaplots

from pandas.tools.plotting import autocorrelation_plot


data = theta_est



def label(ax, string):
    ax.annotate(string, (1, 1), xytext=(-8, -8), ha='right', va='top',
                size=14, xycoords='axes fraction', textcoords='offset points')


fig, axes = plt.subplots(nrows=2, figsize=(10, 4))
fig.tight_layout()

axes[0].plot(data)
label(axes[0], 'MCMC')

#axes[1].acorr(data, maxlags=data.size-1)
#label(axes[1], 'Matplotlib Autocorrelation')


pd.tools.plotting.autocorrelation_plot(data, ax=axes[1])
label(axes[1], 'Autocorrelation')

# Remove some of the titles and labels that were automatically added
for ax in axes.flat:
    ax.set(title='', xlabel='')
plt.show()



    
    

    

def acorr(x, ax=None):
    if ax is None:
        ax = plt.gca()
    n = len(x)
    x = x - x.mean()
    autocorr = np.correlate(x, x, mode='full')
    autocorr = autocorr[x.size:]
    c0 = sum((x - x.mean()) ** 2) / float(n)
    autocorr =autocorr/float(n) / c0
    #autocorr /= autocorr.max()
    return ax.stem(autocorr)

figure(figsize=(12, 2))
#fig, ax = plt.subplots()
acorr(data)
plt.title('Autocorrelation ')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.show()    













figure(figsize=(12, 2))
plt.acorr(data, usevlines=True, normed=True, maxlags=len(data)/4, lw=2)
plt.title('Autocorrelation ')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.show()



def acf(series):
    n = len(series)
    D = asarray(series)
    c0 = sum((D - D.mean()) ** 2) / float(n)

    def r(h):
        acf_lag = ((D[:n - h] - D.mean()) * (D[h:] - D.mean())).sum() / float(n) / c0
        return round(acf_lag, 3)
    x = arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs



def autocorr(wave):
    lags = range(len(wave.ys)//2)
    corrs = [serial_corr(wave, lag) for lag in lags]
    return lags, corrs   

def serial_corr(wave, lag=1):
    n = len(wave)
    y1 = wave.ys[lag:]
    y2 = wave.ys[:n-lag]
    corr = np.corrcoef(y1, y2, ddof=0)[0, 1]
    return corr