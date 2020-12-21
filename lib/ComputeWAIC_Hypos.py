

import numpy as np
import sys, os


def maxlike(log_lik):
    return np.max(log_lik)

def waic(log_lik):
    lpd=np.sum(np.log(np.mean(np.exp(log_lik), axis=0)))
    p_waic = np.sum(np.var(log_lik, axis=0, dtype=np.float64))
    elpd_waic=lpd - p_waic
    waic = -2 * elpd_waic
    return dict(p_waic=p_waic, elpd_waic=elpd_waic,  waic=waic)

# def waic(log_lik):
#     lpd=np.sum(np.log(np.mean((log_lik), axis=0)))
#     p_waic = np.sum(np.var(np.log(log_lik), axis=0, dtype=np.float64))
#     elpd_waic=lpd - p_waic
#     waic = -2 * elpd_waic
#     return dict(p_waic=p_waic, elpd_waic=elpd_waic,  waic=waic)

# input_filenames = sys.argv[1:2]
# for filename in input_filenames:
#         data = np.load(filename)
#         ndata =data['nt'] 


#filenames = sys.argv[3:]

import re
import glob
from operator import itemgetter

from parse_csv import parse_csv


cwd = os.getcwd()
Outfile_dir=os.path.join(cwd, str(sys.argv[3]))
print(Outfile_dir)

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


output_filenames = sorted(glob.glob(Outfile_dir + "/*.csv"), key=numericalSort)


print('............................................................')
#print(output_filenames)
print('............................................................')

results = {}

elpd_waic_values=[]
p_waic_values=[]


for filename in output_filenames:

        output = parse_csv(filename)

        likelihood = -1*output['log_lik']

        result = waic(likelihood)

        results[filename]=(result)
        elpd_waic_values.append(result['elpd_waic'])
        p_waic_values.append(result['p_waic'])


        #output = np.load(filename)
        
        # nparams=int(output['num_params'][0])
        # print (  'num of parameters:' , nparams)
        # ndata =int(output['num_data'][0])
        # print (  'num of data:' , ndata)
        # nsamples = np.shape(likelihood)[0]
        # print (  'num of samples:' , nsamples)


from tempfile import TemporaryFile
outfile = TemporaryFile()

np.savez('Outfile_WAIC_CV_'+str(sys.argv[1])+'_'+str(sys.argv[2])+'.npz', elpd_waic_values=elpd_waic_values, p_waic_values=p_waic_values ) 

print("Writing WAIC completed")


#np.savez('Outfile_WAIC_CV_'+str(sys.argv[2].split(".")[0])+'.npz', elpd_waic_values=elpd_waic_values, p_waic_values=p_waic_values ) 