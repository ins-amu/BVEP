#!/usr/bin/env python3
"""
@author: meysamhashemi INS Marseille


"""

import os
import sys

def converged_hmc(filename):
    pass_hmc=0
    if 'Elapsed Time' in open(filename).read():
        pass_hmc+=1
    return pass_hmc


if __name__ == "__main__":
  
    out_pass=converged_hmc(sys.argv[1])

    if out_pass==1:
       print('hmc converged')   
    sys.exit(0)

