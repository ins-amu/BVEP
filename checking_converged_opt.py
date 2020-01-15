#!/usr/bin/env python3
"""
@author: meysamhashemi INS Marseille

"""

import os
import sys

def converged_opt(filename):
    pass_opt=0
    if 'Optimization terminated normally:' in open(filename).read():
        pass_opt+=1
    return pass_opt


if __name__ == "__main__":
  
    out_pass=converged_opt(sys.argv[1])

    if out_pass==1:
       print('opt converged')   
    sys.exit(0)


