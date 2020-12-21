#!/usr/bin/env python3
"""
@author: meysamhashemi INS Marseille

"""

import os
import sys

def converged_advi(filename):
    pass_advi=0
    if 'COMPLETED.' in open(filename).read():
        pass_advi+=1
    return pass_advi


if __name__ == "__main__":
  
    out_pass=converged_advi(sys.argv[1])

    if out_pass==1:
       print('advi converged')   
    sys.exit(0)

