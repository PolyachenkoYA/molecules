#!/usr/bin/python3 -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import numpy as np
import math
import os
import matplotlib.pyplot as plt

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    argc = len(args)
    if(argc < argc_min):
        print('usage: ./cycle.py start_name')
        sys.exit(1)
    
    base_param_name = args[0];
    params = my.read_params(base_param_name + '_param.dat')
    N = params['Ntot']
    m = params['mu']
    Tmp = params['Tmp']
    
    n0 = 0.3
    Tmp0 = 1
    m0 = 29
    tau_0 = 4
    dt0 = 1024
    dumpDT0 = 32
    
    #n = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    #T = [0.01, 0.02, 0.04, 0.065, 0.1, 0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 2, 3.2, 6.4, 12.8, 25.6]
    #n = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    #T = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]    
    #n = [0.2, 0.3]
    #T = [1, 2]
    
    nT = {}
    alph = [0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.02, 1.04, 1.06, 1.08]    
    nT[0.2] = [_x*0.93 for _x in alph]
    nT[0.15]= [_x*0.9 for _x in alph]
    nT[0.125]= [_x*0.82 for _x in alph]
    nT[0.1]= [_x*0.75 for _x in alph]
    nT[0.075]= [_x*0.64 for _x in alph]
    nT[0.05]= [_x*0.52 for _x in alph]
        
    for ni in nT:
        for Ti in nT[ni]:
            R = math.pow(N/ni, 1/3)/2
            k_t = math.sqrt(m*Tmp0 / (Ti*m0))
            endT = tau_0 * math.pow(n0/ni, 1/3) * k_t
            #dt = dt0 / k_t
            dt = dt0
            dumpDT = dumpDT0 / k_t
            Tmp_str = str(my.str_sgn_round(Ti,5))
            new_name = 'n' + str(my.str_sgn_round(ni, 5)) + '_Tmp' + Tmp_str
            
            if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' Tmp ' + Tmp_str + ' R ' + str(R) + ' dt ' + str(dt) + ' endT ' + str(endT) + ' dumpDT ' + str(dumpDT))):
                return
        
            #if(not my.run_it('./full_cycle.py ' + new_name + ' -energy-diff-pics-gen')):
            if(not my.run_it('./full_cycle.py ' + new_name + ' -all')):
                return
    
            if(not my.run_it('rm ' + new_name + '_param.dat')):
                return
            if(not my.run_it('rm ' + new_name + '_particles.xyz')):
                return

    if(not my.run_it('./clear.sh')):
        return
    
    #os.chdir('../../RES/')
    #if(not my.run_it('./cycle.py > diff_tot.txt')):
    #    return
  
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
