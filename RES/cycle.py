#!/home/ypolyach/anaconda3/bin/python3 -tt
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
    argc_min = 0
    argc = len(args)
    if(argc < argc_min):
        print('usage: ./cycle.py')
        sys.exit(1)
            
    n = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    T = [0.01, 0.02, 0.04, 0.065, 0.1, 0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 2, 3.2, 6.4, 12.8, 25.6]
    
    n = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    T = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]    
    
    n = [0.2, 0.3]
    T = [1, 2]
    
    n = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    T = [0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 2, 3.2]
    
    n = [0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
    #T = [0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4]       
    T = [0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 1.75, 2]
    
    for ni in n:
        for Ti in T:
            new_name = 'N2048/n' + my.str_sgn_round(ni, 3) + '_Tmp' + my.str_sgn_round(Ti, 3)
            
            #os.system('./E.py ' + new_name + ' -')
            
            if(not my.run_it('./diffusion.py ' + new_name + ' -pics-percent')):
                return
            #os.system('./diffusion.py ' + new_name + ' -') # just lmd(n,T)
            #os.system('./E.py ' + new_name + ' -full-pics')
            
            #os.system('cp -rv ./bfr_N100/' + new_name + '/graph ./11/' + new_name)
            
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
