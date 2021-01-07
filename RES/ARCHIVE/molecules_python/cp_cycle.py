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
    argc_ok = 0
    argc = len(args)
    if(argc != argc_ok):
        print('usage: ./cycle.py')
        sys.exit(1)
            
    n = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    T = [0.01, 0.02, 0.04, 0.065, 0.1, 0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 2, 3.2, 6.4, 12.8, 25.6]
    
    n = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
    T = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]    
    
    n = [0.2, 0.3]
    T = [1, 2]
    
    n = [0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
    T = [0.2, 0.4, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4]   
    
    #base_path = os.path.join('/home', 'ypolyach', '\\!molecules', 'RES', 'DATA')
    base_path = os.path.join('~', '\\!molecules', 'RES', 'DATA')
    for ni in n:
        for Ti in T:
            new_name = 'N2048/n' + str(my.str_sgn_round(ni, 5)) + '_Tmp' + str(my.str_sgn_round(Ti,5))
            

            main_path = os.path.join(base_path, new_name)
            if(not my.run_it('mkdir --parents '+ os.path.join(main_path, 'graph'))):
                return
            if(not my.run_it('cp -av ' + os.path.join('.', new_name, '*.dat') + ' ' + main_path)):
                return            
            if(not my.run_it('cp -av ' + os.path.join('.', new_name, '*.log') + ' ' + main_path)):
                return            
            if(not my.run_it('cp -av ' + os.path.join('.', new_name, '*.txt') + ' ' + main_path)):
                return                        
            if(not my.run_it('cp -av ' + os.path.join('.', new_name, 'graph/*') + ' ' + os.path.join(main_path, 'graph'))):
                return
            
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
