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
import re

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc = len(args)
    if(argc > 1):
        print('usage:\n' + sys.argv[0] + '   [output_filename]')
        sys.exit(1)
    output_filename = (args[0] if argc == 1 else 'eff_file.dat')
            
    N_arr = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]
    #N_arr = [512, 1024]
    
    Np_arr = [1, 2, 3, 4, 5, 6]
    #Np_arr = [1, 2, 4]
    
    base_path = os.path.join('/home', 'ypolyach', '!molecules', 'RES', 'DATA')
    f_total = open(output_filename, 'w')
    for Ni in N_arr:
        for Npi in Np_arr:
            model_name = 'p' + str(Npi) + '_N' + str(Ni)
            log_path = os.path.join(base_path, model_name, my.complog_filename)
            log_file = open(log_path, "r")
            log_data = log_file.read()
            log_file.close()
            
            eff_str = re.search(r'efficiency \(e \= endT\/dt\*N\^2\/real\_t\) = \d\.\d+e\+\d+', log_data).group(0)[-11:]
            f_total.write(eff_str + ' ')
        f_total.write('\n')
    f_total.close()
    print('DONE')

            
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
