#!/home/ypolyach/anaconda3/bin/python3 -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import numpy as np
import math as mth
import os
import matplotlib.pyplot as plt

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print('usage: ' + sys.argv[0] + '     model_name     [keys,    N0,     N1]')
        sys.exit(1)
        
    model_name = args[0]
    if(len(args) > 1):
        keys = args[1]
        if(keys == '-all'):
            keys = '-gen' + my.all_proc_flags
    else:
        keys = my.all_proc_flags
    keys = keys.split('-')
    
    my.full_cycle_fnc(model_name, keys, args)
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
