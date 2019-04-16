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
    argc_min = 2
    argc = len(args)
    if((argc < argc_min) or (argc % 2 != 0)):
        print('usage: ./change_params.py      old_name      new_name     [key1   value1   key2   value2 ...]')
        sys.exit(1)
    
    old_name = args[0] + '_param.dat'
    new_name = args[1] + '_param.dat'
    params = my.read_params(old_name)
    
    i = 1
    while(i*2 < argc):
        params[args[2*i]] = args[2*i+1]
        i += 1
    
    my.save_params(new_name, params)
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
