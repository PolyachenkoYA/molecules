#!/usr/bin/python3 -tt
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
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./full_proc.py      model_name     [keys,    N0,     N1]')
        sys.exit(1)
        
    model_name = args[0]    
    if(len(args) > 1):
        keys = args[1]
        if(keys == '-all'):
            keys = '-gen' + my.all_proc_flags
    else:
        keys = my.all_proc_flags
    keys = keys.split('-')
        
    extra_args_str = ''
    for i in range(1,len(args)):
        extra_args_str += (args[i] + ' ')
        
    if(my.find_key(keys, 'gen')):
        if(not my.run_it('./gen ' + model_name)):
            return

    if(not my.run_it('./comp ' + model_name + '            ')):
        return

    my.run_it('mv ' + model_name + '_gen.log ' + os.path.join('./', model_name, 'gen.log') +  '            ')
    
    if(not my.run_it('./post_proc ' + model_name + ' 1' + '                               ')):
        return
    
    if(my.find_key(keys, 'cond') or my.find_key(keys, 'condition')):
        if(not my.run_it('./post_proc ' + model_name + ' 3 1000' + '                                 ')):
            return
    
    if(my.find_key(keys, 'pics') or my.find_key(keys, 'keep') or my.find_key(keys, 'move_res')):
        command_to_run = 'cd ../../' + '            '
        print(command_to_run)
        os.chdir('../../')
        
        my.run_it('mv ' + os.path.join('./', '!COMP', '!go', model_name) + ' ' + os.path.join('./', 'RES', 'DATA') + '            ')
        
        command_to_run = 'cd ./RES' + '            '
        print(command_to_run)
        os.chdir('./RES')   
        
        my.run_it('./full_post_proc.py ' + model_name + ' ' + extra_args_str)
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
