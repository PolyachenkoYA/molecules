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
        print('usage: python ' + sys.argv[0] + '    model_name     [keys,    N0,     N1, ...]')
        sys.exit(1)
        
    model_name, keys, model_dir, graph_dir, time_gaps_str, N0, N1, Nfrm, E, P, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 1, 2, 3)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
    
    extra_args_str = ''
    for i in range(1,len(args)):
        extra_args_str += (args[i] + ' ')
        
    if(my.find_key(keys, 'energy')):
        my.run_it('python E.py ' + model_name + ' ' + extra_args_str)
        
    if(my.find_key(keys, 'maxwell')):
        if(extra_args_str):
            command_to_run = 'python maxwell.py ' + model_name + ' 100 ' + extra_args_str
        else:
            command_to_run = 'python maxwell.py ' + model_name + ' 100 ' + str(stabTind)
        my.run_it(command_to_run)
        
    if(my.find_key(keys, 'diff') or my.find_key(keys, 'diffusion')):
        my.run_it('python diffusion.py ' + model_name + ' '+ extra_args_str)

    if(my.find_key(keys, 'cond') or my.find_key(keys, 'condition')):
        if(extra_args_str):
            command_to_run = 'python condition.py ' + model_name + ' ' + extra_args_str
        else:
            command_to_run = 'python condition.py ' + model_name + ' ' + str(stabTind)
        my.run_it(command_to_run)

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
