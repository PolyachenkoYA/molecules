import sys
import os

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print('usage: ' + sys.argv[0] + '     model_name     input_params_name   [keys,    N0,     N1]')
        sys.exit(1)
        
    model_name = args[0]
    input_params_name = args[1]
    if(len(args) > 2):
        keys = args[2]
        if(keys == '-all'):
            keys = '-gen' + my.all_proc_flags
    else:
        keys = my.all_proc_flags
    keys = keys.split('-')
    
    my.full_cycle_fnc(input_params_name, model_name, keys, args)
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
