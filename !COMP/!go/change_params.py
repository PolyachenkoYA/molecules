#!/home/ypolyach/anaconda3/bin/python3 -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import os

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc = len(args) - 2
    if((argc < 0) or (argc % 2 != 0)):
        print('usage:\npython' + sys.argv[0] + '      old_name      new_name     [key1   value1   key2   value2 ...]')
        sys.exit(1)
    
    my.change_params_fnc(args[0] + '_' + my.param_file_suff, args[1] + '_' + my.param_file_suff, args[2:])
    print('DONE')
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
