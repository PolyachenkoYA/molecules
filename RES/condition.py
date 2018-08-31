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
import re

import mylib_molecules as my

def getGx(x_data, i0, params, dr):
    N_slice = math.ceil(params['R']/dr)
    y = np.zeros(N_slice)
    for i in range(params['Ntot']):
        if(i != i0):
            slice_i = math.floor(my.length(my.shiftCrd(x_data[i], params['R']))/dr)
            if(slice_i < N_slice):
                y[slice_i] += 1
            
    k = 1/(params['Ntot']/math.pow(2*params['R'], 3)*4*math.pi*dr)
    return [y[_i]*k/((_i+1)*dr) for _i in range(N_slice)]

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./condition.py    model_name     [N0,     N1,     R/dr]')
        sys.exit(1)
    
    model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 1, 2, 3)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
    
    if(len(args) < 4):
        path = os.path.join(my.RAW_DATA_PATH, model_name)
        all_filenames = [name_of_smth for name_of_smth in os.listdir(path) if os.path.isfile(os.path.join(path, name_of_smth))]
        Nslice = max([int(_x) for _x in re.findall(r'condition_(\d+).dat', ' '.join(all_filenames))])
    else:
        Nslice = int(args[3])
        
    input_filename = 'condition_' + str(Nslice) + '.dat'
    data = my.load_file(model_name, input_filename)
    x_dr = data[0];
    if(len(x_dr) != Nslice):
        print('wrong condition file\n' + input_filename + '\n' + 'Nslice = ' + str(Nslice) + '; len(x_dr) = ' + str(len(x_dr)))
        exit(1)
    data = data[1:]
    if(N1 > len(data)):
        print('wrong condition file\n' + input_filename + '\n' + 'max_frame_n = ' + str(len(data)-1))
        exit(1)
    draw_on_screen = my.find_key(keys, 'keep')
    if(not draw_on_screen):
        plt.switch_backend('Agg') 
    
    y_gFnc = data[0];
    for i in range(N0,N1):
        for j in range(Nslice):
            y_gFnc[j] += data[i,j]
            
        if(draw_on_screen):
            print('condition progress: ' + str((i+1-N0)/Nfrm*100) + '%                     \r', end='')
    y_gFnc = [_x/(N1-N0) for _x in y_gFnc]
    
    fig_c = -1
    fig = []
    time_gaps_str += ('_Nslice_' + str(Nslice))
    y0 = [math.exp(-my.phi_r(_x)) for _x in x_dr]
    
    path = os.path.join(graph_dir, 'g(t)_' + time_gaps_str + '.png')
    fig_c, fig = my.plot_error(fig_c, fig, x_dr, y_gFnc, y0=1, y_th=y0, y_lbl='n/n0', pic_path=path, show_key=draw_on_screen)
    
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
