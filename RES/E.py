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
from scipy import signal
import copy as cp

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./E.py    model_name    [keys   N0    N1]')
        sys.exit(1)
        
    model_name, keys, model_dir, graph_dir, time_gaps_str, N0, N1, Nfrm, E, P, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 1, 2, 3)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
    
    #for i in range(N0, N1):
    #    for j in range(4):
    #        E[i,j] /= params['Ntot']
    E /= params['Ntot']
        
    fig_c = -1
    fig = []
    draw_on_screen = my.find_key(keys, 'keep')
    if(not draw_on_screen):
        plt.switch_backend('Agg') 
    
    fig_c += 1
    # 0
    fig.append(plt.figure(fig_c))
    plt.plot(t[N0:N1], E[N0:N1,2], '-', label = 'Etot')
    plt.plot(t[N0:N1], E[N0:N1,0], '-', label = 'Ek')
    plt.plot(t[N0:N1], E[N0:N1,1], '-', label = 'Ep')
    plt.xlabel('time')
    plt.ylabel('E')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0.)
    if(draw_on_screen):
        fig[fig_c].show()
    path = os.path.join(graph_dir, 'Energy3_' + time_gaps_str + '.png')
    fig[fig_c].savefig(path)
    
    # 1
    E0 = np.mean(E[N0:N1,0]) + abs(np.mean(E[N0:N1,1]))
    path = os.path.join(graph_dir, 'dE_norm_' + time_gaps_str + '.png')
    E_av = np.mean(E[N0:N1,2])
    #y = [(e_el - E[N0,2])/E0 for e_el in E[N0:N1,2]]
    y = [(e_el/E_av - 1) for e_el in E[N0:N1,2]]
    fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y0=0, y_lbl='E/<E>-1', tit='E/<E> - 1 | std = ' + my.str_sgn_round(np.std(y),3), pic_path=path, show_key=draw_on_screen)
    
    # shadow error
    E0 = np.mean(E[N0:N1,0]) + abs(np.mean(E[N0:N1,1]))
    path = os.path.join(graph_dir, 'Hshd_norm_' + time_gaps_str + '.png')
    H_av = np.mean(E[N0:N1,3])
    y = [(e_el/H_av - 1) for e_el in E[N0:N1,3]]
    #print(E[N0:N1, :])
    fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y0=0, y_lbl='H/<H>-1', tit='H/<H> - 1 | std = ' + my.str_sgn_round(np.std(y), 3), pic_path=path, show_key=draw_on_screen)    
    
    # shadow 
    fig_c += 1
    fig.append(plt.figure(fig_c))
    plt.plot(t[N0:N1], E[N0:N1,2], '-', label = 'Etot')
    plt.plot(t[N0:N1], E[N0:N1,3], '-', label = 'd2H')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    if(draw_on_screen):
        fig[fig_c].show()
    path = os.path.join(graph_dir, 'Energy_H_' + time_gaps_str + '.png')
    fig[fig_c].savefig(path)
    
    # shadow term
    fig_c += 1
    fig.append(plt.figure(fig_c))
    plt.plot(t[N0:N1], E[N0:N1,3] - E[N0:N1,2], '-', label = '$d^2 H$')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.)
    if(draw_on_screen):
        fig[fig_c].show()
    path = os.path.join(graph_dir, 'd2H_' + time_gaps_str + '.png')
    fig[fig_c].savefig(path)    
    
    # Pressure
    path = os.path.join(graph_dir, 'Pressure_' + time_gaps_str + '.png')
    P_th = params['n'] * params['Tmp']
    y = P[N0:N1]
    fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y_th = [_x * P_th for _x in np.ones(len(y))], y0 = np.mean(P[N0:N1]), y_lbl='P', 
                                tit='P(time) | std_rel = ' + my.str_sgn_round(np.std(P[N0:N1]) / np.mean(P[N0:N1]),3), pic_path=path, show_key=draw_on_screen)
        
    # 2
    #path = os.path.join(graph_dir, 'Etot_' + time_gaps_str + '.png')
    #y = [_x*params['Ntot'] for _x in E[N0:N1, 2]]
    #fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y0=np.mean(y), y_lbl='Etot', tit='Etot | std = ' + my.str_sgn_round(np.std(y),3), pic_path=path, show_key=draw_on_screen)
       
    # 3
    #path = os.path.join(graph_dir, 'dE_normSQRT(N)_' + time_gaps_str + '.png')
    ##y = [(e_el - E[N0,2])*mth.sqrt(params['Ntot']) for e_el in E[N0:N1,2]]
    #y = [(e_el - E_av)*mth.sqrt(params['Ntot']) for e_el in E[N0:N1,2]]
    #fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y0=0, y_lbl='(E - <E>)/sqrt(N)', tit='Esqrt | std = ' + my.str_sgn_round(np.std(y),3), pic_path=path, show_key=draw_on_screen)
    
    # 4
    #path = os.path.join(graph_dir, 'dE_' + time_gaps_str + '.png')
    ##y = [(e_el - E[N0,2])*params['Ntot'] for e_el in E[N0:N1,2]]
    #y = [(e_el - E_av)*params['Ntot'] for e_el in E[N0:N1,2]]
    #fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], y, y0=0, y_lbl='E - <E>', tit='E - <E> | std = ' + my.str_sgn_round(np.std(y),3), pic_path=path, show_key=draw_on_screen)
       
    # 5
    path = os.path.join(graph_dir, 'Tmp_' + time_gaps_str + '.png')
    y = Tmp
    fig_c, fig = my.plot_error(fig_c, fig, t[N0:N1], Tmp, y0=Tmp_av, y_lbl='Tmp', tit='Tmp | std_rel = ' + my.str_sgn_round(np.std(Tmp)/Tmp_av,3), pic_path=path, show_key=draw_on_screen)
        
    if(draw_on_screen):
        input()
    
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
