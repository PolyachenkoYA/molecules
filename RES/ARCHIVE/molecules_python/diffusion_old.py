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
import copy as cp

import mylib_molecules as my

def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./diffusion.py model_name [keys, N0, N1, Nprt]')
        sys.exit(1)
        
    model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 1, 2, 3)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
            
    if(len(args) < 4):
        Nprt = params['Ntot']
    else:
        Nprt = int(args[3])
        
    if(len(args) == 1):
        full_mode = 1
    else:
        full_mode = (my.find_key(keys, 'full_mode') or my.find_key(keys, 'full'))
    if(full_mode and (Nfrm < int(Nfrm/2)+5)):
        print('too few frames for full_mode')
        sys.exit(1)
        
    draw_on_screen = my.find_key(keys, 'keep')
    if(not draw_on_screen):
        plt.switch_backend('Agg') 
    make_pics = (my.find_key(keys, 'pics') or draw_on_screen)
    
    
    n0 = params['Ntot'] / math.pow(2*params['R'], 3)
    sgm_th0 = 1
    lmd_th0 = 1/(n0 * math.pi * math.sqrt(2))
    
    sgm_th2 = pow(2, 1/6) * math.sqrt(1 + 2/(3*Tmp_av))
    lmd_th2 = lmd_th0 / (sgm_th2*sgm_th2)
        
    sgm_th1 = math.pow((1+math.sqrt(1+1.5*Tmp_av + ((2*n0)**2)*(n0**2-1)))/2, -1/6)
    lmd_th1 = lmd_th0 / (sgm_th1*sgm_th1)
    k_th = 4*lmd_th1 * math.sqrt(2*Tmp_av / (math.pi*params['mu']))
    
    r2 = []
    dr = [0,0,0]
    cubeSize = 2*params['R']
    cubeSzEff = params['R']
    x_cm = []
    r_cm = []
    for n in range(N0,N1):
        x,v,m = my.read_frame(model_name, n)
        
        if(n == N0):
            # set start conditions
            x0 = cp.deepcopy(x)
            x_real = cp.deepcopy(x)
            x_prev = cp.deepcopy(x)
            r2.append(0)
            
            v_start_av = np.mean([my.length(_x) for _x in v])
        else:        
            r2_sum = 0
            for i in range(Nprt): # find <r^2> considering coords's shifts
                for j in range(3): 
                    dx = x[i][j] - x_prev[i][j]
                    if(dx > cubeSzEff): # re-shift coords
                        dx -= cubeSize
                    elif(dx < -cubeSzEff):
                        dx += cubeSize
                    x_real[i][j] += dx
                    # find real displacement
                    dr[j] = x_real[i][j]-x0[i][j]
                r2_sum += my.dot(dr,dr)
            r2.append(r2_sum/Nprt)
            x_prev = cp.deepcopy(x)
        x_cm.append(my.find_cm(x_real,m))
        r_cm.append(my.length(x_cm[-1]))
                    
        if(draw_on_screen):
            print('diffusion progress: ' + str((n+1-N0)/Nfrm*100) + '%                     \r', end='')
            
    fig_c = -1
    figs = []

    # -------------------------------- r_cm(t) ------------------------------------------
    if(make_pics):
        fig_c += 1
        figs.append(plt.figure(fig_c))
        cm = []
        std_cm = []
        x = t[N0:N1]
        tits = ['x', 'y', 'z']
        for _j in range(3):
            cm.append([x_cm[_i][_j]-x_cm[0][_j] for _i in range(Nfrm)])
            std_cm.append(np.std(cm[_j]))
            plt.plot(x, cm[_j], '-', label = tits[_j] + '; std = ' + str(my.str_sgn_round(std_cm[_j],3))) # experiment data
        x = [t[N0], t[N1-1]]
        plt.plot(x, [0,0], '--')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel('t')
        plt.ylabel('xyz_cm')
        plt.grid(True)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=3, mode="expand", borderaxespad=0.)        
        if(draw_on_screen):
            figs[fig_c].show()
        path = os.path.join(graph_dir, 'xyz_cm(t)_' + time_gaps_str + '.png')
        figs[fig_c].savefig(path)    
            
    c2 = int(Nfrm/2)
    if(full_mode):
        c1 = int(Nfrm/4)        
        c3 = int(Nfrm/100)
        c4 = int(Nfrm/8)
        x = my.arrFnc(t[(N0+c2):(N1-c4)], math.log)
        y = my.arrFnc(r2[c2:-c4], math.log)
        logEnd_approx = np.poly1d(np.polyfit(x, y, 1))
        x = my.arrFnc(t[N0+1:c3], math.log)
        y = my.arrFnc(r2[1:c3], math.log)
        logBeg_approx = np.poly1d(np.polyfit(x, y, 1))
        x = t[N0+1:c3]
        y = my.arrFnc(r2[1:c3], math.sqrt)
        approx_r1 = np.poly1d(np.polyfit(x, y, 1))
    else:
        c1 = 0
    approx_r2 = np.poly1d(np.polyfit(t[(N0+c1):N1], r2[c1:], 1))
    k_exp = approx_r2.c[0]
    bad_k_exp = (k_exp < my.myeps)
    if(bad_k_exp):
        if(k_exp < 0):
            k_exp = -k_exp
        if(k_exp < my.myeps):
            k_exp = my.myeps
            
    # -------------------------------- r2(t) ------------------------------------------
    if(make_pics):
        fig_c += 1
        figs.append(plt.figure(fig_c))
        plt.plot(t[N0:N1], r2, '-') # experiment data
        x = [t[N0], t[N1-1]]
        plt.plot(x, [approx_r2(_x) for _x in x], '--') # line experiment approximation
        
        tau = abs(approx_r2(0))/k_exp
        x = t[N0:N1]
        y = [k_exp*(_x - tau + tau*pow(math.e, - _x/tau)) for _x in x]
        #plt.plot(x, [k_exp*(_x - tau + tau*math.exp(-_x/tau)) for _x in x], '--') # non-line experiment approximation
        plt.plot(x, y, '--') # non-line experiment approximation
        
        #plt.plot([t[c2], t[c2]], [min(r2), max(r2)], '--') # last part boundary
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel('t')
        plt.ylabel('<r^2>')
        plt.grid(True)
        #plt.title('r^2(t) | k_exp = ' + str(my.str_sgn_round(k_exp,3)) + ', k_th = ' + str(my.str_sgn_round(k_th,3)))
        plt.title('r^2(t)')
        if(draw_on_screen):
            figs[fig_c].show()
        path = os.path.join(graph_dir, 'r2(t)_' + time_gaps_str + '.png')
        figs[fig_c].savefig(path)
        
    lmd_1 = k_exp/4 * math.sqrt(math.pi*params['mu'] / (Tmp_av*2))    
    sgm_1 = 1 / math.sqrt(n0 * math.pi * lmd_1 * math.sqrt(2))
    lmd_4 = math.sqrt(abs(approx_r2(0)) * math.pi * 3)/4
        
    if(draw_on_screen):
        print('Tmp_av = ', Tmp_av)
        print('sgm_th1 = ', sgm_th1)
        print('lmd_th1 = ', lmd_th1)
        print('sgm_1 = ', sgm_1)
        print('lmd_1 = ', lmd_1)
        print('lmd_4 = ', lmd_4)
    if(full_mode):        
        i_log = 2
        while(1):
            tg = math.log(r2[i_log+10]/r2[i_log]) / math.log(t[N0+i_log+10]/t[N0+i_log])
            if(abs(2-tg) > 0.2):
                break
            i_log += 1
            if(i_log >= N1-N0-10):
                break
        stab_time_log = t[i_log]
        
        i_r1 = 2
        r1 = my.arrFnc(r2, math.sqrt)
        while(1):
            tg = (r1[i_r1+10]-r1[i_r1]) / (t[N0+i_r1+10] - t[N0+i_r1]) / v_start_av
            if(abs(tg-1) > 0.18):
                break
            i_r1 += 1
            if(i_r1 >= N1-N0-10):
                break
        stab_time_r1 = t[i_r1]
        
        #lmd_2 = 2 * stab_time * math.sqrt(2*Tmp_av / (math.pi*params['mu'])) 
        # it's incorrect to calculate lmd this way here, because in the beggining V doesn't have maxwell distribution                
        lmd_2 = stab_time_log * v_start_av
        lmd_3 = stab_time_r1 * v_start_av
        sgm_2 = 1 / math.sqrt(n0 * math.pi * lmd_2 * math.sqrt(2))
        sgm_3 = 1 / math.sqrt(n0 * math.pi * lmd_3 * math.sqrt(2))
        stab_time_th = lmd_th1/2 * math.sqrt(math.pi*params['mu']/(2*Tmp_av))
        if(draw_on_screen):
            print('stab_time_th = ', stab_time_th)
            print('stab_time_log = ', stab_time_log)
            print('stab_time_r1 = ', stab_time_r1)
            print('sgm_2 = ', sgm_2)
            print('lmd_2 = ', lmd_2)
            print('r2 = ', math.sqrt(r2[i_log]))
            print('sgm_3 = ', sgm_3)
            print('lmd_3 = ', lmd_3)
            print('r1 = ', r1[i_r1])
        else:
            print(n0, Tmp_av, lmd_th0, lmd_th1, lmd_th2, lmd_1, lmd_2, lmd_3, lmd_4)
                
        if(make_pics):
            # -------------------------------- log(r2)(log(t)) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            plt.plot(t[N0:N1], r2, '-') # experiment data
            x = [t[N0+1], t[N1-1]]
            plt.plot(x, [math.exp(logBeg_approx(math.log(_x))) for _x in x], '--') # beginning approximation
            plt.plot(x, [math.exp(logEnd_approx(math.log(_x))) for _x in x], '--') # ending approximation
            #plt.plot(x, [lmd_th1**2, lmd_th1**2], '--')
            #plt.plot(x, [lmd_1**2, lmd_1**2], '--')
            #plt.plot([stab_time_th, stab_time_th], [min(r2), max(r2)], '--')
            #plt.plot([t[c2], t[c2]], [min(r2), max(r2)], '--')
                    
            plt.xlabel('t')
            plt.ylabel('<r^2>')
            plt.grid(True)
            plt.xscale('log')
            plt.yscale('log')
            plt.title('log(r2) | k_begin = ' + str(my.str_sgn_round(logBeg_approx.c[0],3)) + '; k_end = ' + str(my.str_sgn_round(logEnd_approx.c[0],3)))
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r2(t)_loglog_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)
            
            # -------------------------------- r(t) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))        
            plt.plot(t[N0:N1], r1, '-') # experiment data
            x = [t[N0+1], t[c4]]
            plt.plot(x, [approx_r1(_x) for _x in x], '--') # beginning approximation 
            #plt.plot(x, [lmd_th1, lmd_th1], '--')
            #plt.plot([stab_time_th, stab_time_th], [min(r1), max(r1)], '--')
            #plt.plot([stab_time_r1, stab_time_r1], [min(r1), max(r1)], '--')
            #plt.plot([t[c2], t[c2]], [min(r1), max(r1)], '--')
                
            plt.xlabel('t')
            plt.ylabel('sqrt(<r^2>)')
            plt.grid(True)
            #plt.title('sqrt(r2) | l_th = ' + str(my.str_sgn_round(lmd_th1,3)) + '; l_3 = ' + str(my.str_sgn_round(lmd_3,3)) + '; l_2 = ' + str(my.str_sgn_round(lmd_2,3)))
            plt.title('sqrt(<r^2>)')
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r1(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)

                
    if(draw_on_screen):
        input()


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
