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
import copy as cp

import mylib_molecules as my

def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./diffusion.py model_name [keys, N0, N1, Nprt]')
        sys.exit(1)
        
    # --------------------------------- std start ------------------------------------------
    model_name, keys, model_dir, graph_dir, time_gaps_str, N0, N1, Nfrm, E, P, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 1, 2, 3)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
            
    # --------------------------------- input handle ------------------------------------------
    if(Nfrm < 10):
        print('Too few grames (Nfrm < 10)')
        return
    if(len(args) < 5):
        Nprt = params['Ntot']
    else:
        Nprt = int(args[4])
                
    draw_on_screen = my.find_key(keys, 'keep')
    pnrt_progress = (my.find_key(keys, 'percent') or my.find_key(keys, 'prcnt'))
    if(not draw_on_screen):
        plt.switch_backend('Agg') 
    make_pics = (my.find_key(keys, 'pics') or draw_on_screen)
    
    # --------------------------------- theory prediction ------------------------------------------
    n0 = params['n']
    r_th0 = 1
    lmd_th0 = 1/(n0 * math.pi * r_th0**2 * math.sqrt(2))
            
    r_th1 = r_th0 * my.r_th1(Tmp_av, n0)
    lmd_th1 = lmd_th0 * (r_th0/r_th1)**2
    
    #r_th2 = pow(2, 1/6) * math.sqrt(1 + 2/(3*Tmp_av))
    r_th2 = r_th0 * my.r_th2(Tmp_av)
    lmd_th2 = lmd_th0 * (r_th0/r_th2)**2
    
    k_th = 4*lmd_th1 * math.sqrt(2*Tmp_av / math.pi)
    c_th = 1 / (6*Tmp_av)
    
    # --------------------------------- experimant data ------------------------------------------
    diff_file_path = os.path.join(model_dir, my.diff_filename)
    if(os.path.isfile(diff_file_path) and (not my.find_key(keys, 'recomp'))):
        data = np.loadtxt(diff_file_path)
        #data = np.delete(data, (0), axis=0)
        #data = np.delete(data, (0), axis=1)
        r2 = cp.deepcopy(data[0])
        l = len(r2);
        x_cm = np.zeros((l, 3))
        v_cm = np.zeros((l, 3))
        for i in range(l):
            for j in range(3):
                x_cm[i][j] = data[j+1][i]
                v_cm[i][j] = data[j+4][i]
        data = None;
    else:
        r2 = []
        dr = [0,0,0]
        cubeSize = 2*params['R']
        cubeSzEff = params['R']
        x_cm = []
        #r_cm = []
        v_cm = []
        for n in range(N0,N1):
            x,v = my.read_frame(model_name, n)
            m = np.ones(params['Ntot'])
            
            if(n == N0):
                # set start conditions
                x0 = cp.deepcopy(x)
                x_real = cp.deepcopy(x)
                r2.append(0)
                
                v_start_av = np.mean([my.length(_x) for _x in v])
            else:
                r2_sum = 0
                for i in range(Nprt): # find <r^2> considering coords's shifts
                    for j in range(3): 
                        dr[j] = x[i][j]-x0[i][j]
                    r2_sum += my.dot(dr,dr)
                r2.append(r2_sum/Nprt)
            x_cm.append(my.find_cm(x, m))
            v_cm.append(my.find_cm(v,m))
                                
            if(pnrt_progress):
                print('diffusion progress: ' + str((n + 1 - N0) / Nfrm * 100) + '%                     \r', end='')
                
        xx = [_x[0] for _x in x_cm]
        xy = [_x[1] for _x in x_cm]
        xz = [_x[2] for _x in x_cm]
        vx = [_x[0] for _x in v_cm]
        vy = [_x[1] for _x in v_cm]
        vz = [_x[2] for _x in v_cm]        
        np.savetxt(diff_file_path, (r2, xx, xy, xz, vx, vy, vz))
                
    # --------------------------------- experimant analisys ------------------------------------------        
    # find proper tg
    r1 = my.arrFnc(r2, math.sqrt)
    big_time = 3
    small_time = 0.5
    approx_r2 = np.poly1d(np.polyfit(t[N0:N1], r2, 1))
    k_exp = approx_r2.c[0]
    tau = c_th*k_exp
    c2 = int(tau*params['dumpDT'] * big_time)
    if(c2 > Nfrm-10):
        c2 = Nfrm-10
    if(c2 < 10):
        c2 = int(Nfrm/2)
    approx_r2 = np.poly1d(np.polyfit(t[(N0+c2):N1], r2[c2:], 1))
    k_exp = abs(approx_r2.c[0])
    tau = c_th * k_exp
    c2 = int(tau*params['dumpDT'] * big_time)
    if(c2 > Nfrm - 5):
        c2 = Nfrm - 5
    if(c2 < 10):
        c2 = int(Nfrm/2)
    c1 = int(tau*params['dumpDT'] * small_time)
    if(c1 < 3):
        c1 = 3
                
    # in case of strange emergencies
    bad_k_exp = (k_exp < my.myeps)
    if(bad_k_exp):
        if(k_exp < 0):
            k_exp = -k_exp
        if(k_exp < my.myeps):
            k_exp = my.myeps
            
    xa0 = abs(approx_r2(0))
    tau = xa0/k_exp
            
    c_exp = xa0 / k_exp**2
    B = (1 - r2[-2] / (k_exp * t[N1-1])) * math.sqrt(t[N1-1])
    subdiff_c = int(round(B*B*params['dumpDT'] + 1.01))
        
    x = my.arrFnc(t[(N0+c2):N1], math.log)
    y = my.arrFnc(r2[c2:], math.log)
    logEnd_approx = np.poly1d(np.polyfit(x, y, 1))
    x = my.arrFnc(t[N0:c1], math.log)
    y = my.arrFnc(r2[1:c1], math.log)
    logBeg_approx = np.poly1d(np.polyfit(x, y, 1))
    x = t[N0:c1]
    y = my.arrFnc(r2[1:c1], math.sqrt)
    approx_r1 = np.poly1d(np.polyfit(x, y, 1))
    
    #lmd_1 = math.sqrt(3*math.pi*xa0)/4
    lmd_1 = my.lmd_1(xa0)
    #lmd_2 = k_exp/4 * math.sqrt(math.pi / (Tmp_av*2))
    lmd_2 = my.lmd_2(k_exp, 1, Tmp_av)
    D_Einst = lmd_2 * math.sqrt(8 * Tmp_av / math.pi) / 3
    
    selfconsist = my.rel_err(lmd_1, lmd_2)
    
    lmd_3 = my.lmd_3(logBeg_approx.c[0], logBeg_approx.c[1], logEnd_approx.c[0], logEnd_approx.c[1], 1, Tmp_av) # tau -> lmd
    lmd_4 = my.lmd_4(logBeg_approx.c[0], logBeg_approx.c[1], logEnd_approx.c[0], logEnd_approx.c[1]) # straight lmd

    if(make_pics):
        fig_c = -1
        figs = []
        
        if(my.find_key(keys, 'r_cm')):
            # -------------------------------- r_cm(t) -----------------------------------------        
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
                            
        if(my.find_key(keys, 'r2')):
            # -------------------------------- r2(t) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            
            #tau = 1
            x = t[N0:N1]
            x_draw = [_x/tau for _x in x]
            plt.plot(x_draw, r2, '-', label="exp") # experiment data
            
            x = t[N0:N1]
            r2_th = [k_exp*tau*(_x/tau - 1 + math.exp(- _x/tau)) for _x in x]
            x_draw = [_x/tau for _x in x]
            plt.plot(x_draw, r2_th, '--', label="my_th")  # my approximation            
            
            if(my.find_key(keys, 'subdiff')):
                x = t[(N0 + subdiff_c):N1]
                subdiff_app = [k_exp*_x*(1 - B/math.sqrt(_x)) for _x in x]
                x_draw = [_x/tau for _x in x]
                plt.plot(x_draw, subdiff_app, '-.', label="subdiff") # subdiff approximation
                                    
            x = [t[N0], t[N1-1]]
            y = [approx_r2(_x) for _x in x]
            x_draw = [_x/tau for _x in x]
            plt.plot(x_draw, y, '--', label="line asimp") # line approximation            
            
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.xlabel('t/tau')
            plt.ylabel('<r^2>')
            plt.grid(True)
            #plt.title('r^2(t) | k_exp = ' + str(my.str_sgn_round(k_exp,3)) + ', k_th = ' + str(my.str_sgn_round(k_th,3)))
            plt.title('r^2(t) | tau = ' + my.str_sgn_round(tau,3) + ' | B = ' + my.str_sgn_round(B,3))
            plt.legend()
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r2(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)
            
            # -------------------------------- log(r2)(log(t)) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            
            x = t[N0:N1]
            x_draw = [_x/tau for _x in x]
            plt.plot(x_draw, r2, '-', label="exp") # experiment data
                        
            x = [_x/tau for _x in t[N0:N1]]
            plt.plot(x, r2_th, '--', label="my_th") # my approximation         
            
            if(my.find_key(keys, 'subdiff')):
                x = t[(N0 + subdiff_c):N1]
                x_draw = [_x/tau for _x in x]            
                plt.plot(x_draw, subdiff_app, '--', label="subdiff") # subdiff approximation
            
            x = [t[N0], t[N0+c2]]
            y = [math.exp(logBeg_approx(math.log(_x))) for _x in x]
            x_draw = [_x/tau for _x in x]
            plt.plot(x_draw, y, '--', label="ballistic") # beginning approximation        
            
            x = [t[N0+c1], t[N1-1]]
            y = [math.exp(logEnd_approx(math.log(_x))) for _x in x]
            x_draw = [_x/tau for _x in x]            
            plt.plot(x_draw, y, '--', label="linear") # ending approximation            
                    
            plt.xlabel('t/tau')
            plt.ylabel('<r^2>')
            plt.grid(True)
            plt.xscale('log')
            plt.yscale('log')
            plt.title('log(r2) | k_begin = ' + str(my.str_sgn_round(logBeg_approx.c[0],3)) + '; k_end = ' + str(my.str_sgn_round(logEnd_approx.c[0],3)))
            plt.legend()
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r2(t)_loglog_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)
            
            # -------------------------------- r(t) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))        
            plt.plot(t[N0:N1], r1, '-', label="exp") # experiment data
                    
            x = t[N0:N1]
            y = my.arrFnc(r2_th, math.sqrt)
            plt.plot(x, y, '--', label="my_th") # my approximation 
            
            if(my.find_key(keys, 'subdiff')):
                x = t[(N0 + subdiff_c):N1]
                y = my.arrFnc(subdiff_app, math.sqrt)
                plt.plot(x, y, '--', label="subdiff") # sudiff approximation             
            
            x = [t[N0], t[c2]]
            y = [approx_r1(_x) for _x in x]
            plt.plot(x, y, '--', label="ballistic") # beginning approximation             
                
            plt.xlabel('t')
            plt.ylabel('sqrt(<r^2>)')
            plt.grid(True)
            #plt.title('sqrt(r2) | l_th = ' + str(my.str_sgn_round(lmd_th1,3)) + '; l_3 = ' + str(my.str_sgn_round(lmd_3,3)) + '; l_2 = ' + str(my.str_sgn_round(lmd_2,3)))
            plt.title('sqrt(<r^2>)')
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r1(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)    
                            
        if(my.find_key(keys, 'r2_err')):
            # -------------------------------- r2_err(t) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            
            x = t[N0:N1]
            y = [r2[_i] - r2_th[_i] for _i in range(len(r2_th))]
            plt.plot(x, y, '-') # my approximation            

            x = t[(N0 + subdiff_c):N1]
            y = [r2[_i + subdiff_c] - subdiff_app[_i] for _i in range(0, N1 - N0 - subdiff_c)]
            plt.plot(x, y, '--') # subdiff approximation                        
            
            plt.plot([t[N0+1], t[N1-1]], [0, 0], '--')
            
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.xlabel('t')
            plt.ylabel('<r^2> - <r^2>_th')
            plt.grid(True)
            #plt.title('r^2(t) | k_exp = ' + str(my.str_sgn_round(k_exp,3)) + ', k_th = ' + str(my.str_sgn_round(k_th,3)))
            plt.title('<r^2> error')
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r2_err(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)        
            
            # -------------------------------- r2_err_rel(t) ------------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            
            tmp = int(round(tau*5*params['dumpDT']))
            x = t[(N0 + subdiff_c + tmp):N1]
            y = [my.rel_err(r2[_i], r2_th[_i]) for _i in range(tmp, N1 - N0 - subdiff_c)]
            plt.plot(x, y, '-') # my approximation
            
            x = t[(N0 + subdiff_c + tmp):N1]
            y = [my.rel_err(r2[_i + subdiff_c],subdiff_app[_i]) for _i in range(tmp, N1 - N0 - subdiff_c)]
            plt.plot(x, y, '-.') # subdiff approximation                        
            
            plt.plot([t[N0+1], t[N1-1]], [0, 0], '--')
            
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.xlabel('t')
            plt.ylabel('err')
            plt.grid(True)
            #plt.title('r^2(t) | k_exp = ' + str(my.str_sgn_round(k_exp,3)) + ', k_th = ' + str(my.str_sgn_round(k_th,3)))
            plt.title('|<r^2> / <r^2>_th - 1|')
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'r2_err_rel(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)                
    
            
        if(my.find_key(keys, 'v_cm')):
            # -------------------------------- v_cm(t) -----------------------------------------
            fig_c += 1
            figs.append(plt.figure(fig_c))
            x = t[N0:N1]
            dv_cm = []
            std_dvcm = []
            tits = ['vx', 'vy', 'vz']        
            for _j in range(3):
                v_av = np.mean(v_cm[:][_j])
                dv_cm.append([v_cm[_i][_j] - v_av for _i in range(Nfrm)])
                std_dvcm.append(np.std(dv_cm[_j]))
                plt.plot(x, dv_cm[_j], '-', label = tits[_j] + '; std = ' + str(my.str_sgn_round(std_dvcm[_j],3))) # experiment data
            x = [t[N0], t[N1-1]]
            plt.plot(x, [0,0], '--')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.xlabel('t')
            plt.ylabel('v_cm')
            plt.grid(True)
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    ncol=3, mode="expand", borderaxespad=0.)        
            if(draw_on_screen):
                figs[fig_c].show()
            path = os.path.join(graph_dir, 'v_cm(t)_' + time_gaps_str + '.png')
            figs[fig_c].savefig(path)    
            
        if(my.find_key(keys, 'human')):
            print('n = ', n0, '; T = ', Tmp_av)
            print('lmd_0 = ', lmd_th0, '; lmd_T<1 = ', lmd_th1, '; lmd_T>1 = ', lmd_th2)
            print('lmd_a = ', lmd_1, '; lmd_tg = ', lmd_2, '; lmd_exp_intercept = ', lmd_3, '; lmd_tau_intercept = ', lmd_4)
            print('inconsistency = ', selfconsist)
        else:
            print(n0, Tmp_av, lmd_th0, lmd_th1, lmd_th2, lmd_1, lmd_2, lmd_3, lmd_4)
        #print('D_Einst = ' + str(D_Einst))
        
    if(draw_on_screen):
        #print(c_th, c_exp, my.rel_err(c_th, c_exp))
        input()

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
