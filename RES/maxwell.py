import sys
import numpy as np
import math
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 2
    if len(args) < argc_min:
        print('usage: python maxwell.py    model_name    N_V_steps    [keys,    N0,    N1]')
        sys.exit(1)
        
    model_name, keys, model_dir, graph_dir, time_gaps_str, N0, N1, Nfrm, E, P, Tmp, Tmp_av, t, stabTind, params = my.std_start(args, 0, 2, 3, 4)
    # std_start(args, model_i, N0_i, N1_i):
    # model_name, keys, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Tmp, Tmp_av, t, stabTind, params
    
    NV = int(args[1])
    draw_on_screen = my.find_key(keys, 'keep')
    if(not draw_on_screen):
        plt.switch_backend('Agg')
        
    if(len(args) < 4):
        if(draw_on_screen):
            print("N0 isn't set, N0 = ", stabTind, ' will be used')
        N0 = stabTind
        Nfrm = N1-N0
    
    Z_sum = math.pow(1/2/np.pi/Tmp_av, 1.5);
    k0 = 2*np.pi*Z_sum
    v_prob = math.sqrt(2*Tmp_av)
    v2_prob = Tmp_av
    total_v2 = np.zeros(NV-1)
    total_v = np.zeros(NV-1)
    for n in range(N0,N1):
        x,v = my.read_frame(model_name, n)
        m = np.ones(params['Ntot'])
        
        v2 = my.arrDot(v, v)
        v_abs = my.arrFnc(v2, math.sqrt)
        if(n == N0):
            #v2max = max(v2)
            #v2min = min(v2)
            v2max = 3*Tmp_av*4
            v2min = 3*Tmp_av/10
            k_v2 = 1/(params['Ntot']*(v2max-v2min)/NV*Nfrm)
            k_v = 1/(params['Ntot']*(math.sqrt(v2max) - math.sqrt(v2min))/NV*Nfrm)
            v2steps = np.linspace(v2min, v2max, num = NV)
            vsteps = np.linspace(math.sqrt(v2min), math.sqrt(v2max), num = NV)
            x_v2 = [(v2steps[i] + v2steps[i+1])/2 for i in range(NV-1)]
            x_v = [(vsteps[i] + vsteps[i+1])/2 for i in range(NV-1)]            
            y0_v2 = [k0*math.sqrt(_x)*math.exp(-_x/2/Tmp_av) for _x in x_v2]
            y0_v = [2*k0*_x*_x*math.exp(-_x*_x/2/Tmp_av) for _x in x_v]
            #y0_ln = [-*_x/2/Tmp_av for _x in x_v2]
            y0_ln = []
            for _x in x_v2:
                _a = _x/2/Tmp_av
                #y0_ln.append(math.log(1-2*_a/3/params['Ntot']) - _a)
                y0_ln.append(-_a)
            
        peaks, bin_edges = np.histogram(v2, bins=v2steps)        
        for i in range(NV-1):
            total_v2[i] += peaks[i]
        
        peaks, bin_edges = np.histogram(v_abs, bins=vsteps)        
        for i in range(NV-1):
            total_v[i] += peaks[i]        
        
        if(draw_on_screen):
            print('maxwell progress: ' + str((n+1-N0)/(N1-N0)*100) + '%                     \r', end='')

    y_v2 = [total_v2[i]*k_v2 for i in range(NV-1)]
    y_v = [total_v[i]*k_v for i in range(NV-1)]
    y_ln = []
    for i in range(NV-1):
        if(y_v2[i] > 0):
            y_ln.append(math.log(y_v2[i]/k0/math.sqrt(x_v2[i])))
        else:
            y_ln.append(0)
    dy_v = [y_v[i]/y0_v[i]-1 for i in range(len(x_v))]
    dy_v2 = [y_v2[i]/y0_v2[i]-1 for i in range(len(x_v2))]
    dy_ln = [y_ln[i]/y0_ln[i]-1 for i in range(len(x_v2))]
    
    p_lnp = np.poly1d(np.polyfit(x_v2, y_ln, 1));
    k_exp = -p_lnp.c[0]
    k_th = 1/2/Tmp_av
    #print('k_exp = ', k_exp)
    #print('k_th = ', k_th)
    s = 0
    for i in range(NV-1):
        s += (p_lnp(x_v2[i]) - y_ln[i])**2
    s /= (NV-2)
    if(draw_on_screen):
        print('S = ', s)
        
    fig_c = -1
    figs = []
    time_gaps_str += ('_NV_' + str(NV))    
    
    path = os.path.join(graph_dir, 'p(v)_' + time_gaps_str + '.png')
    fig_c, figs = my.plot_error(fig_c, figs, x_v, y_v, y_th=y0_v, x0=v_prob, x_lbl='v', y_lbl='p', pic_path=path, draw_linfit='n', show_key=draw_on_screen)
    
    path = os.path.join(graph_dir, 'p(v2)_' + time_gaps_str + '.png')
    fig_c, figs = my.plot_error(fig_c, figs, x_v2, y_v2, y_th=y0_v2, x0=v2_prob, x_lbl='v^2', y_lbl='p', pic_path=path, draw_linfit='n', show_key=draw_on_screen)

    path = os.path.join(graph_dir, 'ln_p(v2)_' + time_gaps_str + '.png')
    #fig_c, figs = my.plot_error(fig_c, figs, x_v2, y_ln, y_th=y0_ln, x_lbl='v^2', y_lbl='ln(p/v^2)', tit='(ln(p/v^2))(v^2) | k_th = ' + my.str_sgn_round(k_th,3) + ',  k_exp = ' + my.str_sgn_round(k_exp,3) + ' | std = ' + my.str_sgn_round(s,3), pic_path=path, show_key=draw_on_screen)
    fig_c, figs = my.plot_error(fig_c, figs, x_v2, y_ln, y_th=y0_ln, x_lbl='v^2', y_lbl='ln(p/v^2)', tit='k_th = ' + my.str_sgn_round(k_th,3) + ',  k_exp = ' + my.str_sgn_round(k_exp,3) + ' | std = ' + my.str_sgn_round(s,3) + ', k_err = ' + my.str_sgn_round(math.exp(abs(math.log(k_exp/k_th))) - 1, 3) + ', b = ' + my.str_sgn_round(p_lnp.c[1],3), pic_path=path, show_key=draw_on_screen)
    
    path = os.path.join(graph_dir, 'd_p(v)_' + time_gaps_str + '.png')
    fig_c, figs = my.plot_error(fig_c, figs, x_v, dy_v, x0=v_prob, y0=0, x_lbl='v', y_lbl='p/p0-1', pic_path=path, show_key=draw_on_screen)
    
    path = os.path.join(graph_dir, 'd_p(v2)_' + time_gaps_str + '.png')
    fig_c, figs = my.plot_error(fig_c, figs, x_v2, dy_v2, x0=v2_prob, y0=0, x_lbl='v^2', y_lbl='p/p0-1', pic_path=path, show_key=draw_on_screen)
    
    path = os.path.join(graph_dir, 'd_ln_p(v2)_' + time_gaps_str + '.png')
    fig_c, figs = my.plot_error(fig_c, figs, x_v2, dy_ln, y0=0, x_lbl='v^2', y_lbl='dln(p)', pic_path=path, show_key=draw_on_screen)
    
    if(draw_on_screen):
        input()
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
