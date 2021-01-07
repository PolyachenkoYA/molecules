import sys
import numpy as np
import math as math
import os
import matplotlib.pyplot as plt
import re
import matplotlib
import scipy
import scipy.signal
import pathlib
import shutil


RES_SUBPATH = 'graph'
FRAMES_SUBPATH = 'frames'
RAW_DATA_PATH = 'DATA'
E_filename = 'E.txt'
P_filename = 'Pressure.txt'
Time_filename = 'Time.txt'
diff_filename = 'diffusion.dat'
param_file_suff = 'param.dat'
particles_file_suff = 'particles.xyz'
complog_filename = 'comp.log'

head_filename = 'head.txt'
params_names = ['Ntot', 'n', 'n_cr', 'r_cut', 'endT', 'dumpDT', 'dt', 'Tmp', 'dissipK', 'TmpStabEps', 'TmpStabGap', 'compMode', 'binOutF', 'Nthreads']
all_proc_flags = '-energy-maxwell-diff-cond-pics-percent-r_cm-r2-v_cm-r2_err-subdiff'
error_str = '!!ERROR OCCURED!!'
myeps = np.finfo(float).eps * 10

def phi_r(r):
    x = pow(r, -6)
    return 4 * x * (x - 1)

def totalMin(arr):
    while(type(arr) is list):
        arr = min(arr)
    return arr

def totalMax(arr):
    while(type(arr) is list):
        arr = max(arr)
    return arr

def sum(v, n1=0, n2=-1, s0=0):
    if(n2 == -1):
        n2 = len(v)
        
    s = s0
    for i in range(n1,n2):
        s += v[i]
    return s

def dot(v1, v2):
    #if(len(v1) != len(v2)):
    #    return NAN
    #return sum([ v1[i]*v2[i] for i in range(0,len(v1))])
    #return sum([ v1[i]*v2[i] for i in range(3)])
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    
def length(v):
    return math.sqrt(dot(v,v))
    
def arrDot(Varr1, Varr2):
    return [dot(Varr1[i], Varr2[i]) for i in range(min(len(Varr1), len(Varr2)))]
    
def arrFnc(v, fnc):
    return [fnc(_x) for _x in v]
    
def load_file(model_name, file_name):
    # numbers only. output is a matrix
    return np.loadtxt(os.path.join('./', RAW_DATA_PATH, model_name, file_name))
    
def read_file(model_name, file_name):
    # just text. output is lines of text
    #return (open(os.path.join('./' + model_name, file_name), 'rU').read()).split('\n')
    return read_random_file(os.path.join('./', RAW_DATA_PATH, model_name, file_name))
    
def read_random_file(file_name):
    return (open(file_name, 'rU').read()).split('\n')

def read_params(filename):
    #params = read_file(model_name, 'param.dat')
    params = read_random_file(filename)
    params = params[0:2]
    kys = params[0].split()
    vals = params[1].split()
    
    if(len(kys) != len(vals)):
        print('wrong params file:\n' + filename)
        sys.exit(1)
    
    params = {}
    for i in range(len(kys)):
        params[kys[i]] = float(vals[i])
    
    int_params = ['Ntot', 'compMode', 'TmpStabGap', 'binOutF', 'Nthreads']
    for i in int_params:
        params[i] = int(params[i])
    params['Nframes'] = int(params['endT']*params['dumpDT']) + 1
    params['R'] = ((params['Ntot'] / params['n'])**(1.0 / 3)) / 2.0
        
    return params
    
def save_params(filename, params):
    s = ''
    ch_per_pair = 14
    for i in range(len(params_names)):
        s += str(params_names[i]).rjust(ch_per_pair)
    s += '\n'
    for i in range(len(params_names)):
        s += str(params[params_names[i]]).rjust(ch_per_pair)
    s += '\n'
                
    np.savetxt(filename, [s], fmt='%s')
    
def load_E(model_name):
    return load_file(model_name, E_filename)
    
def load_pressure(model_name):
    return load_file(model_name, P_filename)    
    
def load_Time(model_name):
    return load_file(model_name, Time_filename)
    
def read_frame(model_name, index):
    #data = ovito.io.import_file(os.path.join('./' + model_name, str(index) + '.xyz'), columns = ["Position.X", "Position.Y", "Position.Z", "vx", "vy", "vz", "m"])
    
    #data = read_file(model_name, str(index) + '.xyz') # NO_FRAMES_VERSION
    data = read_file(model_name, os.path.join(FRAMES_SUBPATH, str(index) + '.xyz'))
    N = int(data[0])
    data = data[2:]
    
    x = []
    v = []
    for i in range(N):
        particle = data[i].split()
        if(particle):
            x.append([float(particle[0]), float(particle[1]), float(particle[2])])
            v.append([float(particle[3]), float(particle[4]), float(particle[5])])
        
    return x,v

def norm_check(v, v0, k):
    if(v/v0 > k):
        return v0*k
    elif(v/v0 < 1/k):
        return v0/k
    else:
        return v
    
def shiftEl(x, R):
    if(x > R):
        x -= 2*R
    elif(x < -R):
        x += 2*R
    return x

def shiftCrd(x, R):
    x[0] = shiftEl(x[0], R)
    x[1] = shiftEl(x[1], R)
    x[2] = shiftEl(x[2], R)
    return x
    
def printCrd(x):
    if(len(x) != 3):
        print('wrong crd for print')
    else:
        print('(' + str(x[0]) + ';' + str(x[1]) + ';' + str(x[2]) + ')', end='')
        
def findStabTind(T, Tav):
    k = 1
    l = len(T)
        
    while((T[0]<Tav) == (T[k]<Tav)):
        k += 1
        if(k >= l):
            break
        
    return k
    
def read_head(model_name):
    return load_file(model_name, head_filename)
        
def std_start(args, model_i, keys_i, N0_i, N1_i):
    model_name = args[model_i]
    params = read_params(os.path.join(RAW_DATA_PATH, model_name, param_file_suff))
    t = load_Time(model_name) 
    hd = read_head(model_name)

    argc = len(args)
    if(argc <= N0_i):
        N0 = 1
    else:
        N0 = int(args[N0_i])
        
    if(argc <= N1_i):
        #N1 = params['Nframes']
        N1 = int(hd[0])
    else:
        N1 = int(args[N1_i])
    
    if(N1 < N0):
        print('N0 = ', N0, '; N1 = ', N1)
        exit(1)
        
    if(argc <= keys_i):        
        keys = '-'.split('-') 
    else:
        if(args[keys_i] == '-all'):
            args[keys_i] = all_proc_flags
        keys = (args[keys_i]).split('-')
        
    E = load_E(model_name)
    Pressure = load_pressure(model_name)
    Tmp = [E[i,0]*2/3/params['Ntot'] for i in range(N0, N1)]    
    Tmp_av = sum(Tmp)/(N1 - N0)
    stabTind = findStabTind(Tmp,Tmp_av)
    Tmp_av = sum(Tmp[stabTind:])/(N1 - N0 - stabTind)
    
    if((N0 == 0) and (N1 == len(t))):
        time_gaps_str = 'whole'
    elif((N0 == 0) and (N1 < len(t))):
        time_gaps_str = '0_' + str_sgn_round(t[N1-1],3)
    elif((N0 > 0) and (N1 == len(t))):
        time_gaps_str = str_sgn_round(t[N0],3) + '_end'
    else:
        time_gaps_str = str_sgn_round(t[N0],3) + '_' + str_sgn_round(t[N1-1],3)
    model_dir = os.path.join(RAW_DATA_PATH, model_name)
    graph_dir = os.path.join(model_dir, RES_SUBPATH)
    if(not os.path.exists(graph_dir)):
        os.makedirs(graph_dir)
        
    Nfrm = N1-N0
    if(Nfrm <= 0):
        print('error: N1 <= N0')
        sys.exit(1)
    
    return model_name, keys, model_dir, graph_dir, time_gaps_str, N0, N1, Nfrm, E, Pressure, Tmp, Tmp_av, t, stabTind, params
    
def plot_error(fig_c, figs, x, y, y_th=[], x0=float('nan'), y0=float('nan'), x_lbl='time', y_lbl='error', tit = '', pic_path='', show_key=0, filter_rng=15, draw_linfit='y'):
    
    fig_c += 1
    figs.append(plt.figure(fig_c))
    plt.plot(x, y, '-')
    if(y_th):
        plt.plot(x, y_th, '--')
    if(not math.isnan(y0)):
        plt.plot([min(x), max(x)], [y0, y0], '--')
    if(not math.isnan(x0)):
        plt.plot([x0, x0], [min(y + y_th), max(y + y_th)], '--')
    
    if(draw_linfit == 'y'):
        err_fit = np.poly1d(np.polyfit(x, y, 1))
        plt.plot([x[0], x[-1]], [err_fit(x[0]), err_fit(x[-1])], '--')
    
    if(not math.isnan(filter_rng)):
        y_filt = scipy.signal.savgol_filter(y, int(round(filter_rng))*2 + 1, 2)
        plt.plot(x, y_filt, '--') 
        
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))        
    plt.xlabel(x_lbl)
    plt.ylabel(y_lbl)
    plt.grid(True)
    if(not tit):
        tit = '(' + y_lbl + ')(' + x_lbl + ')'
    plt.title(tit)
    if(show_key):
        figs[fig_c].show()
    
    if(pic_path):
        figs[fig_c].savefig(pic_path)
    
    return fig_c, figs
    
def find_key(keys, key0):
    for _key in keys:
        if(_key == key0):
            return 1
    
    return 0
    
def find_cm(x ,m):
    l = len(m)
    m0 = sum(m)
    xc = [0,0,0]
    for i in range(3):
        xc[i] = sum([x[j][i]*m[j] for j in range(l)]) / m0
    return xc


def run_it(command):
    print(command)
    ok = (os.system(command) == 0)
    if(not ok):
        print(error_str)
    return ok    
    
def str_sgn_round(x, n):
    return '%s' % float(('%.' + str(n) + 'g') % x)
    
def r_th1(T, n):
    #return pow((1 + math.sqrt(1 + 3/2 * T - 4 * n**2 * (1 - n**2))) / 2, -1/6)
    return pow((1 + math.sqrt(1 + 3/2 * T)) / 2, -1/6)
    
def r_th2(T):
    return math.sqrt(1 + 2 / (3 * T))
    
def lmd_1(a):
    return 1/4 * math.sqrt(3 * math.pi * a)
    
def lmd_2(k, m, T):
    return k/4 * math.sqrt(math.pi * m / (2 * T))
    
def lmd_3(k1, b1, k2, b2, m, T):
    return 3/4 * math.sqrt(math.pi * T / (2 * m)) * math.exp(- (b2 - b1) / (k2 - k1))
    
def lmd_4(k1, b1, k2, b2):
    return 1/4 * math.sqrt(3 * math.pi / (1 + math.exp(-2))) * math.exp((b1*k2 - k1*b2) / (2 * (k2 - k1)))
    
def rel_err(a, b):
    x_min = min(abs(a),abs(b))
    if(x_min < myeps):
        if(max(abs(a),abs(b)) < myeps):
            return 0
        else:
            return 1/myeps
    else:
        return abs(a - b)/x_min

def change_params_fnc(old_name, new_name, args=[]):
    params = read_params(old_name)
    for i in range(int(round(len(args)/2))):
        params[args[2*i]] = args[2*i+1]
    save_params(new_name, params)

def full_cycle_fnc(input_params_name, model_name, keys, args):
    extra_args_str = ' '.join(args[2:])
        
    if(find_key(keys, 'gen')):
        if(not run_it('./gen ' + input_params_name + ' ' + model_name)):
            return

    if(not run_it('./comp ' + model_name + '            ')):
        return

    if(find_key(keys, 'gen')):
        src = model_name + '_gen.log'
        dst = os.path.join('./', model_name, 'gen.log')
        cmd = 'mv ' + src + ' ' + dst +  '            '
        print(cmd)
        shutil.move(src, dst)
        #run_it(cmd)
        
    if(find_key(keys, 'cond') or find_key(keys, 'condition')):
        if(not run_it('./post_proc ' + model_name + ' 3 150' + '                                 ')):
            return
    
    if(find_key(keys, 'pics') or find_key(keys, 'keep') or find_key(keys, 'move_res')):
        command_to_run = 'cd ../../' + '            '
        print(command_to_run)
        os.chdir('../../')
        
        src = os.path.join('./', '!COMP', '!go', model_name)
        dst = os.path.join('./', 'RES', 'DATA')
        cmd = 'mv ' + src + ' ' + dst + '            '
        print(cmd)
        shutil.move(src, dst)
        #run_it(cmd)
        
        command_to_run = 'cd ./RES' + '            '
        print(command_to_run)
        os.chdir('./RES')   
        
        run_it('python full_post_proc.py ' + model_name + ' ' + extra_args_str)



