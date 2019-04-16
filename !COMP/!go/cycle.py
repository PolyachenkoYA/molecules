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

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    argc = len(args)
    if(argc < argc_min):
        print('usage:\n./cycle.py start_name')
        sys.exit(1)
    
    base_param_name = args[0];
    params = my.read_params(base_param_name + '_' + my.param_file_suff)        
    
    dt_arr = [128,          161,          203,          255,          320,          403,          507,          638,          802,         1009,         1269,         1596,         2008,         2525,         3176,         3995,         5025,         6321,         7950,        10000]
    dt_arr = [128,          161,          203,          255,          320,          403,          507,          638,          802,         1009]
    dt_arr = [403,          507,          638,          802,         1009]
    dt_arr = [128,          161,          203,          255,          320]
    dt_arr = [128, 256, 512, 1024]
    
    N_arr = [100, 200]
    N_arr = [2800, 100,  139,  195,  271,  379,  528,  737, 1028, 1434, 2000]    
    N_arr = [512, 1024, 2048, 4096]
    n = 0.35
        
    rc = [3, 4]
    rc = [2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17, 20, 23, 26, 30, 35, 40, 45, 50]
    rc = [55, 60, 65, 70, 75, 80, 85, 50*math.sqrt(3), 100]
    rc = [0.1, 90]
    
    dt_arr = [128,    256,    512,    640,   1024,   1280,   2048,   2560,   3200,   4096,   5120,   6400,   8192,  10240,  12800,  16000,  16384,  20480,  25600,  32000,  32768,  40960,  51200]
    #dt_arr = [128, 256, 640]
    N_rep = range(4, 12);
    
    #T = [1.5, 1.8, 2, 2.2, 2.5, 3]
    T = [2]
    n = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.64, 0.68, 0.7, 0.72, 0.75, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9]
    n = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.64, 0.68, 0.7, 0.72, 0.75, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9]
    n = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.64, 0.67, 0.7, 0.73, 0.75, 0.77, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85]
    #n = [0.01, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8]
    #n = [0.64, 0.68, 0.7, 0.72, 0.75, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9]
    #n = [0.2, 0.3]
    
    T = [1, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 3]
    n = [0.8, 0.85, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.25, 1.3, 1.35, 1.4]
    #T = [1.5, 2, 2.5]
    #n = [1.2, 1.4, 1]
    
    # ------------------------------ phase --------------------------------
    for Ti in T:
        for ni in n:
            new_name = 'T' + str(Ti) + '_n' + str(ni)
                
            if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' Tmp ' + str(Ti) + ' n ' + str(ni))):
                return        
            if(not my.run_it('./full_cycle.py ' + new_name + ' -gen-pics-energy')):
                return
        
            if(not my.run_it('rm ' + new_name + '_param.dat')):
                return
            if(not my.run_it('rm ' + new_name + '_particles.xyz')):
                return                
        
    
    # ------------------------------ pressure --------------------------------
    #for Ni in N_arr:
        #new_name = 'N' + str(Ni)
            
        #if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' Ntot ' + str(Ni))):
            #return        
        #if(not my.run_it('./full_cycle.py ' + new_name + ' -gen-pics-energy')):
            #return
    
        #if(not my.run_it('rm ' + new_name + '_param.dat')):
            #return
        #if(not my.run_it('rm ' + new_name + '_particles.xyz')):
            #return                
    
    #for Ti in T:
        #for ni in n:
            #new_name = 'T' + str(Ti) + '_n' + str(ni)
                
            #if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' n ' + str(ni) + ' Tmp ' + str(Ti))):
                #return        
            #if(not my.run_it('./full_cycle.py ' + new_name + ' -gen-pics-enegry')):
                #return
        
            #if(not my.run_it('rm ' + new_name + '_param.dat')):
                #return
            #if(not my.run_it('rm ' + new_name + '_particles.xyz')):
                #return                
    
    # ------------------------- dynamic mmemory time -------------------------
    #for i in N_rep:
        #new_base_name = base_param_name + str(i)
        #if(not my.run_it('cp ' + base_param_name + '_' + my.param_file_suff + ' ' + new_base_name + '_' + my.param_file_suff)):
            #return        
        #if(not my.run_it('./gen ' + new_base_name)):
            #return
        
        #for dti in dt_arr:
            #new_name = 'dt' + str(dti) + '_' + str(i)
                
            #if(not my.run_it('./change_params.py ' + new_base_name + ' ' + new_name + ' dt ' + str(dti))):
                #return
            #if(not my.run_it('cp ' + new_base_name + '_' + my.particles_file_suff + ' ' + new_name + '_' + my.particles_file_suff)):
                #return
            #if(not my.run_it('./full_cycle.py ' + new_name + ' -move_res')):
                #return
        
            #if(not my.run_it('rm ' + new_name + '_' + my.param_file_suff)):
                #return
            #if(not my.run_it('rm ' + new_name + '_' + my.particles_file_suff)):
                #return   
                
        #if(not my.run_it('rm ' + new_base_name + '_' + my.param_file_suff)):
            #return                
        #if(not my.run_it('rm ' + new_base_name + '_' + my.particles_file_suff)):
            #return
        #if(not my.run_it('rm ' + new_base_name + '_gen.log')):
            #return                            
    
    # --------------------- rc time --------------------------
    #for ri in rc:
        #new_name = 'rc' + str(ri)
            
        #if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' r_max ' + str(ri))):
            #return        
        #if(not my.run_it('./full_cycle.py ' + new_name + ' -gen-pics-enegry')):
            #return
    
        #if(not my.run_it('rm ' + new_name + '_param.dat')):
            #return
        #if(not my.run_it('rm ' + new_name + '_particles.xyz')):
            #return                
    
    # ---------------------- N time ------------------------
    #for Ni in N_arr:
        #new_name = 'N' + str(Ni)
        #R = ((Ni/n)**(1/3)) / 2
            
        #if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' Ntot ' + str(Ni) + ' R ' + str(R))):
            #return        
        #if(not my.run_it('./full_cycle.py ' + new_name + ' -gen-pics-enegry')):
            #return
    
        #if(not my.run_it('rm ' + new_name + '_param.dat')):
            #return
        #if(not my.run_it('rm ' + new_name + '_particles.xyz')):
            #return            
    
    # ---------------------- dt time ------------------------
    #for dti in dt_arr:
        #new_name = 'dt' + str(dti)
            
        ##if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' dt ' + str(dti) + ' dumpDT ' + str(dmp_i))):
        #if(not my.run_it('./change_params.py ' + base_param_name + ' ' + new_name + ' dt ' + str(dti))):
            #return        
        #if(not my.run_it('./full_cycle.py ' + new_name + ' -energy-pics-gen')):
            #return
    
        #if(not my.run_it('rm ' + new_name + '_param.dat')):
            #return
        #if(not my.run_it('rm ' + new_name + '_particles.xyz')):
            #return   
    
    #if(not my.run_it('./clear.sh')):
    #    return
      
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
