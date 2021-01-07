import sys
import os
import re
import pathlib

import mylib_molecules as my

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc = len(args)
    if(argc > 1):
        print('usage:python \n' + sys.argv[0] + '   [output_filename]')
        sys.exit(1)
    output_filename = (args[0] if argc == 1 else 'eff_file.dat')
            
    N_arr = [16, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]
    N_arr = [32]
    
    Np_arr = [1, 2, 4, 5, 8]
    #Np_arr = [1, 2, 4]
    
    base_path = os.path.join(pathlib.Path.home(), '!molecules', 'RES', 'DATA')
    f_total = open(output_filename, 'w')
    for Ni in N_arr:
        for Npi in Np_arr:
            model_name = 'p' + str(Npi) + '_N' + str(Ni)
            log_path = os.path.join(base_path, model_name, my.complog_filename)
            log_file = open(log_path, "r")
            log_data = log_file.read()
            log_file.close()
            
            found_substr = re.search(r'efficiency \(e \= endT\/dt\*N\^2\/real\_t\) = \d\.\d+e\+\d+', log_data).group(0)
            eff_str = found_substr[found_substr.rfind(' '):]
            f_total.write(eff_str + ' '*(12 - len(eff_str)))
        f_total.write('\n')
    f_total.close()
    print('DONE')

            
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
