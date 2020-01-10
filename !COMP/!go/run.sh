./change_params.py tst tst compMode -5 Nthreads 1
./cycle.py tst gpu
./change_params.py tst tst compMode 5
./cycle.py tst CPU
