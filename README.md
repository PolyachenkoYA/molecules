# molecules
This is an educational project in molecular dynamics

C++, OpenMP, CUDA for numerical integration of Newton's equations

Matlab + Python3 + [ovito](https://ovito.org) for visualisation of results

# installation & usage

## install 
Download the repo:

1) `git clone https://github.com/PolyachenkoYA/molecules.git`

Put the python code to where it will be visible (substitute your python version or path):

2) `cd molecules/\!COMP`

3) `cp ./mylib_molecules.py $HOME/anaconda3/lib/python3.7/site-packages/mylib_molecules.py`

`nvcc` is assumed to be in `/usr/local/cuda-10.1/bin/nvcc`

If it is not, then you need to write a correct way to call `nvcc` on your machine. It 

Compile:

4) `./recomplie_all.sh`

## simple test
Go to the working directory:

5) `cd \!go`

Run the full computation from generating a system to plotting results:

6) `python full_cycle.py tst -all`

Go see the results:

7) `cd ../../RES/DATA/tst/graph`

These are visualized diagnostics of the system.

Parameters of the system are listed in 'RES/DATA/tst/param.dat'.
