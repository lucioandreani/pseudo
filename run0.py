#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import subprocess
import sys
#sys.path.append('../.')               # allows to search for modules in parent directory, if not found in current directory
from read import read_bands,read_vl   # functions that read results from pseudo.out
from plot import plot_bands           # function that plots the photonic bands


########### starts function that writes data of fortran code #########

def write_data(T):
    """Writes the data file of the fortran code."""
    f=open('pseudo.dat','w')
    f.write('1#   5.43      aret=lattice constant \n')
    f.write('2# -0.211 0.04 0.08    v3s,v8s,v11s=simmetric pseudopotential form factors \n')
    f.write('3# 0. 0. 0.            v4a,v8a,v11a=antisimmetric pseudopot. form factors  \n')
    f.write('4# 20                  nk= number of kappa points along each line in BZ \n')
    f.write('5# 5      gmax=maximum modulus of reciprocal lattice vector, in units of 2\pi/a --> determines the number npw of plane waves: min is 2 \n')
    f.write('6# 5      nmax= index number for generating reciprocal lattice vectors, choose high enough to get the bands up to the energy you wish \n')
    f.write('7# 1e-12  upper=numerical tolerance parameter: upper limit for evaluating zero\n')
    f.write('8# -15. 10. 6          ymin,ymax,ny=limiting values and tick number for plot\n')
    f.write('9# 0                   jplot=flag for plot: =0 does not plot, =1 plots\n')
    f.write('####  end of data, do not change or remove this line! ####   \n')
    f.close()



########### starts function that runs the fortran executable #########

def execute_code(T):
    """ Calculates energy bands by running fortran executable """

    name='MAIN'                     # name of working folder
    try:                            # Create folder if not present
        os.stat(name)
    except:
        os.mkdir(name)
    shutil.copy('pseudo.py',name+'/.')  # copy python code pseudo.py inside folder
    os.chdir(name)                    # change to work directory

    write_data(T)                    # write input file pseudo.dat

    #execute fortran code
    print('executing pseudo.py code')
    p1=subprocess.run("./pseudo.py", stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    #print('returncode=',p1)            # useful to inspect if something goes wrong
    #print('stdout=',p1.stdout)         # standard output: contains error messages from the code
    if p1.stdout[0:4] == b'STOP':       # error message from the code! stop and look what is wrong
       print('error message in code execution:',p1.stdout)
       quit()
    if p1.returncode !=0:               # error message from the process, better to stop and have a look
       print("\nerror returned:",p1.stderr)
       quit()

    #copy interesting output files into parent directory
    shutil.copy("pseudo.out","../.")
    shutil.copy("pseudo_vl.out","../.")

    os.chdir('..')                  # cd to parent directory
    #shutil.rmtree(name)             # remove working folder

    return


########### starts function that reads the results from output of fortran code #########

def read_results():
    """ Read energy bands from output of fortran code 
    Files must be in directory MAIN """

    name='MAIN'                     # name of working foloder
    try:                            # Create folder if not present
        os.stat(name)
    except:
        print('from read_bands: directory MAIN not present')
        exit()
    os.chdir(name)                  # change to work directory


    #reads energy bands from file pseudo.out
    f=open('pseudo.out','r')
    nbands,nktot,kappa,energ,ymin,ymax,ny=read_bands(f)
    f.close()

    #reads k-vectors of vertical lines from file pseudo_vl.out
    f=open('pseudo_vl.out','r')
    nkv,kv=read_vl(f)
    f.close()

    os.chdir('..')                  # cd to parent directory
    #shutil.rmtree(name)             # remove working folder

    return nbands, nktot, kappa, energ, nkv, kv, ymin, ymax, ny


########### starts driving code #################à

T=[]

execute_code(T)
nbands, nktot, kappa, energ, nkv, kv, ymin, ymax, ny = read_results()
color='blue'
linewidth=2
line='True'
plot_bands(nbands, nktot, kappa, energ, kv, color, ymin, ymax, ny, linewidth, line)

plt.show()




