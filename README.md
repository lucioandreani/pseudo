# pseudo

## Energy bands of tetrahedral semiconductors by empirical pseudopotential

#### Lucio Andreani, Physics Department, University of Pavia, lucio.andreani@unipv.it, November 2020

To calculate empty lattice bands: empty.py or empty1.py, data and functions are inside the script

To calculate pseudopotential energy bands: pseudo.py with data file pseudo.dat and functions read.py, plot.py, sub.py

Gnuplot files for plotting: empty.gnu and pseudo.gnu<br>
Can be called using scripts ple and pl

Requirements: python3, scipy, numpy, matplotlib<br>
For plotting: gnuplot with pdf viewer

Installation: clone the repository, or [download](https://github.com/lucioandreani/pseudo/pseudo_v2.0.tar) and unpack the tar file. 

Version 1.0, 8 November 2020: reads data and runs single-shot

Version 2.0, 29 November 2020: driving code run0.py, run1.py, run2.py to perform automatic loop over gmax and convergence test
