# pseudo package:

## Energy bands of tetrahedral semiconductors by pseudopotential calculation

### written by Lucio Andreani, lucio.andreani@unipv.it, November 2020

To calculate empty lattice bands: empty.py or empty1.py, data and functions are inside the script

To calculate pseudopotential energy bands: pseudo.py with data file pseudo.dat and functions read.py, plot.py, sub.py

Gnuplot files for plotting: empty.gnu and pseudo.gnu
Can be called using scripts ple and pl

Requirements: python3, scipy, numpy, matplotlib
For plotting: gnuplot with pdf viewer

Version 1.0, 8 November 2020: reads data and runs single-shot

Version 2.0, 29 November 2020: driving code run0.py, run1.py, run2.py to perform automatic loop over gmax and convergence test
