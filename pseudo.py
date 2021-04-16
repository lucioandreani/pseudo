#!/usr/bin/env python3
# this is the shebang syntax, allow to execute simply with the command pseudo.py

import numpy as np
from scipy.linalg import eigh             # scipy function for solving hermitian linear eigenvalue problem
import matplotlib.pyplot as plt           # plotting routine
import sys
sys.path.append('../.')                   # allows to search for modules in parent directory, if not found in current directory
from sub import fcc_points, kline, sort4  # user-defined function
from plot import plot_bands               # user-defined function for plotting
from read import read_data                # user-defined function for reading data

"""Code pseudo.py, written by Lucio Andreani, lucio.andreani@unipv.it
   It calculates the energy bands of tetrahedral semiconductors by the empirical pseudopotential method
   Reference: Yu-Cardona, Fundamentals of Semiconductors, Springer
   Atomic units are used (Bohr radius for length, Hartree for energy)"""

# reads data from file pseudo.dat
f=open('pseudo.dat','r')
aret,v3s,v8s,v11s,v3a,v4a,v11a,nk,gmax,nmax,upper,ymin,ymax,ny,jplot = read_data(f)
print(aret,v3s,v8s,v11s,v3a,v4a,v11a,nk,gmax,nmax,upper,ymin,ymax,ny,jplot)
f.close


#aret=5.43     # lattice constant in Angstrom
#v3s,v8s,v11s= -0.211, 0.04, 0.08    #v3s,v8s,v11s=simmetric pseudopotential form factors in Rydberg
#v3a,v4a,v11a=0.   , 0.  , 0.  ,    #v4a,v4a,v11a=antisimmetric pseudopot. form factors in Rydberg
#nk=50        # number of k-points along each line in BZ
#gmax=5.      # maximum modulus of reciprocal lattice vector, in units of 2\pi/a
#nmax=5       # max index number for reciprocal lattice vectors, choose high enough to get all bands up to max energy of the plot
#upper=1e-12  # numerical tolerance parameter: upper limit for evaluating zero


# fundamental constants
abohr  =0.529177210903   # Bohr radius in Angstrom
hartree=27.211386245988  # Hartree energy in eV

# conversion to atomic units
aret=aret/abohr
v3s,v8s,v11s=v3s/2,v8s/2,v11s/2
v3a,v4a,v11a=v3a/2,v4a/2,v11a/2

# primitive vectors of Bravais lattice
a1=[0., aret/2, aret/2]
a2=[aret/2, 0., aret/2]
a3=[aret/2, aret/2, 0.]

# volume of unit cell
ucvolume=a1[0]*(a2[1]*a3[2]-a2[2]*a3[1]) \
        +a1[1]*(a2[2]*a3[0]-a2[0]*a3[2]) \
        +a1[2]*(a2[0]*a3[1]-a2[1]*a3[2])
print('Primitive vectors of Bravais lattice: a1,a2,a3,ucvolume=',a1,a2,a3,ucvolume)

# primitive vectors of reciprocal lattice

b1,b2,b3=np.empty(3),np.empty(3),np.empty(3)

b1[0]=2*np.pi/ucvolume*(a2[1]*a3[2]-a2[2]*a3[1])
b1[1]=2*np.pi/ucvolume*(a2[2]*a3[0]-a2[0]*a3[2])
b1[2]=2*np.pi/ucvolume*(a2[0]*a3[1]-a2[1]*a3[0])

b2[0]=2*np.pi/ucvolume*(a3[1]*a1[2]-a3[2]*a1[1])
b2[1]=2*np.pi/ucvolume*(a3[2]*a1[0]-a3[0]*a1[2])
b2[2]=2*np.pi/ucvolume*(a3[0]*a1[1]-a3[1]*a1[0])

b3[0]=2*np.pi/ucvolume*(a1[1]*a2[2]-a1[2]*a2[1])
b3[1]=2*np.pi/ucvolume*(a1[2]*a2[0]-a1[0]*a2[2])
b3[2]=2*np.pi/ucvolume*(a1[0]*a2[1]-a1[1]*a2[0])

# volume of Brillouin zone
bzvolume=b1[0]*(b2[1]*b3[2]-b2[2]*b3[1]) \
        +b1[1]*(b2[2]*b3[0]-b2[0]*b3[2]) \
        +b1[2]*(b2[0]*b3[1]-b2[1]*b3[0])

print('Primitive vectors of reciprocal lattice: b1,b2,b3/(4\pi),bzvolume/(4\pi)^3)=',b1/4/np.pi,b2/4/np.pi,b3/4/np.pi,bzvolume/(4*np.pi)**3)

# calculates reciprocal lattice vectors

n=0                         # initialize counter of reciprocal lattice vectors
gx,gy,gz,gmod=[],[],[],[]   # initialize lists of wavevectors
for n1 in range(-nmax,nmax+1):
    for n2 in range(-nmax,nmax+1):
        for n3 in range(-nmax,nmax+1):
            #print(n1,n2,n3,n)
            g=n1*b1+n2*b2+n3*b3    # reciprocal lattice vector (list with 3 elements)
            gx.append(g[0])        # array with x-component
            gy.append(g[1])        # array with y-component
            gz.append(g[2])        # array with z-component
            gmod.append(np.sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]))   # array with modulus
            #print(n,n1,n2,n3,gx[n],gy[n],gz[n],gmod[n])
            n=n+1                  # increment counter
ng=n                               # this is the total number of reciprocal lattice vectors
print('reciprocal lattice vectors have been calculated: ng=',ng)

# sorts according to modulus of reciprocal lattice vectors
#print('\n','before sort: n, gx gy gz gmod in units of 2\pi/aret:')
#for n in range(0,ng):
#    print(n,gx[n]/(2*np.pi/aret),gy[n]/(2*np.pi/aret),gz[n]/(2*np.pi/aret),gmod[n]/(2*np.pi/aret))
gmod,gx,gy,gz = sort4(gmod,gx,gy,gz)  # this routing sorts the four arrays in increasing order of gmod
print('reciprocal lattice vectors have been sorted in order of increasing modulus')
#print('\n','after sort: n, gx gy gz gmod in units of 2\pi:')
#for n in range(0,ng):
#    print(n,gx[n]/(2*np.pi/aret),gy[n]/(2*np.pi/aret),gz[n]/(2*np.pi/aret),gmod[n]/(2*np.pi/aret))
#print()

# writes file with reciprocal lattice vectors
zg=aret/2./np.pi
f=open('pseudo_gvec.out','w')
f.write('# reciprocal lattice vectors: n gx gy gz gmod   \n'  )
for n in range (0,ng):                 # this is the loop over the band index
    f.write('%6i %11.5f' '%11.5f' '%11.5f' '%11.5f' '\n'  %(n,gx[n]*zg,gy[n]*zg,gz[n]*zg,gmod[n]*zg) )
f.close()
print('reciprocal lattice vectors have been written on file pseudo_gvec.out')



# prepares k-vectors for calculating energy bands
fcc_g, fcc_x, fcc_l, fcc_u, fcc_w, fcc_k=fcc_points(aret)    # defines symmetry points in fcc Brillouin zone
#print('symmetry points in BZ zone: G,X,L,K,U,W=', fcc_g,fcc_x,fcc_l,fcc_k,fcc_u,fcc_w)
n_symmetry_points=4
nktot=n_symmetry_points*nk+1  # total number of wavevectors (nk in each segment)

#kappa,ax=[fcc_g],[0]  # starting values for kappa-vectors (kx ky kz), and for array ax used for plotting
#kappa,ax,ncheck=kline(fcc_g,fcc_x,kappa,ax,nk,0)
#kappa,ax,ncheck=kline(fcc_x,fcc_w,kappa,ax,nk,nk)
#kappa,ax,ncheck=kline(fcc_w,fcc_l,kappa,ax,nk,2*nk)
#kappa,ax,ncheck=kline(fcc_l,fcc_g,kappa,ax,nk,3*nk)
#kappa,ax,ncheck=kline(fcc_g,fcc_k,kappa,ax,nk,4*nk)
#kappa,ax,ncheck=kline(fcc_u,fcc_x,kappa,ax,nk,5*nk)
#kappa,ax,kv,ncheck=kline(fcc_k,fcc_x,kappa,ax,kv,nk,5*nk)

# this is another possible path in BZ
kappa,ak=[fcc_l],[0]  # starting values for kappa-vectors (kx ky kz), and for array ax used for plotting,
kappa,ak,ncheck=kline(fcc_l,fcc_g,kappa,ak,nk,0)
kappa,ak,ncheck=kline(fcc_g,fcc_x,kappa,ak,nk,nk)
kappa,ak,ncheck=kline(fcc_x,fcc_u,kappa,ak,nk,2*nk)
kappa,ak,ncheck=kline(fcc_k,fcc_g,kappa,ak,nk,3*nk)
if ncheck+1 != nktot:
    print('ncheck is not equal to nktot, error in code',ncheck,nktot)
    quit()

# normalize arrays ak for later plotting of the bands
#print(kappa,ak)
for jk in range(0,nktot):                 # this is the wavevector axis that has to run from 0 to 1
    ak[jk]=ak[jk]/ak[nktot-1]
#print(kappa,ak)
print('symmetry points and lines have been defined')

# calculates vertical lines for plot along BZ paths
kv=[]
for j in range(1,n_symmetry_points):
    kv.append(ak[j*nk])
#print('kv=',kv)


# find number of reciprocal lattice vectors with modulus up to gmax
for n in range(0,ng):
    if gmod[n]/(2*np.pi/aret) <= gmax:
        npw=n+1
print('gmax,npw=',gmax,npw)

#   calculates array with pseudopotential form factors
form_factor=np.zeros([npw,npw],dtype='complex')
#print(form_factor)
en0=4.*np.pi*np.pi/aret/aret
for n1 in range(0,npw):
  for n2 in range(0,npw):
    #g=[gx[n],gy[n],gz[n]]
    ggx,ggy,ggz=gx[n1]-gx[n2],gy[n1]-gy[n2],gz[n1]-gz[n2]
    ggmod=ggx*ggx+ggy*ggy+ggz*ggz
    g_dot_s=(ggx+ggy+ggz)*aret/8.
    if   np.abs(ggmod-3.*en0 ) < upper:  # (111) reciprocal lattice vector
        form_factor[n1,n2]=complex(v3s*np.cos(g_dot_s) ,-v3a *np.sin(g_dot_s))
    elif np.abs(ggmod-4.*en0 ) < upper:  # (200) reciprocal lattice vector
        form_factor[n1,n2]=complex(0.                  ,-v4a *np.sin(g_dot_s))
    elif np.abs(ggmod-8.*en0 ) < upper:  # (220) reciprocal lattice vector
        form_factor[n1,n2]=complex(v8s*np.cos(g_dot_s) ,0.                  )
    elif np.abs(ggmod-11.*en0) < upper:  # (311) reciprocal lattice vector
        form_factor[n1,n2]=complex(v11s*np.cos(g_dot_s),-v11a*np.sin(g_dot_s))
#   --- end of loop over reciprocal lattice vectors ---

# starts loop over wavevector for calculating energy bands
bands=np.empty([npw,nktot], dtype='float')  # initializes array with band energies
for jk in range(0,nktot):                  # this is the loop over the wavevector
    kk=kappa[jk]
    #print(jk,kk)

#   calculates hamiltonian
    ham=np.zeros([npw,npw], dtype='complex') # initialization
    for n1 in range(0,npw):         # first loop over reciprocal lattice vectors
        g1=[gx[n1],gy[n1],gz[n1]]
        #energ[n1,jk]=((kk[0]+g1[0])**2+(kk[1]+g1[1])**2+(kk[2]+g1[2])**2)/(2.*np.pi/aret)**2  # these are the empty lattice bands in units of Ex
#       second loop over reciprocal lattice vectors
        for n2 in range(0,npw):     # second loop over reciprocal lattice vectors
            if n1 == n2:
                ham[n1,n2]=((kk[0]+g1[0])**2+(kk[1]+g1[1])**2+(kk[2]+g1[2])**2)/2.                # diagonal part: kinetic energy
            else:
                ham[n1,n2]=form_factor[n1,n2]
    #print(ham)

#   solves linear eigenvalue problem to find band energies, and stores them in array "bands"
    if (jk//nk)*nk == jk:   # writes only for the starting value of kappa in each line
      print('diagonalizing for kappa=',kk[0]/(2*np.pi/aret),kk[1]/(2*np.pi/aret),kk[2]/(2*np.pi/aret))
    a = eigh(ham, eigvals_only=True)
    for n in range(0,npw):
        bands[n,jk] = a[n]

#   defines zero of energy for rescaling and minimum energy for plotting vertical lines
    if (abs(kk[0])<upper) and (abs(kk[1])<upper and (abs(kk[2])<upper)):    # this selects the Gamma point
        ezero=a[3]                                                          # this is the energy of the 4th band at the Gamma point
        emin =a[0]                                                          # energy of lowest band, for writing vertical lines in file
#   --- end of loop over wavevector ---

#   rescales energies with zero at topmost valence band
for jk in range(0,nktot):
    for n in range(0,npw):
        bands[n,jk] = (bands[n,jk]-ezero)*hartree

emin = (emin-ezero)*hartree
print('emin=',emin)

# writes energy bands on file pseudo.out, in units of eV
nbands=min(100,npw)
f=open('pseudo.out','w')
f.write('# nbands = %5i   nktot = %5i \n'  %(nbands,nktot))
f.write('# ymin = %10.3f ymax = %10.3f    ny = %3i  \n'  %(ymin,ymax,ny))
for n in range (0,nbands):                 # this is the loop over the band index
    f.write('# band number = %3i    ak kx ky kz E(eV) \n'  %(n+1))
    for jk in range(0,nktot):          # this is the loop over the wavevector, for each band
        kk=kappa[jk]
        z=2*np.pi/aret
        f.write('%12.5f' '%12.5f' '%12.5f' '%12.5f' '%12.5f' '\n'  %(ak[jk],kk[0]/z,kk[1]/z,kk[2]/z,bands[n,jk]) )
    f.write('\n')
f.write('\n')
f.close()
print('energy bands have been written on file pseudo.out')

# writes vertical lines in file pseudo_vl.out
nvl=len(kv)
f=open('pseudo_vl.out','w')
f.write('# number of vertical lines = %2i  \n'  %(nvl))
for k in kv:
   f.write('%11.5f' '%11.5f' '\n'      %(k, -100))
   f.write('%11.5f' '%11.5f' '\n' '\n' %(k, 100))
f.close()
print('vertical lines have been written on file pseudo_vl.out')


# plots energy bands
if jplot == 1:
   plt.figure(figsize=(6,6))
   color='red'
   linewidth='2'    # linewidth or markersize for plot
   line=True
   label='energy bands'
   plot_bands(npw, nktot, ak, bands, kv, color, ymin, ymax, ny, linewidth, line, label)
   plt.legend(loc='upper left', fontsize=16)
   plt.show()
