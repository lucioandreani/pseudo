#!/usr/bin/env python3  
# this is the shebang syntax, allow to execute simply with the command empty.py

import numpy as np
import matplotlib.pyplot as plt


def fcc_points():
    """main symmetry points of fcc Brillouin zone"""
    fcc_g=[0,0,0]
    fcc_x=[2.*np.pi,0,0]
    fcc_l=[np.pi,np.pi,np.pi]
    fcc_u=[2.*np.pi,np.pi/2.,np.pi/2.]
    fcc_w=[2.*np.pi,np.pi,0]
    fcc_k=[3.*np.pi/2.,3.*np.pi/2.,0.]
    return fcc_g,fcc_x,fcc_l,fcc_u,fcc_w,fcc_k


def kline(p,q,kappa,ax,nk,jstart):
    """defines line in Brillouin zone, from point p to q"""
    diff=[q[0]-p[0],q[1]-p[1],q[2]-p[2]]
    dist=np.sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2])
    if nk != 0:
      delta=1./nk
    else:
      delta=0.
    #print('nk,jstart',nk,jstart)
    for jk in range(1,nk+1):
        newk=[p[0]+diff[0]*jk*delta,p[1]+diff[1]*jk*delta,p[2]+diff[2]*jk*delta] # new wavevector point
        kappa.append(newk)                                                       # append new wavevector point to the array
        ax.append(ax[jstart]+dist*jk*delta)                                      # append new point to array ax that contains x-scale
    ncheck=jstart+nk
    return kappa,ax,ncheck

def sort4(gmod,gx,gy,gz):
    """sorts lists in ascending order of gmod"""

    index=np.argsort(gmod)   # this is the index list of the sorted gmod list
    #print(index,gmod)

    n =len(gmod)
    nx=len(gx)
    ny=len(gy)
    nz=len(gz)
    #print('\n','n nx ny nz=',n,nx,ny,nz,'\n')
    if (n != nx) or (n != ny) or (n != nz):
        print('from sort: wrong array dimensions')
        quit()

    gmod1=np.empty(n,dtype='float') # we need to define new arrays: if we use the same for input and output ones, python goes wrong
    gx1  =np.empty(n,dtype='float')
    gy1  =np.empty(n,dtype='float')
    gz1  =np.empty(n,dtype='float')
    for j in range(0,n):            # we reorder all arrays according to sorting of array gmod
        gmod1[j]=gmod[index[j]]
        gx1[j]=gx[index[j]]
        gy1[j]=gy[index[j]]
        gz1[j]=gz[index[j]]

    return gmod1,gx1,gy1,gz1



def plot_bands(nbands, nktot, ax, energ, kv, color, size=10):
    """ plot the energy bands with matplotlib """ 

    #plt.text(0.4, 0.1, 'r/a='+str(T[0]), fontsize=15)
    plt.xlim([0, 1])
    x_ticks=[0.]  # starts defining xticks array corresponding to vertical lines
    for k in kv:
       x_ticks.append(k)
    x_ticks.append(1.)
    ##print(x_ticks)
    #x_labels=['L', '$\Gamma$', 'X', 'K', '$\Gamma$']              # and here are the labels for the BZ
    x_labels=['$\Gamma$', 'X', 'W', 'L', '$\Gamma$', 'K,U', 'X']     # and here are the labels for the BZ
    plt.xlabel('Wavevector k', fontsize=14)
    plt.xticks(ticks=x_ticks, labels=x_labels, fontsize=12)

    ymin,ymax,ny=0.,5.,6        # change here to set range of y axis
    y_ticks = np.linspace(ymin,ymax,ny)
    #print(y_ticks)
    plt.yticks(y_ticks, fontsize=12)
    plt.ylim([0., ymax])
    plt.ylabel('E/E$_x$', fontsize=14)

    # vertical lines for symmetry points Brillouin zone
    en_vl=[ymin,ymax]
    for k in kv:
       x_ticks.append(k)
       kappa_vl=[]
       kappa_vl.append(k)
       kappa_vl.append(k)
       #print(kappa_vl,en_vl)
       plt.plot(kappa_vl,en_vl,linestyle='dotted',color='grey')

    # now plot the bands, eventually!
    for i in range(nbands):
       #plt.plot(ax,energ[i,:],linestyle='dotted',color=color, linewidth=size)  # use this to plot with lines
       plt.plot(ax,energ[i,:],'o',marker='.', color=color, markersize=size)     # use this to plot with points




# MAIN CODE STARTS HERE

aret=1.  # lattice constant, we use dimensionless units
nmax=5   # max index number for reciprocal lattice vectors, choose high enough to get all bands up to max energy of the plot
nk=50    # number of k-points along each line in BZ

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
#print('\n','before sort: n, gx gy gz gmod in units of 2\pi:')
#for n in range(0,ng):
    #print(n,gx[n]/2/np.pi,gy[n]/2/np.pi,gz[n]/2/np.pi,gmod[n]/2/np.pi)
gmod,gx,gy,gz = sort4(gmod,gx,gy,gz)  # this routing sorts the four arrays in increasing order of gmod
print('reciprocal lattice vectors have been sorted in order of increasing modulus')
#print('\n','after sort: n, gx gy gz gmod in units of 2\pi:')
#for n in range(0,ng):
    #print(n,gx[n]/2/np.pi,gy[n]/2/np.pi,gz[n]/2/np.pi,gmod[n]/2/np.pi)
#print()


# prepares k-vectors for calculating energy bands
fcc_g, fcc_x, fcc_l, fcc_u, fcc_w, fcc_k=fcc_points()    # defines symmetry points in fcc Brillouin zone
#print('symmetry points in BZ zone: G,X,L,K,U,W=', fcc_g,fcc_x,fcc_l,fcc_k,fcc_u,fcc_w)
n_symmetry_points=6
nktot=n_symmetry_points*nk+1  # total number of wavevectors (nk in each segment)
kappa,ax=[fcc_g],[0]  # starting values for kappa-vectors (kx ky kz), and for array ax used for plotting
kappa,ax,ncheck=kline(fcc_g,fcc_x,kappa,ax,nk,0)
kappa,ax,ncheck=kline(fcc_x,fcc_w,kappa,ax,nk,nk)
kappa,ax,ncheck=kline(fcc_w,fcc_l,kappa,ax,nk,2*nk)
kappa,ax,ncheck=kline(fcc_l,fcc_g,kappa,ax,nk,3*nk)
kappa,ax,ncheck=kline(fcc_g,fcc_k,kappa,ax,nk,4*nk)
kappa,ax,ncheck=kline(fcc_u,fcc_x,kappa,ax,nk,5*nk)
#kappa,ax,kv,ncheck=kline(fcc_k,fcc_x,kappa,ax,kv,nk,5*nk)
# this is another possible path in BZ
#kappa,ax=[fcc_l],[0]  # starting values for kappa-vectors (kx ky kz), and for array ax used for plotting,
#kappa,ax=kline(fcc_l,fcc_g,kappa,ax,kv,nk,0)
#kappa,ax=kline(fcc_g,fcc_x,kappa,ax,kv,nk,nk)
#kappa,ax=kline(fcc_x,fcc_k,kappa,ax,kv,nk,2*nk)
#kappa,ax=kline(fcc_k,fcc_g,kappa,ax,kv,nk,3*nk)
if ncheck+1 != nktot:
    print('ncheck is not equal to nktot, error in code',ncheck,nktot)
    quit()

# normalize arrays ax for later plotting of the bands
#print(kappa,ax)
for jk in range(0,nktot):                 # this is the vector axis that has to run from 0 to 1
    ax[jk]=ax[jk]/ax[nktot-1]
#print(kappa,ax)
print('symmetry points and lines have been defined')

# calculates vertical lines for plot along BZ paths
kv=[]
for j in range(1,n_symmetry_points):
    kv.append(ax[j*nk])
#print('kv=',kv)



# calculates empty lattice bands
energ=np.empty([ng,nktot], dtype='float')  # initializes array with band energies
for jk in range(0,nktot):                  # this is the loop over the wavevector
    kk=kappa[jk]
    #print(jk,kk)
    for n in range (0,ng):                 # this is the loop over the band index
        g=[gx[n],gy[n],gz[n]]
        #print(n,g,gmod[n])
        energ[n,jk]=((kk[0]+g[0])**2+(kk[1]+g[1])**2+(kk[2]+g[2])**2)/(2.*np.pi)**2  # these are the empty lattice bands in units of Ex
print('empty lattice bands have been calculated')


# writes energy bands on file empty.out
f=open('empty.out','w')
for n in range (0,min(100,ng)):                 # this is the loop over the band index
    f.write('# band number=%3i    ax kx ky kz E/Ex \n'  %(n+1))
    for jk in range(0,nktot):          # this is the loop over the wavevector, for each band
        kk=kappa[jk]
        f.write('%11.5f' '%11.5f' '%11.5f' '%11.5f' '%11.5f' '\n'  %(ax[jk],kk[0],kk[1],kk[2],energ[n,jk]) )
    f.write('\n')
f.write('\n')
f.close()
print('empty lattice bands have been written on file empty.out')

# writes vertical lines in file empty_vl.out
f=open('empty_vl.out','w')
for k in kv:
   f.write('%11.4f' '%9.2e' '\n' %(k,0))
   f.write('%11.4f' '%9.2e' '\n' '\n' %(k,1e6))
f.close()
print('vertical lines have been written on file empty_vl.out')


# plots energy bands
color='red'
size='5'    # markersize for plot
plot_bands(ng, nktot, ax, energ, kv, color, size)

plt.show()
