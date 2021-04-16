import numpy as np

def fcc_points(aret):
    """main symmetry points of fcc Brillouin zone"""
    fcc_g = [0,                0,                0            ]
    fcc_x = [2.*np.pi/aret,    0,                0            ]
    fcc_l = [np.pi/aret,       np.pi/aret,       np.pi/aret   ]
    fcc_u = [2.*np.pi/aret,    np.pi/2./aret,    np.pi/2./aret]
    fcc_w = [2.*np.pi/aret,    np.pi/aret,       0            ]
    fcc_k = [3.*np.pi/2./aret, 3.*np.pi/2./aret, 0.           ]
    return fcc_g,fcc_x,fcc_l,fcc_u,fcc_w,fcc_k


def kline(p,q,kappa,ak,nk,jstart):
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
        ak.append(ak[jstart]+dist*jk*delta)                                      # append new point to array ak that contains x-scale
    ncheck=jstart+nk
    return kappa,ak,ncheck

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
