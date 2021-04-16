#!/usr/bin/env python3

import re
import numpy as np


####### starts function read_data  #########

def read_data(f):
    """reads data from input file
    f should be a filehandle to a previously opened file"""
    for line in f:
        line=re.sub("^\s+","",line).rstrip() # cleans line from initial spaces
        nums=re.split('\s+',line)            # splits the line
        if (nums[1]=='end') and (nums[2]=='of') and (nums[3]=='data,'):        # last line after data, quits and returns
            return aret,v3s,v8s,v11s,v3a,v4a,v11a,nk,gmax,nmax,upper,ymin,ymax,ny,jplot
        elif (nums[0]=='1#'):                  # reading lattice constant
            aret=float(nums[1])
        elif (nums[0]=='2#'):                # reading sym form factors
            v3s,v8s,v11s=float(nums[1]),float(nums[2]),float(nums[3])
        elif (nums[0]=='3#'):                # reading antisym form factors
            v3a,v4a,v11a=float(nums[1]),float(nums[2]),float(nums[3])
        elif (nums[0]=='4#'):                # reading number of k-points
            nk=int(nums[1])
        elif (nums[0]=='5#'):                # reading maximum g
            gmax=float(nums[1])
        elif (nums[0]=='6#'):                # reading maximum n
            nmax=int(nums[1])
        elif (nums[0]=='7#'):                # reading numerical tolerance parameter
            upper=float(nums[1])
        elif (nums[0]=='8#'):                # reading limiting values for plot
            ymin,ymax,ny=float(nums[1]),float(nums[2]),int(nums[3])
        elif (nums[0]=='9#'):                # reading flag for plotting or not
            jplot=int(nums[1])
        else:
            print('wrong starting characters: nums[0]=',nums[0])

####### ends function read_data  #########



############ starts function read_bands ################

def read_bands(f):
    """read bands from an output file in the format of pseudo
    f should be a filehandle to a previously opened file"""
    count=0
    kappa=[]
    jk=0
    for line in f:
       line=re.sub("^\s+","",line).rstrip() # cleans line from initial spaces
       nums=re.split('\s+',line)            # splits the line
       #print('line=',line)
       #print('nums=',nums)
       if (nums[0]=='#') and (nums[1]!='band'):                 # reading comment lines at the beginning of file
           #print('I am going through the first commented lines')
           if (nums[1]=='nbands') and (nums[4]=='nktot'):       # line with number of bands and of wavevectors
             #print('reading line with nbands and nktot')
             #print('I found the line with band and k number')
             nbands=int(nums[3])                                # assign number of bands. Warning! Line in file has to be formatted with spaces
             nktot =int(nums[6])                                # assign number of wavevectors
             #print('nbands=',nbands,'   nktot=',nktot)
             energ = np.zeros([nbands,nktot], dtype="float")    # initializing array with bands
           elif (nums[1]=='ymin') and (nums[4]=='ymax'):        # line with ymin,ymax,ny
             #print('reading line with ymin, ymax,ny')
             ymin,ymax,ny=float(nums[3]),float(nums[6]),int(nums[9])  # assign ymin,ymax,ny
             #print('ymin,ymax,ny=',ymin,ymax,ny)
       elif (nums[0]=='#') and (nums[1]=='band'):               # reading header of a band
           #print('reading header of band number', int(nums[4]))
           count=count+1
           jband=int(nums[4])
           if (count != jband):                                 # check band number, one never knows
             print('something wrong in file=',f)
       elif (line==''):                                         # blank line
           #print('blank line')
           jk=0                                                 # resets the number of wavevectors, as we are about to read another band
       else:                                                    # interesting line, with kappa and energy
           #print(float(nums[0]),float(nums[1]),'hello')
           #print(line)
           #print(nums)
           jk=jk+1
           #print('jband,jk,energ[0,0]=',jband-1,jk,energ[0,0])
           energ[jband-1,jk-1]=float(nums[4])
           if (count == 1):                                     # wavevector array is defined only once, for the first band
             kappa.append(float(nums[0]))
           #print(kappa)
           #print('hello')

    #print(energ)
    return nbands,nktot,kappa,energ,ymin,ymax,ny
############ end function read_bands ################



############ starts function read_vl ################

def read_vl(f):
    """read vertical lines for plotting in Brillouin zone, from an output file in the format of pseudo
    f should be a filehandle to a previously opened file"""
    count=0
    kv=[]
    nkv=0
    first=True
    for line in f:
       line=re.sub("^\s+","",line).rstrip() # cleans line from initial spaces
       nums=re.split('\s+',line)            # splits the line
       #print('line=',line)
       #print('nums=',nums)
       if (nums[0]=='#') and (nums[3]=='vertical'):                 # reading header with number of vertical lines
           #print('reading header of band number', int(nums[6]))
           count=count+1
       elif (line==''):                                         # blank line
           #print('blank line')
           nk=0                                                 # resets the number of wavevectors, as we are about to read another vertical line
       else:                                                    # interesting line, with kappa and energy
           if (first):                                          # this is the first of the two lines, stores the wavevector
             #print('I am in the true')
             nkv=nkv+1
             kv.append(float(nums[0]))
             first=False
           else:                                                # this is the second line of the wavevector, does not need to store
             #print('I am in the false')
             first='True'

    print (kv)
    return nkv,kv

############ ends function read_vl ################



############ starts driving code for test ##############
if __name__ == "__main__":
    import sys

    desc="""read and print pseudo.e output file.
    This code serves just for testing the above library functions."""

    filename='pseudo.out'
    f=open(filename,'r')
    nbands,nktot,kappa,energ,ymin,ymax,ny=read_bands(f)
    #print('nbands,nk=',nbands,nktot)
    #print('k-vectors=',kappa)
    #print('band energies=',energ)
    f.close()


    filename='pseudo_vl.out'
    f=open(filename,'r')
    nkv,kv=read_vl(f)
    #print(nkv,kv)




