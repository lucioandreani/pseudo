import numpy as np
import matplotlib.pyplot as plt
import sub
from matplotlib.offsetbox import OffsetImage,AnnotationBbox
#from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.image as mpimg

def plot_bands(nbands, nktot, ak, energ, kv, color, ymin, ymax, ny, linewidth=2, line=True, label=''):
    """ plot the energy bands with matplotlib """ 

    #plt.text(0.4, 0.1, 'r/a='+str(T[0]), fontsize=15)
    #fig, ax = plt.subplots(figsize=(6,8))

    #plt.figure(figsize=(6,8))
    plt.xlim([0, 1])
    x_ticks=[0.]  # starts defining xticks array corresponding to vertical lines
    for k in kv:
       x_ticks.append(k)
    x_ticks.append(1.)
    ##print(x_ticks)
    x_labels=['L', '$\Gamma$', 'X', 'U,K', '$\Gamma$']              # and here are the labels for the BZ
    #x_labels=['$\Gamma$', 'X', 'W', 'L', '$\Gamma$', 'K,U', 'X']     # and here are the labels for the BZ
    plt.xlabel('Wavevector k', fontsize=14)
    plt.xticks(ticks=x_ticks, labels=x_labels, fontsize=12)

    #ymin,ymax,ny=-15,5,5        # change here to set range of y axis
    y_ticks = np.linspace(ymin,ymax,ny)
    #print(y_ticks)
    plt.yticks(y_ticks, fontsize=12)
    plt.ylim([ymin, ymax])
    plt.ylabel('Energy (eV)', fontsize=14)

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
       if (i!=0):
           label=''
       if line:
           plt.plot(ak,energ[i,:],linestyle='solid',color=color, linewidth=linewidth, label=label)  # use this to plot with lines
       else:
           plt.plot(ak,energ[i,:],'o',marker='.', color=color, markersize=linewidth, label=label)     # use this to plot with points


#   insert fcc Brillouin zone
    #arr_bz = mpimg.imread('bz3.png')
    #imagebox = OffsetImage(arr_bz, zoom=0.25)
    #ab = AnnotationBbox(imagebox, (0.7, 1), bboxprops=dict(edgecolor='none'))
    #ax.add_artist(ab)


