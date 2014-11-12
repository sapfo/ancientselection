#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import scipy
 #./plotting_2d.py ASIP_runMay4_2
import argparse

parser = argparse.ArgumentParser(description='''Takes as input some arrays in npy format and produce a likelihood surface with the maximas (over the grid) indicated. It produces a 2D graph by slicing one Ne at a time''')

parser.add_argument('-i','--likarray', help='The lik array file, expects also the grid used to be in the same place',type=str,required=True)
#parser.add_argument('--datadir', help='directory where the arrays are',type=str,required=True)
parser.add_argument('--dirout', help='directory where write the output data -- default LikSurface',type=str,required=False,default="LikSurface")

args = parser.parse_args()

likarray = args.likarray
project = likarray.rstrip(".npy").split("Lik_array_")[1]
datadir = os.path.split(likarray)[0]
dirout = args.dirout
print "Run id: ",project
print "datadir: ",datadir
print "dirout: ",dirout

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
	
if not os.path.isdir(dirout):
    print "Your outputdir does not exist, created now"
    mkdir_p(dirout)


 
#smaller font on the legend
params = {'legend.fontsize': 6,}
axislabelfontsize = 8
py.rcParams.update(params)

T0=np.load("%s/T0_%s.npy"%(datadir,project))
Gamma=np.load("%s/Gamma_%s.npy"%(datadir,project))
NE=np.load("%s/NE_%s.npy"%(datadir,project))
Tinput=np.load("%s/Tinput_%s.npy"%(datadir,project))

Lik_array=np.load(likarray)

dims=Lik_array.shape
argmax=scipy.argmin(Lik_array,axis=None)
idminall = np.unravel_index(argmax, dims)
idminall_t0,idminall_gamma,idminall_Ne=idminall
minlikall= Lik_array[idminall]
maxlikall= np.max(Lik_array)

PoptAll=T0[idminall_t0],Gamma[idminall_gamma],NE[idminall_Ne]

print "Grid: "
print "T0: ",T0
print "Gamma: ",Gamma
print "NE: ",NE


print "Min lik all: ",minlikall 
print "\'MLEs\':  ",PoptAll
print "index for min lik: ",idminall
print "shape Likarray: ",dims

#raw_input("raw_input")
jump=int(round(len(NE)/3.1))
#jump=int(len(NE)/3.)

print "len(NE) ",len(NE)
print "Slice Lik array every ",jump," NE"
print "NE slices:" 
if len(NE) == 1:
    print NE[0]
    jump = 1
else:
    print NE[0::jump]

for cutI in range(0,len(NE),jump):
    #slice for several NE

    cut=NE[cutI]

#    Lik_array_slice=Lik_array[:,:,idmin_Ne], optimal slice
    Lik_array_slice=Lik_array[:,:,cutI]
    print "log-likelihood surface: ",Lik_array_slice

#    print "SLICING at ",NE[idmin_Ne]
    print "Slicing at NE = ",NE[cutI]

    title=r"log-likelihood surface $\gamma\times t_0$ for $N_e=%i$"%(int(cut))
    if 0:
        if cut==PoptAll[-1]:
            title+=" the global mle"
    #title+="_"+project
        
    dims=Lik_array_slice.shape
    argmax=scipy.argmin(Lik_array_slice,axis=None)
    idmin = np.unravel_index(argmax, dims)
    idmin_t0,idmin_gamma=idmin
    minlik= Lik_array_slice[idmin]
    Popt=T0[idmin_t0],Gamma[idmin_gamma]

    print "Min lik (for that slice) ",minlik 
    print "MLEs (for that slice) ",Popt    
    #print "Local optimum (for that slice of Ne) ",idmin_t0,idmin_gamma    
    #print "Indexidmin ",idmin
    
    if len(Gamma)>1 and len(T0)>1:
        
        fig=py.figure()
        fig.subplots_adjust(left=0.2, bottom=0.15, right=None, top=0.85) 
#    levels = np.arange(np.min(-Lik_array_slice), np.max(-Lik_array_slice),0.2)
        levels = np.arange(-maxlikall,-minlikall,2)

#    CS = plt.contour(T0,Gamma,Lik_array_slice,levels,alpha=0.5)
        CS = plt.contour(Gamma,T0,-Lik_array_slice,levels,alpha=0.5)

        CS2 = plt.contour(Gamma,T0,-Lik_array_slice,[-minlikall-2.996],colors='black',linestyles="dashed",lw=3)

        plt.hlines(Tinput,[min(Gamma)],[max(Gamma)],lw=2,label='sampling times')

        CB = plt.colorbar(CS, shrink=0.8, extend='both',format='%.1f')

        py.xlabel(r'$\gamma$',fontsize=axislabelfontsize)

        py.ylabel(r'$t_0$ [generations] units',fontsize=axislabelfontsize)
 
        py.scatter(Gamma[idmin_gamma],T0[idmin_t0],color='black',s=10,label='argmax=($%.1f$, $%.1f$)'%(Gamma[idmin_gamma],T0[idmin_t0]))
        
        plt.ylim(np.min(T0),np.max(T0))
        plt.xlim(np.min(Gamma),np.max(Gamma))
        titlefontsize = 8
        py.title(title,fontsize=titlefontsize)
        py.legend(loc="lower right")
        print "Likelihood surface saved under: %s/%s_slice%s.pdf"%(dirout,project,cutI)
        py.savefig("%s/%s_slice%s.pdf"%(dirout,project,cutI))
    else:
        py.figure()
        if len(Gamma)==1 and len(T0)>1:
            Xaxis = T0
            py.xlabel(r'$t_0$ [generations] units',fontsize=axislabelfontsize)
            Yaxis = list(-Lik_array[:,0,0])
            plt.vlines(Tinput,min(-Lik_array[:,0,0]),max(-Lik_array[:,0,0]),lw=2,label='sampling times')

            py.title(r"log-likelihood for $\gamma=%s$ and $N_e=%i$"%(Gamma[0],int(cut)))
        elif len(T0)==1 and len(Gamma)>1:
            Xaxis = Gamma
            py.xlabel(r'$\gamma$',fontsize=axislabelfontsize)
            Yaxis = list(-Lik_array[0,:,0])
            py.title(r"log-likelihood for $t_0=%s$ and $N_e=%i$"%(T0[0],int(cut)))
        else:
            print "You array has dimensions ",Lik_array.shape
            print "You need at least one of the dimensions to be bigger than 1"
            sys.exit()
        
        
        py.plot(Xaxis,Yaxis)
        print "Likelihood surface saved under: %s/%s_1dim.pdf"%(dirout,project)              
        py.savefig("%s/%s_1dim.pdf"%(dirout,project))
        #py.show()
