#! /usr/bin/env python


import numpy as np
import sys
from binoms import binoms


def binom(n=100,m=10):
    """ compute binomial coefficient n is bigger than m """
    b=[0]*(n+1)
    b[0]=1
    for i in xrange(1,n+1):
        b[i]=1
        j=i-1
        while j > 0:
            b[j]+=b[j-1]
            j-=1
    return b[m]



def WFmatrix(s=0.1,h=0.5,Pop=100,verboseWF=1):
    """ compute WF matrix with selection:
    AA: 1+s    Aa: 1+sh    aa: 1 
    Ewens p.24
    We have a a precomputed table for Pop<=1000"""

    if verboseWF:
        print "---------"
        print "Pop ",Pop
        print "s ",s
        print "h ",h
        print "---------"
    Binoms = []
    if Pop>1000:
        for j in range(Pop+1):        
            Binoms.append(binom(Pop,j))
    A=np.zeros((Pop+1,Pop+1))    
    for i in range(Pop+1):
        for j in range(Pop+1):
            eta_nom=(1+s)*i**2+(1+s*h)*i*(Pop-i)
            #print 'eta_nom ',eta_nom
            eta_denom=(1+s)*i**2+2*(1+s*h)*i*(Pop-i)+(Pop-i)**2
            #print 'i,j ',i,j
            #print (1+s)*i**2
            #print 2*(1+s*h)*i*(Pop-i)
            #print (Pop-i)**2
            #print 'eta_denom ',eta_denom
            eta=1.0*eta_nom/eta_denom
            eta=float(eta)
            if Pop>1000:      
                A[i,j]=Binoms[j]*eta**j*(1-eta)**(Pop-j)
            else:             
                #if verboseWF:
                #    print "binoms[Pop][j] ",binoms[Pop][j],type(binoms[Pop][j])
                #    print "eta ",eta,type(eta)
                    
                A[i,j]=binoms[Pop][j]*eta**j*(1-eta)**(Pop-j)

    return A

def sampling(n=10,mut=3,f=0.1):
    """ sampling from a Pop with frequ of mutants f , samplesize n, numbuer of mutants in sample mut"""
    if n<1000:
        #print mut,n,f
        #print binoms[n][mut]
        #print 'haha'
        #print 'hallo ',(f)**mut
        #print 'pfffff ',(1-f)**(n-mut)
        return float(binoms[n][mut])*(f)**mut*(1-f)**(n-mut)
    else:
        return float(binom(n,mut))*(f)**mut*(1-f)**(n-mut)

def Apower(mat=np.zeros((3,3)),timestep=3,verbose=0):
    """ matrix to the power timestep. input matrix A """
    if verbose: print "timestep: ",timestep
    Aout=np.linalg.matrix_power(mat, timestep)
    return Aout



def nonneg(A=None,fillvalue=0):
    condition=A<0
    Acond=np.ma.masked_array(A,mask=condition)
    Anoneg=Acond.filled(fill_value=fillvalue)
    return Anoneg

def norm_sum(A=np.eye(10)):    
    return np.sum(abs(A))

def process_expmethod(Ne=None,gamma=None,dominance=None,grid=None,gammathreshold=40,H=None,Hsecond=100,expmethod=None,expmethodsecond='pade',verboseP=0):
    if grid!='default':
        print "the grid is different than default and we have not done the cartography"
        print "WARNING: have not looked into what the gammamax should be"

    if verboseP:        
        print "H before ",H," gamma before ",gamma
    # ------------- step 1 --------: check gamma not too large.
    expmethodnew=expmethod
    if abs(gamma)>gammathreshold:
        print "gamma = ",gamma
        print "abs(Gamma) too high, need more precision "
        expmethodnew=expmethodsecond
        if expmethodnew=='alglin':
            print "!! WARNING: gamma too high but you picked alglin as alternative !!"    
    # -------------- step 2 -------: check H large enough. 
    if expmethodnew=='alglin':
        Hnew=H
    else:
        Hnew=Hsecond;
    
    if dominance==0.5:
        Hmin=int(10.+abs(gamma)*0.84)
        if Hnew<Hmin:
            Hnew=Hmin
    elif dominance==0 or dominance==1:
        Hmin=int(10.+abs(gamma)*1.1)
        if Hnew<Hmin:
            Hnew=Hmin
    else:
        raise ValueError("Only looked at dominance 0, 0.5, or 1")
    if verboseP:        
        print "H after ",H," gamma after ",gamma
   
    return expmethodnew,Hnew

def domains(data=None,smallert=-10,verboseD=1,fixed_time=None):
    """
    Define list of valid intervals
    """
    if fixed_time!=None:
        return [(int(fixed_time),int(fixed_time))]

    T=data[-1];I=data[0];M=data[1]

    if smallert>T[0]:        
        print "Your lower bound on the time is larger than the smaller sampling time....please check the bounds....the first samples wont be used..."
        print "T: ",data[-1]
        print "Presumably lower bound on T: ",smallert
        sys.exit()

    T=[smallert]+data[-1];I=data[0];M=data[1]
    
    if verboseD:
        print "original input"
        print "I ",I,"M ",M,"T ",T
    Domains=[]
    stopat=False
    while not stopat and len(I)>1:
        Domains.append((int(T[0])+1,int(T[1])))
        if I[0]!=0:
            stopat=True
            break
        I=I[1:];M=M[1:];T=T[1:]
    if verboseD:print "Domains: ",Domains
    return Domains


