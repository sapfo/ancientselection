#!/usr/bin/env python


""" 
compute the different matrices and such
"""
import numpy as np
from scipy.linalg import expm as expm
import math
import wrap_lap
import funcs

def Qmatrix(Ne=None,s=None,gamma=None,xx=np.array([0.0,0.1,0.4,0.7,1.0]),verbose=0,h=0.5):
    if gamma==None:
        gamma=2*Ne*s

    H=len(xx)
    G=H-1    

    #print len(xx)
    
    Q=np.zeros((H,H))
    if verbose:
        Delta=[];Beta=[];

    for ele in range(1,H-1):
        #print ele,

        xiplus=xx[ele+1]-xx[ele]
        ximinus=xx[ele-1]-xx[ele]

        xup=xx[ele+1]
        x=xx[ele]
        xdown=xx[ele-1]

        ci=xx[ele]*(1-xx[ele])
        if verbose and 0:
            print "ele ",ele
            print "x ",x
            print "xup ",xup
            print "xdown ",xdown
            print "xiplus ",xiplus
            print "ximinus ",ximinus
            print "diff1",xiplus*ximinus-ximinus**2
            print "h ",h
            raw_input()
        beta=((-1 + x)*x*(-1-(x**2)*gamma + h*(-1 + 2*x)*(x - xdown)*gamma+x*xdown*gamma))/((x - xup)*(xdown -xup))
        delta=-((-1 + x)*x*(-1-(x**2)*gamma  + h*(-1 + 2*x)*(x - xup)*gamma+x*xup*gamma ))/((x - xdown)*(xdown - xup))

        #if verbose:    
        #    if h==0.5:
        #        delta2=ci*(0.5*gamma*xiplus-1)/(xiplus*ximinus-ximinus**2)
        #        beta2=ci*(0.5*gamma*ximinus-1)/(xiplus*ximinus-xiplus**2)
        #        print "delta, delta2 ",delta,delta2
        #        print "beta, beta2 ",beta,beta2
                
        #        raw_input()
        
        if verbose:
            #print 'delta, beta ',delta,beta
            Delta.append(delta)
            Beta.append(beta)

        if delta<=0 or beta<=0:
  #          print 's ',s
#            print 'Ne ',Ne
            print "ele ",ele
            print "x ",x
            print "xup ",xup
            print "xdown ",xdown
            print "xiplus ",xiplus
            print "ximinus ",ximinus
            print "nom dom0",(1+x**2*gamma-x*xup*gamma)
            print "ci ",ci
            print 'delta ',delta
            print 'beta ',beta
            print 'gamma ',gamma
            print 'maxi ',2./max(np.diff(xx))
            #raw_input()
            return None,None,None
            #return None
            
            raw_input("Problem with the definiton of delta or beta")
            #raise   ValueError("Problem with the definiton of delta or beta")
        Q[ele,ele]=-(delta+beta)
        Q[ele,ele+1]=beta
        Q[ele,ele-1]=delta

        #print sum(Q[ele,])
        
    if verbose:
        return Q,np.array(Delta),np.array(Beta)
    else:
        return Q


def diagonalize(Q=np.zeros((3,3)),verbose=0):
    evals, evecs  =np.linalg.eig(Q)
    U = evecs
    L = np.diag(evals)
    Ui = np.linalg.inv(U)
    if verbose:
        return U,L,Ui,evals
    return U,L,Ui


def exponentiate(U=np.ones((3,3)),L=np.zeros((3,3)),Ui=np.ones((3,3)),time=0.1):
    if time<0:
        raise ValueError('time has to be poistive to exponentiate!! there is a bug!!')
    #print 'shape ',L.shape[0]
    Lnew=np.zeros(L.shape)
    for i in range(L.shape[0]):
        value = L[i,i]
        value = math.e**(value*time)
        Lnew[i,i] = value
    
    MatrixExp=np.dot(np.dot(U,Lnew),Ui)

    return MatrixExp



def expmScipy(A,q=7):
    """Compute the matrix exponential using Pade approximation.

    Parameters
    ----------
    A : array, shape(M,M)
        Matrix to be exponentiated
    q : integer
        Order of the Pade approximation

    Returns
    -------
    expA : array, shape(M,M)
        Matrix exponential of A

    """
    A = asarray(A)
    ss = True
    if A.dtype.char in ['f', 'F']:
        pass  ## A.savespace(1)
    else:
        pass  ## A.savespace(0)

    # Scale A so that norm is < 1/2
    nA = np.norm(A,Inf)
    if nA==0:
        return identity(len(A), A.dtype.char)
    from numpy import log2
    val = log2(nA)
    e = int(floor(val))
    j = max(0,e+1)
    A = A / 2.0**j

    # Pade Approximation for exp(A)
    X = A
    c = 1.0/2
    N = eye(*A.shape) + c*A
    D = eye(*A.shape) - c*A
    for k in range(2,q+1):
        c = c * (q-k+1) / (k*(2*q-k+1))
        X = dot(A,X)
        cX = c*X
        N = N + cX
        if not k % 2:
            D = D + cX;
        else:
            D = D - cX;
    F = solve(D,N)
    for k in range(1,j+1):
        F = dot(F,F)
    pass  ## A.savespace(ss)
    return F

def expm2(A):
    """Compute the matrix exponential using eigenvalue decomposition.

    Parameters
    ----------
    A : array, shape(M,M)
        Matrix to be exponentiated

    Returns
    -------
    expA : array, shape(M,M)
        Matrix exponential of A

    """
    A = asarray(A)
    t = A.dtype.char
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'
    s,vr = eig(A)
    vri = inv(vr)
    return dot(dot(vr,diag(exp(s))),vri).astype(t)

def expm3(A,q=20):
    """Compute the matrix exponential using Taylor series.

    Parameters
    ----------
    A : array, shape(M,M)
        Matrix to be exponentiated
    q : integer
        Order of the Taylor series

    Returns
    -------
    expA : array, shape(M,M)
        Matrix exponential of A

    """
    A = asarray(A)
    t = A.dtype.char
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'
    A = mat(A)
    eA = eye(*A.shape,**{'dtype':t})
    trm = mat(eA, copy=True)
    castfunc = cast[t]
    for k in range(1,q):
        trm *= A / castfunc(k)
        eA += trm
    return eA



def tools_expo(Q=None,Beta=None,Delta=None,spteqr=0,threshold=1e-5,verbose=0,verbose2=0):
    N=Q.shape[0]
    

    q=Q[1:-1,1:-1]
    
    delta=Delta[1:]
    beta=Beta[:-1]

    #### compute D and Dinv
    dnext=1.
    
    #ddiag=np.array([1]+[dnext*delta[ele-1]/beta[ele-1] for ele in range(1,N-2)]+[1])
    
    ddiag=[1]
    for ele in range(1,N-2):
        ddiag.append(dnext)
        dnext*=(1.*delta[ele-1]/beta[ele-1])**(1./2)
        

    ddiag.append(dnext)
    exposant=-0.5*np.log10(dnext)
    
    ddiag.append(1)
    ddiag=np.array(ddiag)
    ddiag[1:-1]=(10**exposant)*ddiag[1:-1]
    
    D=np.diag(ddiag)
    Dinv=np.diag(1./ddiag)
    
    #### compute vecpropres valpropres of s
    dinv=Dinv[1:-1,1:-1]
    d=D[1:-1,1:-1]
    s=np.dot(np.dot(dinv,q),d)
    
    if spteqr:
        sdiagonal=-np.diagonal(s)
        soffdiagonal=-np.diagonal(s,1)
        compz = 'i'
        dimension = len(sdiagonal)
        z = np.empty((dimension,dimension))
        info = wrap_lap.spteqr(compz,sdiagonal,soffdiagonal,z)
        print 'info ',info
        oevals = -sdiagonal
        oevecs = z 
        #raw_input()
    else:
        oevals,oevecs=np.linalg.eigh(s)

    o=oevecs
    oT=np.transpose(o)
    
    #### compute O,
    Ot=np.eye(N)
    Ot[1:-1,1:-1]=oT
    O=np.eye(N)
    O[1:-1,1:-1]=o
    
    #--------------------------- check accuracy --------------------------
    flag='ok'
    shouldBeId=np.dot(np.dot(D,O),np.dot(Ot,Dinv))    
    totest=np.abs(funcs.norm_sum(shouldBeId)-np.trace(shouldBeId))
    if verbose:
        print "threshold: ",threshold
        print "totest: ",totest
    
    if np.abs(funcs.norm_sum(shouldBeId)-np.trace(shouldBeId))>=threshold:
        flag='warning'
    if verbose: print "flag: ",flag

    #### compute Lam
    lam=np.zeros(N)
    lam[1:-1]=oevals
    Lam=np.diag(lam)
    #### compute Lamtild
    lamtild=np.zeros(N)
    lamtild[1:-1]=1./oevals
    Lamtild=np.diag(lamtild)

    ### compute T
    R=np.dot(np.dot(Dinv,Q),D)

    T=np.dot(np.dot(Ot,R),O)
    
    ### compute V
    
    V=np.zeros((N,N))

    V[:,0]=T[:,0]
    V[:,-1]=T[:,-1]
    if verbose2:
        return D,Dinv,s,o,oT,O,Ot,V,lam,Lam,Lamtild,totest,flag

    return D,Dinv,s,o,oT,O,Ot,V,lam,Lam,Lamtild,flag


def tools_q_expo(Q=None,Beta=None,Delta=None,spteqr=0,threshold=None,verbose=0):
    N=Q.shape[0]

    q=Q[1:-1,1:-1]
    
    delta=Delta[1:]
    beta=Beta[:-1]

    #### compute D and Dinv
    dnext=1.
    
    #ddiag=np.array([1]+[dnext*delta[ele-1]/beta[ele-1] for ele in range(1,N-2)]+[1])
    
    ddiag=[1]
    for ele in range(1,N-2):
        ddiag.append(dnext)
        dnext*=(1.*delta[ele-1]/beta[ele-1])**(1./2)
        

    ddiag.append(dnext)
    exposant=-0.5*np.log10(dnext)
    
    ddiag.append(1)
    ddiag=np.array(ddiag)
    ddiag[1:-1]=(10**exposant)*ddiag[1:-1]

    D=np.diag(ddiag)
    Dinv=np.diag(1./ddiag)

    #### compute vecpropres valpropres of s
    dinv=Dinv[1:-1,1:-1]
    d=D[1:-1,1:-1]
    s=np.dot(np.dot(dinv,q),d)
    
    if spteqr:
        sdiagonal=-np.diagonal(s)
        soffdiagonal=-np.diagonal(s,1)
        compz = 'i'
        dimension = len(sdiagonal)
        z = np.empty((dimension,dimension))
        info = wrap_lap.spteqr(compz,sdiagonal,soffdiagonal,z)
        print 'info ',info
        oevals = -sdiagonal
        oevecs = z 
        #raw_input()
    else:
        oevals,oevecs=np.linalg.eigh(s)
        #print oevals
        #print len(oevals)

        
    o=oevecs
    #print o
    #raw_input()
    oT=np.transpose(o)
    #--------------------------- check accuracy --------------------------
    flag='ok'
    shouldBeId=np.dot(np.dot(d,o),np.dot(oT,dinv)) 
    totest=np.abs(funcs.norm_sum(shouldBeId)-np.trace(shouldBeId))

    if verbose:
        print "threshold (HAHA) ",threshold
        print "totest: ",totest
    
    if np.abs(funcs.norm_sum(shouldBeId)-np.trace(shouldBeId))>=threshold:
        flag='warning'
    if verbose: print "flag: ",flag

    #### compute Lam 
    lam=oevals  

    return d,dinv,s,o,oT,lam,flag

def expBabikq(d=None,dinv=None,o=None,ot=None,lam=None,time=0,verbose=0,positive=False):
    n=d.shape[0]
    
    if time==0:
        return np.eye(n)

    do=np.dot(d,o)
    doinv=np.dot(ot,dinv)
    factor=1. #desperate attempt to correct the problem with small times.
    expLam=np.diag(factor*np.exp(time*lam))

    if verbose:
        print "lam ",lam
        print "np.exp(time*lam) ",factor*np.exp(time*lam)
        print "np.min(exp(time*lam) ",np.min(factor*np.exp(time*lam))
        print "np.max(exp(time*lam) ",np.max(factor*np.exp(time*lam))

        print 'expLam ',expLam
        #raw_input()
    
    expqt=np.dot(np.dot(do,expLam),doinv)    
    if positive:
      expqt=funcs.nonneg(expqt)  
    return 1/factor*expqt


def expBabik(D=None,Dinv=None,O=None,Ot=None,V=None,Lam=None,Lamtild=None,lam=None,time=0,verbose=0,positive=False):
    N=D.shape[0]
    
    if time==0:
        return np.eye(N)

    DO=np.dot(D,O)
    DOinv=np.dot(Ot,Dinv)
    
    expLam=np.diag(np.exp(time*lam))
    if verbose:
        print 'expLam ',expLam
        raw_input()
    other=np.dot(expLam-np.eye(N)-time*Lam,np.dot(Lamtild,V))
    if verbose:
        print 'other ',other
        raw_input()

    middle=expLam+time*V+other
    if verbose:
        print 'middle ',middle
        raw_input()
    expQt=np.dot(np.dot(DO,middle),DOinv)
    if verbose:
        print 'DO DO-1 ',np.dot(DO,DOinv)
        print np.sum(np.dot(DO,DOinv))
        raw_input()
    if verbose:
        print 'expQt ',np.sum(expQt)
        print expQt
        raw_input()
    if positive:
        expQt=funcs.nonneg(expQt)

    return expQt

def ExpMatrices_q(Qmat=np.eye(10),Betamat=np.eye(10),Deltamat=np.eye(10),Times=[None,None],verbose=0,threshold=1e-5):
    
    Times=np.array(Times)
    d,dinv,s,o,oT,lam,flag=tools_q_expo(Q=Qmat,Beta=Betamat,Delta=Deltamat,threshold=threshold,verbose=verbose)    

    ExpMatrices=[expBabikq(d=d,dinv=dinv,o=o,ot=oT,lam=lam,time=tau,positive=True) for tau in Times]
    
    Times=list(Times)
    return ExpMatrices,flag

def ExpMatrices_Q(Qmat=np.eye(10),Betamat=np.eye(10),Deltamat=np.eye(10),Times=[None,None],verbose=1,threshold=1e-5):
    Times=np.array(Times)
    
    D,Dinv,s,o,oT,O,Ot,V,lam,Lam,Lamtild,flag=tools_expo(Q=Qmat,Beta=Betamat,Delta=Deltamat,threshold=threshold)

    ExpMatrices=[expBabik(D=D,Dinv=Dinv,O=O,Ot=Ot,V=V,Lam=Lam,Lamtild=Lamtild,lam=lam,time=tau,positive=False) for tau in Times]
    
    Times=list(Times)

    return ExpMatrices,flag

if 0:
    xx=np.linspace(0,1,10)
    print xx
    Qmat,Deltamat,Betamat=Qmatrix(gamma=10,xx=xx,verbose=1,h=0.5)
    print Qmat
    raw_input()
    qmat=Qmat[1:-1,1:-1]
    print qmat
    print qmat.shape
    tau=0.01
    d,dinv,s,o,oT,lam=tools_q_expo(Q=Qmat,Beta=Betamat,Delta=Deltamat)

    D,Dinv,s,o,oT,O,Ot,V,lam,Lam,Lamtild=AlgLin.tools_expo(Q=Qmat,Beta=Betamat,Delta=Deltamat)
    

    Pt=expBabikq(d=d,dinv=dinv,o=o,ot=oT,lam=lam,time=tau)
    print "------------------------------------"
    print Pt
    print "-------------------------------------"
    Pt=expm((tau)*qmat)
    print Pt
