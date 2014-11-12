#!/usr/bin/env python
""" 
Computing the likelihood of the data
"""
import numpy as np
import Grids
import AlgLin
import funcs
import os

#######################################################################
def preprocess(time=None,Iori=[1,2],Mori=[10,10],Tori=[-10,1],verbose=1):
    """
    Decide whether the proposed time is 
    1. not feasible
    2. needs to ignore first Iori    
    """
    ## CHECK
    if len(Iori)!=len(Mori) or len(Iori)!=len(Tori):
        print "original input:"
        print "Iori ",Iori,"Mori ", Mori,"Tori ",Tori
        raise ValueError("input I M T should be same size")
    if verbose:
        print "original input:"
        print "Iori ",Iori,"Mori ", Mori,"Tori ",Tori
    T=Tori;M=Mori;I=Iori;
    problem=False
 
    while time>T[0]:##atention: equal should be ok.means frequency is 1/2Ne and first step is 0
        if I[0]!=0:  
            problem=True
            break
        I=I[1:];M=M[1:];T=T[1:]
        problem=False

    if problem:
        if verbose:
            print 'Time intervals ',T
            print 'time ',time
            print "the time %i is not what it should be"%time
        return None,None,None
    if verbose: print "returned values ",I,M,T
    return I,M,T

#preprocess(time=-2,Iori=[0,0,2],Mori=[10,10,10],Tori=[-10,-3,1])


#######################################################################
def ll_q(params=[None,None,None],data=[None,None,None],H=100,verbose=1,expmethod='alglin',expmethodsecond='pade',loga=1,grid='default',positive=1,dominance=0.5,threshold=1e-5,numberprecision=128,thresholdfrerot=1e-5):
   
    """
    Computes the likelihood when given the data, conditioned you never reach fixation until the last sampling time.    
    """

    Hsmall=H-2
    if verbose:    
    #if 1:
        print "data I,M,T",data
        print "params t0,gamma,Ne",params
        print "dominance ",dominance
        print "H ",H
        print "Hsmall ",Hsmall
        raw_input(" Press a key to continue:here,now ")    
    
    #------------- define parameters -------------------------------
    
    t_=params[0];gamma=params[1];Ne=params[2];
    #-------------- redefine expmethod H ----------------------------
    expmethodOri=[expmethod][0]

    expmethod,H=funcs.process_expmethod(Ne=Ne,gamma=gamma,dominance=dominance,grid=grid,gammathreshold=40,H=H,expmethod=expmethod,expmethodsecond=expmethodsecond)

    if expmethodOri!=expmethod: 
        print "expmethod adjusted because high gamma"
        print "expmethod new,H,dominance ",expmethod,H,dominance

    Hsmall=H-2
    if verbose:        
        print "after processing...."
        print "data I,M,T",data
        print "params t0,gamma,Ne",params
        print "H ",H
        print "Hsmall ",Hsmall
        #raw_input(" Press a key to continue ")    
    #------------- define data -------------------------------------

    I=data[0];M=data[1];T=data[2]

    #---- check allele age is appropriate, if not return np.nan ----
    I,M,T=preprocess(time=t_,Iori=I,Mori=M,Tori=T,verbose=0)
    if verbose:
        print "Input data after preprocessing I,M,T",I,M,T
        print "input allele age ",t_
    
    if I==None and M==None and T==None:
        if verbose: 
            print "likelihood is",0
            print "OUTPUT: ",np.nan
        return 5*[np.nan]   #Attention! number of outputs.
    
    # Compute all t[j]-t[j-1]/2Ne (store in DeltaT).
    
    DeltaT=list((np.array(T)-np.array([t_]+T[:-1]))/(2.0*Ne))
    if verbose: print "DeltaT ",DeltaT
    DeltaT.append(sum(DeltaT))

    # -----------------Define the space grids -----------------------
                  
    # xx is for the grid for transition from any frequ
    # oriFI is the index of the 1/2Ne freq for the transition from 1/2Ne
    # xxSmall is for the grid for transition from any frequ, without boundaries
    # oriFI-1 is the index of the 1/2Ne freq for the transition from 1/2Ne

    if grid=='default':
        xx,oriFI=Grids.Ne_grid(H,Ne)        
    elif grid=='symmetric':
        xx=Grids.symmetric_grid(H)
    elif grid=='uniform':
        xx,oriFI=Grids.uniform_grid(H,Ne)         
    elif grid=='expo':
        indice=int(0.10*H)
        xx,oriFI=Grids.exponential_grid_Ne(indice=indice,H=H,Ne=Ne)
        
    xxSmall=xx[1:-1]                        
    maxi=max(np.diff(xx))

    if verbose: 
        print "grid ",grid
        print "xxSmall ",xxSmall[0:10],'...',xxSmall[-10:]        
        print "orginal frequency ",xxSmall[oriFI-1]
        print "what the frq should be ",1./(2*Ne)

    #Check gamma is appopropriate, gamma more than 60 gives unstabilities ---------------------------
    
    #if abs(gamma)>60: 
    #    raise ValueError("gamma is too high or too small for current double precision, max is 60") 
    #elif gamma>=1./maxi or gamma<=-1./maxi:
    #    raise ValueError("gamma is too high or too small for grid space, max is  %f"%(1./maxi))    
    
    Qmat,Deltamat,Betamat=AlgLin.Qmatrix(gamma=gamma,xx=xx,verbose=1,h=dominance)
    
    qmat=Qmat[1:-1,1:-1]
    if verbose:
        print "xx.shape ",xx.shape
        print "Qmat.shape ",Qmat.shape
        print "rate out 1/2Ne of qmat\n ",qmat[oriFI-1,oriFI]
        print "rate out 1/2Ne of Qmat\n ",Qmat[oriFI,oriFI+1]
    
    # --------- expqt for all t --------------------------------

    #---------- matrix diagonalization-babik--------------------
    if expmethod=='alglin':
      #  d,dinv,s,o,oT,lam=AlgLin.tools_q_expo(Q=Qmat,Beta=Betamat,Delta=Deltamat)    
        PtMatrices,flag=AlgLin.ExpMatrices_q(Qmat=Qmat,Betamat=Betamat,Deltamat=Deltamat,Times=DeltaT,verbose=verbose,threshold=threshold)
        if verbose: 
            print "FLAG: ",flag
        if flag == 'warning':
            print "warning, the exponentiation seems problematic"
    #----------- pade approximation ---------------------------
    elif expmethod=="pade":
        PtMatrices=[AlgLin.expm((dt)*qmat) for dt in DeltaT]

    #---------- high precision call ----------------------------
    elif expmethod=='prec':        
        if grid=="default":
            
            argu_frerot=[numberprecision,thresholdfrerot,gamma,dominance,H,1,Ne]+DeltaT
            argu_frerot=" ".join([str(w) for w in argu_frerot])
        else:
            raise ValueError("only grid implemented is the default grid")
        command="./ancientSelectionApaSmallQ %s"%argu_frerot      
        print command
        os.system(command)       
        PtMatrices=[np.loadtxt('exp_Qt_%i.dat'%dtindex, dtype=np.float64,delimiter=" ") for dtindex in range(len(DeltaT))]
        
        [os.system("rm exp_Qt_%i.dat"%dtindex) for dtindex in range(len(DeltaT))]
        
    if len(PtMatrices)!=len(T)+1:
        raise ValueError("something went wrong with matrix computation")
    #PtMatrices=[] #where all matrices are stored
    #for dtindex,dt in enumerate(DeltaT):
    #    if verbose: print 'Interval of time dt ',dt
    #    if expmethod=="babik":
    #        Ptmatrice=AlgLin.expBabikq(d=d,dinv=dinv,o=o,ot=oT,lam=lam,time=dt,positive=True)
    #    elif expmethod=="pade":
    #        Ptmatrice=AlgLin.expm((dt)*qmat)
    #    elif expmethod=="frerot":    
    #        #import them from the right file
    #        matname='exp_Qt_%i.dat'%dtindex
    #        print matname
    #        Ptmatrice=np.eye(Hsmall)
    #        #Ptmatrice=np.loadtxt(matname, dtype=np.float64,delimiter=" ")
    #    if verbose:
    #        print "shape of Ptmatrice ",Ptmatrice.shape
    #    PtMatrices.append(Ptmatrice)

    Ptfinal=PtMatrices[-1] ## "renormalisation"

    PtOri=PtMatrices[0]    ## I.C.

    ExpMat=PtMatrices[1:-1]## all others
    
    if verbose:
        #check starting point, i.e. that grids works well
        print "starting point ",(2.*Ne)*xxSmall[oriFI-1]

        print "PtOri ",PtOri
    #----------- nofixation proba -----------------------------

    prob_nofix=np.sum(Ptfinal[oriFI-1,])     
 
    #----------- initial conditions ---------------------------

    IC=[]
    for j1 in range(Hsmall):  
        #print funcs.sampling(n=M[0],mut=I[0],f=xxSmall[j1])
        #print PtOri[oriFI-1,j1]
        IC.append(funcs.sampling(n=M[0],mut=I[0],f=xxSmall[j1])*PtOri[oriFI-1,j1])
                
    #----------- Loop over next times ---------------------------
        
    L=IC # store in L, elements of the sum
    if verbose:
        print "IC ",L
        print "len(L) ",len(L)

    for ind in range(1,len(T)):
        Pt=ExpMat[ind-1];n=M[ind];mut=I[ind];
        Lnext=[]
        if verbose:
            print "Pt.shape ",Pt.shape
            print "Pt sum rows ",np.sum(Pt,axis=1)
        for j in range(Hsmall):
            somme =sum(Pt[:,j]*L)
            echant=funcs.sampling(mut=mut,n=n,f=xxSmall[j])
            #echant=1
            Lnext.append(echant*somme)
        L=Lnext
    if verbose:
        print "L\n",L
    
    #----------------------------------------- return values -----------------------------------
    if prob_nofix==0 and prob_nofix==sum(L):
        proba_norm=0        
    else:
        proba_norm=sum(L)/prob_nofix
    if loga:
        if verbose:
            print "Returning -log lik................"
            print "sum(L) ",sum(L)
            print "prob_nofix ",prob_nofix
            print "lik ",proba_norm
        print "OUTPUT: -log lik ",-np.log(proba_norm)
        return -np.log(proba_norm),-np.log(sum(L)),PtOri,ExpMat,prob_nofix

    else:
        if verbose:
            print "Returning lik..................."
            print "not loga "
            print "len(L) ",len(L)
            print "sum(L) ",sum(L)
            print "prob_nofix ",prob_nofix
        print  "OUTPUT: sum(L)/prob_nofix ",proba_norm
        return proba_norm,sum(L),PtOri,ExpMat,prob_nofix

#######################################################################

def ll(params=[None,None,None],data=[None,None,None],H=100,verbose=0,expmethod='alglin',expmethodsecond='pade',loga=0,grid='default',dominance=0.5,threshold=1e-5,numberprecision=128,thresholdfrerot=1e-5):
    """
    Computes the (unconditional) likelihood when given the data.    
    """
    if verbose:
        print "data I,M,T",data
        print "params t0,gamma,Ne",params
        print "H ",H
        raw_input(" Press a key to continue ")    

    #------------- define parameters -------------------------------

    t_=params[0];gamma=params[1];Ne=params[2];

    #-------------- redefine expmethod H ----------------------------

    expmethod,H=funcs.process_expmethod(Ne=Ne,gamma=gamma,dominance=dominance,grid=grid,gammathreshold=40,H=H,expmethod=expmethod,expmethodsecond=expmethodsecond)

    if verbose:
        print "expmethod,H,dominance,verbose ",expmethod,H,dominance,verbose

    if verbose:        
        print "after processing...."
        print "data I,M,T",data
        print "params t0,gamma,Ne",params
        print "H ",H
        #raw_input(" Press a key to continue ")    

    #------------- define data -------------------------------------

    I=data[0];M=data[1];T=data[2]

    #---- check allele age is appropriate, if not return np.nan ----

    I,M,T=preprocess(time=t_,Iori=I,Mori=M,Tori=T,verbose=0)
    if verbose:
        print "Input data after preprocessing I,M,T",I,M,T
        print "input allele age ",t_
    
    if I==None and M==None and T==None:
        if verbose: 
            print "likelihood is",0
            print "OUTPUT: ",np.nan
        return 4*[np.nan]    #Attention! number of outputs.

    # --------Compute all t[j]-t[j-1]/2Ne (store in DeltaT).--------
    
    DeltaT=list((np.array(T)-np.array([t_]+T[:-1]))/(2.0*Ne))
    if verbose: print "DeltaT ",DeltaT

    # ------------Define the space grids ---------------------------
    if grid=='default':
        xx,oriFI=Grids.Ne_grid(H,Ne)   
    elif grid=='symmetric':
        xx=Grids.symmetric_grid(H)
    elif grid=='uniform':
        xx,oriFI=Grids.uniform_grid(H,Ne)
    elif grid=='expo':
        indice=int(0.10*H)
        xx,oriFI=Grids.exponential_grid_Ne(indice=indice,H=H,Ne=Ne)
    
    maxi=max(np.diff(xx))

    # ---------Check gamma is appopropriate -------------------------
    
    #if abs(gamma)>60: 
    #    raise ValueError("gamma is too high or too small for current double precision, max is 60") 
    #elif gamma>=1./maxi or gamma<=-1./maxi:
    #    raise ValueError("gamma is too high or too small for grid space, max is  %f"%(1./maxdiff))    

    #Compute the one and unique Q
    Qmat,Deltamat,Betamat=AlgLin.Qmatrix(gamma=gamma,xx=xx,verbose=0,h=dominance)

    # --------- expqt for all t --------------------------------
    #---------- matrix diagonalization-babik--------------------
  
    if expmethod=='alglin':
        #print "HERE", sys.exit()
        PtMatrices,flag=AlgLin.ExpMatrices_Q(Qmat=Qmat,Betamat=Betamat,Deltamat=Deltamat,Times=DeltaT,verbose=verbose,threshold=threshold) 
        if verbose: 
            print "FLAG: ",flag
        if flag == 'warning':
            print "warning, the exponentiation seems problematic"
        
    #----------- pade approximation ---------------------------
    elif expmethod=="pade":
        PtMatrices=[AlgLin.expm((dt)*Qmat) for dt in DeltaT]
    #---------- high precision call ----------------------------
    elif expmethod=='prec':        
        if grid=="default":
            argu_frerot=[numberprecision,thresholdfrerot,gamma,dominance,H,1,Ne]+DeltaT
            argu_frerot=" ".join([str(w) for w in argu_frerot])
        else:
            raise ValueError("only grid implemented is the default grid")
        command="./ancientSelectionApaBigQ %s"%argu_frerot      
        print command
        os.system(command)       

        PtMatrices=[np.loadtxt('exp_Qt_%i.dat'%dtindex, dtype=np.float64,delimiter=" ") for dtindex in range(len(DeltaT))]

        [os.system("rm exp_Qt_%i.dat"%dtindex) for dtindex in range(len(DeltaT))]

    if len(PtMatrices)!=len(T):
        print len(PtMatrices)
        raise ValueError("something went wrong with matrix computation: %s")
    
    PtOri=PtMatrices[0]    ## I.C.
    
    ExpMat=PtMatrices[1:]  ## all others    
    #if verbose and expmethod=='babik':
    #    stepfinal=(T[-1]-t_)/(2.0*Ne))
    #    Ptfinal,flag=AlgLin.ExpMatrices_Q(Qmat=Qmat,Betamat=Betamat,Deltamat=Deltamat,Times=[stepfinal],verbose=1)[0]
    #    print "expected  no fix",1-Ptfinal[oriFI,0]-Ptfinal[oriFI,-1]

    #----------- initial conditions ---------------------------
    IC=[]
    for j1 in range(H):
        #if verbose:
        #print 'sampling1 ',sampling(n=M[0],mut=I[0],f=xxOri[j1])
        IC.append(funcs.sampling(n=M[0],mut=I[0],f=xx[j1])*PtOri[oriFI,j1])
    
    #----------- Loop over next times ---------------------------
        
    L=IC
    for ind in range(1,len(T)):
        Pt=ExpMat[ind-1];n=M[ind];mut=I[ind];
        Lnext=[]
        for j in range(H):
            somme =sum(Pt[:,j]*L)
            echant=funcs.sampling(mut=mut,n=n,f=xx[j])
            Lnext.append(echant*somme)
        L=Lnext
    
    if verbose:
        print "L\n",L
        #------------- return values --------------------------

    if loga:
        print "LOGA.............."
        if verbose:
            print "Returning -log lik................"
            print "-log lik sum(L) ",-np.log(sum(L))
            #print "prob_nofix ",prob_nofi

        return -np.log(sum(L)),oriFI,PtOri,ExpMat
    else:
        if verbose: 

            print "Returning lik..................."
            print "not loga "
            print "len(L) ",len(L)
            print "lik sum(L) ",sum(L)
            print "lik ",sum(L)

        return sum(L),oriFI,PtOri,ExpMat

#poba_Q,orIF,PtOri,ExpMat=ll(data=data,params=params,H=H,verbose=0,expmethod='frerot',loga=1,grid='default',threshold=1e-5,numberprecision=128,thresholdfrerot=1e-5)
#Proba_q_cond,Proba_q,PtOri,ExpMat,prob_nofix=ll_q(data=data,params=params,H=H,verbose=0,expmethod='frerot',loga=1,grid='default',threshold=1e-5,numberprecision=300,thresholdfrerot=1e-5)
#ll_q(data=data,params=params,H=H,verbose=1,babik=1,loga=0,positive=1)

#ll(data=data,params=params,H=H,verbose=1,babik=1,loga=0)
