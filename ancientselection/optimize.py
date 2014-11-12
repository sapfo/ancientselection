#!/usr/bin/env python

import sys
import scipy.optimize
import numpy as np
#from Inference_new import ll_q 
import inference
#from time import time
import time


_out_of_bounds_val=np.nan

#: Storage for times at which each stream was flushed
__times_last_flushed = {}

def delayed_flush(stream=sys.stdout, delay=1):
    """
    COPIED FROM Ryan Gutengkunst' Misc.py function

    Flush a stream, ensuring that it is only flushed every 'delay' *minutes*.
    Note that upon the first call to this method, the stream is not flushed.

    stream: The stream to flush. For this to work with simple 'print'
            statements, the stream should be sys.stdout.
    delay: Minimum time *in minutes* between flushes.

    This function is useful to prevent I/O overload on the cluster.
    """
    global __times_last_flushed

    curr_time = time.time()
    # If this is the first time this method has been called with this stream,
    # we need to fill in the times_last_flushed dict. setdefault will do this
    # without overwriting any entry that may be there already.
    if stream not in __times_last_flushed:
        __times_last_flushed[stream] = curr_time
    last_flushed = __times_last_flushed[stream]

    # Note that time.time() returns values in seconds, hence the factor of 60.
    if (curr_time - last_flushed) >= delay*60:
        stream.flush()
        __times_last_flushed[stream] = curr_time

def _project_params_down(pin, fixed_params):
    """
    COPIED FROM Ryan Gutengkunst'Optimize function
    Eliminate fixed parameters from paramater input
    """
    if fixed_params is None:
        return pin

    if len(pin) != len(fixed_params):
        print pin,fixed_params
        raise ValueError('fixed_params list must have same length as input '
                         'parameter array.')

    pout = []
    for ii, (curr_val,fixed_val) in enumerate(zip(pin, fixed_params)):
        if fixed_val is None:
            pout.append(curr_val)

    return np.array(pout)


def _project_params_up(pin, fixed_params):
    """    
    COPIED FROM Ryan Gutengkunst'Optimize function
    Fold fixed parameters into pin.
    """
    if fixed_params is None:
        return pin

    pout = np.zeros(len(fixed_params))
    orig_ii = 0
    for out_ii, val in enumerate(fixed_params):
        if val is None:
            pout[out_ii] = pin[orig_ii]
            orig_ii += 1
        else:
            pout[out_ii] = fixed_params[out_ii]
    return pout

#: Counts calls to object_func
_counter = 0

def _object_func(params, model_func,  
                 lower_bound_obj=None, upper_bound_obj=None, 
                 verbose_obj=1,flush_delay=0, fixed_params=None, ll_scale=1,                 
                 args_ll={}
                 ):
    """
    Objective function for optimization.
    params has to be the first argument, that is the key for the optimization. 
    Moreover uses a clever way to increase or decrease the number of fixed paramters, 
    modified from Ryan Gutenkunst...do something about that...
    """
    global _counter
    _counter += 1

    # Deal with fixed parameters
    params = _project_params_up(params, fixed_params)

    # Returns an _out_of_bounds value/ll_scale if algo tries something out of the specified bounds 
    # This is for algo that are not dealing explicitely with bounds, e.g. Nelder Mead
    # I have an _out_of_bounds value of np.nan and it seems to work

    if (lower_bound_obj is not None) and (lower_bound_obj is not (None,None,None)):
        for pval,bound in zip(params, lower_bound_obj):
            if bound is not None and pval < bound:
                print "Exit because out of lower bound "
                return _out_of_bounds_val/ll_scale

    if (upper_bound_obj is not None) and (upper_bound_obj is not (None,None,None)):
        for pval,bound in zip(params, upper_bound_obj):
            if bound is not None and pval > bound:
                print "Exit because out of upper bound "
                return _out_of_bounds_val/ll_scale

#    ns = data.sample_sizes 
    #print 'params before ',params


    #print all_args
    #print "params ",params
    
    #if verbose_obj:
    print "params: ",params
    print "Fixed_params: ",fixed_params
    
    #print "all args ",all_args
    #all_args = [params, H] + list(func_args) + [H]
    #print "all args ",all_args
    #print "all arguments ",all_args

    result = model_func(params,**args_ll)
    result = result[0]
    if verbose_obj: print "result -----------------", result
    # Bad result
    if np.isnan(result):
        result = _out_of_bounds_val

    if (verbose_obj > 0) and (_counter % verbose_obj == 0):
        param_str = 'array([%s])' % (', '.join(['%- 12g'%v for v in params]))
        print 'intermediate result: %-8i, %-12g, %s' % (_counter, result, param_str)
        delayed_flush(delay=flush_delay)

    return result/ll_scale

                 #params, data, model_func, H=None, 
                 #lower_bound=None, upper_bound=None, 
                 #verbosefunc=0, babik=True,loga=True,grid='ori',positive=True, 
                 #flush_delay=0, fixed_params=None, ll_scale=1

def timethis(func=None,data=None,params=None,H=None,dominance=None,nrep=1,verbose=1,verboselik=0):
    starttime=time.time()
    if verbose:
        print "### you are in time this ###"
    func(data=data,params=params,H=H,dominance=dominance,loga=1,grid='default',positive=1,verbose=verboselik)        
    endtime=time.time()
    exec_time=endtime-starttime
    if verbose:
        print "startime, endtime, exec_time",starttime, endtime, exec_time    
    return exec_time/(1.*nrep)


def exhaustive_search(model_func,data=(None,None,None),H=None,Grid=None,running_time=None,Bounds=None,fixed_params=(None,None,None),dominance=0.5,verboseexhaustive=0,verboselik=0,**otherargs):
    
#p0, data, model_func, H, lower_bound=None, upper_bound=None,
#           verbose=0, flush_delay=0.5, epsilon=1e-3, 
#           pgtol=1e-5, babik=True, maxiter=None, full_output=False,
#           func_args=[], fixed_params=None, ll_scale=1
    """
    INPUT:
    - func     - data    - either Grid and runningtime==None
                         - or Grid==None and runningtime (in seconds)
    (Check Grid==None or runningtime=None )
    - if Grid!=None:
        need NED=Grid
             TOD=Grid
             GammaD=Grid
             Grid=NE,T0,Gamma, so a cube.
    - else:
         need upper bound lower bound NE, Gamma, 
         have two outside def, one for running time, one for intervals

    OUTPUT:
    - Popt 
    - Whole grid of values
    """
    #model_func,Grid=None,running_time=None,lower_bound=None,upper_bound=None,fixed_params=None,verbose=0
    if (Grid,running_time)==(None,None):
        raise ValueError("need to input either grid or running time for exhaustive search")
    if (Grid!=None and running_time!=None):
        raise ValueError("need either Grid or running time to be None")
    if len(fixed_params)!=3:
        raise ValueError("need 3 values for the fixed param even if None,None,None")

    I,M,T=data
    #time one run:
    params_timing=(T[0]-100,10,400)
    print "timing one call of the function....for params (hard coded): ",params_timing
    timeforonerun=timethis(nrep=1,func=model_func,data=data,params=params_timing,H=H,dominance=dominance,verboselik=verboselik,**otherargs)
    
    for ii,fixed_val in enumerate(fixed_params):   
            if fixed_val is not None:
                Bounds[ii]=(fixed_val,fixed_val)

    BoundsT0,BoundsGamma,BoundsNE=Bounds   
    if verboseexhaustive:
        print "BoundsT0 ",BoundsT0
        print "BoundsGamma ",BoundsGamma
        print "BoundsNE ",BoundsNE

    if verboseexhaustive:
        print "running time for one function call",timeforonerun
    if Grid:
    
        T0D,GammaD,NED=Grid

        T0=np.unique(np.linspace(BoundsT0[0],BoundsT0[1],T0D))
        Gamma=np.unique(np.linspace(BoundsGamma[0],BoundsGamma[1],GammaD))
        NE=np.unique(np.linspace(BoundsNE[0],BoundsNE[1],NED)
)
        dims=(T0.shape[0],Gamma.shape[0],NE.shape[0])

        
        #if verboseexhaustive:
        print "Expect a total time of (rough approximate): ",timeforonerun*dims[0]*dims[1]*dims[2],"[sec]"
        ncalls = dims[0]*dims[1]*dims[2]

            #raw_input("press enter of ok")
    #define the function, the params for timing can be anython
    if running_time:
        division=1.*fixed_params.count(None)
        if verboseexhaustive: print "division ",division
        
        pts=int((1.*running_time/timeforonerun)**(1/division))
        print "(1.*running_time/timeforonerun) ",(1.*running_time/timeforonerun)**(1/3.)
        if verboseexhaustive:
            print "total allowed running time ",running_time
            print "number of points per direction ",pts
       
        #### deal with the fixed params
        for ii,fixed_val in enumerate(fixed_params):   
            if fixed_val is not None:
                Bounds[ii]=(fixed_val,fixed_val)

        BoundsT0,BoundsGamma,BoundsNE=Bounds   
        if verboseexhaustive:
            print "BoundsT0 ",BoundsT0
            print "BoundsGamma ",BoundsGamma
            print "BoundsNE ",BoundsNE


        ####
        T0=np.unique(np.linspace(BoundsT0[0],BoundsT0[1],pts))
        Gamma=np.unique(np.linspace(BoundsGamma[0],BoundsGamma[1],pts))
        NE=np.unique(np.linspace(BoundsNE[0],BoundsNE[1],pts))
        dims=(T0.shape[0],Gamma.shape[0],NE.shape[0])
        if verboseexhaustive:
            print "T0: ",T0,"\nGamma: ",Gamma,"\nNE: ",NE
            print "dims ",dims
        print "ncalls to be run: ",dims[0]*dims[1]*dims[2]
        #raw_input("check dimensions and press key")

        #Grid=np.array(T0,Gamma,NE)
    
    Lik_array=np.zeros(dims)
    where=-1
    for j_t_ind,j_t in enumerate(T0):
        for j_g_ind,j_g in enumerate(Gamma):
            for j_N_ind,j_N in enumerate(NE):
                where+=1
                print "------------------------------------"
                if verboseexhaustive:
                    
                    if int(100.*where/ncalls)%10==0:
                        print "We are at: ",100.*were/ncalls,"% of run ----------------------------- "
                        #raw_input()
                t0=j_t;gamma=j_g;Ne=j_N;  
                params=t0,gamma,Ne
                #if verboseexhaustive:

                print "Current params: ",params
#               print "data ",data

                #lik,oriFI,PtOri,ExpMat=model_func(params,data,H=H,dominance=dominance)
#               print "just before dominance",dominance
                lik,summation,PtOri,ExpMat,prob_nofixlik=model_func(params,data,H=H,dominance=dominance,verbose=verboselik)
                                
                #raw_input("again ")
                Lik_array[j_t_ind,j_g_ind,j_N_ind]=lik

    #find the maximum **** some kind of argmax stuff
    argmax=scipy.argmin(Lik_array,axis=None)
    if verboseexhaustive: print 'argmin ',argmax
    idmin = np.unravel_index(argmax, dims)
    if verboseexhaustive: print 'idmin ',idmin    
    idmin_t0,idmin_gamma,idmin_Ne=idmin
    minlik=Lik_array[idmin]
    if verboseexhaustive: print 'min lik',minlik

    Popt=T0[idmin_t0],Gamma[idmin_gamma],NE[idmin_Ne]

    return Popt,minlik,Lik_array,T0,Gamma,NE,T 

def saveoutput(Lik_array=None,T0=None,Gamma=None,NE=None,Tinput=None,project='',directory='/home/sapfo/Projects/AncientSelection/Output_Arrays',verbose=1):
    if Lik_array!=None:        
        name="Lik_array_%s"%project
        if verbose:
            print "Saving Lik array under the dir ",directory,"with name ",name
        np.save("%s/%s"%(directory,name),Lik_array)

    if T0!=None:
        name="T0_%s"%project
        if verbose:
            print "Saving T0 under the dir ",directory,"with name ",name
        np.save("%s/%s"%(directory,name),T0)

    if Gamma!=None:
        name="Gamma_%s"%project
        if verbose:
            print "Saving Gamma under the dir ",directory,"with name ",name
        np.save("%s/%s"%(directory,name),Gamma)

    if NE!=None:
        name="NE_%s"%project
        if verbose:
            print "Saving NE under the dir ",directory,"with name ",name
        np.save("%s/%s"%(directory,name),NE)
    if Tinput!=None:
        name="Tinput_%s"%project
        if verbose:
            print "Saving Tinput under the dir ",directory,"with name ",name
        np.save("%s/%s"%(directory,name),Tinput)

 
    return None


def neldermead(model_func=None,p0=(None,None,None),
               fixed_params=(None,None,None),
               verbose_obj=1,flush_delay=0.5,full_output_here=1,
               ll_scale=1,lower_bound_obj=None,upper_bound_obj=None,
               args_nelder=dict(xtol=0.0001,ftol=0.0001, 
               maxiter=100, maxfun=20, disp=1, retall=1, 
               callback=None),         
               args_ll={}):


    '''
    Nelder-Mead Simplex algorithm
    scipy.optimize.fmin
    '''

    args = (model_func, lower_bound_obj, upper_bound_obj, verbose_obj, flush_delay, fixed_params, ll_scale,args_ll)

    #print args
    #raw_input()

    #print args
    #print verbose_obj

    p0 = _project_params_down(p0, fixed_params)
    
    outputs=scipy.optimize.fmin(_object_func, p0, args=args,**args_nelder)
    print "OUTPUTS: ", outputs
    #need full_output and retall of fmin set to 1

    xopt,fopt,iteration,funcalls,warnflag,allvecs=outputs

    xopt = _project_params_up(xopt, fixed_params)

    if not full_output_here:
        return xopt
    else:
        return xopt,fopt,iteration,funcalls,warnflag,allvecs
    
#func_ex=inference.ll_q

#dico_ll=dict(data=(I,M,T),H=H,verbose=0, expmethod='babik',
#          expmethodsecond='pade',loga=1,grid='default',
#          positive=1,dominance=0.5,threshold=1e-5,
#          numberprecision=128,thresholdfrerot=1e-5)

#dico_scipy_bfgs=dict(approx_grad=True,
#                     m=10,
#                     factr=10000000.0,
#                     pgtol=1e-05,
#                     epsilon=1e-08,
#                     iprint=3, 
#                     maxfun=15000)

#dico_nelder=dict(xtol=0.0001, ftol=0.0001, 
#                 maxiter=maxiterbase, maxfun=maxfunbase, full_output=1, 
#                 disp=1, retall=1, callback=None)

#outputs=neldermead(model_func=func_ex,p0=params,fixed_params=[-1000,None,500],             
#                   full_output_here=True, flush_delay=0.5,verbose_obj=1,
#                   lower_bound_obj=(-1000,0,500), upper_bound_obj=(-1000,100,500),
#                   args_nelder=dico_nelder,
#                   args_ll=dico_ll)

#popt = optimize(model_func=func_ex,p0=params, 
#             fixed_params=(-1000,None,500),  
#             lower_bound=(-1000,0,500), upper_bound=(-1000,100,500),
#             flush_delay=0.5,verbose_obj=2,
#             ll_scale=1,             
#             full_output=True,
#             args_scipy_bfgs=dico_scipy_bfgs,   
#             args_ll=dico_ll)

