#!/usr/bin/env python

import sys,os,random,errno
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

# ---------------- Arguments
parser = argparse.ArgumentParser(description='''Take as input time serial data (one locus). The output can be either mles or likelihood surface over a grid or the mles'
An example of a run: ./ancientselection.py -i TestData_sel0_Jeff.py --dirout Out --codedir ../ancientselection/ --run run1 --exhaust''',formatter_class=RawTextHelpFormatter)
#./ancientselection.py -i TestData_sel0_Jeff.py --dirout Out --codedir ../ancientselection/ --run run1 --exhaust cube

#parser = argparse.ArgumentParser(description='''Take as input time serial data. The output can be either a likelihood values for a grid or the mles'
#An example of a run: ./main_optimize_nelder.py TestData_Jeff.py 400 Oct4run1''',formatter_class=argparse.ArgumentDefaultsHelpFormatter)


#parser.add_argument('--mikefile', help='mikeinfile',type=argparse.FileType('r'),required=True)
#parser.add_argument('--toselect', help='list of indivs to select from mikefile, one per line',type=str,required=True)
#parser.add_argument('--toorder', help='list of pops to reorder indivs from mikefile, one per line',type=str,required=True)
#parser.add_argument('--mikenew', help='mike formatted output file',type=argparse.FileType('wa'),required=True)

helpdatafile = '''datafile is a python script that contains
the data in the following format

For the data itself, here is an example:
M_ = [10,10,10] #for 10 chrom at 3 time points
I_ = [3,3,3] #for 3 derived alleles at each time point
T_ = [-100,-50,0] #the time points in generations 
dominance_ = 0.5 #if codominance 

where; 
 M_: python list with the total number of chromosomes 
 I_: python list wiyth the number of derived alleles
 T_: python list with the sampling times generations
 dominance_: float with the  dominance coefficient for the data (usually 0, 0.5 or 1)


For the parameters:
e.g. 
Upper_bounds_ = [0,10,1000] #(t0_up,gamma_up,Ne_up)
Lower_bounds_ = [-150,-10,500] #(t0_low,gamma_low,Ne_low)
fixed_params_ = [None,None,1000] #for t0 and gamma to be free while the pop size is set to 1000.

where:
 Upper_bounds_ = python list (t0_up,gamma_up,Ne_up) upper bounds for 
                 t0, gamma and Ne (in this order!!)
 Lower_bounds_ = python list (t0_low,gamma_low,Ne_low) lower bounds for 
                 t0, gamma and Ne (in this order!!)
 fixed_params = python list indicating which parameters should be fixed 
                (same order: t0,gamma,Ne). The value 
                is set to None if the parameters is not to be fixed or to the value it 
                should be fixed at. The fixed values should be compatible with the bounds.

'''

parser.add_argument('--version','-v', action='version', version='%(prog)s 0.0')


parser.add_argument('--datafile','-i', help=helpdatafile,type=str,required=True)

parser.add_argument('--run','-r', help='''added string to the project name,
only used to label output files (by default it is the datafile name minus '.py' extension)''',type=str,required=False,default='')

parser.add_argument('--dirout', help='directory out (if does not exist, will be created) --default Out',type=str,required=False, default='Out')

parser.add_argument('--codedir', help='directory where the code lives (to be added to your path)',type=str,required=True,default="../bin")

parser.add_argument('--gridsize', help='size of the grid (H) -- default 400',type=int,required=False,default=400)

parser.add_argument('--gridtype', help='type of grid, either of (default,symmetric,uniform,expo) --default default',type=str,required=False,default='default')

parser.add_argument('--expmethod1', help='''exponential method 1, used always if gamma small enough,\neither of (alglin,pade,prec)
alglin: in detail in the paper
pade: implemented in scipy
prec: arbitrary precision, the grid has to be the default grid: !!not ready yet!!!
 --default alglin''',type=str,required=False,default='alglin')

parser.add_argument('--expmethod2', help='''exponential method 2 (see above)
used for large abs(gamma),\neither of (pade,prec)
prec: not implemented yet
 -- default pade''',type=str,required=False,default='pade')

parser.add_argument('--exhaust', help='''computes the likelihood on a grid: 
either cube or predefinite running time 
(usage --exhaust cube or --exhaust time)
Note: if --exhaust not specified will try 
to find the maximum likelihood using a 
nelder-mead algorithm -- default cube)''',required=False,default=False)

parser.add_argument('--T0dim', help='number of evaluations for the age (default 5), only in use if --exhaust cube',type=int,required=False,default=5)
parser.add_argument('--Gammadim', help='number of evaluations for gamma  (default 5), only in use if --exhaust cube',type=int,required=False,default=5)
parser.add_argument('--NEdim', help='number of evaluations for Ne (default 5), only in use if --exhaust cube',type=int,required=False,default=5)

parser.add_argument('--runningtime', help='''only in use if --cube time, 
you can specify how long you want it to run. 
The number of points per paramaters will be the same (if not fixed) -- default 300''',type=int,required=False,default=5*60)

parser.add_argument('--nonconditional', help='''likelihood either conditional (default) on allele segragating at the last sampling time 
-- default is to condition, i.e. without --nonconditional flag''',action='store_true',required=False,default=False)

parser.add_argument('--verbose', help='increase standard out',action='store_true',required=False,default=False)
parser.add_argument('--debug', help='debug, lots of standard out',action='store_true',required=False,default=False)

args = parser.parse_args()

# parse all arguments
verbose = args.verbose
debug = args.debug

datafile = args.datafile
run = args.run
project = datafile.split('.py')[0]+'_'+run
dirout = args.dirout
codedir = args.codedir #(to add to the path)
H = args.gridsize
gridtype = args.gridtype
expmethod1 = args.expmethod1
expmethod2 = args.expmethod2
exhaust = args.exhaust
runningtime = args.runningtime
T0dim = args.T0dim
Gammadim = args.Gammadim
NEdim = args.NEdim


nonconditional = args.nonconditional

#parse datafile
execfile(datafile)

#append path etc.
sys.path.append(codedir)
import inference ## ll function
import optimize  ## nelder mead and exhaustive
import funcs ##for the domain definition

#create output directory 
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
	
if not os.path.isdir(dirout):
    if verbose:
        print "your outputdir does not exist, created now"
    mkdir_p(dirout)

#"check" the datafile
try:
    data=(I_,M_,T_)
    dominance_
except:
    print "Something is missing in your datafile!"
    print "Check that you have all variables\nI_,M_,T_,dominance_\nUpper_bounds_,Lower_bounds_,fixed_params_ defined in your datafile"    
    sys.exit()

#check the mode
if exhaust and exhaust not in ['cube','time']:
    print "--exhaust can only be followed by cube or time"
    sys.exit()
# print output

if verbose:
    print "Datafile:\t",datafile
    print "project name:\t",project

    print "Data:"
    print "M:\t",M_
    print "I:\t",I_
    print "T:\t",T_
    print "dominance:\t",dominance_

    print "Parameters:"
    print "Upper Bounds:\t",Upper_bounds_
    print "Lower Bounds:\t",Lower_bounds_
    print "Fixed params:\t",fixed_params_

    print "Grid size (H):\t",H    
    print "Grid type:\t",gridtype
    
    print "expmethod1: ",expmethod1
    print "expmethod2: ",expmethod2
    
    if nonconditional == True:
        print "Will use the unconditional process (Q matrix)..."
    else:
        print "Will use the conditional matrix (q matrix)..."
    if exhaust:
        print "exhaustive mode"
    else:
        print "optimization mode (nelder-mead)"
    if debug: print "Debugging ..."
#further checks: same size of M, I, T

if len(set([len(M_),len(T_),len(I_)]))!=1:
    print "M_,I_,T_: all have to be the same length!!"
    print "ckeck: "
    print M_,I_,T_
    sys.exit()

if len(set([len(Upper_bounds_),len(Lower_bounds_),len(fixed_params_)]))!=1:
    print "Upper_bounds_,Lower_bounds_,fixed_params_: all have to be the same length!!"
    print "ckeck: "
    print Upper_bounds_,Lower_bounds_,fixed_params_
    sys.exit()


t0_low,gamma_low,Ne_low=Lower_bounds_
t0_up,gamma_up,Ne_up=Upper_bounds_

if nonconditional:
    ll_func = inference.ll
else:
    ll_func = inference.ll_q
    

dico_ll=dict(data=data,H=H,verbose=debug, expmethod=expmethod1,
          expmethodsecond=expmethod2,loga=1,grid='default',
          positive=1,dominance=dominance_,threshold=1e-4,
          numberprecision=128,thresholdfrerot=1e-5)

### compute the likelihood on one set of values (make sure dictionaries have all they need)


if debug:
    
    p0gamma=0.01
    p0Ne=random.uniform(Ne_low,Ne_up)
    p0t=random.uniform(t0_low,T_[0])

    pinput = [p0t,p0gamma,p0Ne]
    print "Testing that the ll function is well defined..."
    print "pinput: ",pinput
    trial_ll=ll_func(pinput,**dico_ll)
    #print "trial_ll ",trial_ll

if debug:
    print "debug ",debug

####
if not exhaust:
    Popt=[]
    P0=[]
    Likopt=[]
    Warnflags=[]

    Domains=funcs.domains(data = data,smallert = t0_low,fixed_time = fixed_params_[0])
    dico_nelder=dict(xtol=0.0001, ftol=0.0001, 
                 maxiter=1000, maxfun=2000, full_output=1, 
                 disp=1, retall=1, callback=None)
    
    
    if verbose:
        print "Trying to find mles for each domain"

    for count,domt in enumerate(Domains):
        count+=1
        
        #raw_input()
        p0t=random.uniform(domt[0],domt[1])
        p0gamma=random.uniform(gamma_low,gamma_up)
        p0Ne=random.uniform(Ne_low,Ne_up)
        
        p0=[p0t,p0gamma,p0Ne]

#        print "initial values: ",p0
        P0.append(p0)
        
        lower_bound=domt[0],gamma_low,Ne_low
        upper_bound=domt[1],gamma_up,Ne_up

        if verbose:
            print "Domain number: ",count,", domt: ",domt
            print "Starting values: ",p0
            print "current lower bound ",lower_bound
            print "current upper bound ",upper_bound
            print "fixed_params ",fixed_params_        
        

        func_ex=ll_func

        #currently default maxiter=100, maxfun=20
        xopt,fopt,iteration,funcalls,warnflag,allvecs = optimize.neldermead(            
            model_func=func_ex,p0=p0,fixed_params=fixed_params_,             
            flush_delay=0.5,verbose_obj=debug,
            lower_bound_obj=lower_bound, upper_bound_obj=upper_bound,
            args_nelder=dico_nelder,
            args_ll=dico_ll)
    
        if verbose: print 'popt: ',xopt

        Popt.append(xopt)
        Likopt.append(fopt)
        Warnflags.append(warnflag)

    #print results to the screen
    print "Run %s finished!!"%run
    print "Starting parameters:\n",P0
    print "MLEs:\n",Popt    
    print "Maximum likelihood:\n",Likopt
    print "Warnflags (warnflag!=0, neldermead did not converge):\n",Warnflags
    Likopt=np.array(Likopt)
    
    np.save("%s/Pstart_%s"%(dirout,project),P0)
    np.save("%s/Popt_%s"%(dirout,project),Popt)
    np.save("%s/Likopt_%s"%(dirout,project),Likopt)
    np.save("%s/Warnflags_%s"%(dirout,project),Warnflags)
 
elif exhaust:
    if verbose:
        print "Starting exhaustive computation ..."
    
    if exhaust=='cube':
        
        #Define the cube
        Grid = T0dim,Gammadim,NEdim
        print "Computation on a grid, where each T0 Gamma and NE have dims: ",Grid

        Popt,maxlik,Lik_array,T0,Gamma,NE,Tinput=optimize.exhaustive_search(ll_func,data=data,H=H,Grid=Grid,running_time=None,Bounds=[(t0_low,t0_up),(gamma_low,gamma_up),(Ne_low,Ne_up)],fixed_params=fixed_params_,verbose=debug,dominance=dominance_,verboselik=debug)
    
    elif exhaust=='time':
       
        Popt,maxlik,Lik_array,T0,Gamma,NE,Tinput=optimize.exhaustive_search(ll_func,data=data,H=H,Grid=None,running_time=runningtime,Bounds=[(t0_low,t0_up),(gamma_low,gamma_up),(Ne_low,Ne_up)],fixed_params=fixed_params_,verbose=debug,dominance=dominance_,verboselik=debug)
    
    optimize.saveoutput(Lik_array=Lik_array,T0=T0,Gamma=Gamma,NE=NE,Tinput=Tinput,project=project,directory=dirout,verbose=1)
    
    print "-------------------------------"
    print "Run %s finished!!"%run
    print "MLEs (over exhaustive search):\n",Popt    
    print "Maximum likelihood (exhaustive search):\n",maxlik


    np.save("%s/Popt_%s"%(dirout,project),Popt)
    np.save("%s/Likopt_%s"%(dirout,project),maxlik)
 
    if verbose:
        print "Grid dimensions: ",Lik_array.shape
        print "Popt: ",Popt
        print "maxlik: ",maxlik

    

    
