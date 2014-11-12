#! /usr/bin/env python

"""
Various necessary grids
"""

import numpy as np
import scipy.optimize

def Ne_grid(num_pts,Ne,verbose=0):
    """
    Modified from Ryan Gutenkunst dadi, v.1.2.3 http://code.google.com/p/dadi/
    to include a point at 1/2Ne

    It is a nonuniform grid of points on [0,1].
    The grid is weighted to be denser near 0 and 1.

    Between 0 and 2/2Ne it is a uniform grid with 1/10 of the total number of points.
    
    Then it smoothly increases like a quadratic function for the remaining 9/10.

    """
    # Rounds down...
    small_pts = int(num_pts / 10)+((int(num_pts / 10)+1)%2)
    if verbose:
        print 1./(2*Ne)
        print "small_pts ",small_pts
    large_pts = num_pts - small_pts+1
    
    #grid = np.linspace(0, 0.05, small_pts)
    grid = np.linspace(0, 2.0/(2*Ne), small_pts)
    startFreq=(small_pts-1)/2
    
    if verbose:
        print "small grid ",grid
        print "startFreq ",startFreq
        print "1/2Ne ",1./(2*Ne),grid[startFreq]
    q = np.linspace(0, 1, large_pts)
    dq = q[1] - q[0]
    x_start = grid[-1]
    dx_start = grid[-1] - grid[-2]

    d = x_start
    c = dx_start/dq
    b = -3*(-dq + dx_start + dq*x_start)/dq
    a = -2 * b/3
    grid = np.concatenate((grid[:-1], a*q**3 + b*q**2 + c*q + d))
    
    if verbose:
        print 'a,b,c ',a,b,c
    #starts=1/(2.0*Ne)*np.ones(num_pts)
    #mini=min(abs(grid-starts))
    #startFreq=list(abs(grid-starts)).index(mini)
    #np.where(grid==1/(2.*Ne))
    #startFreq=np.where(np.around(grid,10)==np.around(1/(2.*Ne),10))[0][0]
    return grid,startFreq

def findbeta(beta,indice,H,Ne,verbose=0):
    if (H-1)%indice!=0:
        raise ValueError("pick another index "+str(H-1)+"mod"+str(indice)+" = "+str((H-1)%indice))

    alpha=-1+2.*indice/(H-1)
    if verbose: print "alpha ",alpha
    return ((1/(1+np.exp(-beta*alpha))-(1/(1+np.exp(beta))))/((1/(1+np.exp(-beta)))-(1/(1+np.exp(beta)))))-1./(2*Ne)


def exponential_grid_Ne(indice=40,H=401,Ne=10000):
    """

    Modified from Ryan Gutenkunst dadi, v.1.5.1 http://code.google.com/p/dadi/
    modified to include a point 1/2Ne.

    An exponentially spaced grid with parameter beta. Denser around 0 and 1.
    
    indice: is the index at wich xx[indice]=1/2Ne, where xx is the grid
    H: grid size (grids includes 0 and 1, so H points total)
    Ne: pop size

    Uses scipy's bentq to find for which beta  xx[indice]=1/2Ne, where
    beta "controls the degree to which grid points crowd against x=0 or x=1."    
    In dadi beta is called crwd.  
  
    """
    args=(indice,H,Ne)
    leftbound=0.000001
    rightbound=1000
    beta,r=scipy.optimize.brentq(findbeta,leftbound,rightbound,args=args,full_output=True)

    if r.flag!='converged':
        raise ValueError('brentq has not converged')

    unif = np.linspace(-1,1,H)
    grid = 1./(1. + np.exp(-beta*unif))

    # Normalize
    grid = (grid-grid[0])/(grid[-1]-grid[0])

    return grid,indice


################## OLD or BIBLIOGRAPHY  ##########################################

def exponential_grid(pts, crwd=8.):
    """
    Taken from dadi 1.5.1
    An exponentially spaced grid. This is now the default grid.

    crwd controls the degree to which grid points crowd against x=0 or x=1.
    The value of crwd=8 seems to be a good default for integration with large
    systems. See estimate_best_exp_grid_crwd for some empirical optimizations
    beyond this.

    This grid was contributed by Simon Gravel.
    """
    unif = np.linspace(-1,1,pts)
    grid = 1./(1. + np.exp(-crwd*unif))

    # Normalize
    grid = (grid-grid[0])/(grid[-1]-grid[0])

    return grid

def uniform_grid(num_points,Ne):
    """
    This is a uniform grid with an added point at 1/2Ne
    Return the grid and the index of that 1/2Ne    
    """
    startF=1/(2.*Ne)
    
    if num_points==(2*Ne+1):
        grid=np.linspace(0,1,num_points)
        startFI=1
        return grid,startFI

    else:
        grid=np.linspace(0,1,num_points-1)
        if startF in grid:
            raise StandardError("you (i mean I) did something stupid...")
        grid=list(grid)
        grid.append(startF)    
        grid.sort()
        startFI=grid.index(startF)
        grid=np.array(grid)

        return grid,startFI
    
   
    
def default_grid(num_pts):
    """
    Taken from dadi, as is.
    A nonuniform grid of points on [0,1].

    The grid is weighted to be denser near 0 and 1, which is useful for
    population genetic simulations. In between, it smoothly increases and
    then decreases the step size.
    """
    # Rounds down...
    small_pts = int(num_pts / 10)
    #small_pts=2
    large_pts = num_pts - small_pts+1
    #Ne=10000.0
    grid = np.linspace(0, 0.05, small_pts)
    #grid = np.linspace(0, 2.0/(2*Ne), small_pts)
    print grid
    #raw_input()
    # The interior grid spacings are a quadratic function of x, being
    # approximately x1 at the boundaries. To achieve this, our actual grid
    # positions must come from a cubic function.
    # Here I calculate and apply the appropriate conditions:
    #   x[q = 0] = x_start  =>  d = x_start
    #   dx/dq [q = 0] = dx_start => c = dx_start/dq
    #   x[q = 1] = 1
    #   dx/dq [q = 1] = dx_start
    
    q = np.linspace(0, 1, large_pts)
    dq = q[1] - q[0]
    x_start = grid[-1]
    dx_start = grid[-1] - grid[-2]

    d = x_start
    c = dx_start/dq
    b = -3*(-dq + dx_start + dq*x_start)/dq
    a = -2 * b/3
    grid = np.concatenate((grid[:-1], a*q**3 + b*q**2 + c*q + d))

    return grid





def symmetric_grid(num_pts):
    """
    A nonuniform grid of points on [0,1].

    The grid is weighted to be denser near 0 and 1, which is useful for
    population genetic simulations. In between, it smoothly increases and
    then decreases the step size.
    """
    # Rounds down...
    
    # The interior grid spacings are a quadratic function of x, being
    # approximately x1 at the boundaries. To achieve this, our actual grid
    # positions must come from a cubic function.
    # Here I calculate and apply the appropriate conditions:
    #   x[q = 0] = x_start  =>  d = x_start
    #   dx/dq [q = 0] = dx_start => c = dx_start/dq
    #   x[q = 1] = 1
    #   dx/dq [q = 1] = dx_start
    
    q = np.linspace(0, 1, num_pts)
    
    d=0
    
    c = 1/(10.0*num_pts)
    b = 3-(3.0/(10*num_pts))
    a = -2+1/(5.0*num_pts)    
    
    
    grid = (a*q**3 + b*q**2 + c*q + d)

    return grid


def selection_grid(num_pts):
    """
    A nonuniform grid of points on [0,1].

    The grid is weighted to be denser near 0 and 1, which is useful for
    population genetic simulations. In between, it smoothly increases and
    then decreases the step size.
    """
    ### for -0.1 to 0.1
    if num_pts%2==0:
        num_pts+=1

    q = np.linspace(0, 1, num_pts)
    
    a=0.02;b=-0.03;c=0.03;d=-0.01;

    #a=-2   
    #b = 3
    #c = 1
    #d = -1
       
    
    grid = (a*q**3 + b*q**2 + c*q + d)

    return grid

def gamma_grid(num_pts):
    """
    A nonuniform grid of points on [0,1].

    The grid is weighted to be denser near 0 and 1, which is useful for
    population genetic simulations. In between, it smoothly increases and
    then decreases the step size.
    """
    ### for -0.1 to 0.1
    if num_pts%2==0:
        num_pts+=1

    q = np.linspace(0, 1, num_pts)
    
    a=20.;b=-30.;c=30.;d=-10.;
    a=100.;b=-150.;c=150.;d=-50.;
    #a=-2   
    #b = 3
    #c = 1
    #d = -1       
    
    grid = (a*q**3 + b*q**2 + c*q + d)

    return grid

def effective_grid(mini=100,maxi=10000000):
    effgrid=[mini]
    while 10*effgrid[-1]<=maxi:
        effgrid.append(10*effgrid[-1])

    #print effgrid
    return effgrid
