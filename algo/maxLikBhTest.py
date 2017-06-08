 # -*- coding: utf-8 -*-
import numpy as np
from mpi4py import MPI
from scipy.linalg import expm
try:
    import algo.BurstSearch as BurstSearch
    import algo.BGrate as BGrate
    import algo.fretAndS as fretAndS
except ImportError:
    import BurstSearch
    import BGrate
    import fretAndS
import datetime,sys
from scipy.optimize import basinhopping
from array import array
from scipy.linalg import expm
#from mpiBurstLikelihood import calcBurstLikelihood



class bhSteps(object):
    def __init__(self, bounds, stepsize=1):
        self.bounds = bounds
        self.xs=len(bounds)
        if stepsize<0:
            stepsize=stepsize*-1
        if stepsize>1:
            stepsize=1/stepsize
        self.stepsizex=stepsize
    def __call__(self, x):
        print(self.stepsizex)
        print('=============')
        #print('old x:'+str(x))
        for idx in range(self.xs):
            xmax=(self.bounds[idx][1]-x[idx])*self.stepsizex
            xmin=(self.bounds[idx][0]-x[idx])*self.stepsizex
            print(x[idx]+xmin,x[idx]+xmax)
            x[idx]=x[idx]+np.random.uniform(xmin,xmax)
        #print('new x:'+str(x))
        return x



def func2d(x):
    f = np.cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + \
                                                            0.2) * x[0]
    df = np.zeros(2)
    df[0] = -14.5 * np.sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2
    df[1] = 2. * x[1] + 0.2
    return f, df

minimizer_kwargs = {"method":"L-BFGS-B", "jac":True}
x0 = [1.0, 1.0]

n_states=2
boundE=[(-8.15,8.999)]*n_states
boundK=[(0.1,100)]*(n_states*(n_states-1))
bound=boundE#+boundK
bhStep=bhSteps(bound,1)
print(bhStep([0.5,0.6,1,99]))
ret = basinhopping(func2d, x0, minimizer_kwargs=minimizer_kwargs,niter=200,take_step=bhStep)
print("global minimum: x = [%.4f, %.4f], f(x0) = %.4f" % (ret.x[0], ret.x[1],ret.fun))
