# from lmfit.minimizer import Minimizer
# def mylorentzian(x, amp, cen, wid):
#     "lorentzian function: wid = half-width at half-max"
#     return (amp / (1 + ((x-cen) / wid)**2))

# fitter = Minimizer()
# fitter.asteval.symtable['lorentzian'] = mylorentzian

#!/usr/bin/env python

# <examples/doc_parameters_basic.py>


def myexpr(s,d):
    return s/d
import numpy as np

from lmfit import Minimizer, Parameters, report_fit

# create data to be fitted
x = np.linspace(0, 15, 301)
data = (5. * np.sin(2*x - 0.1) * np.exp(-x*x*0.01) +
        np.random.normal(size=len(x), scale=0.2))


# define objective function: returns the array to be minimized
def fcn2min(params, x, data):
    """Model a decaying sine wave and subtract data."""
    amp = params['amp']
    shift = params['shift']
    omega = params['omega']
    decay = params['decay']
    model = amp * np.sin(x*omega + shift) * np.exp(-x*x*decay)
    return model - data


# create a set of Parameters
params = Parameters()
params.add('amp', value=10, min=0)
params.add('shift', value=0.0, min=-np.pi/2., max=np.pi/2)
params.add('div',  value=-10, vary=True, min=-20., max=20)


params.add('omega', value=3.0)


# do fit, here with leastsq model
minner = Minimizer(fcn2min, params, fcn_args=(x, data))
minner.asteval.symtable['exprd'] = myexpr
params.add('decay', expr='myexpr')
result = minner.minimize()

# calculate final result
final = data + result.residual

# write error report
report_fit(result)

# try to plot results
try:
    import matplotlib.pyplot as plt
    plt.plot(x, data, 'k+')
    plt.plot(x, final, 'r')
    plt.show()
except ImportError:
    pass
# <end of examples/doc_parameters_basic.py>