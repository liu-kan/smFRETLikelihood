# my first python program
#iterative convolution of IRF with flourescenece decay, measured by TCSPC
# thank you Dr. Jon M. Pikal ( my adviser, University of Wyoming)
# thank you Dr. Matt Newville for pointing me in right direction

import numpy as np
from lmfit import Model
import matplotlib.pyplot as plt
plt.close('all')
# read data from file
x,decay1,irf=np.loadtxt(r"data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')
# plot the raw data file ( irf and decay1)
plt.figure(1)
plt.semilogy(x,decay1,x,irf)
plt.show()
# Calculate the weighting factor for tcspc
wWeights=1/np.sqrt(decay1+1)# check for divide by zero,  have used +1 to avoid dived by zero

# define the single exponential model
def jumpexpmodel(x,tau1,ampl1,tau2,ampl2,y0,x0,args=(irf)):# Lifetime decay fit Author: Antonino Ingargiola - Date: 07/20/2014
    ymodel=np.zeros(x.size)
    t=x
    c=x0
    n=len(irf)
    irf_s1=np.remainder(np.remainder(t-np.floor(c)-1, n)+n,n)
    irf_s11=(1-c+np.floor(c))*irf[irf_s1.astype(int)]
    irf_s2=np.remainder(np.remainder(t-np.ceil(c)-1,n)+n,n)
    irf_s22=(c-np.floor(c))*irf[irf_s2.astype(int)]
    irf_shift=irf_s11+irf_s22
    irf_reshaped_norm=irf_shift/sum(irf_shift)
    ymodel = ampl1*np.exp(-(x)/tau1)
    ymodel+= ampl2*np.exp(-(x)/tau2)
    z=Convol(ymodel,irf_reshaped_norm)
    z+=y0
    return z

def Convol(x,h): # change in convolution calcualataion
    #need same length of x and h
    X=np.fft.fft(x)
    H=np.fft.fft(h)
    xch=np.real(np.fft.ifft(X*H))
    return xch

# this is just for testing propse, test the exponential decay function
#plt.figure(2)
#def jumpexpmodel(x,tau1,ampl1,tau2,ampl2,y0,x0,args=(irf)):
#y=jumpexpmodel(x,50,2632.85,220.36,543.9,21.17,8.56280)
#plt.semilogy(x,y,'ro',x,decay1,'bo',x,irf)
#plt.title("test the model")

# assign the model for fitting and initialize the parameters
mod=Model(jumpexpmodel)
pars = mod.make_params(tau1=10,ampl1=1000,tau2=10,ampl2=1000,y0=0,x0=10,args=irf)
pars['x0'].vary =True
pars['y0'].vary =True

# fit this model with weighting , parameters as initilized and print the result
result = mod.fit(decay1,params=pars,weights=wWeights,method='leastsq',x=x)
print(result.fit_report())
plt.figure(5)
plt.subplot(2,1,1)
plt.semilogy(x,decay1,'r-',x,result.best_fit,'b')
plt.subplot(2,1,2)
plt.plot(x,result.residual)
plt.show()

''' references and helpful notes
1.http://www.igorexchange.com/node/4201
  TCSPC deconvolution Igor
2. http://nbviewer.ipython.org/github/tritemio/notebooks/blob/master/Lifetime_decay_fit.ipynb
   lifetime decay fit
3. Tcspcfit - A Matlab package for fitting multiexponential fluorescence decay curves
4. lmfit nonlinear curve fitting manual
5. I have used ideas from these and other reference
'''
