from scipy import stats
import numpy as np
from scipy.optimize import minimize
import pylab as py

ydata = np.array([0.1,0.15,0.2,0.3,0.7,0.8,0.9, 0.9, 0.95])
xdata = np.array(range(0,len(ydata),1))

class mle:
    
    def sigmoid(self,params):
        k = params[0]
        x0 = params[1]   
        sd = params[2]
    
        yPred = 1 / (1+ np.exp(-k*(xdata-x0)))
    
        # Calculate negative log likelihood
        LL = -np.sum( stats.norm.logpdf(ydata, loc=yPred, scale=sd ) )
    
        return(LL)
    def calcMle(self):    
    
        initParams = [1, 1, 1]

        results = minimize(self.sigmoid, initParams, method='Nelder-Mead')
        print (results)
        return(results.x)
    
mm=mle()
estParms = mm.calcMle()
yOut = yPred = 1 / (1+ np.exp(-estParms[0]*(xdata-estParms[1])))

py.clf()
py.plot(xdata,ydata, 'go')
py.plot(xdata, yOut)
py.show()

