import pickle
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from array import array


def func(x, a, b):
    return a * np.exp(-1/b * x) 

TauD = pickle.load( open( "../data/TauDraw.pickle", "rb" ) )
burstTauD = TauD[0]
ad=np.dot(burstTauD,1e-9)
ad=burstTauD
thebin=200
[hist,bin]=np.histogram(ad, bins=thebin)
x=array("d")
y=array("d")
sidx=np.argmax(hist) #只拟合最大值之后的数据    
lenhist=len(hist)
for vh in range(sidx,lenhist):
    if hist[vh]>0 and bin[vh]>bin[sidx]+0.2:        
        y.append(hist[vh])
        x.append(bin[vh]-bin[sidx]-0.2)     
        #print(bin[vh])           
xx=np.array(x)
yy=np.array(y)
p0 = [np.min(ad),np.max(ad)]
#plsq = least_squares(fun, p0, jac=jac, gtol=0,bounds=([0, p0[0]], [9999,p0[1]]), args=(xx,yy), verbose=1)
popt, pcov = curve_fit(func, xx, yy, \
                verbose=1,method='trf',loss='cauchy')
#,bounds=([-np.inf, p0[0]/20,-np.inf], [np.inf,p0[1]*20,np.inf]))#,\
print(popt)
        
fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = plt.hist(ad, thebin, facecolor='g', alpha=0.75)


l = plt.plot(xx+bin[sidx]+0.2, func(xx,*popt), 'b--', linewidth=2)

plt.legend(loc='best')
plt.show()    
                    