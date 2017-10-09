import pickle
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import numpy as np
from array import array
from sklearn.metrics import r2_score

def fun(p, x, y):
    return p[0] * np.exp(-p[1] * x) - y

#detaT=0.68

def fitExpTau(burstTauD,thebin=200):
    detaT=max(burstTauD)/thebin*4.35
    
    [hist,bin]=np.histogram(burstTauD, bins=thebin)
    x=array("d")
    y=array("d")
    sidx=np.argmax(hist) #只拟合最大值之后的数据    
    lenhist=len(hist)
    startT=0
    lastH=hist[0]
    for vh in range(0,sidx):
        if hist[vh]>lastH*2.5:
            startT=bin[vh]
            break
        lastH=hist[vh]
        
    for vh in range(sidx,lenhist):
        if hist[vh]>0 and bin[vh]>bin[sidx]+detaT:        
            y.append(hist[vh])
            x.append(bin[vh]-bin[sidx]-detaT)     
            #print(bin[vh])           
    xx=np.array(x)
    yy=np.array(y)
    p0 = [xx[0],(2*np.log(2)/xx[0])]

    plsq = least_squares(fun, p0, args=(xx,yy))#bounds=([0, p0[0]], [9999,p0[1]]),
    
    coefficient_of_dermination = r2_score(yy, fun(plsq.x,xx,yy)+yy)
    print(coefficient_of_dermination)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # n, bins, patches = plt.hist(burstTauD, thebin, facecolor='g', alpha=0.75)
    # l = plt.plot(xx+bin[sidx]+detaT, fun(plsq.x,xx,yy)+yy, 'b--', linewidth=2)
    # #plt.legend(loc='best')
    # plt.show()    
    return [bin[sidx]-startT+1.0/plsq.x[1]+detaT,startT]

if __name__=='__main__':
    TauD = pickle.load( open( "../data/TauDraw.pickle", "rb" ) )
    [Tau,T0]=fitExpTau(TauD[0])
    print([Tau,T0])