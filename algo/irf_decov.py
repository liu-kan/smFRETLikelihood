# -*- coding: utf-8 -*-
import binRawData
import BGrate
import numpy as np
from lmfit import Model,Parameters
import matplotlib.pyplot as plt

def prepareData(irfdbname,binData,sampleNum=1000):
    '''
    TOFIX: sa.searchsorted(bins[-1], 'right')]
    两组数据的 sampleNum,b 需要兼容
    '''
    h,b=binRawData.statsDelayTime(binData,sampleNum,"D")#,bin0=100,binLen=2)
    irfbr=BGrate.calcBGrate(irfdbname,20,400)#,T0=0.0,Tlen=600)
    irfbinData=binRawData.binRawData(br,irfdbname,binData['binMs'])
    hi,bi=binRawData.statsDelayTime(irfbinData,sampleNum,"D")#,bin0=100,binLen=2)
    histlen=int(sampleNum/2)
    startt=0
    binRawData.statsBins(irfbinData)
    print("IRF photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    print(bi[0][0],b[0][0],bi[0][10],b[0][10])        
    return h[0][startt:histlen],hi[0][startt:histlen],np.linspace(startt,histlen-1,histlen-startt)


# define the single exponential model
def jumpexpmodel(t,irf,y0,x0,tau1,ampl1,tau2,ampl2,tau3,ampl3,tau4,ampl4,tau5,ampl5):# Lifetime decay fit Author: Antonino Ingargiola - Date: 07/20/2014
    ymodel=np.zeros(t.size)
    c=x0
    n=len(irf)
    irf_s1=np.remainder(np.remainder(t-np.floor(c)-1, n)+n,n)
    irf_s11=(1-c+np.floor(c))*irf[irf_s1.astype(int)]
    irf_s2=np.remainder(np.remainder(t-np.ceil(c)-1,n)+n,n)
    irf_s22=(c-np.floor(c))*irf[irf_s2.astype(int)]
    irf_shift=irf_s11+irf_s22
    irf_reshaped_norm=irf_shift/sum(irf_shift)
    ymodel = ampl1*np.exp(-(t)/tau1)
    ymodel+= ampl2*np.exp(-(t)/tau2)
    ymodel+= ampl3*np.exp(-(t)/tau3)
    ymodel+= ampl4*np.exp(-(t)/tau4)
    ymodel+= ampl5*np.exp(-(t)/tau5)
    z=Convol(ymodel,irf_reshaped_norm)
    z+=y0
    return z

def plot_fit(time,irf, ydata, params, yscale='log', zoom_origin=False):
    """
    Function to plot data and model function.
    """
    plt.semilogy(time, ydata, marker='.')
    plt.semilogy(time, jumpexpmodel(time,irf, **{k: v.value for k, v in params.items()}),
         color='r', lw=2.5, alpha=0.7)
    if zoom_origin:
        dt = time[1] - time[0]
        t0 = params['offset'] - 50*dt
        plt.xlim(t0 - 50*dt, t0 + 50*dt)

def Convol(x,h): # change in convolution calcualataion
    #need same length of x and h
    X=np.fft.fft(x)
    H=np.fft.fft(h)
    xch=np.real(np.fft.ifft(X*H))
    return xch

def fitIRF(decayData,irf,delayT,pars):
    wWeights=1/np.sqrt(decayData+1)# check for divide by zero,  have used +1 to avoid dived by zero
    mod=Model(jumpexpmodel,independent_vars=['t','irf'])    
    # fit this model with weighting , parameters as initilized and print the result
    result = mod.fit(decayData,params=pars,weights=wWeights,method='leastsq',t=delayT,irf=irf)
    print(result.fit_report())
    return result

if __name__=='__main__':
    import pickle
    irfdbname="/dataB/smfretData/irf/alexa488_IRF_32MHz_PIE_3KCPS.sqlite"
    dbname="/dataB/smfretData/rsv86c224c.sqlite"
    savefn='/dataB/tmp/'+\
        dbname.split('/')[-1].split('.')[-2]+'_gd'+".pickle"

    # br=BGrate.calcBGrate(dbname,20,400,T0=6,Tlen=615.56797)
    # binData=binRawData.binRawData(br,dbname,1)
    # binRawData.statsBins(binData)
    # print("Data photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    # decay1,irf1,x=prepareData(irfdbname,binData,sampleNum=1000)    
    # with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([decay1,irf1,x], f,protocol=-1)

    decay1,irf1,x=pickle.load(open(savefn,'rb'))    
    # x,decay1,irf1=np.loadtxt(r"data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')
    params = Parameters()
    params.add('x0', value=0)#, vary=False)
    params.add('y0', value=120)#, vary=False)
    params.add('tau1', value=20)
    params.add('tau2', value=50)
    params.add('ampl1', value=4500)#, vary=False)
    params.add('ampl2', value=4500)
    params.add('tau3', value=50)
    params.add('ampl3', value=4500)#
    params.add('tau4', value=50)
    params.add('ampl4', value=4500)#
    params.add('tau5', value=50)
    params.add('ampl5', value=4500)#
    # plt.close('all')
    # plt.figure(1)
    # plot_fit(x,irf1, decay1, params)
    # # plt.figure(2)
    # # plt.semilogy(x,irf1)
    # plt.show()

    result=fitIRF(decay1,irf1,x,params)    
    plt.figure(5)
    plt.subplot(2,1,1)
    plt.semilogy(x,decay1,'r-',x,result.best_fit,'b')
    plt.subplot(2,1,2)
    plt.plot(x,result.residual)
    plt.show()