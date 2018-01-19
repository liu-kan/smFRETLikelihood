# -*- coding: utf-8 -*-
import binRawData
import BGrate
import numpy as np
from lmfit import Model

def prepareData(irfdbname,binData,sampleNum=1000):
    '''
    TOFIX: sa.searchsorted(bins[-1], 'right')]
    两组数据的 sampleNum,b 需要兼容
    '''
    h,b=binRawData.statsDelayTime(binData,sampleNum,"D")#,bin0=100,binLen=2)
    irfbr=BGrate.calcBGrate(irfdbname,20,400)#,T0=0.0,Tlen=600)
    irfbinData=binRawData.binRawData(br,irfdbname,binData['binMs'])
    hi,bi=binRawData.statsDelayTime(binData,b,"D")#,bin0=100,binLen=2)
    binRawData.statsBins(irfbinData)
    print("IRF photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    print(bi[0][0],b[0][0],bi[0][10],b[0][10])
    return h[0],hi[0],b[0]



#x,decay1,irf=np.loadtxt(r"../data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')

irf=()
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

def fitIRF(decayData,irf,delayT):
    # this is just for testing propse, test the exponential decay function
    #plt.figure(2)
    #def jumpexpmodel(x,tau1,ampl1,tau2,ampl2,y0,x0,args=(irf)):
    #y=jumpexpmodel(x,50,2632.85,220.36,543.9,21.17,8.56280)
    #plt.semilogy(x,y,'ro',x,decay1,'bo',x,irf)
    #plt.title("test the model")
    # assign the model for fitting and initialize the parameters
    wWeights=1/np.sqrt(decayData+1)# check for divide by zero,  have used +1 to avoid dived by zero
    mod=Model(jumpexpmodel)
    pars = mod.make_params(tau1=10,ampl1=1000,tau2=10,ampl2=1000,y0=0,x0=10,args=irf)
    pars['x0'].vary =True
    pars['y0'].vary =True

    # fit this model with weighting , parameters as initilized and print the result
    result = mod.fit(decayData,params=pars,weights=wWeights,method='leastsq',x=delayT)
    print(result.fit_report())

if __name__=='__main__':
    irfdbname="/dataB/smfretData/irf/alexa488_IRF_32MHz_PIE_3KCPS.sqlite"
    dbname="/dataB/smfretData/rsv86c224c.sqlite"
    br=BGrate.calcBGrate(dbname,20,400,T0=6,Tlen=615.56797)
    binData=binRawData.binRawData(br,dbname,1)
    binRawData.statsBins(binData)
    print("Data photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    prepareData(irfdbname,binData,sampleNum=1000)
