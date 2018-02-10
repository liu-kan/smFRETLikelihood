# -*- coding: utf-8 -*-
import binRawData
import BGrate
import numpy as np
from array import array
import irf_decov
import matplotlib.pyplot as plt

def prepareData(binData,sampleNum=1000,startt=0,fft=True):

    h,b=binRawData.statsDelayTime(binData,sampleNum,"D",byBin=True)#,bin0=100,binLen=2)    
    histlen=int(sampleNum/2)    
    if fft:
        return h[0][startt:histlen],np.linspace(0,histlen-startt-1,histlen-startt)
    else:
        return h[0][startt:histlen],b[0][startt:histlen]

def prepareCalcData(binData,histIRF,sampleNum=1000,fft=True,startt=0):
    lenbin=len(binData['chs']["All"]['chl'])
    lf=array('d')
    for i in range(0,lenbin):
        if binData['chs']["All"]['ntag'][i]>20:
            l,chi=irf_decov.calcTauOf1Bin(histIRF,binData,i,20,6.8695,'leastsq')
            if chi>=1 and chi<5000:
                lf.append(l)
    
    h,b=np.histogram(lf, sampleNum)            
    histlen=int(sampleNum/2)    
    if fft:
        return h[startt:histlen],np.linspace(0,histlen-startt-1,histlen-startt)
    else:
        return h[startt:histlen],b[startt:histlen]


if __name__=='__main__':
    import pickle,sys,getopt
    irfdbname="/dataB/smfretData/irf/alexa488_IRF_32MHz_PIE_3KCPS.sqlite"
    dbname="/dataB/smfretData/21c_224c.sqlite"
    savefn='/dataB/tmp/'+\
        dbname.split('/')[-1].split('.')[-2]+'_delayhist'+".pickle"
    fromdb=True
    binMs=1
    decay1=None;irf1=None;x=None;
    sampleNum=1000
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "fs:b:", ["fromfile", "samplenum=", "binms="])  
        for o, v in opts: 
            if o in ("-f", "--fromfile"):
                fromdb=False
            if o in ("-s", "--samplenum"):
                sampleNum = int(v)
            if o in ("-b", "--binms"):
                binMs = v                
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")        
        sys.exit(1); 
    if fromdb:
        br=BGrate.calcBGrate(dbname,0,400,T0=6,Tlen=167.97)
        binData=binRawData.binRawData(br,dbname,binMs)
        irfbr=BGrate.calcBGrate(irfdbname,20,400)#,T0=0.0,Tlen=600)
        irfbinData=binRawData.binRawData(irfbr,irfdbname,binMs)                
        with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([binData,irfbinData], f,protocol=-1)
    else:
        binData,irfbinData=pickle.load(open(savefn,'rb'))    
        
    # x,decay1,irf1=np.lo oadtxt(r"data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')
    hi,bi,_dt=binRawData.statsDelayTime(irfbinData,sampleNum,"D")#,bin0=100,binLen=2)    
  
    plt.figure()
    # decay1,x=prepareData(binData,sampleNum,fft=False)
    decay1,x=prepareCalcData(binData,hi,sampleNum)
    x=x*64/sampleNum
    plt.plot(x,decay1)
    plt.show()    




