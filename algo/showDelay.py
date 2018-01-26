# -*- coding: utf-8 -*-
import binRawData
import BGrate
import numpy as np

import matplotlib.pyplot as plt

def prepareData(binData,sampleNum=1000,startt=0,fft=True):

    h,b=binRawData.statsDelayTime(binData,sampleNum,"D")#,bin0=100,binLen=2)    
    histlen=int(sampleNum/2)    
    if fft:
        return h[0][startt:histlen],np.linspace(0,histlen-startt-1,histlen-startt)
    else:
        return h[0][startt:histlen],b[0][startt:histlen]

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

        with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump(binData, f,protocol=-1)
    else:
        binData=pickle.load(open(savefn,'rb'))    
    # x,decay1,irf1=np.lo oadtxt(r"data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')
    decay1,x=prepareData(binData,sampleNum)    
  
    plt.figure()
    decay1,x=prepareData(binData,sampleNum,fft=False)
    plt.semilogy(x,decay1)
    plt.show()    




