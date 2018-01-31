# -*- coding: utf-8 -*-
import BurstSearch
import BGrate
import numpy as np
from array import array

import matplotlib.pyplot as plt


if __name__=='__main__':
    import pickle,sys,getopt

    dbname="/dataB/smfretData/21c_224c.sqlite"
    savefn='/dataB/tmp/'+\
        dbname.split('/')[-1].split('.')[-2]+'_burstwhist'+".pickle"
    fromdb=True
    sampleNum=1000
    burst=None
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
        br=BGrate.calcBGrate(dbname,20,400,chs=["All"])#,T0=10.0,Tlen=160.0)    
        burst=BurstSearch.findBurst(br,dbname,["All"],30,6)
        with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump(burst, f,protocol=-1)
    else:
        burst=pickle.load(open(savefn,'rb'))    
        
    lenbin=len(burst['chs']["All"]['burstW'])
    lf=array('d')
    for i in range(lenbin):
        lf.append(burst['chs']["All"]['burstW'][i])
    h,b=np.histogram(lf, sampleNum)                

    plt.figure()
    size=sampleNum/10
    plt.plot(b[0:size],h[0:size])
    plt.show()    




