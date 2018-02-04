# -*- coding: utf-8 -*-
import BurstSearch
import BGrate
import numpy as np
from array import array

import matplotlib.pyplot as plt
from array import array
from scipy.interpolate import spline 

if __name__=='__main__':
    import pickle,sys,getopt

    dbname="/dataB/smfretData/21c_224c.sqlite"
    savefn='data/'+\
        "21c_224c_0.5_[2.3509349643187933, 1.1838989238159361, 3.3364004127243341, 3.5348338881347292].pickle"
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
        burstTau, burstFRET,burst=pickle.load(open(savefn,'rb'))   
    binw=27
    bins=(binw,binw)    
    Hp, xedgesp, yedgesp = np.histogram2d(burstFRET,burstTau, bins=(27,27))              
    H, xedges, yedges = np.histogram2d(burstFRET,burstTau, bins=bins)              
    H=H.transpose()[::-1]
    #H[1]=np.zeros(binw27)
    gzw=1#/binw
    lfy=array('d')
    lfx=array('d')
    for idx in range(binw):
        y=H[:,idx]
        toty=0
        ynow=0
        yv=0
        for idy in range(binw):
            toty+=H[idy,idx]*gzw
        if toty>4000:
            for idy in range(binw):
                ynow+=H[idy,idx]*gzw
                if ynow>=toty/2:
                    lfy.append(yedges[idy])
                    lfx.append(xedges[idx])
                    break
        lfy[idx]=yedges[4]
    xnew = np.linspace(min(lfx),max(lfx),300) 
    _smooth = spline(lfx,lfy,xnew)
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,1)
    im=ax.imshow(Hp.transpose()[::-1], interpolation='sinc', \
                       cmap=cm.jet,extent=[xedgesp[0],xedgesp[-1],yedgesp[0],yedgesp[-1]])
    # ax[1].set_title(title)
    # ax.plot(xnew,_smooth)  
    ax.plot(lfx,lfy)  
    fig.colorbar(im)                       
    

    plt.show()  




