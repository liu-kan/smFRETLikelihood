#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 10:31:42 2016

@author: liuk
"""

import numpy as np
import sys

import BurstSearch
import BGrate
import binRawData
sys.path.append('./ui')
from qtPlot import ControlMainWindow 

from array import array
from PyQt5 import QtWidgets



def FretAndS(burst,bins=(25,25),bgrate=None,bgfilter=True,ESm='k',chl='All'):
    #conn = sqlite3.connect(dbname)
    #c = conn.cursor()
    rESm=0
    if ESm=='K' or ESm=='k':
        rESm=0
    else:
        rESm=1
    lenburst=len(burst['chs'][chl]['chl'])
    print("lenburst:",lenburst)
    burstFRET = array("d")#np.zeros(lenburst)
    burstSeff = array("d")#np.zeros(lenburst)
    wei = array("d")#np.zeros(lenburst)
    #fw = np.zeros(lenburst)
    markDel=False
    if 'markDel' in burst['chs']["All"]:
            markDel=True    
    isBurst=False
    if 'burstW' in burst['chs'][chl]:
        isBurst=True
    for i in range(lenburst):
#        c.execute("select Dtime,ch from fretData_All where TimeTag>=? and TimeTag<= ?",
#                  (burst['chs']['All'].stag[i],burst['chs']['All'].etag[i]))
#        data=c.fetchall()
        if markDel:
            if burst['chs']['All']['markDel'][i]:
                continue
        data=burst['chs'][chl]['chl'][i]
        w=len(data)
        if type(w)!=type(1):
            continue
        if len(data)<1:
            continue      
        bgAA=0
        bgDD=0
        bgDA=0  
        if bgfilter:
            if bgrate!=None:
                tt=burst['chs'][chl]['timetag'][i]
                backgT=0
                if isBurst:
                    backgT=burst['chs'][chl]['burstW'][i]/2+tt[0]*\
                        bgrate["SyncResolution"] #中点时刻
                else:
                    backgT=burst['chs'][chl]['binMs']*0.5e-3+tt[0]*\
                        bgrate["SyncResolution"]                        
                bgAA=BurstSearch.getBGrateAtT(bgrate,"AexAem",backgT)
                bgDD=BurstSearch.getBGrateAtT(bgrate,"DexDem",backgT)
                bgDA=BurstSearch.getBGrateAtT(bgrate,"DexAem",backgT)            
            elif not isBurst:
                bgAA= burst['chs'][chl]['AAmean'] + burst['chs'][chl]['AAstd']#每个bin中的光子数
                bgDD=burst['chs'][chl]['DDmean']
                bgDA=burst['chs'][chl]['DAmean']

        nda=0#ch1
        ndd=0;#ch2
        naa=0;#ch3
        nad=0#ch4
        for d in data:
            if d==1:
                nda+=1
            elif d==2:
                ndd+=1
            elif d==3:
                naa+=1
            elif d==4:
                nad+=1
        if bgfilter:
            if bgrate!=None:
                if isBurst:
                    naa=naa-bgAA*burst['chs'][chl]['burstW'][i]
                    ndd=ndd-bgDD*burst['chs'][chl]['burstW'][i]
                    nda=nda-bgDA*burst['chs'][chl]['burstW'][i]
                    if naa< bgAA*burst['chs'][chl]['burstW'][i] or ndd<0 or nda<0:
                        continue
                else:
                    naa=naa-bgAA*burst['chs'][chl]['binMs']*1e-3
                    ndd=ndd-bgDD*burst['chs'][chl]['binMs']*1e-3
                    nda=nda-bgDA*burst['chs'][chl]['binMs']*1e-3
                    if naa< bgAA*burst['chs'][chl]['binMs']*1e-3 or ndd<0 or nda<0:
                        continue 
            elif not isBurst:
                naa=naa-bgAA
                ndd=ndd-bgDD
                nda=nda-bgDA
                if naa< bgAA or ndd<0 or nda<0:
                    continue                                                           
        wei.append(w)
        gamma=0.34
        beta=1.42
        DexDirAem=0.08
        Dch2Ach=0.07
        e=0;s=0
        if (nda+ndd)==0:
            burstFRET.append(1)
            burst['chs'][chl]['e'][i]=1
        else:
            if rESm==0:
                e=(nda)/(nda+gamma*ndd)
            else:                
                e=(nda*(1-DexDirAem)-Dch2Ach*ndd)/((1-DexDirAem)*nda+(gamma-Dch2Ach)*ndd)
            burstFRET.append(e)
            burst['chs'][chl]['e'][i]=e
        if (nda+ndd+naa)==0:
            burstSeff.append(1)
            burst['chs'][chl]['s'][i]=1
        else:
            if rESm==0:
                s=(nda+gamma*ndd)/(nda+gamma*ndd+naa/beta)
            else:
                s=((1-DexDirAem)*nda+(gamma-Dch2Ach)*ndd)/ \
                    ((1-DexDirAem)*nda+(gamma-Dch2Ach)*ndd+naa/beta)
            burstSeff.append(s)
            burst['chs'][chl]['s'][i]=s
    H, xedges, yedges = np.histogram2d(burstFRET,burstSeff, bins=bins)#, weights=wei)
    #conn.close()
    #fig, ax = plt.subplots()
    #plt.subplots_adjust(bottom=0.15)

    #im=plt.imshow(H.transpose()[::-1], interpolation='bessel',
    #              cmap=cm.jet,extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    #plt.colorbar(im)



    #rs=RectSect(ax,xedges,yedges)
    #toggle_selector(rs.toggle_selectorRS,plt)
    #RectBuilder(ax,xedges,yedges,rs.toggle_selectorRS)

    #callback = MaxLikehood(ax,burstSeff, burstFRET,burst)
    #axMLH = plt.axes([0.1, 0.02, 0.1, 0.05])

    #bnMLH = Button(axMLH, 'calc ML')
    #bnMLH.on_clicked(callback.calc)
    #import seaborn as sns
    #g = sns.JointGrid(x=burstFRET, y=burstSeff)

    #g.plot_marginals(sns.distplot)
    #g.plot_joint(plt.hist2d)
    return burstSeff, burstFRET,wei,H,xedges, yedges

if __name__ == '__main__':
    import pickle
    dbname='/home/liuk/sf/oc/data/38.sqlite'

    dbname='E:/liuk/proj/ptu/data/55.sqlite'
    #dbname='E:/sf/oc/data/38.sqlite'
    dbname="/Users/lp1/liuk/proj/data/LS9_150pM_poslineardiUb25c101c_alex488cy5_32MHz.sqlite"
    br=BGrate.calcBGrate(dbname,20,400,30,60)
    # burst=BurstSearch.findBurst(br,dbname,["All"],30,6)
    binTime=1
    dddaaaT=[7.1,4.1,3.1,9.1]
    if len(sys.argv)>1:
        binTime=float(sys.argv[1])      
    burst=binRawData.binRawData(br,dbname,binTime)
    binRawData.statsBins(burst,['All'])
    savefn='/dataZ1/tmp/lineardiub/'+\
        dbname.split('/')[-1].split('.')[-2]+'_'+str(binTime)+'_'+\
        str(dddaaaT)+".pickle"
    # with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([burst], f,protocol=-1)
    
    binRawData.burstFilterByBin(burst,dddaaaT)
    # binRawData.statsBins(burst,['AllBurst'])
    burstSeff, burstFRET,wei,H,xedges, yedges=FretAndS(burst,(27,27),None,True,'z'\
                ,"All")

    # app = QtWidgets.QApplication(sys.argv)
    # mySW = ControlMainWindow(H,xedges, yedges)
    # mySW.show()
    # sys.exit(app.exec_())
    title= "bin:"+str(binTime)+"ms,E-S"
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(2,1)
    im=ax[1].imshow(H.transpose()[::-1], interpolation='sinc', \
                       cmap=cm.jet,extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    # ax[1].set_title(title)
    fig.colorbar(im)                       
    
    # import matplotlib.pyplot as plt
    ax[0].hist(burstFRET, bins=40)#,weights=wei) 
    ax[0].set_title(title)
    plt.show()

    # import seaborn as sns
    # g = sns.JointGrid(x=burstFRET, y=burstSeff)

    # g.plot_marginals(sns.distplot)
    # g.plot_joint(plt.hist2d)    