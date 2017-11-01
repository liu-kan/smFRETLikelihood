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



def FretAndS(burst,bins=(25,25),bgrate=None,filter=True):
    #conn = sqlite3.connect(dbname)
    #c = conn.cursor()
    lenburst=len(burst['chs']['All']['chl'])
    print("lenburst:",lenburst)
    burstFRET = array("d")#np.zeros(lenburst)
    burstSeff = array("d")#np.zeros(lenburst)
    wei = array("d")#np.zeros(lenburst)
    #fw = np.zeros(lenburst)
    isBurst=False
    if 'burstW' in burst['chs']['All']:
        isBurst=True
    for i in range(lenburst):
#        c.execute("select Dtime,ch from fretData_All where TimeTag>=? and TimeTag<= ?",
#                  (burst['chs']['All'].stag[i],burst['chs']['All'].etag[i]))
#        data=c.fetchall()
        data=burst['chs']['All']['chl'][i]
        w=len(data)
        if type(w)!=type(1):
            continue
        if len(data)<1:
            continue      
        bgAA=0
        bgDD=0
        bgDA=0  
        if filter:
            if bgrate!=None:
                tt=burst['chs']['All']['timetag'][i]
                backgT=0
                if isBurst:
                    backgT=burst['chs']['All']['burstW'][i]/2+tt[0]*bgrate["SyncResolution"] #中点时刻
                else:
                    backgT=burst['chs']['All']['binMs']*0.5e-3+tt[0]*bgrate["SyncResolution"]                        
                bgAA=BurstSearch.getBGrateAtT(bgrate,"AexAem",backgT)
                bgDD=BurstSearch.getBGrateAtT(bgrate,"DexDem",backgT)
                bgDA=BurstSearch.getBGrateAtT(bgrate,"DexAem",backgT)            
            elif not isBurst:
                bgAA= burst['chs']['AexAem']['mean'] + burst['chs']['AexAem']['std']#每个bin中的光子数
                bgDD=burst['chs']['DexDem']['mean']
                bgDA=burst['chs']['DexAem']['mean']

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
        if filter:
            if bgrate!=None:
                if isBurst:
                    naa=naa-bgAA*burst['chs']['All']['burstW'][i]
                    ndd=ndd-bgDD*burst['chs']['All']['burstW'][i]
                    nda=nda-bgDA*burst['chs']['All']['burstW'][i]
                    if naa< bgAA*burst['chs']['All']['burstW'][i] or ndd<0 or nda<0:
                        continue
                else:
                    naa=naa-bgAA*burst['chs']['All']['binMs']*1e-3
                    ndd=ndd-bgDD*burst['chs']['All']['binMs']*1e-3
                    nda=nda-bgDA*burst['chs']['All']['binMs']*1e-3
                    if naa< bgAA*burst['chs']['All']['binMs']*1e-3 or ndd<0 or nda<0:
                        continue 
            elif not isBurst:
                naa=naa-bgAA
                ndd=ndd-bgDD
                nda=nda-bgDA
                if naa< bgAA or ndd<0 or nda<0:
                    continue                                                           
        wei.append(w)
        gamma=0.31        
        beta=1.42
        if (nda+ndd)==0:
            burstFRET.append(1)
            burst['chs']['All']['e'][i]=1
        else:
            burstFRET.append((nda)/(nda+gamma*ndd))
            burst['chs']['All']['e'][i]=(nda)/(nda+gamma*ndd)
        if (nda+ndd+naa)==0:
            burstSeff.append(1)
            burst['chs']['All']['s'][i]=1
        else:
            burstSeff.append((nda+gamma*ndd)/(nda+gamma*ndd+naa/beta))
            burst['chs']['All']['s'][i]=(nda+gamma*ndd)/(nda+gamma*ndd+naa/beta)

    H, xedges, yedges = np.histogram2d(burstFRET,burstSeff, bins=bins, weights=wei)
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
    dbname="../data/rsv86c224c.sqlite"
    br=BGrate.calcBGrate(dbname,20,400)
    burst=BurstSearch.findBurst(br,dbname,["All"],15,3.5)
    #burst=binRawData.binRawData(br,dbname,2)
    #binRawData.statsBins(burst)
    burstSeff, burstFRET,wei,H,xedges, yedges=FretAndS(burst,(27,27),None)

    # with open('E:/tmp/objs.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([burstSeff, burstFRET,wei,H,xedges], f)

    # Getting back the objects:
    #with open('objs.pickle') as f:  # Python 3: open(..., 'rb')
        #obj0, obj1, obj2 = pickle.load(f)

    app = QtWidgets.QApplication(sys.argv)
    mySW = ControlMainWindow(H,xedges, yedges)
    mySW.show()
    sys.exit(app.exec_())
