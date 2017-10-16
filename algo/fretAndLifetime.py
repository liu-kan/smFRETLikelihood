#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 10:31:42 2016

@author: liuk
"""

import numpy as np
try:
    import algo.BurstSearch as BurstSearch
    import algo.BGrate as BGrate
except ImportError:
    import BurstSearch
    import BGrate

import matplotlib
from array import array
from scipy.optimize import leastsq

import matplotlib.cm as cm
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QVBoxLayout
#from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import sys
import matplotlib as mpl
#mpl.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib.pyplot as plt
import pickle

class ControlMainWindow(QtWidgets.QMainWindow):
    def __init__(self,H,xedges,yedges, parent=None):
        super(ControlMainWindow, self).__init__(parent)
        self.H=H
        self.xedges=xedges
        self.yedges=yedges
        self.main_frame = QtWidgets.QWidget()
        self.setupUi()

    def setupUi(self):
        vbox = QVBoxLayout()
        figure = mpl.figure.Figure(figsize=(10, 10))
        leftcanvas = MatplotlibWidget(figure,self.H,self.xedges,self.yedges)
        #vbox.addWidget(self.mpl_toolbar)
        leftcanvas.setParent(self.main_frame)
        vbox.addWidget(leftcanvas)
        self.mpl_toolbar = NavigationToolbar(leftcanvas, self.main_frame)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

class MatplotlibWidget(FigureCanvasQTAgg):
    def __init__(self, fig,H,xedges,yedges):
        super(MatplotlibWidget, self).__init__(fig)

        #-- set up an axe artist --

        self.dax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        im=self.dax.imshow(H.transpose()[::-1], interpolation='bessel',
                      cmap=cm.jet,extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        fig.colorbar(im)
        self.rs=RectSect(self.dax,xedges,yedges)
        self.rb=RectBuilder(self.dax,xedges,yedges,self.rs.toggle_selectorRS)
        self.draw()
        self.updatingLine=True

        #---- setup event ----
        #self.mpl_connect('key_press_event', self.onkeypress)
        #self.mpl_connect('motion_notify_event', self.onMouseMove)
        #self.mpl_connect('button_press_event', self.onclick)
    def onkeypress(self, event):
        print(' Key pressed.')
        if event.key in ['D', 'd'] and self.rs.active:
            print(' RectangleSelector deactivated.')
            self.rs.set_active(False)
        if event.key in ['A', 'a'] and not self.rs.active:
            print(' RectangleSelector activated.')
            self.rs.set_active(True)
    def onMouseMove(self,event):
        if not event.inaxes:
            return
        if self.updatingLine:
            self.dax.lines = [self.dax.lines[0]] # keep the first two lines
            self.dax.axhline(y=event.ydata, color="k")
            self.draw()

    def onclick(self, event):
        if not event.inaxes:
            return
        self.updatingLine=not self.updatingLine



#burstD如果是浮点数，则为Donor only的TauD，否则为Donor only的提取的burst
def FretAndLifetime(burst,bins=(25,25),bgrate=None,burstD=4.1,bgrateD=None,T0=6.8695):
    #conn = sqlite3.connect(dbname)
    #c = conn.cursor()
    Tau_D=1e-9
    if type(burstD) is float:
        Tau_D=burstD*1e-9
        print('已经制定TauD')
    else:
        lenburstD=len(burstD["All"]['chl'])
        burstTauD = array("d")#np.zeros(lenburst)
        weiD = array("d")#np.zeros(lenburst)
            
        for i in range(lenburstD):
            data=burstD["All"]['chl'][i]
            if bgrateD==None:
                tt=burstD["All"]['timetag'][i]
                #print(tt)
                backgT=burstD["All"]['burstW'][i]/2+tt[0]*bgrate["SyncResolution"] #中点时刻
                bgAA=BurstSearch.getBGrateAtT(bgrateD,"AexAem",backgT)
                bgDD=BurstSearch.getBGrateAtT(bgrateD,"DexDem",backgT)
                bgDA=BurstSearch.getBGrateAtT(bgrateD,"DexAem",backgT)            
            w=len(data)       
            nda=0#ch1
            ndd=0;#ch2
            naa=0;#ch3
            nad=0#ch4        
            sumdtimed=array("d")
            for idxd in range(w):
                d=data[idxd]
                if d==1:
                    nda+=1
                elif d==2:
                    ndd+=1
                    sumdtimed.append(burstD["All"]['dtime'][i][idxd]*burstD["DelayResolution"]*1e9)
                elif d==3:
                    naa+=1
                elif d==4:
                    nad+=1                
            if bgrateD==None:
                naa=naa-bgAA*burstD["All"]['burstW'][i]
                ndd=ndd-bgDD*burstD["All"]['burstW'][i]
                nda=nda-bgDA*burstD["All"]['burstW'][i]
                if naa< bgAA*burstD["All"]['burstW'][i] or ndd<0 or nda<0:
                    continue
            weiD.append(w)
            Tau=np.mean(sumdtimed)
            burstTauD.append(Tau)
        [hist,bin]=np.histogram(burstTauD, bins=200)
        x=array("d")
        y=array("d")
        sidx=np.argmax(hist) #只拟合最大值之后的数据
        lenhist=len(hist)
        for vh in range(sidx,lenhist):
            if hist[vh]>0:        
                y.append(np.log(hist[vh]))
                x.append(bin[vh])            
        m = 2  #多项式的次数
        #先随机产生一组多项式分布的参数
        p0 = np.random.randn(m)
        plsq = leastsq(BGrate.residuals, p0, args=(y, x))
        print(plsq)
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # n, bins, patches = plt.hist(burstTauD, 200, facecolor='g', alpha=0.75)
        # plt.legend(loc='best')
        # plt.show()
        Tau_D=-1*plsq[0][0]*1e-9#        
        with open('../data/TauD.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([burstTauD], f)
    print('Tau_D:',Tau_D)
    lenburst=len(burst["All"]['chl'])
    print("lenburst:",lenburst)
    burstFRET = array("d")#np.zeros(lenburst)
    burstSeff = array("d")#np.zeros(lenburst)
    burstTau = array("d")#np.zeros(lenburst)
    wei = array("d")#np.zeros(lenburst)
    #fw = np.zeros(lenburst)

    '''     
    i=1000
    data=burst["All"]['chl'][i]
    w=len(data)
    print(w)
    nda=0#ch1
    ndd=0;#ch2
    naa=0;#ch3
    nad=0#ch4
    histdtime=array("d")
    for idxd in range(w):
        d=data[idxd]
        dtime=burst["All"]['dtime'][i][idxd]*burst["DelayResolution"]
        if d==1:
            nda+=1
        elif d==2:
            ndd+=1
            histdtime.append(dtime)
        elif d==3:
            naa+=1
        elif d==4:
            nad+=1
            #sumdtimed+=dtime
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = plt.hist(histdtime, 100,  facecolor='g', alpha=0.75)
    print(len(histdtime))
    plt.legend(loc='best')
    plt.show() 
    '''
    

    for i in range(lenburst):
#        c.execute("select Dtime,ch from fretData_All where TimeTag>=? and TimeTag<= ?",
#                  (burst["All"].stag[i],burst["All"].etag[i]))
#        data=c.fetchall()
        data=burst["All"]['chl'][i]
        if bgrate!=None:
            tt=burst["All"]['timetag'][i]
            #print(tt)
            backgT=burst["All"]['burstW'][i]/2+tt[0]*bgrate["SyncResolution"] #中点时刻
            bgAA=BurstSearch.getBGrateAtT(bgrate,"AexAem",backgT)
            bgDD=BurstSearch.getBGrateAtT(bgrate,"DexDem",backgT)
            bgDA=BurstSearch.getBGrateAtT(bgrate,"DexAem",backgT)            
        w=len(data)
        #print(w)
        nda=0#ch1
        ndd=0;#ch2
        naa=0;#ch3
        nad=0#ch4
        sumdtimed=array("d")
        for idxd in range(w):
            d=data[idxd]            
            if d==1:
                nda+=1
            elif d==2:
                ndd+=1
                detime=burst["All"]['dtime'][i][idxd]*burst["DelayResolution"]-T0*1e-9
                if detime>=0:
                    sumdtimed.append(detime)
            elif d==3:
                naa+=1
            elif d==4:
                nad+=1
                #sumdtimed+=dtime
        
        if bgrate!=None:
            naa=naa-bgAA*burst["All"]['burstW'][i]
            ndd=ndd-bgDD*burst["All"]['burstW'][i]
            nda=nda-bgDA*burst["All"]['burstW'][i]
            if naa< bgAA*burst["All"]['burstW'][i] or ndd<0 or nda<0:
                continue        
        Tau=np.mean(sumdtimed)/(Tau_D)
        if True:#Tau<=1:
            wei.append(w)
            burstTau.append(Tau)
            burst["All"]['lifetime'][i]=Tau        
            gamma=0.31        
            beta=1.42
            if (nda+ndd)==0:
                burstFRET.append(1)
                burst["All"]['e'][i]=1
            else:
                theFret=(nda)/(nda+gamma*ndd)
                burstFRET.append(theFret)
                burst["All"]['e'][i]=theFret
            if (nda+ndd+naa)==0:
                burstSeff.append(1)
                burst["All"]['s'][i]=1
            else:
                theFret=(nda+gamma*ndd)/(nda+gamma*ndd+naa/beta)
                burstSeff.append(theFret)
                burst["All"]['s'][i]=theFret

    H, xedges, yedges = np.histogram2d(burstFRET,burstTau, bins=bins, weights=wei)
    #print(burstTau[0:100])
    #conn.close()
    fig, ax = plt.subplots()
    #plt.subplots_adjust(bottom=0.15)

    im=plt.imshow(H.transpose()[::-1], interpolation='bessel',
                  cmap=cm.jet,#extent=[0,1,0,1])
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plt.colorbar(im)
    plt.show()


    #rs=RectSect(ax,xedges,yedges)
    #toggle_selector(rs.toggle_selectorRS,plt)
    #RectBuilder(ax,xedges,yedges,rs.toggle_selectorRS)

    #callback = MaxLikehood(ax,burstSeff, burstFRET,burst)
    #axMLH = plt.axes([0.1, 0.02, 0.1, 0.05])

    #bnMLH = Button(axMLH, 'calc ML')
    #bnMLH.on_clicked(callback.calc)
    #import seaborn as sns
    #g = sns.JointGrid(x=burstFRET, y=burstTau)

    #g.plot_marginals(sns.distplot)
    #g.plot_joint(plt.hist2d)
    return burstTau, burstFRET,wei,H,xedges, yedges

def betweenXY(v,arr):
    si=-1
    ei=-1
    if v==None:
        return si,ei
    for i, av in enumerate(arr):
        if ei==-1 and av<=v:
            si=i
        if si!=-1 and av>=v:
            ei=i
            break
    return si,ei

class toggle_selector:
    def __init__(self, rs,pltc):

        self.rs=rs
        rs.set_active(False)
        pltc.connect('key_press_event', self)
    def __call__(self, event):
        print(' Key pressed.')
        if event.key in ['D', 'd'] and self.rs.active:
            print(' RectangleSelector deactivated.')
            self.rs.set_active(False)
        if event.key in ['A', 'a'] and not self.rs.active:
            print(' RectangleSelector activated.')
            self.rs.set_active(True)

class RectSect:
    def __init__(self, ax,xedges,yedges):

        self.xedges=xedges
        self.yedges=yedges
        self.rects = ax
        self.toggle_selectorRS = matplotlib.widgets.RectangleSelector(ax, self,
                                       drawtype='box', useblit=False,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=False)
        self.toggle_selectorRS.set_active(False)
        self.n_m1=lambda x: True if x!=-1 else False

    def __call__(self, eclick,erelease):
        'eclick and erelease are the press and release events'
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        #print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        #print(" The button you used were: %s %s" % (eclick.button, erelease.button))

        ax=betweenXY(x1,self.xedges)
        ay=betweenXY(y1,self.yedges)

        ax2=betweenXY(x2,self.xedges)
        ay2=betweenXY(y2,self.yedges)

        if all(map(self.n_m1,ax)) and all(map(self.n_m1,ay)) and all(map(self.n_m1,ax2)) and all(map(self.n_m1,ay2)):
            rx=self.xedges[ax2[1]]-self.xedges[ax[0]]
            ry=self.yedges[ay2[1]]-self.yedges[ay[0]]
            print(rx,ry)
            b=matplotlib.patches.Rectangle((self.xedges[ax[0]],
                                        self.yedges[ay[0]]), rx,
                                       ry,color='black')
            self.rects.add_patch(b)
            self.rects.figure.canvas.draw()
        print(len(self.rects.patches))

class RectBuilder:
    def __init__(self, rect,xedges,yedges,rs):
        self.rects = rect
        self.rx=xedges[1]-xedges[0]
        self.ry=yedges[1]-yedges[0]
        self.rs=rs
        self.xedges=xedges
        self.yedges=yedges
        self.cid = rect.figure.canvas.mpl_connect('button_press_event', self)
        self.n_m1=lambda x: True if x!=-1 else False

    def __call__(self, event):
        #print('click', event)
        if self.rs.active: return
        if event.inaxes!=self.rects.axes: return

        ax=betweenXY(event.xdata,self.xedges)
        ay=betweenXY(event.ydata,self.yedges)

        if all(map(self.n_m1,ax)) and all(map(self.n_m1,ay)):
            b=matplotlib.patches.Rectangle((self.xedges[ax[0]],
                                        self.yedges[ay[0]]), self.rx,
                                       self.ry,color='black')
            #print((self.xedges[ax[0]],self.yedges[ay[0]]))
            self.rects.add_patch(b)
            self.rects.figure.canvas.draw()
        print(len(self.rects.patches))

if __name__ == '__main__':
    import pickle
    dbname="/home/liuk/proj/data/LS35_RSV86C224C.sqlite"
    dbname="/home/liuk/proj/data/LS1_150pM_48diUb_NC_488_cy5_32MHz_1.sqlite"
    dbTau_D="/home/liuk/proj/data/Tau_D.sqlite"
    br=BGrate.calcBGrate(dbname,20,400)
    if type(br)==type(1):
        exit(-1)

    burst=BurstSearch.findBurst(br,dbname,["All"])
    #brD=BGrate.calcBGrate(dbTau_D,20,400)
    #burstD=BurstSearch.findBurst(br,dbTau_D,["All"])

    burstSeff, burstFRET,wei,H,xedges, yedges=FretAndLifetime(burst,(25,25),br,4.1)

    # with open('E:/tmp/objs.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([burstSeff, burstFRET,wei,H,xedges], f)

    # Getting back the objects:
    #with open('objs.pickle') as f:  # Python 3: open(..., 'rb')
        #obj0, obj1, obj2 = pickle.load(f)

    #app = QtWidgets.QApplication(sys.argv)
    #mySW = ControlMainWindow(H,xedges, yedges)
    #mySW.show()
    #sys.exit(app.exec_())
