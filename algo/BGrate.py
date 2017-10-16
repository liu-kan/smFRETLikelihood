# -*- coding: utf-8 -*-
#作者:liuk
#博客:liuk.io
#邮箱:liukan@126.com

import numpy as np
#import scipy as sp
from scipy.optimize import leastsq #这里就是我们要使用的最小二乘的函数

import sqlite3
from array import array
from cycler import cycler
#import collections
from sklearn.metrics import r2_score

def fake_func(p, x):
    f = np.poly1d(p) #多项式分布的函数
    return f(x)
    #return p*np.exp(-p*x)

#残差函数
def residuals(p, y, x):
    return y - fake_func(p, x)
#残差函数

def getBG(ch,timeSp,lenBin,c,MeasDesc_GlobalResolution,T0=0,Tlen=-1,axBG=None, axR2=None):
    c.execute("select TimeTag from fretData_"+ch+" ORDER BY TimeTag limit 1")
    #c.execute('SELECT * FROM stocks WHERE symbol=?', t)
    t1= c.fetchone()[0]

    hasData=True
    buf = array("d")
    r2buf = array("d")
    timeline=array("d")
    #idxbuf=0
    t0=int(T0/MeasDesc_GlobalResolution)
    tt1=max([t1,t0])
    if Tlen<=0:
        Tlen=-T0-100
    tEnd=int((T0+Tlen)/MeasDesc_GlobalResolution)  
    while hasData:
        t2=1+(int(timeSp/MeasDesc_GlobalResolution)+tt1)
        tt2=min([t2,tEnd+1])
        if tt2>=tt1:
            sql="select TimeTag from fretData_"+ch+" where TimeTag >= ? and TimeTag < ? ORDER BY TimeTag"
            c.execute(sql,(tt1,tt2))
        elif tEnd>0:
            hasData=False
            break;
        else:
            #print("to end")
            sql="select TimeTag from fretData_"+ch+\
                " where TimeTag >= ? and TimeTag < ? ORDER BY TimeTag"
            c.execute(sql,(tt1,t2))
        data=c.fetchall()
        if len(data)<5:
            hasData=False
            break;
        TimeTag=np.zeros(len(data)-1)
        for idxd in range(len(data)-1):
            TimeTag[idxd]=data[idxd+1][0]-data[idxd][0]
        #print(TimeTag[1:10])
        #print((t1,t2))
        TimeTag=TimeTag*MeasDesc_GlobalResolution
        [hist,bin]=np.histogram(TimeTag, bins=lenBin)
        #hist = np.log(hist)
        x=array("d")
        y=array("d")
        iy=0
        for vh in hist:
            if vh>0:
                y.append(np.log(vh))
                x.append(bin[iy])
            iy=iy+1
        #x=bin[0:lenbin]
        #y=np.log(hist)#/(timeSp/lenbin)

        m = 2  #多项式的次数
        #先随机产生一组多项式分布的参数
        p0 = np.random.randn(m)
        plsq = leastsq(residuals, p0, args=(y, x))
        #print(plsq)
        coefficient_of_dermination = r2_score(fake_func(plsq[0],x),y)
        r2buf.append(coefficient_of_dermination)
        buf.append(-1*plsq[0][0])
        timeline.append(tt1*MeasDesc_GlobalResolution)
        tt1=t2
        #hasData=False
    if axBG!=None:
        axBG.plot(np.frombuffer(timeline, dtype=np.double),np.frombuffer(buf, dtype=np.double), label='BG_'+ch+"(cps)")
    if axR2!=None:
        axR2.plot(np.frombuffer(timeline, dtype=np.double), np.frombuffer(r2buf, dtype=np.double), label='r2_BG_'+ch)
    #显示R2 拟合结果
    #bgRate = collections.namedtuple('bgRate', ['time', 'bgrate'])
    #br=bgRate(timeline,buf)
    #return br
    bgRate=dict({'time':timeline,'bgrate':buf})
    return bgRate


def calcBGrate(dbname,timeSp=20.0,lenbin=300,T0=0,Tlen=-1,axBG=None, axR2=None, chs=[]):    
    if (Tlen>=0 and Tlen<timeSp) or T0<0:
        print('Tlen must bigger than timeSp, T0 must >=0')
        return -1
    conn = sqlite3.connect(dbname)
    #timeSp=20
    #lenbin=300
    c = conn.cursor()
    c.execute("SELECT value FROM fretData_Var where varName='MeasDesc_GlobalResolution'")
    MeasDesc_GlobalResolution= c.fetchone()[0] #6.2497500099996e-08
    c.execute("SELECT value FROM fretData_Var where varName='MeasDesc_Resolution'")
    DelayResolution= c.fetchone()[0] #6.2497500099996e-08

    if len(chs)==0:
        chs=["DexAem","DexDem","AexAem","All"]

    #bgRateCh = collections.namedtuple('bgRateCh', ['ch', 'bgRate'])
    brcd=dict()
    for ch in chs:
        brcd[ch]=getBG(ch,timeSp,lenbin,c,MeasDesc_GlobalResolution,T0,Tlen,axBG, axR2)
        #.append(bgRateCh(ch,getBG(ch,timeSp,lenbin,c)))    
    brcd["SyncResolution"]=MeasDesc_GlobalResolution
    brcd["DelayResolution"]=DelayResolution
    brcd["T0"]=T0
    brcd["Tlen"]=Tlen    
    conn.close()
    return brcd
if __name__ == '__main__':
    #calcBGrate('/home/liuk/sf/oc/data/38.sqlite')
    import matplotlib.pyplot as pl
    f, (axBG, axR2) = pl.subplots(1,2)
    axBG.set_prop_cycle(cycler('color', ['red', 'black', 'yellow','blue']))
    axR2.set_prop_cycle(cycler('color', ['red', 'black', 'yellow','blue']))
    br=calcBGrate('/home/liuk/proj/data/LS35_RSV86C224C.sqlite',timeSp=30,lenbin=100,T0=10.0,Tlen=600.0,axBG=axBG,axR2=axR2)
    axBG.legend(loc='best')
    axR2.legend(loc='best')
    pl.show()