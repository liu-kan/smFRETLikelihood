# -*- coding: utf-8 -*-
#作者:liuk
#博客:liuk.io
#邮箱:liukan@126.com

import numpy as np
#import scipy as sp
from scipy.optimize import leastsq #这里就是我们要使用的最小二乘的函数
#import matplotlib.pyplot as pl
import sqlite3
from array import array
from cycler import cycler
#import collections

def fake_func(p, x):
    f = np.poly1d(p) #多项式分布的函数
    return f(x)
    #return p*np.exp(-p*x)

#残差函数
def residuals(p, y, x):
    return y - fake_func(p, x)
#残差函数

def getBG(ch,timeSp,lenBin,c,MeasDesc_GlobalResolution):
    c.execute("select TimeTag from fretData_"+ch+" ORDER BY TimeTag limit 1")
    #c.execute('SELECT * FROM stocks WHERE symbol=?', t)
    t1= c.fetchone()[0]

    hasData=True
    buf = array("d")
    timeline=array("d")
    #idxbuf=0
    while hasData:
        t2=1+(int(timeSp/MeasDesc_GlobalResolution)+t1)
        sql="select TimeTag from fretData_"+ch+" where TimeTag >= ? and TimeTag < ? ORDER BY TimeTag"
        c.execute(sql,(t1,t2))
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
        m = 2  #多项式的次数
        #先随机产生一组多项式分布的参数
        p0 = np.random.randn(m)
        plsq = leastsq(residuals, p0, args=(y, x))
        #print(plsq)
        buf.append(-1*plsq[0][0])
        timeline.append(t1*MeasDesc_GlobalResolution)
        t1=t2
        #hasData=False
    #pl.plot(np.frombuffer(timeline, dtype=np.double), np.frombuffer(buf, dtype=np.double), label='BG_'+ch+"(cps)")    
    bgRate=dict({'time':timeline,'bgrate':buf})
    return bgRate


def calcBGrate(dbname,timeSp=20,lenbin=300):    
    conn = sqlite3.connect(dbname)
    #timeSp=20
    #lenbin=300
    c = conn.cursor()
    c.execute("SELECT value FROM fretData_Var where varName='MeasDesc_GlobalResolution'")
    MeasDesc_GlobalResolution= c.fetchone()[0] #6.2497500099996e-08
    c.execute("SELECT value FROM fretData_Var where varName='MeasDesc_Resolution'")
    DelayResolution= c.fetchone()[0] #2.50000003337858e-11

    #print (MeasDesc_GlobalResolution)
    chs=["DexAem","DexDem","AexAem","All"]
    #fig = pl.figure()
    #ax = fig.add_subplot(111)
    #ax.set_prop_cycle(cycler('color', ['red', 'black', 'yellow','blue']))#set_color_cycle(['red', 'black', 'blue'])    
    brcd=dict()
    for ch in chs:
        brcd[ch]=getBG(ch,timeSp,lenbin,c,MeasDesc_GlobalResolution)        
    #pl.legend(loc='best')
    #pl.show()
    brcd["SyncResolution"]=MeasDesc_GlobalResolution
    brcd["DelayResolution"]=DelayResolution
    conn.close()
    return brcd
if __name__ == '__main__':
    #calcBGrate('/home/liuk/sf/oc/data/38.sqlite')
    br=calcBGrate('/home/liuk/proj/data/RSV89C224C.sqlite')
