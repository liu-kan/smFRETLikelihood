#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 10:31:42 2016

@author: liuk
"""

import BGrate
import sqlite3
from array import array
#import collections
#import matplotlib.pyplot as pl


def getBGrateAtT(bgra, ch, timing):
    lent = len(bgra[ch]['time'])
    ridxt=-1
    if timing<bgra[ch]['time'][0] or timing>bgra[ch]['time'][lent-1]:
        return 0
    for idxt in range(lent-1):
        if timing<bgra[ch]['time'][idxt+1] and timing>=bgra[ch]['time'][idxt]:
            ridxt=idxt
            break
    if ridxt>=0:
        return bgra[ch]['bgrate'][ridxt]
    else:
        return 0
def data2Dcol(data,row0,row1,col):
    rows=slice(row0,row1)
    if row1==0:
        rows=slice(row0,None)
    d=data[rows]
    r=[]
    for line in d:
        r.append(line[col])
    return r


def findBurst(br,dbname,chs,continuousPhoton=30,F=5):
#dbname='/home/liuk/prog/ptu/data/30.sqlite'
#dbname='E:/doc/proj/ptu/data/37.sqlite'

    T0=br["T0"]
    Tlen=br["Tlen"]    
    conn = sqlite3.connect(dbname)
    #timeSp=20
    #lenbin=300
    c = conn.cursor()
    #chs=["DexAem","DexDem","AexAem"]

    ##滑动窗口寻找相邻的光子
    #chs=["All"]
    burst=dict()
    burst['chs']=dict()
    blockNum=1000000
    if Tlen<=0:
        Tlen=-T0-100    
    tEnd=int((T0+Tlen)/br["SyncResolution"])
    for ch in chs:
        c.execute("select TimeTag from fretData_"+ch+" ORDER BY TimeTag limit 1")
        #c.execute('SELECT * FROM stocks WHERE symbol=?', t)
        t1= c.fetchone()[0]-1
        
        T1=int(T0/br["SyncResolution"])
        tt1=max(t1,T1-1)
        hasData=False
        if t1>=0:
            hasData=True
        buf = array("d")
        timeline=array("d")
        burstW=array("d")
        ntag=array("l")
        fretE=array("d")
        fretS=array("d")
        fretZ=array("d")
        lifetime=array("d")
        timetag=[]
        dtime=[]
        chl=[]
        #idxbuf=0
        while hasData:
            if tEnd>0:
                sql="select TimeTag,Dtime,ch from fretData_"+ch+" where TimeTag >= ? and TimeTag < ? ORDER BY TimeTag limit ?"
                c.execute(sql,(tt1,tEnd,blockNum))
            else:
                sql="select TimeTag,Dtime,ch from fretData_"+ch+" where TimeTag >= ?  ORDER BY TimeTag limit ?"
                c.execute(sql,(tt1,blockNum))                
            data=c.fetchall()
            lendata=len(data)

            if lendata<continuousPhoton*2:
                hasData=False
                break;
            #timetag 如果有相同的可能掉数据
            i=0
            while i <(lendata-continuousPhoton):
                jugeD = 0
                etiming=0
                aLocalRate=0
                bigSNR=False
                #print(jugeD)
                while jugeD<lendata-continuousPhoton-i:
                    ts=(data[i+jugeD+continuousPhoton-1][0]-data[i+jugeD][0])*br["SyncResolution"]
                    localRate=continuousPhoton/ts
                    timing=data[i+jugeD+(continuousPhoton>>1)][0]*br["SyncResolution"]
                    if localRate<F*getBGrateAtT(br,ch,timing):
                        break;
                    else:
                        bigSNR=True
                        jugeD=jugeD+1

                if bigSNR:
                    #print(("i",i))
                    ntag.append(continuousPhoton+jugeD)
                    #print(data[i:i+continuousPhoton+jugeD][0])
                    #hasData=False
                    
                    timetag.append(data2Dcol(data,i,i+continuousPhoton+jugeD,0))
                    dtime.append(data2Dcol(data,i,i+continuousPhoton+jugeD,1))
                    chl.append(data2Dcol(data,i,i+continuousPhoton+jugeD,2))
                    # timetag.append(df[i:i+continuousPhoton+jugeD][0])
                    # dtime.append(df[i:i+continuousPhoton+jugeD][1])
                    # chl.append(df[i:i+continuousPhoton+jugeD][2])
                    etiming=data[i+continuousPhoton+jugeD-1][0]*br["SyncResolution"]
                    stiming=data[i][0]*br["SyncResolution"]
                    aLocalRate=(jugeD+continuousPhoton)/(etiming-stiming)
                    buf.append(aLocalRate)
                    timeline.append(stiming)
                    burstW.append(etiming-stiming)
                    fretS.append(-1)
                    fretE.append(-1)
                    fretZ.append(-1)
                    lifetime.append(-1)
                    i=i+continuousPhoton+jugeD
                else:
                    i=i+1
            tt1=data[-continuousPhoton][0]  #todo 检查数组内是否有重复
            #print(("t1",t1))
        #fig = pl.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(timeline,buf)
        #cburst = collections.namedtuple('burst', ['ntag', 'burstW','timetag','dtime','chl','e','s'])
        #burst[ch]=cburst(ntag,burstW,timetag,dtime,chl,fretE,fretS)
        cburst=dict({'ntag':ntag, 'burstW':burstW,'timetag':timetag,'dtime':dtime,\
                    'chl':chl,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime})
        burst['chs'][ch]=cburst
    conn.close()
    burst["SyncResolution"]=br["SyncResolution"]
    burst["DelayResolution"]=br["DelayResolution"]
    return burst

if __name__ == '__main__':
    import time

    dbname='E:/liuk/proj/ptu/data/46.sqlite'
    dbname='/home/liuk/proj/data/86c_224c.sqlite'
    br=BGrate.calcBGrate(dbname,20,400,chs=["All"],T0=10.0,Tlen=60.0)
    start = time.time()
    findBurst(br,dbname,["All"])
    end = time.time()
    print(end - start)
