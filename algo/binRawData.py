
import BurstSearch
import BGrate
import sqlite3
import numpy as np
import time
from array import array

def findTagLastLessThanIdx(data,lessThan):
    #注意！！！timetag 如果有相同的可能掉数据
    lendata=len(data)
    r=-1
    for i in range(lendata):
        if data[i][0]==lessThan:
            r=i
        elif data[i][0]>lessThan:
            r=i-1
            break
    return r

def binRawData(  bgrate, dbname, binMs = 1,chs=["DexAem","DexDem","AexAem","All"] ):
    binTag=(binMs*1e-3)/bgrate["SyncResolution"] #每个bin对应的timetag值
    print(binTag)
    span=int(10/binMs)+1
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    #暂时只处理All通道
    binData=dict()
    binData['chs']=dict()
    nowpercent=0
    T0=bgrate["T0"]
    Tlen=bgrate["Tlen"]
    c.execute("select min(TimeTag) from fretData_All")
    idx= c.fetchone()[0] # Start TimeTag
    if idx!=None:
        hasData=True
    c.execute("select max(TimeTag) from fretData_All")
    TimeLast= c.fetchone()[0]
    tt0=T0/bgrate["SyncResolution"]
    ttEnd=-1
    if Tlen>0:
        ttEnd=(T0+Tlen)/bgrate["SyncResolution"]
    if ttEnd>0:
        TimeLast=min(TimeLast,ttEnd)
    TimeStart=max(idx,tt0)    
    for ch in chs:
        timetag=[]
        dtime=[]
        chl=[]
        fretE=array("d")
        fretS=array("d")
        fretZ=array("d")
        lifetime=array("d")      
        ntag=array("l")  
        if ch=='All':
            ndaTag=array('l')
            nddTag=array('l')
            naaTag=array('l')
            nadTag=array('l')
        hasData=False
        c.execute("select min(TimeTag) from fretData_"+ch)
        idx= c.fetchone()[0] # Start TimeTag
        if idx!=None:
            hasData=True
            idx=TimeStart
        percent=0
        start = time.time()
        timeIdx=0
        while hasData:        
            nowpercent= (idx-TimeStart)/(TimeLast-TimeStart)
            if nowpercent>=percent+0.1:
                now = time.time()
                print(nowpercent,'% speed =',(idx-TimeStart)*bgrate["SyncResolution"]/(now-start),"s/s")
                percent=nowpercent
                                     
            if (idx+binTag*span)>=TimeLast:
                sql="select TimeTag,dtime,ch from fretData_"+ch+" where TimeTag>=? and TimeTag<=? ORDER BY TimeTag" 
                c.execute(sql,(idx,TimeLast))
                hasData=False
            else:
                sql="select TimeTag,dtime,ch from fretData_"+ch+" where TimeTag>=? and TimeTag<? ORDER BY TimeTag" 
                c.execute(sql,(idx,idx+binTag*span))
            data=c.fetchall()
            r0=0
            r1=0
            # r00=0
            for i in range(span):  
                # if i==0:
                #     r00=r0   
                r1=findTagLastLessThanIdx(data,TimeStart+timeIdx*binTag*span+(i+1)*binTag)    
                if r1<0 and i!=span-1:
                    continue   
                elif i==span:
                    r1=0
                tt=BurstSearch.data2Dcol(data,r0,r1,0)
                numtag=len(tt)
                if numtag>0:
                    timetag.append(tt)
                    dtime.append(BurstSearch.data2Dcol(data,r0,r1,1))
                    cha=BurstSearch.data2Dcol(data,r0,r1,2)
                    chl.append(cha)   
                    lifetime.append(0)
                    fretE.append(-1)
                    fretS.append(-1)
                    fretZ.append(-1)             
                    ntag.append(numtag)
                    if ch=='All':
                        ndaTag.append(cha.count(1))
                        ddf=cha.count(2)
                        aaf=cha.count(3)
                        adf=cha.count(4)
                r0=r1
            idx=idx+binTag*span
            # print(data[r00][0],data[r0][0])
            # print("=============")
            timeIdx+=1
        binData['chs'][ch]=dict({'timetag':timetag,'dtime':dtime,'chl':chl,'binMs':binMs \
                        ,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime,'ntag':ntag})
    if nowpercent<0.95:
        print("!!!!!span的合理范围 >= 7 !!!!!!!!")
    conn.close()
    binData["SyncResolution"]=bgrate["SyncResolution"]
    binData["DelayResolution"]=bgrate["DelayResolution"]    
    return binData
def countBin(binData,chl):
    if chl in binData['chs']:
        binDataCh=binData['chs'][chl]
        lenbin=len(binDataCh['timetag'])
        print("bin#:",lenbin)
        sumBin=np.empty(lenbin)
        for i in range(lenbin):
            sumBin[i]=len(binDataCh['timetag'][i])
        return sumBin
    else:
        return None
def statsBins(binData):
    for chl in binData['chs'].keys():
        binDataCh=binData['chs'][chl]
        lenbin=len(binDataCh['timetag'])
        stag=binDataCh['timetag'][0][0]
        etag=binDataCh['timetag'][lenbin-1][0]
        sumBin=np.zeros(lenbin)
        for i in range(lenbin):
            sumBin[i]=binDataCh['ntag'][i]
            if i==lenbin-1 and sumBin[i]>1:
                etag=binDataCh['timetag'][lenbin-1][binDataCh['ntag'][i]-1]
        binDataCh['mean']=np.mean(sumBin)
        binDataCh['std']=np.std(sumBin)
        binDataCh['photondiffMean']=(etag-stag)* \
            binData["SyncResolution"]/np.sum(sumBin)
def statsDelayTime(binData):
    lenBin=len(binData['chs']["All"]['chl'])
    a=array('d')
    d=array('d')
    for i in range(lenBin):
        w=binData['chs']["All"]['ntag'][i]
        for idxd in range(w): 
            if binData['chs']["All"]['chl'][i][idxd]==1 \
             or binData['chs']["All"]['chl'][i][idxd]==3: #DA or AA
                a.append(1e9*binData['chs']["All"]['dtime'][i][idxd]\
                *binData["DelayResolution"])
            else:
                d.append(1e9*binData['chs']["All"]['dtime'][i][idxd]\
                *binData["DelayResolution"])
    histA,binA= np.histogram(a, 200)                 
    histD,binD= np.histogram(d, 200)
    return [histA,histD],[binA,binD]
if __name__=='__main__':
    dbname="../data/rsv86c224c.sqlite"
    br=BGrate.calcBGrate(dbname,20,400)#,T0=10.0,Tlen=200)
    binData=binRawData(br,dbname,2,["DexDem"])
    # h,b=statsDelayTime(binData)
    # import sys
    # sys.path.append('./ui')
    # from histBar_stacked import plotStackedHist
    # plotStackedHist(h,b)
    statsBins(binData)
    print("DexDem photondiffMean(ms)",binData['chs']['DexDem']['photondiffMean']* \
        1e3)

    
    