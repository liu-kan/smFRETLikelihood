# -*- coding: utf-8 -*-
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
    print("binTag:",binTag)
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
        # if ch=='All':
        #     ndaTag=array('l')
        #     nddTag=array('l')
        #     naaTag=array('l')
        #     nadTag=array('l')
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
                speed=(idx-TimeStart)*bgrate["SyncResolution"]/(now-start)
                print('\r>>',"%2.5f"%(nowpercent*100),'% speed =',"%3.5f"%speed \
                    ,"s/s", end='',flush=True)
                percent=nowpercent
                                     
            if (idx+binTag*span)>=TimeLast:
                sql="select TimeTag,dtime,ch from fretData_"+ch+\
                    " where TimeTag>=? and TimeTag<=? ORDER BY TimeTag" 
                c.execute(sql,(idx,TimeLast))
                hasData=False
            else:
                sql="select TimeTag,dtime,ch from fretData_"+ch+\
                    " where TimeTag>=? and TimeTag<? ORDER BY TimeTag" 
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
                    # if ch=='All':
                    #     ndaTag.append(cha.count(1))
                    #     ddf=cha.count(2)
                    #     aaf=cha.count(3)
                    #     adf=cha.count(4)
                r0=r1
            idx=idx+binTag*span
            # print(data[r00][0],data[r0][0])
            # print("=============")
            timeIdx+=1
        binData['chs'][ch]=dict({'timetag':timetag,'dtime':dtime,'chl':chl,\
                    'binMs':binMs ,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime,\
                    'ntag':ntag})
        print("\t readed ch ",ch)
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
def burstFilter(binData,dddaaaT):
    fDD=dddaaaT[0];fDA=dddaaaT[1];fAA=dddaaaT[2]
    T0=0
    daT0=0;winStDA=[]
    ddT0=0;winStDD=[]
    aaT0=0;winStAA=[]
    adT0=0;winStAD=[]
    toBeDel=[]
    chAll=binData['chs']['All']
    lenbin=len(chAll['timetag'])
    for i in range(lenbin):
        chs=chAll['chl'][i]
        leninbin=len(chs)
        for j in range(leninbin):
            if chs[j]==1:
                daDetaT=chAll['timetag'][i][j]-daT0
                '''
                间隔比均值大 肯定属于新的burst 开新窗口
                '''
                if daDetaT>=binData['chs']['All']['photonDAdiffMean']:
                    lenWin=len(winStDA)
                    '''
                    判断旧窗口的内容是否删除：
                    判断之前已经存储的间隔和要求的点的数目，如果低于阈值就标记为删除
                    如果高于阈值表明前面的点已经符合要求，把现在的点加入判别窗口，这样
                    如果下一个点的间隔也大于均值这个点就会被删除，如果不大于就继续判断
                    '''
                    if lenWin>0 and lenWin<fDA:
                        for k in winStDA:
                            toBeDel.append(k)
                    winStDA=[(i,j)]
                else:
                    winStDA.append((i,j))
                daT0=chAll['timetag'][i][j]
            elif chs[j]==2:
                ddDetaT=chAll['timetag'][i][j]-ddT0
                if ddDetaT>=binData['chs']['All']['photonDDdiffMean']:
                    lenWin=len(winStDD)
                    if lenWin>0 and lenWin<fDD:
                        for k in winStDD:
                            toBeDel.append(k)
                    winStDD=[(i,j)]
                else:
                    winStDD.append((i,j))
                ddT0=chAll['timetag'][i][j]
            elif chs[j]==3:
                aaDetaT=chAll['timetag'][i][j]-aaT0
                if aaDetaT>=binData['chs']['All']['photonAAdiffMean']:
                    lenWin=len(winStAA)
                    if lenWin>0 and lenWin<fAA:
                        for k in winStAA:
                            toBeDel.append(k)
                    winStAA=[(i,j)]
                else:
                    winStAA.append((i,j))
                aaT0=chAll['timetag'][i][j]
            # else chs[j]==2:
            #     adDetaT=chAll['timetag'][i][j]-adT0
            #     if adDetaT>=binData['AexDem']['photondiffMean']:
            #         lenWin=len(winStAD)
            #         if lenWin>0 and lenWin<fAD:
            #             for k in winStAD:
            #                 toBeDel.append(k)
            #         winStAD=[(i,j)]
            #     else:
            #         winStAD.append((i,j))
            #     adT0=chAll['timetag'][i][j]    
    addChAllf(binData,toBeDel)      
    binData['chs'].pop('All',None)
    binData['chs']['All']=binData['chs'].pop('Allf',None)

def addChAllf(binData,toBeDel):
    toBeDel=sorted(toBeDel)    
    allCh=binData['chs']['All']
    binMs=allCh['binMs']
    timetag=[]
    dtime=[]
    chl=[]
    fretE=array("d")
    fretS=array("d")
    fretZ=array("d")
    lifetime=array("d")      
    ntag=array("l")  
    lenBins=len(allCh['timetag'])
    ib=0
    p100=1
    lendel=len(toBeDel)
    print(lendel)
    for i in range(lenBins):
        lenbin=allCh['ntag'][i]   
        wtimetag=[]
        wdtime=[]
        wchl=[]     
        if i!=toBeDel[ib][0]:
            wtimetag=allCh['timetag'][i]
            wdtime=allCh['dtime'][i]
            wchl=allCh['chl'][i]    
        else:        
            for j in range(lenbin):
                if j!=toBeDel[ib][1]:
                    wtimetag.append(allCh['timetag'][i][j])
                    wdtime.append(allCh['dtime'][i][j])
                    wchl.append(allCh['chl'][i][j])
                else:
                    if ib<lendel-1:
                        ib=ib+1                
        nntag=len(wchl)
        if nntag>0:
            timetag.append(wtimetag)
            dtime.append(wdtime)
            chl.append(wchl)
            fretE.append(allCh['e'][i])
            fretZ.append(allCh['z'][i])
            fretS.append(allCh['s'][i])
            ntag.append(nntag)
            lifetime.append(allCh['lifetime'][i])
        now=100.0*i/lenBins
        if now>p100:
            print("\r>> add finished:","%2.5f" % p100,'%', end='',flush=True)
            p100=now+5
    binData['chs']['Allf']=dict({'timetag':timetag,'dtime':dtime,'chl':chl,\
                    'binMs':binMs ,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime,\
                    'ntag':ntag})
    print()


def delRec(binData,toBeDel):
    toBeDel=sorted(toBeDel,reverse=True)
    lasti=int(toBeDel[0][0])
    len2beDel=len(toBeDel)
    i=0.0
    p100=1;
    for k in toBeDel:
        if len(binData['chs']['All']['timetag'][lasti])==0:
            del binData['chs']['All']['timetag'][lasti]
            del binData['chs']['All']['dtime'][lasti]
            del binData['chs']['All']['chl'][lasti]
            del binData['chs']['All']['e'][lasti]
            del binData['chs']['All']['s'][lasti]
            del binData['chs']['All']['z'][lasti]
            del binData['chs']['All']['lifetime'][lasti]
            del binData['chs']['All']['ntag'][lasti]
        del binData['chs']['All']['timetag'][k[0]][k[1]]
        del binData['chs']['All']['dtime'][k[0]][k[1]]
        del binData['chs']['All']['chl'][k[0]][k[1]]
        binData['chs']['All']['ntag'][k[0]]=binData['chs']['All']['ntag'][k[0]]-1
        lasti=k[0]
        i=i+1
        if i/len2beDel*100>p100:
            p100=i/len2beDel*100
            print("\r>> del finished:",p100, end='',flush=True)
            p100=p100+2
    print("filtered binData")   
    
def statsBins(binData,chs=[]):
    if len(chs)<1:
        for chl in binData['chs'].keys():
            realStatsBins(binData,chl)
    else:
        for chl in chs:
            realStatsBins(binData,chl)

def realStatsBins(binData,chl):
    binDataCh=binData['chs'][chl]
    lenbin=len(binDataCh['timetag'])
    stag=binDataCh['timetag'][0][0]
    etag=binDataCh['timetag'][lenbin-1][0]
    sumBin=np.ndarray(1)
    print("Ch:",chl,'has',lenbin,'bins')
    if chl=="All":
        sumBin=np.zeros((4,lenbin))
    else:
        sumBin=np.zeros(lenbin)
    for i in range(lenbin):
        if chl=="All":
            unique, counts = np.unique(binDataCh['chl'][i], return_counts=True)
            chdict=dict(zip(unique, counts))
            for chk in chdict.keys():
                sumBin[chk-1,i]=chdict[chk]
                if i==lenbin-1 and sumBin[chk-1,i]>1:
                    etag=binDataCh['timetag'][lenbin-1][binDataCh['ntag'][i]-1]
        else:
            sumBin[i]=binDataCh['ntag'][i]
            if i==lenbin-1 and sumBin[i]>1:
                etag=binDataCh['timetag'][lenbin-1][binDataCh['ntag'][i]-1]
    binDataCh['mean']=np.mean(sumBin)
    binDataCh['std']=np.std(sumBin)
    if chl=='All':        
        m=np.mean(sumBin,1)
        binDataCh['DAmean']=m[0]
        binDataCh['DDmean']=m[1]
        binDataCh['AAmean']=m[2]
        binDataCh['ADmean']=m[3]
        std=np.std(sumBin,1)
        binDataCh['DAstd']=std[0]
        binDataCh['DDstd']=std[1]
        binDataCh['AAstd']=std[2]
        binDataCh['ADstd']=std[3]                
        dm=(etag-stag)/np.sum(sumBin,1)
        binDataCh['photonDAdiffMean']=dm[0]
        binDataCh['photonDDdiffMean']=dm[1]
        binDataCh['photonAAdiffMean']=dm[2]
        binDataCh['photonADdiffMean']=dm[3]        
    binDataCh['photondiffMean']=(etag-stag)/np.sum(sumBin)
    
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
def statsPhotonDiff(binData,chs=[]):
    if len(chs)<1:
        for chl in binData['chs'].keys():
            realStatsBins(binData,chl)
    else:
        for chl in chs:
            realStatsBins(binData,chl)

if __name__=='__main__':
    dbname="../data/rsv86c224c.sqlite"
    br=BGrate.calcBGrate(dbname,20,400,T0=10.0,Tlen=200)
    binData=binRawData(br,dbname,2)
    # h,b=statsDelayTime(binData)
    # import sys
    # sys.path.append('./ui')
    # from histBar_stacked import plotStackedHist
    # plotStackedHist(h,b)
    statsBins(binData)
    burstFilter(binData,5.1,4.1,3.1)
    print("DexDem photondiffMean(us)",\
    binData['chs']['DexDem']['photondiffMean']*\
            binData["SyncResolution"]*1e6)

    
    