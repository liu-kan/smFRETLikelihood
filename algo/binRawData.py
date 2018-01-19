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
                    ,"s/s\t", end='',flush=True)
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
        print(" readed ch ",ch)
    if nowpercent<0.95:
        print("!!!!!span的合理范围 >= 7 !!!!!!!!")
    conn.close()
    binData["SyncResolution"]=bgrate["SyncResolution"]
    binData["DelayResolution"]=bgrate["DelayResolution"]
    binData['binMs']=binMs
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
def burstFilterByBin(binData,dddaaaT,burstPhotonNum=30):
    fDD=dddaaaT[0];fDA=dddaaaT[1];fAA=dddaaaT[2]
    T0=0
    fDex=dddaaaT[3]
    winStDA=[]
    winStDD=[]
    winStAA=[]
    chAll=binData['chs']['All']
    lenbin=len(chAll['ntag'])
    toBeDel=np.ndarray(lenbin,dtype=np.int32)
    for i in range(lenbin):
        chs=chAll['chl'][i]
        leninbin=len(chs)
        for j in range(leninbin):
            if chs[j]==1:
                winStDA.append((i,j))
            elif chs[j]==2:
                winStDD.append((i,j))
            elif chs[j]==3:
                winStAA.append((i,j))
        if len(winStAA)< fAA:
            toBeDel[i]=True
        elif len(winStDA)<fDA:
            toBeDel[i]=True
        elif len(winStDD)<fDD:
            toBeDel[i]=True
        elif len(winStDA)+len(winStDD)<fDex:
            toBeDel[i]=True
        else:
            toBeDel[i]=False
        winStDA=[]
        winStDD=[]
        winStAA=[]
    burst=[]
    burstData=[]
    for i in range(lenbin):
        if not toBeDel[i]:
            burst.append(i)
        else:
            if len(burst)>0:
                if numPhotonInBurst(binData,burst)<burstPhotonNum:
                    toBeDel[slice(burst[0],burst[-1]+1)]=True
                else:
                    burstData.append((burst[0],burst[-1]))
            burst=[]
    print("toBeDel",np.sum(toBeDel))
    chAll['markDel']=toBeDel
    chAll['burst']=burstData

def numPhotonInBurst(binData,burst):
    chAll=binData['chs']['All']
    ntags=chAll['ntag'][slice(burst[0],burst[-1]+1)]
    return(np.sum(ntags))

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
                        toBeDel.extend(winStDA)
                    winStDA=[(i,j)]
                else:
                    winStDA.append((i,j))
                daT0=chAll['timetag'][i][j]
            elif chs[j]==2:
                ddDetaT=chAll['timetag'][i][j]-ddT0
                if ddDetaT>=binData['chs']['All']['photonDDdiffMean']:
                    lenWin=len(winStDD)
                    if lenWin>0 and lenWin<fDD:
                        toBeDel.extend(winStDD)
                    winStDD=[(i,j)]
                else:
                    winStDD.append((i,j))
                ddT0=chAll['timetag'][i][j]
            elif chs[j]==3:
                aaDetaT=chAll['timetag'][i][j]-aaT0
                if aaDetaT>=binData['chs']['All']['photonAAdiffMean']:
                    lenWin=len(winStAA)
                    if lenWin>0 and lenWin<fAA:
                        toBeDel.extend(winStAA)
                    winStAA=[(i,j)]
                else:
                    winStAA.append((i,j))
                aaT0=chAll['timetag'][i][j]

    formBurst(binData,toBeDel)
    addChAllf(binData,toBeDel)
    binData['chs'].pop('All',None)
    binData['chs']['All']=binData['chs'].pop('Allf',None)
def formBurst(binData,toBeDel):
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
    wtimetag=[]
    wdtime=[]
    wchl=[]
    burstTag=[]
    i=0
    for i in range(lenBins):
        lenbin=allCh['ntag'][i]
        if i!=toBeDel[ib][0]:
            wtimetag.extend( allCh['timetag'][i])
            wdtime.extend(allCh['dtime'][i])
            wchl.extend(allCh['chl'][i])
        else:
            for j in range(lenbin):
                if j!=toBeDel[ib][1]:
                    wtimetag.append(allCh['timetag'][i][j])
                    wdtime.append(allCh['dtime'][i][j])
                    wchl.append(allCh['chl'][i][j])
                else:
                    if ib<lendel-1:
                        ib=ib+1
                        nntag=len(wtimetag)
                        if nntag>0:
                            timetag.append(wtimetag)
                            dtime.append(wdtime)
                            chl.append(wchl)
                            fretE.append(allCh['e'][i])
                            fretZ.append(allCh['z'][i])
                            fretS.append(allCh['s'][i])
                            ntag.append(nntag)
                            lifetime.append(allCh['lifetime'][i])
                            burstTag.append((wtimetag[0],wtimetag[-1]))
                            wtimetag=[]
                            wdtime=[]
                            wchl=[]

        now=100.0*i/lenBins
        if now>p100:
            print("\r>> add finished:","%2.5f" % p100,'%', end='',flush=True)
            p100=now+5
    binData['chs']['AllBurst']=dict({'timetag':timetag,'dtime':dtime,'chl':chl,\
                    'binMs':binMs ,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime,\
                    'ntag':ntag,'burstTag':burstTag})
    print()


def groupDel(binData,toBeDel):
    toBeDel=sorted(toBeDel)
    allCh=binData['chs']['All']
    lenBins=len(allCh['timetag'])
    ib=0;idx=[0,0]
    group=[]
    groupDel=[]
    lendel=len(toBeDel)
    groupNum=[]
    for i in range(lenBins):
        lenbin=allCh['ntag'][i]
        if i!=toBeDel[ib][0]:
            idx[0]=i
            idx[1]=lenbin-1
            continue
        for j in range(lenbin):
            if j!=toBeDel[ib][1]:
                idx[0]=i
                idx[1]=j
            else:
                if ib<lendel-1:
                    ib=ib+1
                if len(groupDel)>0:
                    if (groupDel[-1][0]==i and groupDel[-1][1]+1==j) or (groupDel[-1][0] \
                        ==i-1 and groupDel[-1][1]==allCh['ntag'][i-1]-1 and j==0):
                        groupDel.append((i,j))
                    else:
                        group.append(groupDel[0])
                        group.append(groupDel[-1])
                        groupNum.append(len(groupDel))
                        groupDel=[(i,j)]
                else:
                    groupDel=[(i,j)]
    print("groupDel#",len(group),"max g#",max(groupNum))

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
        print("AAmean",binDataCh['AAmean'])
    binDataCh['photondiffMean']=(etag-stag)/np.sum(sumBin)
    binDataCh['photonEffTime']=(etag-stag)*binData["SyncResolution"]

def statsDelayTime(binData,binNum=200,chl="Both",bin0=0,binLen=-1):
    lenBin=len(binData['chs']["All"]['chl'])
    a=array('d')
    d=array('d')
    if bin0<0:
        bin0=0
    if binLen<0 or binLen+bin0>=lenBin:
        binEnd=lenBin
    else:
        binEnd=binLen+bin0
    print((bin0,binEnd),np.percentile(binData['chs']["All"]['ntag'],99))
    print((bin0,binEnd),np.percentile(binData['chs']["All"]['ntag'],90))
    for i in range(bin0,binEnd):
        w=binData['chs']["All"]['ntag'][i]
        for idxd in range(w):
            if binData['chs']["All"]['chl'][i][idxd]==1 \
             or binData['chs']["All"]['chl'][i][idxd]==3: #DA or AA
                a.append(1e9*binData['chs']["All"]['dtime'][i][idxd]\
                *binData["DelayResolution"])
            else:
                d.append(1e9*binData['chs']["All"]['dtime'][i][idxd]\
                *binData["DelayResolution"])
    histA,binA= np.histogram(a, binNum)
    histD,binD= np.histogram(d, binNum)
    if chl=="Both":
        return [histA,histD],[binA,binD]
    elif chl=="A":
        return [histA],[binA]
    else:
        return [histD],[binD]
def statsPhotonDiff(binData,chs=[]):
    hist=[];bin=[]
    if len(chs)<1:
        for chl in binData['chs'].keys():
            realStatsPhotonTagDiff(binData,chl)
            histA,binA= np.histogram(binData['chs'][chl]['photonTagDiff'], 2000)
            hist.append(histA);bin.append(binA)
    else:
        for chl in chs:
            realStatsPhotonTagDiff(binData,chl)
            histA,binA= np.histogram(binData['chs'][chl]['photonTagDiff'], 2000)
            hist.append(histA);bin.append(binA)
    return hist,bin
def realStatsPhotonTagDiff(binData,chl):
    lenBins=len(binData['chs'][chl]['timetag'])
    lastpt=0;lastA=np.zeros(1)#;rA=np.zeros(1)
    rA=[]
    for i in range(lenBins):
        cpd=np.diff(binData['chs'][chl]['timetag'][i])
        if len(binData['chs'][chl]['timetag'][i])>0:
            if len(cpd)>0:
                if lastpt>0 and i>0:
                    lastA[0]=binData['chs'][chl]['timetag'][i][0]-lastpt
                    cpd=np.concatenate((lastA,cpd),axis=0)
            else:
                cpd=np.zeros(1);cpd[0]=binData['chs'][chl]['timetag'][i][0]-lastpt
            lastpt=binData['chs'][chl]['timetag'][i][-1]
            # rA=np.concatenate(( rA,cpd),axis=0)
            rA.extend(cpd.tolist() )
    binData['chs'][chl]['photonTagDiff']=rA

def exportBin2PQdotDat(binData,path):
    lenBins=len(binData['chs']['All']['timetag']);f=[]
    if lenBins>0:
        for i in range(1,5):
            f.append(open(path+str(binData['chs']['All']['binMs'])+"_"+str(i)+".dat",'w'))
    for i in range(lenBins):
        time=i*(1e-3)*binData['chs']['All']['binMs']
        unique, counts = np.unique(binData['chs']['All']['chl'][i], return_counts=True)
        rchdict=dict(zip(unique, counts))
        chdict=dict({1:0,2:0,3:0,4:0})
        for k in rchdict.keys():
            chdict[k]=rchdict[k]
        for chk in range(1,4):
            f[chk-1].write('%.3f\t%d'%(time,chdict[chk])+'\r\n')
        f[3].write('%.3f\t%d'%(time,chdict[1]+chdict[2])+'\r\n')
    for fi in f:
        fi.close()



if __name__=='__main__':
    dbname="/dataB/smfretData/irf/alexa488_IRF_32MHz_PIE_3KCPS.sqlite"
    br=BGrate.calcBGrate(dbname,20,400)#,T0=0.0,Tlen=600)
    binData=binRawData(br,dbname,1)
    hi,bi=statsDelayTime(binData,1000,"D")#,bin0=100,binLen=2)
    statsBins(binData)
    print("IRF photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    dbname="/dataB/smfretData/rsv86c224c.sqlite"
    br=BGrate.calcBGrate(dbname,20,400,T0=6,Tlen=615.56797)
    binData=binRawData(br,dbname,1)
    h,b=statsDelayTime(binData,1000,"D")#,bin0=100,binLen=2)
    import sys
    sys.path.append('./ui')
    from histBar_stacked import plotStackedHist
    plotStackedHist([h[0],hi[0]],[b[0],bi[0]],log=True)
    #plotStackedHist((hi),(bi),log=True)
    statsBins(binData)
    print("Data photonEffTime:",binData['chs']['DexDem']["photonEffTime"])
    # burstFilter(binData,[5.1,4.1,3.1])
    print("DexDem photondiffMean(us)",\
    binData['chs']['DexDem']['photondiffMean']*\
            binData["SyncResolution"]*1e6)
    savefn='/dataZ1/tmp/lineardiub/'+\
        dbname.split('/')[-1].split('.')[-2]+'_gd'+".pickle"
    # with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([burst], f,protocol=-1)

    # h,b=statsPhotonDiff(binData,['All','DexDem','DexAem','AexAem'])
    # import sys
    # sys.path.append('./ui')
    # from histBar_stacked import plotSubHist
    # plotSubHist(h,b)
    # savePath='/dataZ1/tmp/lineardiub/'+\
    #     dbname.split('/')[-1].split('.')[-2]+'_'
    # exportBin2PQdotDat(binData,savePath)
