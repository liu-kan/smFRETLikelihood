try:
    import algo.BurstSearch as BurstSearch
    import algo.BGrate as BGrate
except ImportError:
    import BurstSearch
    import BGrate
import sqlite3

def findLastLessThanIdx(data,lessThan):
    #注意！！！timetag 如果有相同的可能掉数据
    pass
def binRawData( dbname, bgrate, binMs = 1 ):
    binTag=binMs*1e-3/bgrate["SyncResolution"]
    span=1000
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    chs=["All"] #暂时只处理All通道
    for ch in chs:
        timetag=[]
        dtime=[]
        chl=[]
        hasData=False
        c.execute("select TimeTag from fretData_"+ch+" ORDER BY TimeTag limit 1")
        idx= c.fetchone()[0]
        if idx!=None:
            hasData=True
        while hasData:
            sql="select TimeTag from fretData_"+ch+" where TimeTag>=? and TimeTag<? ORDER BY TimeTag"
            c.execute(sql,(idx,idx+binTag*span))
            data=c.fetchall()
            lendata=len(data)
            if lendata<=binTag:
                hasData=False
            mvidx=idx

            for i in range(span):     
                r1=findLastLessThanIdx(data,idx+(i+1)*binTag)         
                timetag.append(BurstSearch.data2Dcol(data,r0,r1,0))
            idx=idx+binTag*span                