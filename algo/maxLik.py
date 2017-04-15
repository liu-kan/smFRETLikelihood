import numpy as np
from scipy.linalg import expm
try:
    import algo.BurstSearch as BurstSearch
    import algo.BGrate as BGrate
    import algo.fretAndS as fretAndS
    from algo.mpBurstLikehood import calcBurstLikehood
except ImportError:
    import BurstSearch
    import BGrate
    import fretAndS
    from mpBurstLikehood import calcBurstLikehood
import datetime
from scipy.optimize import minimize
from array import array
import ctypes
import multiprocessing,time

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

class GS_MLE():
    def __init__(self, burst,Sth=0.88):
        self.timemes=datetime.datetime.now()
        self.burst=burst
        self.n_states=1 #调用likelihood前更新n
        self.n_burst=len(burst["All"].chl)
        self.Sth=Sth
        self.minIter=0
        self.cpu_count=multiprocessing.cpu_count()
        self.chunks=list(chunks(range(self.n_burst), int(self.n_burst/(self.cpu_count-1))))
        self.procNum=len(self.chunks)
    def MaxLikehood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        """
        self.params=params
        self.sharedArrP=multiprocessing.Array("d",params[:self.n_states*self.n_states])
        self.procRec=[]
        #lockRec=[]
        self.queueOut = multiprocessing.Queue()
        self.numB = multiprocessing.Value('l', 0)
        self.numWorkingProc = multiprocessing.Value(ctypes.c_int16, 0)
        #同一时间正在进行计算的进程，如果为负，结束所有进程
        self.sumCanStartEvent=multiprocessing.Event()
        self.lkhCanStartEvent=multiprocessing.Event()
        for idx_proc in range(self.procNum):
            #lk=multiprocessing.Lock()
            cp=calcBurstLikehood(self.chunks[idx_proc],self.queueOut,self.burst,self.n_states\
                                ,self.Sth,self.procNum,self.sharedArrP,self.numB,\
                                self.numWorkingProc,self.sumCanStartEvent,self.lkhCanStartEvent)
            self.procRec.append(multiprocessing.Process(target=cp))
            #lockRec.append(lk)
        for pid in self.procRec:
            pid.daemon=True
            pid.start()
        starttime = datetime.datetime.now()
        results = minimize(self.lnLikelihood, params, method='Nelder-Mead')
        print (results)
        for pid in self.procRec:
            pid.terminate()
        endtime = datetime.datetime.now()
        print ("Total Max likehood time:",endtime - starttime," s")

    def lnLikelihood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        调用likelihood前更新self.n_states
        ln(L)=sum(ln(L_j))
        \[L = {1^T}\prod\limits_{k = 2}^{N_j} [F(c_k)exp(K\tau _k)]F(\c_1)p_{eq} \]
        """
        sumLnL=0
        self.minIter+=1
        E=genMatE(self.n_states,params[:self.n_states])
        K=genMatK(self.n_states,params[self.n_states:self.n_states*self.n_states])
        p=genMatP(K)

        if self.minIter%10==0:
            print("==================================")
            print(p)
            print(K)
            print(E)
            if self.minIter%10==0:
                oldtime=self.timemes
                self.timemes=datetime.datetime.now()
                print("The speed of analysis is %f burst/s" % ((10*self.n_burst)/float((self.timemes-oldtime).seconds)))
        queueIn = multiprocessing.Queue()
        queueOut = multiprocessing.Queue()
        procRec=[]
        lockRec=[]
        numB = multiprocessing.Value('l', 0)
        for idx_burst in range(self.n_burst):
            queueIn.put(idx_burst)
        for idx_proc in range(self.cpu_count-1):
            queueIn.put(-1)
            lk=multiprocessing.Lock()
            cp=calcBurstLikehood(queueIn,queueOut,self.burst,self.n_states,self.Sth,self.E,K,p,lk,numB)
            procRec.append(multiprocessing.Process(target=cp))
            lockRec.append(lk)

        for pid in procRec:
            pid.daemon=True
            pid.start()

        # while queueIn.empty()==False:
        #     print("queueIn no empty,",queueIn.qsize())
        #     time.sleep(2)
        #
        for idx_proc in range(self.cpu_count-1):
            lockRec[idx_proc].acquire()
            #procRec[idx_proc].terminate()
            #lockRec[idx_proc].release()
        # for pidt in procRec:
        #     #pid.daemon=True
        #     pidt.terminate()

        #print("Joined")
        #print(queueOut.qsize(),numB.value)
        qs=numB.value
        for idx_burst in range(qs):
            #resBLH=0
            #try:
            resBLH=queueOut.get(True,3)
            #except:
            #    continue
            #print(sumLnL[0][0]+5000)
            sumLnL+=resBLH
        for idx_proc in range(self.cpu_count-1):
            #lockRec[idx_proc].acquire()
            procRec[idx_proc].terminate()
            lockRec[idx_proc].release()
        queueIn.close()
        queueOut.close()

        #print(-sumLnL[0][0])
        return -sumLnL[0][0]

def mdotl(*args):
    if len(args)<2:
        return None
    r=args[0]
    for i in range(len(args)-1):
        r=np.dot(r,args[i+1])
    return r

if __name__ == '__main__':
    import matplotlib,os
    if os.name!='posix':
        print("no fork,freeze_support need.")
        multiprocessing.freeze_support()
    starttime = datetime.datetime.now()


    dbname='E:/liuk/proj/ptu/data/55.sqlite'
    #dbname='E:/sf/oc/data/38.sqlite'
    dbname="E:/dbox/sf/oc/data/1min.sqlite"
    dbname='/home/liuk/sf/oc/data/1min.sqlite'
    br=BGrate.calcBGrate(dbname,20,400)
    burst=BurstSearch.findBurst(br,dbname,["All"])
    burstSeff, burstFRET,wei,H,xedges, yedges=fretAndS.FretAndS(dbname,burst,(27,27),br)
    #matplotlib.pyplot.plot(burst["All"].s)
    gsml=GS_MLE(burst,0.891)
    gsml.n_states=2
    gsml.MaxLikehood([0.3,0.7,0.2, 3,3,3, 3,3,3])

    endtime = datetime.datetime.now()
    print (endtime - starttime)
