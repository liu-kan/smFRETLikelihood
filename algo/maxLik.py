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
        self.cpu_count=4#multiprocessing.cpu_count()
        self.chunks=list(chunks(range(self.n_burst), int(self.n_burst/(self.cpu_count-1))))
        self.procNum=len(self.chunks)
        self.allLikhDone=multiprocessing.Semaphore(self.procNum-1)
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
                                self.numWorkingProc,self.sumCanStartEvent,self.lkhCanStartEvent\
                                ,self.allLikhDone)
            self.procRec.append(multiprocessing.Process(target=cp))
            #lockRec.append(lk)

        for pid in self.procRec:
            pid.daemon=True
            pid.start()
        starttime = datetime.datetime.now()
        results = minimize(self.lnLikelihood, params[:self.n_states*self.n_states], method='Nelder-Mead',options=dict(maxiter=1000,disp=True),tol=1e-10)
        print (results)
        self.numWorkingProc=-1
        endtime = datetime.datetime.now()
        print ("Total Max likehood time:",endtime - starttime," s")
        for pid in self.procRec:
            pid.join(1)
            pid.terminate()

    def lnLikelihood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        调用likelihood前更新self.n_states
        ln(L)=sum(ln(L_j))
        \[L = {1^T}\prod\limits_{k = 2}^{N_j} [F(c_k)exp(K\tau _k)]F(\c_1)p_{eq} \]
        """

        self.lkhCanStartEvent.clear()
        sumLnL=0
        self.sharedArrP[:self.n_states*self.n_states]=params[:self.n_states*self.n_states]

        self.lkhCanStartEvent.set()

        self.minIter+=1
        #print(self.minIter)
        if self.minIter%10==0:

            E=genMatE(self.n_states,params[:self.n_states])
            K=genMatK(self.n_states,params[self.n_states:self.n_states*self.n_states])
            p=genMatP(K)
            print("==================================")
            print(p)
            print(K)
            print(E)
            if self.minIter%10==0:
                oldtime=self.timemes
                self.timemes=datetime.datetime.now()
                timesp=float((self.timemes-oldtime).seconds)
                if timesp>0:
                    print("The speed of analysis is %f burst/s" % ((10*self.n_burst)/timesp))

        self.sumCanStartEvent.wait()
        qs=self.numB.value
        for idx_burst in range(qs):
            #resBLH=0
            try:
                resBLH=self.queueOut.get(True,2)
            except:
                print("queueOut except")
                continue
            #print(sumLnL[0][0]+5000)
            sumLnL=sumLnL+resBLH[0][0]
        self.numB.value=0
        #print(sumLnL)
        return -sumLnL

def mdotl(*args):
    if len(args)<2:
        return None
    r=args[0]
    for i in range(len(args)-1):
        r=np.dot(r,args[i+1])
    return r

def genMatE(n,args):
    if len(args)<1:
        return None
    if len(args)!=n:
        return None
    matE=np.zeros([n,n])
    for i in range(n):
        matE[i,i]=args[i]
    return matE
def genMatK(n,args):
    if len(args)<1:
        return None
    if len(args)!=n*n-n:
        return None
    matK=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if i<j:
                matK[i,j]=args[i*(n-1)+j-1]
            elif i>j:
                matK[i,j]=args[i*(n-1)+j]
    for i in range(n):
        for j in range(n):
            if i==j:
                matK[i,j]=-np.sum(matK[:,j])
    return matK
def genMatP(matK):
    n=matK.shape[0]
    if n<1:
        return None
    matP=np.empty([n,1])
    ap=0
    for i in range(n):
        ap+=matK[i,i]
    for i in range(n):
        matP[i,0]=matK[i,i]/ap
    return matP

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
