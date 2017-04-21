from scipy.linalg import expm
import numpy as np
import multiprocessing
class calcBurstLikehood():#multiprocessing.Process):
    def __init__(self,chunkRange,queueOut,burst,n_states,\
            Sth,procNum,paramsArr,numB,numWorkingProc,sumCanStartEvent\
            ,lkhCanStartEvent,allLikhDone):
        #multiprocessing.Process.__init__(self)
        self.burst=burst
        self.allLikhDone=allLikhDone #Semaphore (procNum)
        self.sumCanStartEvent=sumCanStartEvent
        self.lkhCanStartEvent=lkhCanStartEvent
        self.n_states=n_states
        self.n_burst=len(burst["All"].chl)
        self.Sth=Sth
        self.params=paramsArr
        E=genMatE(self.n_states,self.params[:self.n_states])
        K=genMatK(self.n_states,self.params[self.n_states:self.n_states*self.n_states])
        p=genMatP(K)
        self.E=E
        self.K=K
        self.p=p
        self.procNum=procNum
        self.queueOut=queueOut
        self.running=True
        self.numB=numB
        self.chunkRange=chunkRange
        self.numWorkingProc=numWorkingProc
    def matF(self,c_k):
        if c_k==1:
            return self.E
        elif c_k==2:
            return np.eye(self.n_states)-self.E
        else:
            return None
    def __call__(self):
        #idx_burst=0
        while self.numWorkingProc.value!=-1:
            startLKH=self.lkhCanStartEvent.wait(1)
            #print("proc:",self.params[:])
            if not startLKH:
                continue
            #self.numWorkingProc.value+=1
            #print("proc:",self.params[:])
            E=genMatE(self.n_states,self.params[:self.n_states])
            K=genMatK(self.n_states,self.params[self.n_states:self.n_states*self.n_states])
            p=genMatP(K)
            self.E=E
            self.K=K
            self.p=p
            for idx_burst in self.chunkRange:
                #print("idx_burst",idx_burst)
                if self.burst["All"].s[idx_burst]>=self.Sth:
                    #self.queueOut.put(0)
                    continue
                lenPhoton=self.burst["All"].ntag[idx_burst]
                if lenPhoton<2:
                    #self.queueOut.put(0)
                    continue
                lnL_j=np.linspace(1,1,self.n_states)
                t_k_0=-1
                prod=np.eye(self.n_states)
                for idx_photon in range(lenPhoton):
                    F=self.matF(self.burst["All"].chl[idx_burst].iloc[idx_photon])
                    if F is not None:
                        if t_k_0<0:
                            t_k_0=self.burst["All"].timetag[idx_burst].iloc[idx_photon]*self.burst["SyncResolution"] \
                                +self.burst["All"].dtime[idx_burst].iloc[idx_photon]*self.burst["DelayResolution"]
                            lnL_j=np.dot(F,p)
                            continue
                        t_k_1=self.burst["All"].timetag[idx_burst].iloc[idx_photon]*self.burst["SyncResolution"] \
                            +self.burst["All"].dtime[idx_burst].iloc[idx_photon]*self.burst["DelayResolution"]
                        tau=t_k_1-t_k_0
                        t_k_0=t_k_1
                        FdotExp=np.dot(F,expm(K)*tau)
                        prod=np.dot(FdotExp,prod)
                if t_k_0<0:
                    #self.queueOut.put(0)
                    continue
                lnL_j=np.dot(FdotExp,lnL_j)
                T1=np.linspace(1,1,self.n_states)
                T1.shape=(self.n_states,1)
                T1=np.transpose(T1)
                L_burst=np.dot(T1,lnL_j)
                lnL_burst=0
                if L_burst<1e-300:
                    print("L_burst is too small:",L_burst)
                    lnL_burst=np.log(L_burst*1e300)-np.log(1e300)
                else:
                    lnL_burst=np.log(L_burst)
                self.queueOut.put(lnL_burst)
                self.numB.value+=1

            #还需保证计算sum前 lkhCanStartEvent被clear,但是如果有进程还没有开始计算（小概率）lkhCanStartEvent不能被clear
            #所以要先等待所有进程开始计算
            if not self.allLikhDone.acquire(False):
                self.lkhCanStartEvent.clear()
                self.sumCanStartEvent.set()
                #print("numWorkingProc allLikhDone.acquire",self.numWorkingProc.value)
                #self.numWorkingProc.value-=1
                #self.sumCanStartEvent.wait()
            else:
                self.sumCanStartEvent.wait()
                #self.numWorkingProc.value-=1
                #print("numWorkingProc ",self.numWorkingProc.value)
                self.allLikhDone.release()

        #print("calc end")
        #self.lock.release()
        #self.terminate()


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
