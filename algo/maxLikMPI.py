import numpy as np
from mpi4py import MPI
from scipy.linalg import expm
try:
    import algo.BurstSearch as BurstSearch
    import algo.BGrate as BGrate
    import algo.fretAndS as fretAndS
except ImportError:
    import BurstSearch
    import BGrate
    import fretAndS
import datetime
from scipy.optimize import minimize
from array import array
from scipy.linalg import expm
from mpiBurstLikelihood import calcBurstLikelihood
class GS_MLE():
    def __init__(self, burst,comm,Sth=0.88):
        self.timemes=datetime.datetime.now()
        self.burst=burst
        self.n_states=1 #调用likelihood前更新self.n_states
        self.n_burst=len(burst["All"].chl)
        self.Sth=Sth
        self.minIter=0
        self.comm=comm
    def MaxLikehood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        """
        self.params=params

        results = minimize(self.lnLikelihood, params, method='Nelder-Mead')


    def lnLikelihood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        调用likelihood前更新self.n_states
        ln(L)=sum(ln(L_j))
        \[L = {1^T}\prod\limits_{k = 2}^{N_j} [F(c_k)exp(K\tau _k)]F(\c_1)p_{eq} \]
        """
        sumLnL=0
        stop=params[self.n_states*self.n_states]
        stop=self.comm.bcast(stop,root=0)
        self.E=genMatE(self.n_states,params[:self.n_states])
        K=genMatK(self.n_states,params[self.n_states:self.n_states*self.n_states])
        p=genMatP(K)

        if self.minIter%10==0:
            print("==================================")
            print(p)
            print(K)
            print(self.E)
            if self.minIter%10==0:
                oldtime=self.timemes
                self.timemes=datetime.datetime.now()
                print("The speed of analysis is %f burst/s" % ((10*self.n_burst)/float((self.timemes-oldtime).seconds)))

        for idx_burst in range(self.n_burst):
            if self.burst["All"].s[idx_burst]>=self.Sth:
                continue
            lenPhoton=self.burst["All"].ntag[idx_burst]
            if lenPhoton<2:
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
            sumLnL+=lnL_burst
        return -sumLnL
    def matF(self,c_k):
        if c_k==1:
            return self.E
        elif c_k==2:
            return np.eye(self.n_states)-self.E
        else:
            return None
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

comm=MPI.COMM_WORLD
rank=comm.Get_rank()

def chunks(l, n):
    """Yield successive nth chunks from l."""
    lenl=len(l)
    stack=[]
    if(lenl%n==0):
        for i in range(0, lenl, lenl/n):
            stack.append(l[i:i + lenl/n])
        for i in range(0, lenl, lenl/n):
            yield stack.pop()
    else:
        for i in range(0, lenl, int(lenl/n)+1):
            stack.append(l[i:i + int(lenl/n)+1])
        for i in range(0, lenl, int(lenl/n)+1):
            yield stack.pop()

if __name__ == '__main__':

    clsize=comm.Get_size()
    if rank==0:
        #starttime = datetime.datetime.now()
        dbname='/home/liuk/sf/oc/data/38.sqlite'
        dbname='E:/liuk/proj/ptu/data/55.sqlite'
        #dbname='E:/sf/oc/data/38.sqlite'
        dbname='/home/liuk/sf/oc/data/1min.sqlite'
        #dbname='/prog/data/1min.sqlite'
        br=BGrate.calcBGrate(dbname,20,400)
        burst=BurstSearch.findBurst(br,dbname,["All"])
        burstSeff, burstFRET,wei,H,xedges, yedges=fretAndS.FretAndS(dbname,burst,(27,27),br)
        n_burst=len(burst["All"].chl)
        if n_burst<clsize:
            clsize=n_burst
        chunkLists=list(chunks(range(n_burst), int(n_burst/clsize)))
        #gsml=GS_MLE(burst,0.891)
        #gsml.n_states=2
        #gsml.MaxLikehood([0.3,0.7,0.2, 3,3,3, 3,3,3])

        #endtime = datetime.datetime.now()
        #print (endtime - starttime)
        n_states=2
    else:
        burst=dict()
        n_states=-1
    clsize=comm.bcast(clsize,root=0)
    burst=comm.bcast(burst,root=0)
    n_states=comm.bcast(n_states,root=0)
    gsml=GS_MLE(burst,comm,0.891)
    gsml.n_states=n_states
    params=[0.3,0.7,0.2, 3,3,3, 3,3,3]
    params=params[:n_states]
    params.append(0)
    #stop=0
    if rank==0:
        gsml.MaxLikehood(params)
        stop=1
        gsml.lnLikelihood(params)
    else:
        gsml.lnLikelihood(params)