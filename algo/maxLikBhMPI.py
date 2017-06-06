 # -*- coding: utf-8 -*-
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
import datetime,sys
from scipy.optimize import basinhopping
from array import array
from scipy.linalg import expm
#from mpiBurstLikelihood import calcBurstLikelihood

def appendResult(fn,results,n,timesp,Sth,dbname):
    fo = open(fn, "a")
    fo.write("n_states:"+str(n)+" Sth:"+str(Sth)+' : ======== '+str(datetime.datetime.now())+'\n')
    fo.write(dbname+'\n')
    fo.write( str (results)+'\n')
    fo.write("E:==========="+'\n')
    fo.write(str(genMatE(n,results.x))+'\n')
    fo.write("K:==========="+'\n')
    k=genMatK(n,results.x+'\n')
    fo.write(str(k)+'\n')
    fo.write("P==========="+'\n')
    fo.write(str(genMatP(k))+'\n')
    fo.write('time spend : '+str(timesp)+'\n')
    fo.close()

class bhSteps(object):
    def __init__(self, bounds, stepsize=1):
        self.bounds = bounds
        if stepsize<0:
            stepsize=stepsize*-1
        if stepsize>1:
            stepsize=1/stepsize
        self.stepsize=stepsize
    def __call__(self, x):
        idx=0
        for bound in self.bounds:
            xmax=(bound[1]-x[idx])*self.stepsize
            xmin=(x[idx]-bound[0])*self.stepsize
            x[idx]=x[idx]+np.random.uniform(xmin,xmax)
            idx=idx+1

class GS_MLE():
    def __init__(self, burst,comm,burstIdxRange,Sth=0.88,dbname=''):
        self.timemes=MPI.Wtime()
        self.burst=burst
        self.burstIdxRange=burstIdxRange
        self.n_states=1 #调用likelihood前更新self.n_states
        self.n_burst=len(burst["All"]['chl'])
        self.Sth=Sth
        self.minIter=0
        self.comm=comm
        self.stop=[0]
        self.params=[]
        self.dbname=dbname
        self.oldIter=0
    def MaxLikehood(self,params):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        """

        self.params=params
        startTime=datetime.datetime.now()
        boundE=[(-0.15,0.999)]*self.n_states
        boundK=[(0.1,100000)]*(self.n_states*(self.n_states-1))
        bound=boundE+boundK
        bhStep=bhSteps(bound,0.5)
        minimizer_kwargs = {"method":"L-BFGS-B","args":(self.stop,)}
        results = basinhopping(self.lnLikelihood, params, \
                           minimizer_kwargs=minimizer_kwargs,\
                           take_step= bhStep)
        stopTime=datetime.datetime.now()
        print(results)
        appendResult('results.txt',results,self.n_states,stopTime-startTime,self.Sth,self.dbname)
        print("Time spend:",stopTime-startTime)
        self.stop=[1]
        self.lnLikelihood(params,self.stop)

    def lnLikelihood(self,params,stopp):
        """calc ln likehood of TCSCP stream.

        params[:n_states] E_{states_i}
        params[n_states:n_states^2] matrix of K
        调用likelihood前更新self.n_states
        ln(L)=sum(ln(L_j))
        \[L = {1^T}\prod\limits_{k = 2}^{N_j} [F(c_k)exp(K\tau _k)]F(\c_1)p_{eq} \]
        """
        sumLnL=0

        stopp[0]=self.comm.bcast(stopp[0], root=0)
        if stopp[0]!=0:
            return 0
        T1=np.linspace(1,1,self.n_states)
        T1.shape=(self.n_states,1)
        T1=np.transpose(T1)
        self.params=self.comm.bcast(self.params, root=0)
        self.E=genMatE(self.n_states,params[:self.n_states])
        #print(self.params)
        K=genMatK(self.n_states,params[self.n_states:self.n_states*self.n_states])
        #print(params[self.n_states:self.n_states*self.n_states])
        p=genMatP(K)
        self.minIter=self.minIter+1
        rank = self.comm.Get_rank()
        if rank ==0:
            if (self.minIter-self.oldIter)*self.n_burst>10000:
                print("==================================")
                print(p)
                print(K)
                print(self.E)
                oldtime=self.timemes
                self.timemes=MPI.Wtime()
                timesp=self.timemes-oldtime
                if timesp<1e-100:
                    timesp=1.0
                print("The speed of analysis is %f burst/s" % (((self.minIter-self.oldIter)*self.n_burst)/timesp))
                self.oldIter=self.minIter
                sys.stdout.flush()
        sumlnAlpha=np.zeros([1,1])
        for burstIdxRangeItem in self.burstIdxRange:
            for idx_burst in burstIdxRangeItem:
                if self.burst["All"]['s'][idx_burst]>=self.Sth:
                    continue
                lenPhoton=self.burst["All"]['ntag'][idx_burst]
                if lenPhoton<2:
                    continue
                lnL_j=np.ones([self.n_states,1])
                t_k_0=-1
                #prod=np.eye(self.n_states)

                for idx_photon in range(lenPhoton):
                    F=self.matF(self.burst["All"]['chl'][idx_burst].iloc[idx_photon])
                    if F is not None:
                        if t_k_0<0:
                            t_k_0=self.burst["All"]['timetag'][idx_burst].iloc[idx_photon]*self.burst["SyncResolution"] \
                                +self.burst["All"]['dtime'][idx_burst].iloc[idx_photon]*self.burst["DelayResolution"]
                            lnL_j=np.dot(F,p)
                            continue
                        t_k_1=self.burst["All"]['timetag'][idx_burst].iloc[idx_photon]*self.burst["SyncResolution"] \
                            +self.burst["All"]['dtime'][idx_burst].iloc[idx_photon]*self.burst["DelayResolution"]
                        tau=t_k_1-t_k_0
                        t_k_0=t_k_1
                        FdotExp=np.dot(F,expm(K*tau))
                        lnL_j=np.dot(FdotExp,lnL_j)
                        alpha=1.0/np.dot(T1,lnL_j)
                        lnL_j=lnL_j*alpha
                        sumlnAlpha+=np.log(alpha)
                if t_k_0<0:
                    continue
                #lnL_j=np.dot(prod,lnL_j)
                L_burst=np.dot(T1,lnL_j)
                lnL_burst=np.log(L_burst)-sumlnAlpha
                sumLnL+=lnL_burst[0,0]
        summ=self.comm.reduce(sumLnL,op=MPI.SUM, root=0)
        if rank ==0:
            return -summ
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


def chunks(l, n):
    """Yield successive nth chunks from l."""
    lenl=len(l)
    stack=[]
    stackList=[]
    for i in range(0, lenl, lenl//n):
        stack.append(l[i:i + lenl//n])
    for i in range(n):
        stackList.append([stack[i]])
    if len(stack)>n:
        for i in range(n,len(stack)):
            stackList[i-n].append(stack[i])
    for i in range(0, n):
        yield stackList.pop()


def main(comm,dbname,n_states,Sth):
    rank=comm.Get_rank()
    clsize=comm.Get_size()
    #print("============size=====",clsize)
    if rank==0:
        br=BGrate.calcBGrate(dbname,20,400)
        burst=BurstSearch.findBurst(br,dbname,["All"])
        burstSeff, burstFRET,wei,H,xedges, yedges=fretAndS.FretAndS(dbname,burst,(27,27),br)
        n_burst=len(burst["All"]['chl'])
        if n_burst<clsize:
            clsize=n_burst
        chunkLists=list(chunks(range(n_burst), clsize))

        #gsml=GS_MLE(burst,0.891)
        #gsml.n_states=2
        #gsml.MaxLikehood([0.3,0.7,0.2, 3,3,3, 3,3,3])

        #endtime = datetime.datetime.now()
        #print (endtime - starttime)

    else:
        burst=dict()

        chunkLists=list()
    clsize=comm.bcast(clsize,root=0)
    burst=comm.bcast(burst,root=0)

    burstIdxRange=comm.scatter( chunkLists, root=0)
    #print(burstIdxRange,rank)
    gsml=GS_MLE(burst,comm,burstIdxRange,Sth,dbname)
    gsml.n_states=n_states

    params=[0.38,0.6,0.5,675.0, 325.0,3,3, 3,3,3]
    params=params[:n_states*n_states]
    #print(params)
    stop=[0]
    if rank==0:
        gsml.MaxLikehood(params)
        #gsml.lnLikelihood(params,stop)
    else:
        while stop[0]==0:
            gsml.lnLikelihood(params,stop)

if __name__ == '__main__':
    #starttime = datetime.datetime.now()
    dbname='/home/liuk/sf/oc/data/38.sqlite'
    dbname='E:/liuk/proj/ptu/data/55.sqlite'
    #dbname='E:/sf/oc/data/38.sqlite'
    dbname='/smfret/1min.sqlite'
    dbname='/home/liuk/prog/smFRETLikelihood/data/1min.sqlite'
    n_states=2
    comm=MPI.COMM_WORLD
    rank=comm.Get_rank()
    Sth=0.90
    if rank==0:
        print("=========running mpi in ",comm.Get_size()," nodes~=============")
        sys.stdout.flush()
        if len(sys.argv)>1:
            dbname=sys.argv[1]
        if len(sys.argv)>2:
            n_states=int(sys.argv[2])
        if len(sys.argv)>3:
            Sth=float(sys.argv[3])
    dbname=comm.bcast(dbname,root=0)
    n_states=comm.bcast(n_states,root=0)
    Sth=comm.bcast(Sth,root=0)
    main(comm,dbname,n_states,Sth)
