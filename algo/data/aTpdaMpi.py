from math import *
import math
from scipy.optimize import differential_evolution
from scipy.optimize import basinhopping
from scipy.optimize import minimize
import pickle,pathlib
import copy,sys
import time,os
import numpy as np
# import numba as nb
import logging
import logging.handlers
from scipy import stats
# import mpi4py
# # mpi4py.rc.initialize = False
# mpi4py.rc.finalize = False
from mpi4py import MPI
import datetime
class mburst():
    def __init__(self, istart, istop,start, stop):
            self.istart=istart
            self.istop=istop
            self.start=start
            self.stop=stop
        
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
def mpiloaddata(comm,logger,filename="pdampi.dat"):
    if comm==None:
        p=pathlib.Path(filename)
        fn=str(p)        
        sbuf=os.path.getsize(fn)                
        rbuf=bytearray(sbuf)
        sdata=None
        with open(fn, 'rb') as f:
            try:
                rbuf=f.read()
                sdata=pickle.loads(rbuf)
            except : # whatever reader errors you care about
                print ("Could not write file:", fn)                
        return sdata        
    else:
        p=pathlib.Path(filename)
        fn=str(p)      
        sbuf=os.path.getsize(fn)     
        rbuf=bytearray(sbuf)
        fhn = MPI.File.Open(comm, fn, amode= MPI.MODE_RDONLY )
        status = MPI.Status()
        fhn.Read(rbuf,status)    
        sdata=pickle.loads(rbuf)
        fhn.Close()
        return sdata

def savedata(sdata,filename="pdampi.dat"):
    # p=pathlib.Path(pathlib.Path.home(),"tmp",filename)    
    pathlib.Path('data').mkdir(parents=True, exist_ok=True)
    p=pathlib.Path("data",filename)
    fn=str(p)    
    buf=pickle.dumps(sdata,protocol=4)
    sbuf=len(buf)
    fh = open(fn, 'wb')
    fh.write(buf)
    fh.close()
    # loginfo(comm,logger,"save dict data")
    # print ("save data at ",p.resolve())
    return sbuf


def loginfo(comm,log,msg,rank=0):
    if comm!=None:
        r=comm.Get_rank()
        name=MPI.Get_processor_name()
        if comm.Get_rank()==rank or rank<0:
            log.info("{} rank{},{}".format(name,r,msg))
    else:
        print(msg)
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
    # print(matK)
    n=matK.shape[0]
    if n<1:
        return None
    matP=np.empty([n,1])
    ap=0
    for i in range(n):
        ap+=matK[i,i]
    for i in range(n):
        matP[i,0]=matK[i,i]/ap
    # print(matP)        
    return matP
# @nb.jit
def drawDisIdx(idx,p,size=1):
    # hight=max(idx)
    # low=min(idx)
    # while True:
    #     X = np.random.random_integers(low,hight)
    #     Y = np.random.uniform(0.0, 1)        
    #     if Y < p[X]:
    #         return X    
    disc = stats.rv_discrete(name='disc', values=(idx, p))
    return disc.rvs(size=size)

# @nb.jit(nb.double(nb.double,nb.int32))
def _drawTau(k,size=1):
    # hight=k*1.5
    # low=0
    # while True:
    #     X = np.random.uniform(low,hight)
    #     Y = np.random.uniform(0.0, k)              
    #     if Y < k*exp(-k*X): 
    #         return X  
    sc=1/k
    b=stats.expon.rvs(scale=sc,size=8)
    ib=np.where( (b>0))[0]
    if len(ib)>0:
        # print (b[ib[0]])
        return b[ib[0]]    
    return stats.expon.rvs(scale=sc,size=size)
#@nb.jit(nb.double(nb.int64,nb.double,nb.double),cache=True)
def _draw_P_B_Tr(pN,T,bg_rate):    
    # while True:
    #     X = np.random.random_integers(0,pN)
    #     Y = np.random.uniform(0.0, 1)
    #     bg_bar=bg_rate*T
    #     fx=bg_bar**X*exp(-bg_bar)/math.factorial(X)
    #     if Y < fx:
    #         return X
    ba=np.empty(T.shape[0])
    for i in range(T.shape[0]):
        b=stats.poisson.rvs(T[i]*bg_rate,size=8)
        ie=np.where((b<pN[i]))[0]
        if len(ie)>0:
            ba[i]=b[ie[0]]
        else:
            ba[i]=np.random.random_integers(0,pN[i])
    return ba
# @nb.jit#(nb.double(nb.double,nb.double),cache=True)
def _drawE(e,v):    
    E=stats.norm.rvs(loc=e,scale=v,size=8)
    ie=np.where((E>0.01) & (E<0.99))[0]
    if len(ie)>0:
        return E[ie[0]]
    return e

# @nb.jit(nb.double(nb.int32,nb.double,nb.int32))
def _drawA_fi_e(fi,e,size=1):
    # while True:
    #     X = np.random.random_integers(0,fi)
    #     Y = np.random.uniform(0.0, 1)        
    #     if Y < math.factorial(fi)/(math.factorial(fi-X)*math.factorial(X))* \
    #         e**X*(1-e)**(fi-X):
    #         return X              
    return stats.binom.rvs(fi, e, size=size)
# @nb.jit
def _drawJ_Si2Sj(matP,n_states,i):
    P_i2j=copy.deepcopy(matP)
    P_i2j[i]=0
    P_i2j=P_i2j/sum(P_i2j)
    j=drawDisIdx(np.arange(n_states),P_i2j)
    return j

class pdaPy:
    def set_nstates(self,n,e,k,v):
        self.n_states=n
        if type(k) is not list:
            k=k.tolist()
            # print(k)
        if type(e) is not list:
            e=e.tolist()            
            # print(e)
        self.k=k
        self.E=e        
        self.Evar=v
        p=copy.deepcopy(e)
        p.extend(k)
        p.extend(v)
        self.params=p
        # print("=============================")
        # loginfo(self.comm,self.logger,"self.params {}".format(self.params))
        # loginfo(self.comm,self.logger,"self.params {}".format(self.params))
        rk=10**np.asarray(k)*1000
        self.matK=genMatK(n,rk)
        self.matP=genMatP(self.matK)
        # print("matP",self.matP)
        # print("matK",self.matK)
    # def __init__(self,bursts_list,timesMs,maskAD,pN,T,SgDivSr,bg_ad_rate,bg_dd_rate,clk_p,histBinNum=50,reSampleTimes=10,burstWindows=(0)):
    def __init__(self,comm,bursts_list,timesMs,maskAD,maskDD,T,SgDivSr,bg_ad_rate,bg_dd_rate,\
            clk_p,logger,histBinNum=50,\
            reSampleTimes=10,burstWindows=(0),maxiter=5000):
        self.comm=comm
        if comm!=None:
            self.rank=comm.Get_rank()
        self.histBinNum=histBinNum
        self.logger=logger
        self.maxiter=maxiter
        self.clk_p=clk_p
        self.reSampleTimes=reSampleTimes
        self.mask_ad=maskAD
        self.mask_dd=maskDD
        self.bursts=bursts_list
        self.burstWindowsStart=0
        self.burstWindowsEnd=T.shape[0]
        self.burstWindowsRa=burstWindows
        self.numWindows=self.burstWindowsEnd-self.burstWindowsStart
        # assert self.numWindows==bg_ad_rate.shape[0]
        # assert self.numWindows==bg_dd_rate.shape[0]
        # self.pN=pN
        self.T=T
        self.SgDivSr=SgDivSr
        self.timesMs=timesMs
        self.bg_ad_rate=bg_ad_rate
        self.bg_dd_rate=bg_dd_rate        
        # self.calcP_N()
        self.stop=[0]
        self.gamma=0.34     
        self.beta=1.42
        self.DexDirAem=0.08
        self.Dch2Ach=0.07         
    def drawJ_Si2Sj(self,i):
        P_i2j=copy.deepcopy(self.matP)
        P_i2j[i]=0
        P_i2j=P_i2j/sum(P_i2j)
        j=drawDisIdx(np.arange(self.n_states),P_i2j)
        return j

    def drawTau(self,si,sj):
        return _drawTau(self.matK[si][sj])

    # def draw_P_B_Tr(self,idx,bg_rate):
    #     return _draw_P_B_Tr(self.pN(idx),self.T[idx],bg_rate) 

    def drawA_fi_e(self,fi,e):
        return _drawA_fi_e(fi,e)

    def monteCarlo(self,idx):
        # np.set_printoptions(precision=13,suppress=True)
        bSgDivSr=np.zeros(self.reSampleTimes)
        stateIdx=np.arange(self.n_states)
        siAr=drawDisIdx(stateIdx,self.matP,size=self.reSampleTimes)
        # drawDisIdx 获取开始态 
        for sampleTime in range(self.reSampleTimes):
            si=siAr[sampleTime]
            mcSpendTime=0
            bins=[self.bursts[idx].start]
            sidx=[si]
            # print("self.T[idx]",self.T[idx])
            while self.T[idx]>mcSpendTime:        
                sj=_drawJ_Si2Sj(self.matP,self.n_states,si)[0]        
                # sj=self.drawJ_Si2Sj(si)
                sidx.append(sj)   
                # print("sidx",sidx)             
                Tau=_drawTau(self.matK[si][sj])
                # Tau=self.drawTau(si,sj)            
                si=sj
                mcSpendTime+=Tau
                # print("Tau",Tau)
                if mcSpendTime>=self.T[idx]:
                    bins.append(self.bursts[idx].stop)
                    
                else:
                    bins.append(bins[0]+mcSpendTime/self.clk_p)
            bins=np.asarray(bins,dtype=np.int64)#/self.clk_p
            # if len(bins)>2:
            #     print(bins)        
            m_ad=self.mask_ad[self.bursts[idx].istart:self.bursts[idx].istop+1]
            b_times=self.timesMs[self.bursts[idx].istart:self.bursts[idx].istop+1]
            m_dd=self.mask_dd[self.bursts[idx].istart:self.bursts[idx].istop+1]
            # mask_F=np.logical_or(m_dd,m_ad)
            f_ia, _ = np.histogram(b_times[m_ad], bins)
            f_id, _ = np.histogram(b_times[m_dd], bins)
            f_if=(self.gamma-self.Dch2Ach)*f_id + (1-self.DexDirAem)*f_ia
            f_i=np.empty_like(f_if)
            if self.bg_dd_rate==0:
                f_i=np.around(f_if).astype(int)
            else:
                # 计算背景噪声
                t_a=np.diff(bins)*self.clk_p
                bg_a=_draw_P_B_Tr((1-self.DexDirAem)*f_ia,t_a,self.bg_ad_rate)
                bg_d=_draw_P_B_Tr((self.gamma-self.Dch2Ach)*f_id,t_a,self.bg_dd_rate)
                f_i=np.around(f_if - bg_d - bg_a).astype(int)
            ############
            stransNum=len(f_i)
            F=sum(f_if)
            for i in range(stransNum):
                de=_drawE(self.E[sidx[i]],self.Evar[sidx[i]])                
                ai=_drawA_fi_e(f_i[i],de)
                bSgDivSr[sampleTime]=bSgDivSr[sampleTime]+ai
            bSgDivSr[sampleTime]=bSgDivSr[sampleTime]/F
        return bSgDivSr.tolist()

    def aTpda(self,EK=False):
        pdaBurstE=[]
        self.pdaBurstE=[]        
        if EK:
            from multiprocessing import Pool
            ra=range(self.burstWindowsStart,self.burstWindowsEnd)
            with Pool(processes=int(os.cpu_count()/2)) as pool:  
                pdaBurstE = pool.map(self.monteCarlo, ra)
            for be in pdaBurstE:
                self.pdaBurstE.extend(be)                
        else:
            for r in self.burstWindowsRa:
                for i in r:
                    pdaBurstE.extend(self.monteCarlo(i))        
            pdaBurstEga=self.comm.gather(pdaBurstE,root=0)            
            if self.rank==0:
                for be in pdaBurstEga:
                    self.pdaBurstE.extend(be)

    def aTpdaEK(self):
        self.set_nstates(self.n_states,self.fitP[:self.n_states],self.fitP[self.n_states:\
            self.n_states*self.n_states],self.fitP[self.n_states*self.n_states:])
        loginfo(self.comm,self.logger,"drawPic matP:{}".format(self.matP))
        loginfo(self.comm,self.logger,"drawPic matK:{}".format(self.matK))
        loginfo(self.comm,self.logger,"drawPic E:{}".format(self.E))
        self.aTpda(True)
        rSgDivSr=self.SgDivSr[self.burstWindowsStart:self.burstWindowsEnd]
        SgDivSrR=np.histogram(rSgDivSr,bins=self.histBinNum) 
        tSgDivSr=np.zeros_like(SgDivSrR[0])
        for SgDivSrV in self.pdaBurstE:            
            idxHistBin=np.searchsorted(SgDivSrR[1], SgDivSrV)
            if idxHistBin<self.histBinNum:  
                tSgDivSr[idxHistBin]+= 1
        tSgDivSr=tSgDivSr/float(self.reSampleTimes)
        chisqr=0
        n=copy.deepcopy(self.histBinNum)
        # sO=sum(SgDivSrR[0])
        # sE=sum(tSgDivSr)
        for O,E in zip(SgDivSrR[0],tSgDivSr):
            if O>0:
                # chisqr+=(O/sO-E/sE)**2/(O/sO)
                chisqr+=(float(O-E))**2/float(O)
            else:
                n=n-1
        chisqr=chisqr/(n-self.n_states*(self.n_states+1))            
        loginfo(self.comm,self.logger,"sum(SgDivSrR) {}".format(sum(SgDivSrR[0])))
        loginfo(self.comm,self.logger,"sum(tSgDivSr) {}".format(sum(tSgDivSr)))
        loginfo(self.comm,self.logger,"chisqr: {}".format(chisqr))
        return SgDivSrR,tSgDivSr

    def chiSqr(self):
        rSgDivSr=self.SgDivSr[self.burstWindowsStart:self.burstWindowsEnd]
        SgDivSrR=np.histogram(rSgDivSr,bins=self.histBinNum)  #,range=(0,5)
        tSgDivSr=np.zeros_like(SgDivSrR[0])
        for SgDivSrV in self.pdaBurstE:            
            idxHistBin=np.searchsorted(SgDivSrR[1], SgDivSrV)
            if idxHistBin<self.histBinNum:  
                tSgDivSr[idxHistBin]+= 1
        tSgDivSr=tSgDivSr/float(self.reSampleTimes)
        chisqr=0
        n=copy.deepcopy(self.histBinNum)
        # sO=sum(SgDivSrR[0])
        # sE=sum(tSgDivSr)
        for O,E in zip(SgDivSrR[0],tSgDivSr):
            if O>0:
                # chisqr+=(O/sO-E/sE)**2/(O/sO)
                chisqr+=(float(O-E))**2/float(O)
            else:
                n=n-1
        chisqr=chisqr/(n-self.n_states*(self.n_states+1))                
        loginfo(self.comm,self.logger,"chiSqr:{}".format(chisqr))
        loginfo(self.comm,self.logger,"=======================================")        
        return chisqr 

    def chiSqrArrT(self,params,stopp):
        # loginfo(self.comm,self.logger,"before bcast stopp{}".format(self.params) )
        stopp[0]=self.comm.bcast(stopp[0], root=0)
        # loginfo(self.comm,self.logger,"after bcast stopp{}".format(self.params) )
        if stopp[0]!=0:
            return 0
        if self.rank==0:
            start = MPI.Wtime()
            loginfo(self.comm,self.logger,"start" )
        n=self.n_states
        
        
        paramsN=np.asarray(self.params,dtype=np.float64)
        # loginfo(self.comm,self.logger,"paramsOld{}".format(paramsN) )    
        # loginfo(self.comm,self.logger,"paramsNew{}".format(params) )    
        lenp=self.n_states*(self.n_states+1)    
        # self.comm.Barrier()
        if self.rank==0:        
            paramsN=np.asarray(params,dtype=np.float64)            
        # for idxp in range(lenp):
        self.comm.Bcast([paramsN,lenp,MPI.DOUBLE], root=0)
        self.params=paramsN.tolist()
        loginfo(self.comm,self.logger,"paramsBca{}".format(self.params) )
        e=self.params[:self.n_states]
        k=self.params[self.n_states:self.n_states*self.n_states]
        evar=self.params[self.n_states*self.n_states:]
        self.set_nstates(n,e,k,evar)
        self.aTpda()           
        if self.rank==0:
            chi=self.chiSqr()            
            end = MPI.Wtime()
            loginfo(self.comm,self.logger,"speed is :{} burst/s at {}".format(self.numWindows/(end-start),datetime.datetime.now()) )
            # self.comm.Barrier()
            return chi
        # else:
            # self.comm.Barrier()
    # def calcP_N(self):
    #     hist, bin_edges = np.histogram(self.pN,bins=self.pNbins)
    #     pNhist=hist/self.numWindows
    #     center = (bin_edges[:-1] + bin_edges[1:]) / 2    
    #     print("pN range",min(self.pN),max(self.pN))
    #     photonN=np.arange(max(0,min(self.pN)),max(self.pN)+1,1)
    #     pNcount = prepData.interpSpl(center,pNhist, photonN)#,'linear')    
    #     pNcount[np.where(pNcount<0)]=0
    #     self.PN = dict(zip(photonN, pNcount))
    def P_N(self,N):
        return self.PN[N]

    def P_F_RT_F(self,F_RT,F_G,e):
        F=round(F_RT+F_G)
        F_RT=round(F_RT)
        F_G=round(F_G)
        return factorial(F)/(factorial(F_RT)*factorial(F_G))* \
        e**F_RT*(1-e)**(F_G)

    def pdaP_B(self,B_bar,B):
        fm=(B_bar**B)*exp(-B_bar)
        # B_bar=round(B_bar)
        B=round(B)
        return fm/factorial(B)

    def drawPdaP_Br(self,B_bar,B):
        fm=(B_bar**B)*exp(-B_bar)
        # B_bar=round(B_bar)
        B=round(B)
        return fm/factorial(B)

    def pSgDivSr(self,e,p_N,F_RT,F_G,bg_dd,bg_ad,b_g_bar,b_r_bar):
        # print(e,p_N,F_RT,F_G,bg_dd,bg_ad,b_g_bar,b_r_bar)
        return p_N*self.P_F_RT_F(F_RT,F_G,e)*self.pdaP_B(b_g_bar,bg_dd)*\
        self.pdaP_B(b_r_bar,bg_ad)
    def setStateNum(self,state):
        self.n_states=state
        self.set_nstates(state, np.random.rand(state).tolist(),[0]*(state-1)*state,[0.1]*state)
    def findP(self):
        # bnds = (0,1)
        xmin = [0.03]*self.n_states
        xmink=[np.log10(0.1)]*(self.n_states*self.n_states-self.n_states)
        xminv = [0.001]*self.n_states
        xmin.extend(xmink)
        xmin.extend(xminv)
        xmax = [0.97]*self.n_states
        xmaxk=[np.log10(10)]*(self.n_states*self.n_states-self.n_states)
        xmaxv = [0.6]*self.n_states
        xmax.extend(xmaxk)
        xmax.extend(xmaxv)
        bounds = [(low, high) for low, high in zip(xmin, xmax)]

        minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds,args=(self.stop,))
        # self.params=[0.31841749529573427, 0.66238514811983, 0.7750358211881881, 1.481928598181776, 0.021950315293281453, 0.34853313837435074]
        # res=minimize(self.chiSqrArrT,self.params, args=(self.stop,),method='L-BFGS-B', bounds=bounds)
        # res = basinhopping(self.chiSqrArrT,self.params, minimizer_kwargs=minimizer_kwargs, \
        #     disp=True)#,interval=5,stepsize=0.1,niter=10
        res = differential_evolution(self.chiSqrArrT, args=(self.stop,),strategy ='currenttobest1bin',bounds=bounds,\
                           disp=True,polish=True,maxiter=self.maxiter)
        loginfo(self.comm,self.logger,"{}".format(res))
        loginfo(self.comm,self.logger,"{}".format(res.x))
        # print (self.numWindows)
        self.fitP=res.x
        loginfo(self.comm,self.logger,'''
        ===========================================
        {}  
        ==========================================='''.format(datetime.datetime.now()))        
        self.stop=[1]
        self.chiSqrArrT(self.params,self.stop)