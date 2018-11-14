from math import *
import math
from scipy.optimize import differential_evolution
from scipy.optimize import basinhopping
from data import prepData
import copy,sys
from progressbar import *
from multiprocessing import Pool
import time,os
import numpy as np
import numba as nb
import logging
import logging.handlers
from scipy import stats

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
@nb.jit(cache=True)
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

@nb.jit(nb.double(nb.double,nb.int32),cache=True)
def _drawTau(k,size=1):
    # hight=k*1.5
    # low=0
    # while True:
    #     X = np.random.uniform(low,hight)
    #     Y = np.random.uniform(0.0, k)              
    #     if Y < k*exp(-k*X): 
    #         return X  
    sc=1/k
    return stats.expon.rvs(scale=sc,size=size)
@nb.jit(nb.double(nb.int64,nb.double,nb.double),cache=True)
def _draw_P_B_Tr(pN,T,bg_rate):    
    while True:
        X = np.random.random_integers(0,pN)
        Y = np.random.uniform(0.0, 1)
        bg_bar=bg_rate*T
        fx=bg_bar**X*exp(-bg_bar)/math.factorial(X)
        if Y < fx:
            return X
@nb.jit(nb.double(nb.double,nb.double),cache=True)
def _drawE(e,v):
    while True:
        E=stats.norm.rvs(loc=e,scale=v,size=4)
        for e in E:
            if e>0.01 and e<0.99:
                return e
    
@nb.jit(nb.double(nb.int32,nb.double,nb.int32),cache=True)
def _drawA_fi_e(fi,e,size=1):
    # while True:
    #     X = np.random.random_integers(0,fi)
    #     Y = np.random.uniform(0.0, 1)        
    #     if Y < math.factorial(fi)/(math.factorial(fi-X)*math.factorial(X))* \
    #         e**X*(1-e)**(fi-X):
    #         return X              
    return stats.binom.rvs(fi, e, size=size)
@nb.jit(cache=True)
def _drawJ_Si2Sj(matP,n_states,i):
    P_i2j=copy.deepcopy(matP)
    P_i2j[i]=0
    P_i2j=P_i2j/sum(P_i2j)
    j=drawDisIdx(np.arange(n_states),P_i2j)
    return j

class ekTakeStep(object):
    def __init__(self, n,bounds,stepsize=0.1):
        self.stepsize = stepsize
        self.n=n
        self.bounds=bounds
        handler = logging.handlers.WatchedFileHandler(os.environ.get("LOGFILE", "./yourapp.log"))
        formatter = logging.Formatter(logging.BASIC_FORMAT)
        handler.setFormatter(formatter)
        self.root = logging.getLogger()
        self.root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
        self.root.addHandler(handler)
        logging.info(bounds) 

    def __call__(self, x):
        # s = self.stepsize 
        logging.info(x) 
        low=np.zeros(self.n)
        hi=np.zeros(self.n)
        for i in range(self.n*self.n):
            low[i]=x[i]-self.bounds[i][0]
            hi[i]=self.bounds[i][1]-x[i]
        # print(low,hi)        
        # x[0:self.n]=x[0:self.n]+ np.random.uniform(-2.*s, 2.*s,self.n)
        # x[self.n:self.n*self.n] = x[self.n:self.n*self.n]+np.random.uniform(-100, 100, self.n*self.n-self.n)
        x=x+np.random.uniform(low,hi)
        logging.info(x) 
        return x.tolist()

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
        print("=============================")
        print("self.params",self.params)
        rk=10**np.asarray(k)*1000
        self.matK=genMatK(n,rk)
        self.matP=genMatP(self.matK)
        print("matP",self.matP)
        print("matK",self.matK)

    def __init__(self,bursts_list,timesMs,maskAD,pN,T,SgDivSr,bg_ad_rate,bg_dd_rate,clk_p,histBinNum=50,\
            reSampleTimes=10,burstWindows=(0)):
        self.histBinNum=histBinNum
        self.clk_p=clk_p
        self.reSampleTimes=reSampleTimes
        self.mask_ad=maskAD
        self.bursts=bursts_list
        self.burstWindowsStart=0
        self.burstWindowsEnd=T.shape[0]
        if burstWindows==(0):
            pass
        elif len(burstWindows)>1:
            self.burstWindowsStart=burstWindows[0]
            self.burstWindowsEnd=burstWindows[1]            
        else:
            self.burstWindowsEnd=burstWindows
        self.burstWindowsRa=range(self.burstWindowsStart,self.burstWindowsEnd)
        self.numWindows=self.burstWindowsEnd-self.burstWindowsStart
        # assert self.numWindows==bg_ad_rate.shape[0]
        # assert self.numWindows==bg_dd_rate.shape[0]
        self.pN=pN
        self.T=T
        self.SgDivSr=SgDivSr
        self.timesMs=timesMs
        self.bg_ad_rate=bg_ad_rate
        self.bg_dd_rate=bg_dd_rate        
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

    def draw_P_B_Tr(self,idx,bg_rate):
        return _draw_P_B_Tr(self.pN(idx),self.T[idx],bg_rate) 

    def drawA_fi_e(self,fi,e):
        return _drawA_fi_e(fi,e)

    def monteCarlo(self,idx):
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
                    bins.append(bins[0]+mcSpendTime)
            bins=np.asarray(bins)#/self.clk_p
            # print("bins",bins)            
            m_ad=self.mask_ad[self.bursts[idx].istart:self.bursts[idx].istop+1]
            b_times=self.timesMs[self.bursts[idx].istart:self.bursts[idx].istop+1]
            m_dd=self.mask_dd[self.bursts[idx].istart:self.bursts[idx].istop+1]
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
            stransNum=len(f_i)
            F=sum(f_if)
            for i in range(stransNum):
                de=_drawE(self.E[sidx[i]],self.Evar[sidx[i]])                
                ai=_drawA_fi_e(f_i[i],de)
                bSgDivSr[sampleTime]=bSgDivSr[sampleTime]+ai
            bSgDivSr[sampleTime]=bSgDivSr[sampleTime]/F
        return bSgDivSr.tolist()

    def aTpda(self):
        # self.pdaBurstE=np.zeros(self.numWindows)
        # progress = ProgressBar()
        # for i in progress(self.burstWindowsRa):              
        #     self.pdaBurstE[i-self.burstWindowsStart]=self.monteCarlo(i)
        pdaBurstE=[]
        self.pdaBurstE=[]
        with Pool(processes=int(os.cpu_count()/2)) as pool:  
            pdaBurstE = pool.map(self.monteCarlo, self.burstWindowsRa)
        for be in pdaBurstE:
            self.pdaBurstE.extend(be)

    def aTpdaEK(self):
        self.set_nstates(self.n_states,self.fitP[:self.n_states],self.fitP[self.n_states:\
            self.n_states*self.n_states],self.fitP[self.n_states*self.n_states:])
        print("drawPic matP",self.matP)
        print("drawPic matK",self.matK)
        print("drawPic E",self.E)
        self.aTpda()
        rSgDivSr=self.SgDivSr[self.burstWindowsStart:self.burstWindowsEnd]
        SgDivSrR=np.histogram(rSgDivSr,bins=self.histBinNum,range=(0,5)) 
        tSgDivSr=np.zeros_like(SgDivSrR[0])
        for SgDivSrV in self.pdaBurstE:            
            idxHistBin=np.searchsorted(SgDivSrR[1], SgDivSrV)
            if idxHistBin<self.histBinNum:  
                tSgDivSr[idxHistBin]+= 1
        tSgDivSr=tSgDivSr/float(self.reSampleTimes)
        print("sum(SgDivSrR)",sum(SgDivSrR[0]))
        print("sum(tSgDivSr)",sum(tSgDivSr))
        return SgDivSrR,tSgDivSr

    def chiSqr(self):
        rSgDivSr=self.SgDivSr[self.burstWindowsStart:self.burstWindowsEnd]
        SgDivSrR=np.histogram(rSgDivSr,bins=self.histBinNum,range=(0,5)) 
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
                # chisqr+=(O/sO-E/sE)**2/(E/sE)
                chisqr+=(float(O-E))**2/float(O)
            else:
                n=n-1
        chisqr=chisqr/(n-self.n_states*(self.n_states+1))                
        print("chiSqr: ",chisqr)
        print("=======================================")        
        return chisqr 

    def chiSqrArrT(self,params):
        # print("params",params)
        start = time.time()
        n=self.n_states
        e=params[:self.n_states]
        k=params[self.n_states:self.n_states*self.n_states]
        evar=params[self.n_states*self.n_states:]
        self.set_nstates(n,e,k,evar)
        self.aTpda()        
        chi=self.chiSqr()
        end = time.time()
        print("speed is ",self.numWindows/(end-start),"burst/s" )
        return chi

    def calcP_N(self):
        hist, bin_edges = np.histogram(self.pN,bins=self.pNbins)
        pNhist=hist/self.numWindows
        center = (bin_edges[:-1] + bin_edges[1:]) / 2    
        print("pN range",min(self.pN),max(self.pN))
        photonN=np.arange(max(0,min(self.pN)),max(self.pN)+1,1)
        pNcount = prepData.interpSpl(center,pNhist, photonN)#,'linear')    
        pNcount[np.where(pNcount<0)]=0
        self.PN = dict(zip(photonN, pNcount))
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

    def findP(self):
        # bnds = (0,1)
        xmin = [0.02]*self.n_states
        xmink=[np.log10(0.01)]*(self.n_states*self.n_states-self.n_states)
        xminv = [0.001]*self.n_states
        xmin.extend(xmink)
        xmin.extend(xminv)
        xmax = [0.98]*self.n_states
        xmaxk=[np.log10(100)]*(self.n_states*self.n_states-self.n_states)
        xmaxv = [1.5]*self.n_states
        xmax.extend(xmaxk)
        xmax.extend(xmaxv)
        bounds = [(low, high) for low, high in zip(xmin, xmax)]
        minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
        # p=copy.deepcopy(self.params)
        # res = basinhopping(self.chiSqrArrT,self.params, minimizer_kwargs=minimizer_kwargs, \
        #     disp=True,stepsize=0.3,niter=20)#,interval=5,stepsize=0.1,niter=10
        res = differential_evolution(self.chiSqrArrT, strategy ='randtobest1exp',maxiter=20,bounds=bounds,\
                           disp=True)
        print (res)
        print (res.x)
        # print (self.numWindows)
        self.fitP=res.x
        return res.x