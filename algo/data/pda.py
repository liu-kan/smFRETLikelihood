import numpy as np
from math import *
from scipy.optimize import minimize
from data import prepData
class pdaPy:

    def __init__(self,pN,SgDivSr,bg_ad_bins,bg_dd_bins,F_RT,F_G,bg_dd_bar,\
    bg_ad_bar,pNbins=79,histBinNum=50):
        self.histBinNum=histBinNum
        self.pNbins=pNbins
        self.numWindows=pN.shape[0]
        assert self.numWindows==SgDivSr.shape[0]
        assert self.numWindows==bg_ad_bins.shape[0]
        assert self.numWindows==bg_dd_bins.shape[0]
        assert self.numWindows==F_G.shape[0]
        assert self.numWindows==F_RT.shape[0]
        self.pN=pN
        self.pNbins=pNbins
        self.SgDivSr=SgDivSr
        self.bg_ad_bins=bg_ad_bins
        self.bg_dd_bins=bg_dd_bins
        self.F_RT=F_RT
        self.F_G=F_G
        self.bg_dd_bar=bg_dd_bar
        self.bg_ad_bar=bg_ad_bar
        print("bg_dd_bar",bg_dd_bar)
        print("bg_ad_bar",bg_ad_bar)
        self.sumN=sum(pN)
        self.calcP_N()
    
    def set_burst_duration(self,T):
        self.T=T
    def set_bg_rate(self,dd,ad):
        self.bgRateDD=dd
        self.bgRateAD=ad
    def draw_P_B_Tr(self,T,r):

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

    def pdaBin(self,i):
        pass
    
    def pdaE(self,e,SgDivSrRange=(0,5)):
        mSgDivSr=self.pdaBinE(e)
        SgDivSrR=np.histogram(self.SgDivSr,bins=self.histBinNum,range=(0,5)) 
        tSgDivSr=np.zeros_like(SgDivSrR[0])
        for i in range(0,self.numWindows):
            SgDivSrV=self.SgDivSr[i]            
            idxHistBin=np.searchsorted(SgDivSrR[1], SgDivSrV)
            if idxHistBin<self.histBinNum:    
                tSgDivSr[idxHistBin]+= mSgDivSr[i]
        return SgDivSrR,tSgDivSr
        
    def pdaBinE(self,e):        
        tSgDivSr=np.zeros_like(self.SgDivSr)
        for i in range(0,self.numWindows):
            SgDivSrV=self.SgDivSr[i]
            tSgDivSr[i]=self.pSgDivSr(e,self.P_N(self.pN[i]),self.F_RT[i],self.F_G[i],\
            self.bg_dd_bins[i],self.bg_ad_bins[i],self.bg_dd_bar,self.bg_ad_bar)        
        return tSgDivSr*self.numWindows
    
    def chiSqr(self,e):
        mSgDivSr=self.pdaBinE(e)
        deltaX=(self.pN-mSgDivSr)/self.numWindows
        return 1/(self.numWindows-1)*sum(np.inner(deltaX,deltaX)/\
        (self.pN/self.numWindows))
    def rSqr(self,e):
        mSgDivSr=self.pdaBinE(e)
        SgDivSrR=np.histogram(self.SgDivSr,bins=self.histBinNum,range=(0,5)) 
        tSgDivSr=np.zeros_like(SgDivSrR[0])
        for i in range(0,self.numWindows):
            SgDivSrV=self.SgDivSr[i]            
            idxHistBin=np.searchsorted(SgDivSrR[1], SgDivSrV)
            if idxHistBin<self.histBinNum:    
                tSgDivSr[idxHistBin]+= mSgDivSr[i]
        deltaX=SgDivSrR[0]-tSgDivSr        
        SStot=SgDivSrR[0]-np.mean(SgDivSrR[0])
        return np.inner(deltaX,deltaX)/np.inner(SStot,SStot)

    def findP(self):
        bnds = (0,1)
        res = minimize(self.chiSqr, 0.5, method='Nelder-Mead')
        print (res)
        print (res.x)
        print (self.numWindows)
        return res.x