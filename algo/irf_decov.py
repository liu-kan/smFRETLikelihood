# -*- coding: utf-8 -*-
import binRawData
import BGrate
import numpy as np
from lmfit import Model,Parameters
import matplotlib.pyplot as plt

def prepareData(histIRF,binData,binIdx=-1,sampleNum=1000,T0=0,fft=True):    
    h=None
    b=None
    if binIdx==-1:
        h,b,sdt=binRawData.statsDelayTime(binData,sampleNum,"D")#,bin0=100,binLen=2)    
    else:
        h,b,sdt=binRawData.statsDelayTime(binData,sampleNum,"D",binIdx,1)
    histlen=int(sampleNum/2)    
    # binRawData.statsBins(irfbinData)
    # print("IRF photonEffTime:",irfbinData['chs']['DexDem']["photonEffTime"])
    # print("Data Max PIE delay(ns):",b[0][-1])       
    startt=max(0,int(T0/64*sampleNum)-1)
    if fft:
        return h[0][startt:histlen],histIRF[0][startt:histlen],np.linspace(0,histlen-startt-1,histlen-startt),sdt
    else:
        return h[0][startt:histlen],histIRF[0][startt:histlen],b[0][startt:histlen],sdt


# define the single exponential model
#jumpexpmodel(t,irf,paras):
def jumpexpmodel(t,irf,**paras):#,y0,x0,tau1,ampl1,tau2,ampl2,tau3,ampl3,tau4,ampl4,tau5,ampl5):#,tau6,ampl6):# Lifetime decay fit Author: Antonino Ingargiola - Date: 07/20/2014
    ymodel=np.zeros(t.size)
    # paras= params.valuesdict()
    c=paras['x0']
    y0=paras['y0']
    n=len(irf)
    irf_s1=np.remainder(np.remainder(t-np.floor(c)-1, n)+n,n)
    irf_s11=(1-c+np.floor(c))*irf[irf_s1.astype(int)]
    irf_s2=np.remainder(np.remainder(t-np.ceil(c)-1,n)+n,n)
    irf_s22=(c-np.floor(c))*irf[irf_s2.astype(int)]
    irf_shift=irf_s11+irf_s22
    irf_reshaped_norm=irf_shift/sum(irf_shift)
    numcom=int((len(paras)-2)/3)
    for i in range(numcom):
        idx=str(i+1)
        ymodel += paras['ampl'+idx]*np.exp(-(t)/paras['tau'+idx])    
    # ymodel+= ampl1*np.exp(-(t)/tau1)        
    # ymodel+= ampl2*np.exp(-(t)/tau2)        
    # ymodel+= ampl3*np.exp(-(t)/tau3)        
    # ymodel+= ampl4*np.exp(-(t)/tau4)        
    # ymodel+= ampl5*np.exp(-(t)/tau5)        
    # ymodel+= ampl6*np.exp(-(t)/tau6)        
    z=Convol(ymodel,irf_reshaped_norm)
    z+=y0
    return z


def exp1model(t,irf,y0,x0,tau1,ampl1):
    ymodel=np.zeros(t.size)
    # c=paras['x0']
    # y0=paras['y0']
    c=x0    
    n=len(irf)
    irf_s1=np.remainder(np.remainder(t-np.floor(c)-1, n)+n,n)
    irf_s11=(1-c+np.floor(c))*irf[irf_s1.astype(int)]
    irf_s2=np.remainder(np.remainder(t-np.ceil(c)-1,n)+n,n)
    irf_s22=(c-np.floor(c))*irf[irf_s2.astype(int)]
    irf_shift=irf_s11+irf_s22
    irf_reshaped_norm=irf_shift/sum(irf_shift)
    # numcom=int((len(paras)-2)/2)
    # for i in range(numcom):
    #     idx=str(i+1)
    #     ymodel+ =paras['ampl'+idx]*np.exp(-(t)/paras['tau'+idx])    
    ymodel+= ampl1*np.exp(-(t)/tau1)        
    # ymodel+= ampl2*np.exp(-(t)/tau2)        
    # ymodel+= ampl3*np.exp(-(t)/tau3)        
    z=Convol(ymodel,irf_reshaped_norm)
    z+=y0
    return z

def plot_fit(time,irf, ydata, params, yscale='log', zoom_origin=False):
    """
    Function to plot data and model function.
    """
    plt.semilogy(time, ydata, marker='.')
    #plt.semilogy(time, jumpexpmodel(time,irf, params),
    plt.semilogy(time, jumpexpmodel(time,irf, **{k: v.value for k, v in params.items()}),    
         color='r', lw=2.5, alpha=0.7)
    if zoom_origin:
        dt = time[1] - time[0]
        t0 = params['offset'] - 50*dt
        plt.xlim(t0 - 50*dt, t0 + 50*dt)

def Convol(x,h): # change in convolution calcualataion
    #need same length of x and h
    X=np.fft.fft(x)
    H=np.fft.fft(h)
    xch=np.real(np.fft.ifft(X*H))
    return xch

def fitIRF(decayData,irf,delayT,pars,model,method='differential_evolution'):
    wWeights=1/np.sqrt(decayData+1)# check for divide by zero,  have used +1 to avoid dived by zero
    mod=Model(model,independent_vars=['t','irf'])    
    # fit this model with weighting , parameters as initilized and print the result  differential_evolution' leastsq'
    result = mod.fit(decayData,pars,weights=wWeights,method=method,t=delayT,irf=irf,iter_cb=iter_cbX)    
    return result

def calcTauOf1Bin(histIRF,binData,binIdx,sampleNum,T0,method='differential_evolution'):    
    decay1,irf1,x,_dt=prepareData(histIRF,binData,binIdx,sampleNum=sampleNum,T0=T0)
    maxI=max(decay1)*max(irf1)
    minI=min(min(decay1),min(irf1))
    params = Parameters()
    params.add('x0', value=0,min=-sampleNum/2,max=sampleNum/2)#, vary=False)
    params.add('y0', value=minI)#,min=-minI,max=maxI)#, vary=False)
    params.add('tau1', value=sampleNum/8,min=sampleNum/2000,max=sampleNum)    
    params.add('ampl1', value=(maxI**0.5)/3)#,min=minI,max=maxI)#, vary=False)
    # params.add('tau2', value=sampleNum/8,min=sampleNum/36,max=sampleNum/2)    
    # params.add('ampl2', value=(minI+maxI)/2,min=minI,max=1.5*maxI)#, vary=False)
    # params.add('tau3', value=sampleNum/8,min=sampleNum/36,max=sampleNum/2)    
    # params.add('ampl3', value=(minI+maxI)/2,min=minI,max=1.5*maxI)#, vary=False)
    result=fitIRF(decay1,irf1,x,params,exp1model,method)   
    return result.params['tau1'].value*64/sampleNum,result.redchi

def percentage(idx,x,irf1,*params):
    numcom=int((len(params)-2)/2)    
    x0=params[0]
    # print("x0",x0)
    y0=params[1]
    # print("y0",y0)
    eall=np.zeros(len(x))
    ev=np.zeros(numcom)        
    for i in range(numcom):
        e1=exp1model(x,irf1,y0,x0,params[2+i*2],params[3+i*2])-y0
        eall+=e1
        ev[i]=sum(e1)
        # print(i,numcom)        
    # for i in range(numcom):
        # idx=str(i+1)
        # print("Tau"+idx+": ",result.params['tau'+idx]*64/1000," ns")
        # print("p:",ev[i]/sum(eall))    
    # print('p'+str(idx),ev[idx-1]/sum(eall))
    # print('tau'+str(idx),params[2+(idx-1)*2])
    # print('ampl'+str(idx),params[3+(idx-1)*2])    
    return ev[idx-1]/sum(eall)


def iter_cbX(params, iter, resid, *args, **kws):
    numcom=int((len(params)-2)/3)
    # print("resid:",resid)
    # print("resid10:",resid*10)
    # print(" ITER ", iter, ["%s" % p for p in kws.keys()])
    x=kws['t']
    irf=kws['irf']
    par=(params['x0'],params['y0'],)
    for i in range(numcom):
        idx=str(i+1)
        par+=(params['tau'+idx],params['ampl'+idx])
    p=np.zeros(numcom)
    for i in range(numcom):             
        p[i]=percentage(i+1,x,irf,*par)
        idx=str(i+1)
        if p[i]<params['p'+idx].min or p[i]>params['p'+idx].max:
            # print('OFp'+str(idx),p)  
            return
    print("no of")  
    dw=[]
    for v in par:
        dw.append(v.value)
    dw.extend(p.tolist())
    appendCSV('data/csv.txt',dw)   
def genParams(pdata,i):
    sizep=pdata.shape[1]
    params=dict()
    numcom=int((sizep-2)/3)
    params['x0']=pdata[i,0]
    params['y0']=pdata[i,1]
    for ip in range(numcom):
        idx=str(ip+1)
        params['tau'+idx]=pdata[i,2+ip*2]
        params['ampl'+idx]=pdata[i,3+ip*2]
    return params

def genParamsList(numcom):
    pstr='x0,y0'
    for i in range(numcom):
        idx=str(i+1)
        pstr+=',tau'+idx
        pstr+=',ampl'+idx
    return pstr

def createCSV(fn):
    fo = open(fn, "w")
    fo.close()
import csv
def appendCSV(fn,line):
    with open(fn, 'a') as f:
        csvwriter = csv.writer(f)
        # for line in lines:
        csvwriter.writerow(line)

def findBetter(time,irf,odata,pdata,plot=False):    
    sizen=pdata.shape[0]
    res=np.zeros(sizen)
    resi=np.zeros([sizen,len(odata)])
    for i in range(sizen):
        params=genParams(pdata,i)
        mdata= jumpexpmodel(time,irf, **{k: v for k, v in params.items()})
        resi[i,]=np.asarray(odata)-np.asarray(mdata)
        res[i]=np.dot(resi[i,],resi[i,])
    idx=np.argmin(res)        
    if not plot:
        return idx
    params=genParams(pdata,idx)
    best_fit=jumpexpmodel(time,irf, **{k: v for k, v in params.items()})
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.semilogy(time,odata,'r-',time,best_fit,'b')
    plt.subplot(2,1,2)
    plt.plot(time,resi[idx,].tolist())
    plt.show() 

    

if __name__=='__main__':
    import pickle,sys,getopt
    irfdbname="data/alexa488_IRF_32MHz_PIE_3KCPS.sqlite"
    dbname="data/21c_224c.sqlite"
    savefn='data/'+\
        dbname.split('/')[-1].split('.')[-2]+'_gd'+".pickle"
    createCSV('data/csv.txt')
    fromdb=True
    binMs=1
    decay1=None;irf1=None;x=None;
    sampleNum=1000
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "fs:b:", ["fromfile", "samplenum=", "binms="])  
        for o, v in opts: 
            if o in ("-f", "--fromfile"):
                fromdb=False
            if o in ("-s", "--samplenum"):
                sampleNum = int(v)
            if o in ("-b", "--binms"):
                binMs = v                
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")        
        sys.exit(1); 
    if fromdb:
        br=BGrate.calcBGrate(dbname,0,400,T0=6,Tlen=267.97)
        binData=binRawData.binRawData(br,dbname,binMs)
        irfbr=BGrate.calcBGrate(irfdbname,20,400)#,T0=0.0,Tlen=600)
        irfbinData=binRawData.binRawData(irfbr,irfdbname,binMs)    
        with open(savefn, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([irfbinData,binData], f,protocol=-1)
    else:
        irfbinData,binData=pickle.load(open(savefn,'rb'))    
    # x,decay1,irf1=np.lo oadtxt(r"data/tcspcdatashifted.csv",delimiter=',',unpack=True,dtype='float')
    hi,bi,meandt=binRawData.statsDelayTime(irfbinData,sampleNum,"D")#,bin0=100,binLen=2)
    decay1,irf1,x,sdt=prepareData(hi,binData,sampleNum=sampleNum,T0=6)
    maxI=max(decay1)*max(irf1)
    minI=min(min(decay1),min(irf1))

    E=[0.3736,0.6739,0.8314]
    upb=11.5
    lowb=0.05
    TauD0=4.1
    npE=np.asarray(E)
    npTau=(1-npE)*TauD0*1000/64
    iniTauD0=TauD0*1000/64
    P=[0.2562,0.2995,0.4443]
    npP=np.asarray(P)
    deltaP=3e-2
    params = Parameters()
    params._asteval.symtable['tt'] =x
    params._asteval.symtable['irf'] =irf1
    params._asteval.symtable['func'] =percentage
    params.add('x0', value=0,min=-sampleNum,max=sampleNum)#, vary=False)
    params.add('y0', value=minI,min=-minI,max=maxI)#, vary=False)
    params.add('tau1', value=npTau[0],min=npTau[0]*lowb,max=npTau[0]*upb)    
    params.add('ampl1', value=(maxI**0.5)/3,min=minI/12,max=maxI*2)#, vary=False)
    params.add('tau2', value=npTau[1],min=npTau[1]*lowb,max=npTau[1]*upb)    
    params.add('ampl2', value=(maxI**0.5)/3,min=minI/12,max=maxI*2)#, vary=False)
    params.add('tau3', value=npTau[2],min=npTau[2]*lowb,max=npTau[2]*upb)    
    params.add('ampl3', value=(maxI**0.5)/3,min=minI/12,max=maxI*2)#, vary=False)
    params.add('p1', expr='func(1,tt,irf,'+genParamsList(3)+')', value=npP[0], min=npP[0]-deltaP, max=npP[0]+deltaP)
    params.add('p2', expr='func(2,tt,irf,'+genParamsList(3)+')', value=npP[1], min=npP[1]-deltaP, max=npP[1]+deltaP)
    params.add('p3', expr='func(3,tt,irf,'+genParamsList(3)+')', value=npP[2], min=npP[2]-deltaP, max=npP[2]+deltaP)

    # params.add('tau6', value=sampleNum/8,min=sampleNum/2000,max=sampleNum)    
    # params.add('ampl6', value=(maxI**0.5)/3,min=minI,max=maxI)#, vary=False)    
    

    result=fitIRF(decay1,irf1,x,params,jumpexpmodel)#,'leastsq')  
    print(result.fit_report()) 
    print(result.params['tau1'].value*64/sampleNum) 
    print(result.redchi) 
    # print("params['tau2']",result.params['tau2'].value) 


    numcom=int((len(result.params)-2)/3)
    y0=result.params['y0'].value
    x0=result.params['x0'].value

    eall=np.zeros(len(x))
    ev=np.zeros(numcom)    
    for i in range(numcom):
        idx=str(i+1)
        e1=exp1model(x,irf1,y0,x0,result.params['tau'+idx].value,result.params['ampl'+idx].value)-y0
        eall+=e1
        ev[i]=sum(e1)
    # plt.subplot(2,1,1)        
    # plt.semilogy(x,eall+y0,'y-')
    print("Tau^{hat}",npTau*64/1000)
    for i in range(numcom):
        idx=str(i+1)
        print("Tau"+idx+": ",result.params['tau'+idx].value,result.params['tau'+idx].value*64/1000," ns")
        print("p"+idx+": ",ev[i]/sum(eall))
    my_data = np.genfromtxt('data/csv.txt', delimiter=',')
    findBetter(x,irf1,decay1,my_data,True)



    
    # plt.figure(2)
    # decay1,irf1,x=prepareData(hi,binData,sampleNum=sampleNum,T0=6,fft=False)
    # # plot_fit(x,irf1, decay1, result.params)    
    # plt.semilogy(x,irf1)
    # plt.semilogy(x,decay1)
   




