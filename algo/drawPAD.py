# coding: utf-8
import copy,os
import logging
from sklearn.metrics import r2_score
import numpy as np
from scipy import interpolate 
import datetime

from data import  aTpdaMpi
from mpi4py import MPI
import pickle,pathlib

def burstBin(full_fname,comm,pick,logname="pdampilogger.log"):  
    comm=None
    logger=None
    serialdict=aTpdaMpi.mpiloaddata(comm,logger,full_fname)
    #dict(bursts=sub_bursts_l,times=times,mask_ad=mask_ad,\
    #    T_burst_duration=T_burst_duration,SgDivSr=SgDivSr,clk_p=clk_p) 
    sub_bursts_l=serialdict['bursts']
    times=serialdict['times']
    mask_ad=serialdict['mask_ad']    
    T_burst_duration=serialdict['T_burst_duration']    
    SgDivSr=serialdict['SgDivSr']    
    clk_p=serialdict['clk_p']    
    mask_dd=serialdict['mask_dd']        
    bg_dd_rate=serialdict['bg_dd_rate']        
    bg_ad_rate=serialdict['bg_ad_rate']        
    # pdaIns=aTpdaMpi.pdaPy(comm,sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,0,0,clk_p,logger,\
    #     50,5,None)
    pdaIns=aTpdaMpi.pdaPy(comm,sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,\
        bg_ad_rate,bg_dd_rate,clk_p,logger,\
        85,5,None)        
    state=2
    pdaIns.setStateNum(state)
    pdaIns.fitP=[0.41204293319149393, 0.3514845133212948, 0.752816275833226, -0.31291654103739597, 0.7933347931616084, 0.8560774685539163, 0.2439534857827601, -0.8744333819896276, 0.22234927666601756, 0.08593357944187233, 0.17929881711426596, 0.005103389253265755]
    pdaIns.fitP=[0.7419905253260204, 0.44440661302580986, -0.9544647400402456, -0.5081247809890648, 0.0028045654508258577, 0.1210317387363557]
    SgDivSrR,tSgDivSr=pdaIns.aTpdaEK()
    r2=r2_score(SgDivSrR[0], tSgDivSr)  
    print("R^2:",r2)
        #############
        #rank=1
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    # plt.plot(SgDivSrR[1][1:],SgDivSrR[0])

    center = (SgDivSrR[1][:-1] + SgDivSrR[1][1:]) / 2 
    width = np.diff(SgDivSrR[1])        
    plt.bar(center,SgDivSrR[0] , align='center', width=width)
    plt.plot(SgDivSrR[1][1:],tSgDivSr,color="y")
    width = np.diff(SgDivSrR[1])
    print(np.sum(SgDivSrR[0]*width/(sum(SgDivSrR[0])*width)))        
    plt.title('PDA E fit')
    plt.savefig('temp.png')

def usage():  
    print("Usage:%s -i|--dat preprocessed.dat -o|--log logfile.log" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    size=0
    pick='out.pickle'
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "l:i:s:o", ["dat=", "hdf=","size=","pickle="])  
        for o, v in opts: 
            if o in ("-l", "--log"):
                savefn = v
            if o in ("-i", "--dat"):
                dbname=v
            if o in ("-o", "--pickle"):
                pick=v                
            if o in ("-s", "--size"):
                size = int(v.strip())
                print(size)

    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    comm=None
    if len(savefn)>1:     
        burstBin(dbname,comm,pick,savefn)
    else:
        burstBin(dbname,comm,pick)

