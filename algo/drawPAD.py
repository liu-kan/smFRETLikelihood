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
        50,5,None)        
    state=2
    pdaIns.setStateNum(state)
    pdaIns.fitP=[0.6681530537400076,0.4784057184204892,0.962157335443822, 0.9454653214094496, 0.31922221822241714, 0.40184797098402764]
    pdaIns.fitP=[0.24847515799718217, 0.5182234987325258, 0.6704834891885445, 0.3846255599778262, 1.1311325208824514, 1.1030106563348268, 1.6278606238336475, 1.4575266144868801, 0.17129647270017628, 0.3601388637318751, 0.4939339490171665, 0.009230574939612601]
    pdaIns.fitP=[0.4822779436560594, 0.6775371214614038, -0.22855479513013477, 0.7752800770484245, 0.12264960511418932, 0.08244340848886908]

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

