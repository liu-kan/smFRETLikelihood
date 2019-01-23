# coding: utf-8
import copy,os
import logging
from sklearn.metrics import r2_score
import numpy as np
from scipy import interpolate 
import datetime

from data import aTpdaMpiGA as aTpdaMpi
from mpi4py import MPI
import pickle,pathlib
import scipy.io as sio
from array import array

def burstBin(full_fname,comm,pick,logname="pdampilogger.log"):  
    comm=None
    logger=None
    serialdict=aTpdaMpi.mpiloaddata(comm,logger,full_fname)
    sub_bursts_l=serialdict['bursts']
    times=serialdict['times']
    mask_ad=serialdict['mask_ad']    
    T_burst_duration=serialdict['T_burst_duration']    
    SgDivSr=serialdict['SgDivSr']    
    clk_p=serialdict['clk_p']    
    mask_dd=serialdict['mask_dd']        
    bg_dd_rate=serialdict['bg_dd_rate']        
    bg_ad_rate=serialdict['bg_ad_rate']        
    pdaIns=aTpdaMpi.pdaPy(comm,sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,\
        bg_ad_rate,bg_dd_rate,clk_p,logger,\
        80,5,None)        
    state=2
    pdaIns.setStateNum(state)

    pdaIns.fitP=[0.48031496062992124, 0.55905511811023623, 4305.0733310748947, 7497.4068475248732, 0.39096480938416417, 0.321871945259042]
    SgDivSrR,tSgDivSr,matp,matk=pdaIns.aTpdaEK()


    center = (SgDivSrR[1][:-1] + SgDivSrR[1][1:]) / 2 
    width = np.diff(SgDivSrR[1])        
    SgDivSrR0=SgDivSrR[0]/(sum(SgDivSrR[0]))
    tSgDivSr=tSgDivSr/(sum(tSgDivSr))
    # plt.bar(center,SgDivSrR0 , align='center', width=width)
    # plt.plot(SgDivSrR[1][1:],tSgDivSr,color="y")    
    sio.savemat(logname, {'fret':SgDivSr*100,'histCounFit':SgDivSrR0,'xcol':SgDivSrR[1][1:]*100,\
        'p':matp,'pdacenter':center*100})    
    print(matp)

def usage():  
    print("Usage:%s -i|--dat preprocessed.dat -o|--log logfile.log" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    size=0
    pick='out.pickle'
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "m:i:", ["mat=", "dat="])  
        for o, v in opts: 
            if o in ("-m", "--mat"):
                savefn = v
            if o in ("-i", "--dat"):
                dbname=v
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    comm=None
    burstBin(dbname,comm,pick,savefn)
    