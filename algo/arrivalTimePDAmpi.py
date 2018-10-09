
# coding: utf-8
import copy,time,os
import logging

import numpy as np
from scipy import interpolate 
import datetime
from logging.handlers import RotatingFileHandler
from data import  aTpdaMpi
# import mpi4py
# # mpi4py.rc.initialize = False
# mpi4py.rc.finalize = False
from mpi4py import MPI
import pickle,pathlib

def burstBin(full_fname,state,comm,pick,logname="pdampilogger.log"):  
    rank=comm.Get_rank()
    testrank=33
    clsize=comm.Get_size()
    # logging.basicConfig(filename=logname,level=0)
    my_handler = RotatingFileHandler(logname, mode='a', maxBytes=3*1024*1024, \
                                    backupCount=5, encoding=None, delay=0)
    my_handler.setLevel(logging.INFO)                                    
    logger = logging.getLogger('pda')
    logger.setLevel(logging.INFO)
    logger.addHandler(my_handler)
    aTpdaMpi.loginfo(comm,logger,'''
    ===========================================
    {}
    ==========================================='''.format(datetime.datetime.now()))

    serialdict=aTpdaMpi.mpiloaddata(comm,logger,full_fname)
    #dict(bursts=sub_bursts_l,times=times,mask_ad=mask_ad,\
    #    T_burst_duration=T_burst_duration,SgDivSr=SgDivSr,clk_p=clk_p) 
    sub_bursts_l=serialdict['bursts']
    times=serialdict['times']
    mask_ad=serialdict['mask_ad']    
    mask_dd=serialdict['mask_dd']    
    T_burst_duration=serialdict['T_burst_duration']    
    SgDivSr=serialdict['SgDivSr']    
    clk_p=serialdict['clk_p']    
    bg_ad_rate=serialdict['bg_ad_rate']    
    bg_dd_rate=serialdict['bg_dd_rate']    
    if rank==0:
        burstlen=len(sub_bursts_l)
        chunkLists=list(aTpdaMpi.chunks(range(burstlen), clsize))
    else:        
        chunkLists=None        
    chunkLists=comm.scatter( chunkLists, root=0)    
    
    aTpdaMpi.loginfo(comm,logger,"T_burst_duration[10] {} after b".format(T_burst_duration[10]),testrank)     
    aTpdaMpi.loginfo(comm,logger,"SgDivSr[10] {} after b".format(SgDivSr[10]),testrank)         
    aTpdaMpi.loginfo(comm,logger,"chunkLists {} after b".format(chunkLists),-1) 
    aTpdaMpi.loginfo(comm,logger,"sub_bursts_l[10] {} after b".format(sub_bursts_l[10]),testrank)     
    aTpdaMpi.loginfo(comm,logger,"times[-1] {} after b".format(times[-1]),testrank)     
    aTpdaMpi.loginfo(comm,logger,"mask_ad[4356] {} after b".format(mask_ad[4356]),testrank)  
    aTpdaMpi.loginfo(comm,logger,"mask_ad[-1] {} after b".format(mask_ad[-10:-1]),testrank)  
    
    pdaIns=aTpdaMpi.pdaPy(comm,sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,\
        bg_ad_rate,bg_dd_rate,clk_p,logger,\
        70,5,chunkLists)
    pdaIns.set_nstates(state, np.random.rand(state).tolist(),[0]*(state-1)*state,[0.1]*state)
    mpistop=[0]
    fitE=None
    if rank==0:
        fitE=pdaIns.findP()
    else:
        while mpistop[0]==0:
            p=copy.deepcopy(pdaIns.params)
            pdaIns.chiSqrArrT(p,mpistop)    
    
    # MPI.Finalize()
    # comm.Barrier()
    # if rank!=0:
    #     MPI.Finalize()
    # else:
    #     aTpdaMpi.loginfo(comm,logger,"Finalizing other than rank 0")  
    # # comm.Barrier()        
    # if rank==0:
    #     p=pathlib.Path(pick)
    #     fn=str(p)    
    #     time.sleep(9)
    #     SgDivSrR,tSgDivSr=pdaIns.aTpdaEK()
    #     with open(fn, 'wb') as f:
    #         try:
    #             dictdata=dict(SgDivSrR=SgDivSrR,tSgDivSr=tSgDivSr)
    #             buf=pickle.dump(dictdata,protocol=4)
    #         except : # whatever reader errors you care about
    #             print ("Could not write file:", fn)
    #     MPI.Finalize()
        #############
        #rank=1
        # import matplotlib as mpl
        # mpl.use('Agg')
        # import matplotlib.pyplot as plt
        # plt.plot(SgDivSrR[1][1:],SgDivSrR[0])
        # plt.plot(SgDivSrR[1][1:],tSgDivSr,color="y")
        # width = np.diff(SgDivSrR[1])
        # print(np.sum(SgDivSrR[0]*width/(sum(SgDivSrR[0])*width)))        
        # plt.title('P(S_G/S_R) fit')        
        # plt.savefig('temp.png')

def usage():  
    print("Usage:%s -i|--dat preprocessed.dat -o|--log logfile.log" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    state=2
    pick='out.pickle'
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "l:i:s:o", ["log=", "dat=","state=","pickle="])  
        for o, v in opts: 
            if o in ("-l", "--log"):
                savefn = v
            if o in ("-i", "--dat"):
                dbname=v
            if o in ("-o", "--pickle"):
                pick=v                
            if o in ("-s", "--state"):
                state = int(v.strip())
                # print(state)

    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    comm=MPI.COMM_WORLD 
    if len(savefn)>1:     
        burstBin(dbname,state,comm,pick,savefn)
    else:
        burstBin(dbname,state,comm,pick)

