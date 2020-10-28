
# coding: utf-8
import copy
import logging

import numpy as np
from scipy import interpolate 

# sns = init_notebook()
from data import  prepData

def burstBin(full_fname,savefn='pdampi.dat',logname="pdampilogger.log"):      
    times_s=0
    mask_ad_s=0
    logging.basicConfig(filename=logname,level=0)
    logger = logging.getLogger('pda')
    comm=None
    sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,bg_ad_rate,bg_dd_rate,clk_p=prepData.prepHdf5(full_fname,logger)
    # burstlen=len(sub_bursts_l)
    # chunkLists=list(prepData.chunks(range(burstlen), clsize))
    # testrank=3

    prepData.loginfo(comm,logger,"T_burst_duration[10] {} before b".format(T_burst_duration[10]),0)     
    prepData.loginfo(comm,logger,"SgDivSr[10] {} before b".format(SgDivSr[10]),0) 
    prepData.loginfo(comm,logger,"sub_bursts_l[10] {} before b".format(sub_bursts_l[10]),0) 
    prepData.loginfo(comm,logger,"times[-1] {} before b".format(times[-1]))    
    prepData.loginfo(comm,logger,"mask_ad[4356] {} before b".format(mask_ad[4356]))  
    prepData.loginfo(comm,logger,"mask_ad[-1] {} before b".format(mask_ad[-100:-1]))          
    dictdata=dict(bursts=sub_bursts_l,times=times,mask_ad=mask_ad,mask_dd=mask_dd,\
        T_burst_duration=T_burst_duration,SgDivSr=SgDivSr,clk_p=clk_p,\
        bg_ad_rate=bg_ad_rate,bg_dd_rate=bg_dd_rate) 
    # sbuf=prepData.savedata(comm,logger,dictdata,savefn)
    # print("Data size is ",sbuf," bytes!")
    # prepData.loginfo(comm,logger,"Data size is {} bytes!".format(sbuf))    
    
    prepData.saveHDF5( savefn, \
        sub_bursts_l,times,mask_ad,mask_dd,\
        T_burst_duration,SgDivSr,clk_p,\
        bg_ad_rate,bg_dd_rate )

def usage():  
    print("Usage:%s -i inputfilename.hf5 -o outputfilename.hdf5" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "o:i:")  
        for o, v in opts: 
            if o in ("-o"):
                savefn = v
            if o in ("-i"):
                dbname=v
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    if len(dbname)<1:
        usage()    
        sys.exit(1)    
    if len(savefn)>1:     
        burstBin(dbname,savefn)
    else:
        burstBin(dbname)

