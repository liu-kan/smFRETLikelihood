# *- encoding: utf-8 -*-
import pickle,pathlib
import numpy as np
import datetime
from scipy import interpolate 
from math import *
from collections import Counter

class mburst():
    def __init__(self, istart, istop,start, stop):
            self.istart=istart
            self.istop=istop
            self.start=start
            self.stop=stop
 

from fretbursts import *
import h5py

# from mpi4py import MPI
def interpSpl(x_orig, y_orig, x_interp,kind="spline"):
    if type(x_orig) is list:
        x_orig=np.asarray(x_orig)
    if type(y_orig) is list:
        y_orig=np.asarray(y_orig)
    if type(x_interp) is list:
        x_interp=np.asarray(x_interp)         
    if kind=='linear':
        f = interpolate.interp1d(x_orig, y_orig)
        return f(x_interp)
    splrep = interpolate.splrep(x_orig, y_orig)
    return interpolate.splev(x_interp, splrep ) 

def loginfo(comm,log,msg,rank=0):
    # r=comm.Get_rank()
    # if comm.Get_rank()==rank or rank<0:
    log.info("{}".format(msg))

def duplicateValArray(a,b):
    # print (a[0])
    # print (b[0])
    s=np.concatenate((a,b),0)
    so=np.sort(s)
    # print(so)
    # print(Counter(so))
    return [item for item, count in Counter(so).items() if count > 1]



def prepHdf5(full_fname,logger,comm=None):
    time_bin = 1e-3  # 1 ms
    loginfo(comm,logger,'''
    ===========================================
    {}
    ==========================================='''.format(datetime.datetime.now()))    
    d = loader.photon_hdf5(full_fname)
    # bpl.plot_alternation_hist(d)
    loader.alex_apply_period(d)
    time_bin_clk = time_bin / d.clk_p   
    bg_time_s=30
    d.calc_bg(fun=bg.exp_fit, time_s=bg_time_s, tail_min_us='auto', F_bg=2.5)
    bg_dd=d.bg[Ph_sel(Dex='Dem')][0]
    bg_ad=d.bg[Ph_sel(Dex='Aem')][0]
    bg_ad_rate=np.mean(bg_ad)*time_bin
    bg_dd_rate=np.mean(bg_dd)*time_bin
    bg_time=np.arange(d.time_min+bg_time_s/2,d.time_max,bg_time_s)
    # bg_dd_splrep = interpolate.splrep(bg_time, bg_dd)
    # bg_ad_splrep = interpolate.splrep(bg_time, bg_ad)
    bg_time_bins_l=[]
    

    # print("print(len(d.bg[Dex='Dem']))",len(d.bg[Ph_sel(Dex='Dem')][0]))
    # print((d.bg[Ph_sel(Dex='Dem')]))
    # dplot(d, timetrace_bg)
    # dplot(d, timetrace)
    # xlim(10, 20)
    # ylim(-50, 50)
    d.burst_search()
    ds = d.select_bursts(select_bursts.naa, th1=3, computefret=False)
    ds = ds.select_bursts(select_bursts.size, th1=25, computefret=False)
    dsfuse = ds.fuse_bursts(ms=0)
    # dsfuse.leakage = 0.07
    # alex_jointplot(dsfuse)
    # dplot(dsfuse, hist_fret, show_kde=True)

    nanotimes = d.nanotimes[0]
    # nanotimes_d = nanotimes[d.get_D_em()]
    # nanotimes_a = nanotimes[d.get_A_em()]
    nanotimes_dd = nanotimes[d.get_D_em_D_ex()]


    # hist_params = dict(bins=range(3000), histtype='step', alpha=0.6, lw=1.5)
    # #hist(nanotimes, color='k', label='Total ph.', **hist_params)
    # hist(nanotimes_d, color='g', label='D. em. ph.', **hist_params)
    # hist(nanotimes_a, color='r', label='A. em. ph.', **hist_params)
    # plt.legend()
    # plt.yscale('log')

    roi = dict(E1=0, E2=1, S1=0.0, S2=0.97, rect=False)
    d_fret_mix = dsfuse.select_bursts(select_bursts.ES, **roi)
    # g = alex_jointplot(d_fret_mix)
    # bpl.plot_ES_selection(g.ax_joint, **roi);

    # dplot(d_fret_mix, hist_fret, show_kde=True)
    from fretbursts.phtools.burstsearch import Burst, Bursts
    times = d.ph_times_m[0]  # timestamps array
    
    bursts = d_fret_mix.mburst[0]
    # logger.info("mes time {}s".format(d.time_max-d.time_min))
    loginfo(comm,logger,"mes time {}s".format(d.time_max-d.time_min))
    loginfo(comm,logger,"Number of bursts: #{}".format(bursts.num_bursts))
    # logger.info("Number of bursts: #{}".format(bursts.num_bursts))
    mask_dd = d.get_ph_mask(ph_sel=Ph_sel(Dex='Dem'))   # donor excitation, donor emission
    mask_ad = d.get_ph_mask(ph_sel=Ph_sel(Dex='Aem'))   # donor excitation, acceptor emission
    mask_aa = d.get_ph_mask(ph_sel=Ph_sel(Aex='Aem'))   # acceptor excitation, acceptor emission
    timesad=times[mask_ad]
    sub_bursts_l = []
    bleachingBurst=0
    # bext.burst_data(d_fret_mix)
    for burst in bursts:
        msburst=times[burst.istart:burst.istop+1]
        m_dd=mask_dd[burst.istart:burst.istop+1]
        m_ad=mask_ad[burst.istart:burst.istop+1]
        m_aa=mask_aa[burst.istart:burst.istop+1]
        
        avgdd=np.mean(msburst[m_dd])*d.clk_p*1e3
        avgda=np.mean(msburst[m_ad])*d.clk_p*1e3
        avgaa=np.mean(msburst[m_aa])*d.clk_p*1e3
        if abs(avgdd-avgda)>0.5 or abs(avgaa-avgda)>0.5 or (len(msburst[m_dd])+len(msburst[m_ad]))<24:
            bleachingBurst=bleachingBurst+1
            continue
        # Compute binning of current bursts
        sub_bursts_l.append(mburst(istart=burst.istart, istop=burst.istop,
                                    start=burst.start, stop=burst.stop))
    # fBursts=Bursts.from_list(sub_bursts_l)    
    loginfo(comm,logger,"bleachingBurst:{}".format(bleachingBurst))
    # logger.info("bleachingBurst:{}".format(bleachingBurst))
    from fretbursts.phtools.burstsearch import count_ph_in_bursts
    PR=[]
    Tau=[]
    SgDivSr=[]
    pN=[]
    pNa=[] 
    T_burst_duration=[]
    dd_count=[]
    ad_count=[]
    aa_count=[]
    gamma=0.34
    beta=1.42
    DexDirAem=0.08
    Dch2Ach=0.07 
    i=-1
    count0bins=0
    for burst in sub_bursts_l:
        msburst = times[burst.istart:burst.istop + 1]
        m_dd = mask_dd[burst.istart:burst.istop + 1]
        m_ad = mask_ad[burst.istart:burst.istop + 1]
        m_aa = mask_aa[burst.istart:burst.istop + 1]
        # print (len(msburst),msburst)
        # print(len(m_ad),m_ad)
        # print(len(msburst[m_dd]))
        counts_dd = len(msburst[m_dd])
        counts_ad = len(msburst[m_ad])
        counts_aa = len(msburst[m_aa])
        T=(msburst[-1]-msburst[0])*d.clk_p
        T_burst_duration.append(T)
        # pr=(counts_ad / (counts_dd*gamma + counts_ad)).tolist()
        fe=(counts_ad *(1-DexDirAem)-Dch2Ach*counts_dd)/ ((gamma-\
                                Dch2Ach)*counts_dd + (1-DexDirAem)*counts_ad)
        N=counts_dd+counts_ad
        Na=counts_ad
        # SgSr=counts_ad/(counts_dd+counts_ad)
        dd=counts_dd
        ad=counts_ad
        i=i+1
        j=-1        
        # fe = SgSr
        if np.isnan(fe) or np.isinf(fe) or fe==0:
            continue
        SgDivSr.append(fe)
        pN.append(N)
        pNa.append(Na)
        bg_time_bins_l.append((msburst[0]+msburst[-1])/2*d.clk_p)
        dd_count.append(dd)
        ad_count.append(ad)
        aa_count.append(counts_aa)


    bg_time_bins=np.asarray(bg_time_bins_l)
    print("len(bg_time_bins):{}".format(len(bg_time_bins)))    
    
    if len(bg_dd)==len(bg_time):
        bg_dd_bins = interpSpl(bg_time,bg_dd, bg_time_bins) *time_bin
        bg_ad_bins = interpSpl(bg_time,bg_ad, bg_time_bins)  *time_bin
    else:
        bg_dd_bins = interpSpl(bg_time,bg_dd[:-1], bg_time_bins) *time_bin
        bg_ad_bins = interpSpl(bg_time,bg_ad[:-1], bg_time_bins)  *time_bin        
    F_RT=np.asarray(ad_count)-bg_ad_bins*time_bin
    F_G=np.asarray(dd_count)-bg_dd_bins*time_bin
    F_RTraw=np.asarray(ad_count)
    F_Graw=np.asarray(dd_count)
    F_AA=np.asarray(aa_count)
    loginfo(comm,logger,"F_RT min:{},num F_RT<=0:{}".format(min(F_RT),np.where(F_RT<=0)[0].shape[0]))
    loginfo(comm,logger,"F_G min:{},num F_G<=0:{}".format(min(F_G),np.where(F_G<=0)[0].shape[0]))
    
    pN=np.asarray(pN)
    bins_idx=np.where((pN>=25 )& (F_RT>0) &(F_G>0))   #& (F_RTraw>=5) &(F_Graw>=8)&(F_AA>=4)
    # =prepData.duplicateValArray(F_idx,N_idx[0])
    pN=pN[bins_idx]
    SgDivSr=np.asarray(SgDivSr)
    SgDivSr=SgDivSr[bins_idx]
    bg_ad_bins=bg_ad_bins[bins_idx]
    bg_dd_bins=bg_dd_bins[bins_idx]
    F_RT=F_RT[bins_idx]
    F_G=F_G[bins_idx]
    # plt.plot(bg_time, bg_dd, "o",  label=u"原始数据")    
    # plt.plot(bg_time_bins, bg_dd_bins, label=u"B-spline插值")
    # plt.legend()
    # plt.show()    
    numWindows=len(pN)
    loginfo(comm,logger,"numWindows:{}".format(numWindows))    
    loginfo(comm,logger,"bg_dd_bins#:{}".format(len(bg_dd_bins)))    
                
    T_burst_duration=np.asarray(T_burst_duration,dtype=np.float64)
    loginfo(comm,logger,"mask_ad.dtype{}".format(mask_ad.dtype))
    print("bg_ad_rate",bg_ad_rate)
    return sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,bg_ad_rate,bg_dd_rate,d.clk_p

def savedata(comm,logger,dictdata,filename="pdampi.dat"):
    # p=pathlib.Path(pathlib.Path.home(),"tmp",filename)    
    # pathlib.Path('data').mkdir(parents=True, exist_ok=True)
    p=pathlib.Path(filename)
    fn=str(p)    
    buf=pickle.dumps(dictdata,protocol=4)
    sbuf=len(buf)
    fh = open(fn, 'wb')
    fh.write(buf)
    fh.close()
    loginfo(comm,logger,"save dict data")
    print ("save data at ",p.resolve())
    return sbuf

def saveHDF5(savefn,sub_bursts_l,times,mask_ad,mask_dd,T_burst_duration,SgDivSr,clk_p,bg_ad_rate,bg_dd_rate):
    print("clk_p",clk_p)
    f=h5py.File(savefn,"w")
    dt=h5py.special_dtype(vlen=str)
    filename = f.create_dataset('filename', (1,), dtype=dt)
    filename[0] = savefn
    bursts = f.create_group('sub_bursts_l')
    sburst = len(sub_bursts_l)
    istart = bursts.create_dataset('istart', (sburst,), dtype=np.int64)
    istop = bursts.create_dataset('istop', (sburst,), dtype=np.int64)
    stop = bursts.create_dataset('stop', (sburst,), dtype=np.int64)
    start = bursts.create_dataset('start', (sburst,), dtype=np.int64)
    for i in range(sburst):
        istart[i] = sub_bursts_l[i].istart
        start[i] = sub_bursts_l[i].start
        istop[i] = sub_bursts_l[i].istop
        stop[i] = sub_bursts_l[i].stop
    f.create_dataset('times', data=times, dtype=np.int64)
    f.create_dataset('mask_ad', data=mask_ad, dtype=np.int8)
    f.create_dataset('mask_dd', data=mask_dd, dtype=np.int8)
    f.create_dataset('T_burst_duration',
                     data=T_burst_duration, dtype=np.float32)
    f.create_dataset('SgDivSr', data=SgDivSr, dtype=np.float32)
    f.create_dataset('clk_p',(1,), data=clk_p, dtype=np.float32)
    f.create_dataset('bg_ad_rate',(1,), data=bg_ad_rate, dtype=np.float32)
    f.create_dataset('bg_dd_rate',  (1,), data=bg_dd_rate,dtype=np.float32)
    f.flush()
    v=f['/clk_p']
    print('delta clk_p', clk_p-v[0])
    v=f['T_burst_duration']
    print('delta T_burst_duration', T_burst_duration[1000]-v[1000])
    f.close()

if __name__ == '__main__':
    a=np.array([1, 2, 3, 5, 7, 9, 10])
    b=[]
    for aa in a:
        b.append(pdaP_B(1,aa))
    print(b)
