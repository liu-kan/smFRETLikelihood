
# coding: utf-8

import numpy as np
from scipy import interpolate 
from fretbursts import *
# sns = init_notebook()
from data import  prepData
from data import  aTpda

def burstBin(full_fname,mat_name):  
    time_bin = 1e-3  # 1 ms
     
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

    roi = dict(E1=0, E2=1, S1=0.0, S2=0.95, rect=False)
    d_fret_mix = dsfuse.select_bursts(select_bursts.ES, **roi)
    # g = alex_jointplot(d_fret_mix)
    # bpl.plot_ES_selection(g.ax_joint, **roi);

    # dplot(d_fret_mix, hist_fret, show_kde=True)
    from fretbursts.phtools.burstsearch import Burst, Bursts
    times = d.ph_times_m[0]  # timestamps array
    bursts = d_fret_mix.mburst[0]
    print("time",(d.time_max-d.time_min))
    print('\nNumber of bursts:', bursts.num_bursts)
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
        if abs(avgdd-avgda)>0.5 or abs(avgaa-avgda)>0.5 or (len(msburst[m_dd])+len(msburst[m_ad]))<20:
            bleachingBurst=bleachingBurst+1
            continue
        # Compute binning of current bursts
        sub_bursts_l.append(Burst(istart=burst.istart, istop=burst.istop,
                                    start=burst.start, stop=burst.stop))
    fBursts=Bursts.from_list(sub_bursts_l)
    print("bleachingBurst",bleachingBurst)
    from fretbursts.phtools.burstsearch import count_ph_in_bursts
    PR=[]
    Tau=[]
    SgDivSr=[]
    pN=[]
    pNa=[] 
    T_burst_duration=[]
    dd_count=[]
    ad_count=[]
    gamma=0.31        
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
        pr=(counts_ad *(1-DexDirAem)-Dch2Ach*counts_dd)/ ((gamma-\
                                Dch2Ach)*counts_dd + (1-DexDirAem)*counts_ad)
        N=counts_dd+counts_ad
        Na=counts_ad
        SgSr=counts_dd/counts_ad
        dd=counts_dd
        ad=counts_ad
        i=i+1
        j=-1        
        fe = SgSr
        if np.isnan(fe) or np.isinf(fe) or fe==0:
            continue
        SgDivSr.append(fe)
        pN.append(N)
        pNa.append(Na)
        bg_time_bins_l.append((msburst[0]+msburst[-1])/2*d.clk_p)
        dd_count.append(dd)
        ad_count.append(ad)


    bg_time_bins=np.asarray(bg_time_bins_l)
    import matplotlib.pyplot as plt
    print("len(bg_time_bins_l)",len(bg_time_bins_l))
    bg_dd_bins = prepData.interpSpl(bg_time,bg_dd, bg_time_bins) *time_bin
    bg_ad_bins = prepData.interpSpl(bg_time,bg_ad, bg_time_bins)  *time_bin
    F_RT=np.asarray(ad_count)-bg_ad_bins*time_bin
    F_G=np.asarray(dd_count)-bg_dd_bins*time_bin
    print("F_RT min ",min(F_RT),"num F_RT<=0",np.where(F_RT<=0)[0].shape[0])
    print("F_G min ",min(F_G),"num F_G<=0",np.where(F_G<=0)[0].shape[0])

    pN=np.asarray(pN)
    bins_idx=np.where((pN>=20 )& (F_RT>0) &(F_G>0))
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
    print("numWindows",numWindows)
    print("bg_dd_bins#",len(bg_dd_bins))
    assert len(PR)==len(Tau)
    import pandas as pd
    # df=pd.DataFrame(dict(x = d_fret_mix.E[0], y = burstsw))
    df=pd.DataFrame(dict(x = PR, y = Tau))
    # import seaborn as sns
    # g = sns.jointplot(x='x',y='y',data=df, kind="kde", size=7, space=0,ylim=(0,1),xlim=(0,1))
    # f, ax = plt.subplots(figsize=(6, 6))
    # sns.kdeplot(df.x, df.y, ax=ax)
    # sns.rugplot(df.x, color="g", ax=ax)
    # sns.rugplot(df.y, vertical=True, ax=ax)
    # plt.hist(nanotimes_dd*25e-12*1e9,120)
    # import matplotlib.pyplot as plt
    import scipy.io as sio
    from array import array

    Fsel=array("d")
    for i in PR:
        if i<0.99 and i >0.01:
            Fsel.append(i)
    sio.savemat(mat_name, {'fret':Fsel})
    #plt.show()

    # import matplotlib.cm as cm
    # import matplotlib.pyplot as plt
    # fig,ax=plt.subplots(2,1)
    # H, xedges, yedges = np.histogram2d(PR,Tau, bins=(27,27))
    # im=ax[1].imshow(H.transpose()[::-1], interpolation='sinc', \
    #                     cmap=cm.jet,extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    # # ax[1].set_title(title)
    # fig.colorbar(im)                       
    T_burst_duration=np.asarray(T_burst_duration,dtype=np.float64)
    
    # plt.hist(T_burst_duration)
    # plt.show()
    pdaIns=aTpda.pdaPy(sub_bursts_l,times,mask_ad,pN,T_burst_duration,SgDivSr,bg_ad_rate,bg_dd_rate,d.clk_p,\
        70,5)#,(10,500)
    pdaIns.set_nstates(2, [0.3,0.7],[0,0],[0.2,0.2])
    # fitE=pdaIns.findP()
    
    pdaIns.fitP=[0.03251902,  0.30312247, -0.07488153 , .43928702, \
        0.86239877,0.001 ]
    
    SgDivSrR,tSgDivSr=pdaIns.aTpdaEK()
    
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt    
    
    # plt.subplot(1,4,2)
    plt.plot(SgDivSrR[1][1:],SgDivSrR[0])
    plt.plot(SgDivSrR[1][1:],tSgDivSr,color="y")
    width = np.diff(SgDivSrR[1])
    print(np.sum(SgDivSrR[0]*width/(sum(SgDivSrR[0])*width)))
    
    plt.title('P(S_G/S_R) fit')
    # plt.subplot(1,4,3)
    # hist, bin_edges = np.histogram(pN,bins=79)
    
    # hista, bin_edgesa = np.histogram(pNa,bins=50)
    # pNhistA=hista/numWindows
    # pNhist=hist/numWindows
    # centerA = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2 
    # width = np.diff(bin_edges)
    # widthA = np.diff(bin_edgesa)
    # center = (bin_edges[:-1] + bin_edges[1:]) / 2    
    # # plt.bar(center,pNhist , align='center', width=width)
    # plt.plot(center,pNhist)
    # photonN=np.arange(max(0,min(pN)),max(pN),1)
    # # P(N)
    # pNcount = prepData.interpSpl(center,pNhist, photonN)#,'linear')    
    # plt.plot (photonN,pNcount,color='y')
    # plt.title('P(N)')
    # plt.subplot(1,4,4)
    # plt.bar(centerA,pNhistA , align='center', width=widthA)
    # # ax.set_xticks(bin_edges)   
    # plt.title('P_r(N)')    
    # print(np.sum(pNhist*width))
    plt.savefig('temp.png')    
    # plt.show()

def usage():  
    print("Usage:%s -i|--hdf inputfilename.hf5 -o|--mat outputfilename.mat" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "o:i:", ["mat=", "hdf="])  
        for o, v in opts: 
            if o in ("-o", "--mat"):
                dbname=v
            if o in ("-i", "--hdf"):
                savefn = v
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    if len(dbname)<1 or len(savefn)<1:
        usage()    
        sys.exit(1)
    burstBin(savefn,dbname)

