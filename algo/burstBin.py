
# coding: utf-8

import numpy as np
from fretbursts import *
# sns = init_notebook()
def burstBin(full_fname,mat_name):
    
    d = loader.photon_hdf5(full_fname)
    # bpl.plot_alternation_hist(d)
    loader.alex_apply_period(d)
    d.calc_bg(fun=bg.exp_fit, time_s=30, tail_min_us='auto', F_bg=1.7)
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

    roi = dict(E1=0, E2=1, S1=0.0, S2=1, rect=False)
    d_fret_mix = dsfuse.select_bursts(select_bursts.ES, **roi)
    # g = alex_jointplot(d_fret_mix)
    # bpl.plot_ES_selection(g.ax_joint, **roi);

    # dplot(d_fret_mix, hist_fret, show_kde=True)
    from fretbursts.phtools.burstsearch import Burst, Bursts
    times = d.ph_times_m[0]  # timestamps array
    bursts = d_fret_mix.mburst[0]
    print('\nNumber of bursts:', bursts.num_bursts)

    time_bin = 0.7e-3  # 0.5 ms
    time_bin_clk = time_bin / ds.clk_p

    sub_bursts_list = []
    fbin=[]
    burstsw=[]
    for burst in bursts:
        # Compute binning of current bursts
        bins = np.arange(burst.start, burst.stop + time_bin_clk, time_bin_clk)
        counts, _ = np.histogram(times[burst.istart:burst.istop+1], bins)
        
        # From `counts` in each bin, find start-stop times and indexes (sub-burst).
        # Note that start and stop are the min and max timestamps in the bin,
        # therefore they are not on the bin edges. Also the burst width is not
        # exactly equal to the bin width.
        sub_bursts_l = []
        sub_start = burst.start
        sub_istart = burst.istart
        burstsw.append((burst.stop-burst.start)*d.clk_p*1e3)
        sub_fbin=[]
        for count in counts:
            # Let's skip bins with 0 photons
            if count == 0:
                continue            
            sub_istop = sub_istart + count - 1
            sub_bursts_l.append(Burst(istart=sub_istart, istop=sub_istop,
                                    start=sub_start, stop=times[sub_istop]))
            # if len(nanotimes_dd[sub_istart:sub_istop])>0:
            sub_fbin.append(dict(ns=nanotimes_dd[sub_istart:sub_istop],
                ms=times[sub_istart:sub_istop],istart=sub_istart,istop=sub_istop))
            sub_istart += count 
            sub_start = times[sub_istart]
        fbin.append(sub_fbin)
        sub_bursts = Bursts.from_list(sub_bursts_l)
        assert sub_bursts.num_bursts > 0
        assert sub_bursts.width.max() < time_bin_clk
        sub_bursts_list.append(sub_bursts)


    mask_dd = d.get_ph_mask(ph_sel=Ph_sel(Dex='Dem'))   # donor excitation, donor emission
    mask_ad = d.get_ph_mask(ph_sel=Ph_sel(Dex='Aem'))   # donor excitation, acceptor emission
    mask_aa = d.get_ph_mask(ph_sel=Ph_sel(Aex='Aem'))   # acceptor excitation, acceptor emission
    from fretbursts.phtools.burstsearch import count_ph_in_bursts
    PR=[]
    Tau=[]
    gamma=0.31        
    beta=1.42
    DexDirAem=0.08
    Dch2Ach=0.07 
    i=-1
    for bursts in sub_bursts_list:
        counts_dd = count_ph_in_bursts(bursts, mask_dd)
        counts_ad = count_ph_in_bursts(bursts, mask_ad)
        counts_aa = count_ph_in_bursts(bursts, mask_aa)
        # print(bursts,times[0:43])
        counts_ddm=0
        counts_adm=0
        for burst in bursts:
            msburst = times[burst.istart:burst.istop + 1]
            m_dd = mask_dd[burst.istart:burst.istop + 1]
            m_ad = mask_ad[burst.istart:burst.istop + 1]        
            counts_ddm+= len(msburst[m_dd])
            counts_adm += len(msburst[m_ad])
        assert counts_adm==sum(counts_ad)
        assert counts_ddm==sum(counts_dd)
        # pr=(counts_ad / (counts_dd*gamma + counts_ad)).tolist()
        pr=((counts_ad *(1-DexDirAem)-Dch2Ach*counts_dd)/ ((gamma-\
                                Dch2Ach)*counts_dd + (1-DexDirAem)*counts_ad)).tolist()    
        i=i+1
        j=-1
        for fe in pr:
            j=j+1
            w=len(fbin[i][j]['ns'])
            if np.isnan(fe):

                continue
            elif w>0:
                # liftt=np.mean(fbin[i][j]['ns'])*2.5e-2
                liftt=((np.mean(fbin[i][j]['ns'])*2.5e-2)-7.9)/4.2
                if True or liftt<=1 and liftt>=0:
                    w=1
                    wfe=[fe]*w
                    wft=[liftt]*w
                    PR.extend( wfe)
                    Tau.extend(wft)

            
    assert len(PR)==len(Tau)
    import pandas as pd
    # df=pd.DataFrame(dict(x = d_fret_mix.E[0], y = burstsw))
    df=pd.DataFrame(dict(x = PR, y = Tau))
    import seaborn as sns
    import matplotlib.pyplot as plt

    # g = sns.jointplot(x='x',y='y',data=df, kind="kde", size=7, space=0,ylim=(0,1),xlim=(0,1))
    # f, ax = plt.subplots(figsize=(6, 6))
    # sns.kdeplot(df.x, df.y, ax=ax)
    # sns.rugplot(df.x, color="g", ax=ax)
    # sns.rugplot(df.y, vertical=True, ax=ax)
    # plt.hist(nanotimes_dd*25e-12*1e9,120)

    import scipy.io as sio
    from array import array

    Fsel=array("d")
    for i in PR:
        if i<0.99 and i >0.01:
            Fsel.append(i)
    sio.savemat(mat_name, {'fret':Fsel})
    plt.show()
    # import matplotlib.cm as cm
    # import matplotlib.pyplot as plt
    # fig,ax=plt.subplots(2,1)
    # H, xedges, yedges = np.histogram2d(PR,Tau, bins=(27,27))
    # im=ax[1].imshow(H.transpose()[::-1], interpolation='sinc', \
    #                     cmap=cm.jet,extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    # # ax[1].set_title(title)
    # fig.colorbar(im)                       

    import matplotlib.pyplot as plt
    plt.hist(burstsw, bins=100) 
    # plt.title(title)
    plt.show()

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
