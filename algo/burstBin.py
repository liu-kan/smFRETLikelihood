
# coding: utf-8

import numpy as np
from fretbursts import *
# sns = init_notebook()
full_fname='LS3.h5'
d = loader.photon_hdf5(full_fname)
# bpl.plot_alternation_hist(d)
loader.alex_apply_period(d)
d.calc_bg(fun=bg.exp_fit, time_s=10, tail_min_us='auto', F_bg=2)
# dplot(d, timetrace_bg)
# dplot(d, timetrace)
# xlim(10, 20)
# ylim(-50, 50)
d.burst_search()

ds = d.select_bursts(select_bursts.size, th1=30)
dsfuse = ds.fuse_bursts(ms=0)
dsfuse.leakage = 0.07
# alex_jointplot(dsfuse)
# dplot(dsfuse, hist_fret, show_kde=True)

nanotimes = d.nanotimes[0]
nanotimes_d = nanotimes[d.get_D_em()]
nanotimes_a = nanotimes[d.get_A_em()]

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

time_bin = 0.5e-3  # 0.5 ms
time_bin_clk = time_bin / ds.clk_p

sub_bursts_list = []
fbin=[]
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
    for count in counts:
        # Let's skip bins with 0 photons
        if count == 0:
            continue            
        sub_istop = sub_istart + count - 1
        sub_bursts_l.append(Burst(istart=sub_istart, istop=sub_istop,
                                  start=sub_start, stop=times[sub_istop]))
        fbin.append(dict(ns=nanotimes[sub_istart:sub_istop],
            ms=times[sub_istart:sub_istop],istart=sub_istart,istop=sub_istop))
        sub_istart += count 
        sub_start = times[sub_istart]
    
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
fNan=np.zeros(len(fbin))
i=-1
for bursts in sub_bursts_list:
    counts_dd = count_ph_in_bursts(bursts, mask_dd)
    counts_ad = count_ph_in_bursts(bursts, mask_ad)
    counts_aa = count_ph_in_bursts(bursts, mask_aa)
    pr=(counts_ad / (counts_dd + counts_ad)).tolist()
    i=i+1
    for fe in pr:
        if np.isnan(fe):
            fNan[i]=1
            continue
        else:
            PR.append( fe)
i=-1
for bin in fbin:
    i=i+1
    if fNan[i]==0:
        Tau.append(((np.mean(bin['ns'])*2.5e-2)-2.5)/7.1)
assert len(PR)==len(Tau)
import pandas as pd
df=pd.DataFrame(dict(PR = PR, Tau = Tau))
# import seaborn as sns
g = sns.JointGrid(x='PR', y='Tau',data=df)
g.plot_marginals(sns.distplot)
g.plot_joint(plt.hexbin,gridsize=20)
#fretbursts.burstlib_ext.calc_mean_lifetime