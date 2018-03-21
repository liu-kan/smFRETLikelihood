
# coding: utf-8

import numpy as np
from fretbursts import *
# sns = init_notebook()
full_fname='data/n25c.h5'
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

roi = dict(E1=0, E2=1, S1=0.0, S2=1, rect=False)
d_fret_mix = dsfuse.select_bursts(select_bursts.ES, **roi)


times = d.ph_times_m[0]  # timestamps array
bursts = d_fret_mix.mburst[0]
print('\nNumber of bursts:', bursts.num_bursts)
time_bin = 5e-3  # 5 ms
time_bin_clk = time_bin / ds.clk_p
stime=0.11
etime=0.17

mask_dd = d.get_ph_mask(ph_sel=Ph_sel(Dex='Dem'))   # donor excitation, donor emission
mask_ad = d.get_ph_mask(ph_sel=Ph_sel(Dex='Aem'))   # donor excitation, acceptor emission
mask_aa = d.get_ph_mask(ph_sel=Ph_sel(Aex='Aem'))   # acceptor excitation, acceptor emission
tgstart=stime/ds.clk_p
tgend=etime/ds.clk_p
itgend=0
itstart=0
for burstistart in burst.istart:
    if times[burstistart]>tgstart:
        itgstart=burstistart
        tgstart=times[burstistart]
        break
for burstistop in burst.istop:
    if burstistop>tgend:
        tgend=burstistop
        break        
    bins = np.arange(tgstart, tgend + time_bin_clk, time_bin_clk)
    counts, _ = np.histogram(times[burst.istart:burst.istop+1], bins)

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
    pr=(counts_ad / (counts_dd*gamma + counts_ad)).tolist()
    
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


plt.show()
