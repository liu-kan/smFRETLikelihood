
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
time_bin = 1e-3  # 5 ms
time_bin_clk = time_bin / ds.clk_p
stime=3.6
etime=3.9

from fretbursts.phtools.burstsearch import Burst, Bursts
mask_dd = d.get_ph_mask(ph_sel=Ph_sel(Dex='Dem'))   # donor excitation, donor emission
mask_ad = d.get_ph_mask(ph_sel=Ph_sel(Dex='Aem'))   # donor excitation, acceptor emission
mask_aa = d.get_ph_mask(ph_sel=Ph_sel(Aex='Aem'))   # acceptor excitation, acceptor emission
tgstart=stime/ds.clk_p
tgend=etime/ds.clk_p
itgend=0
itgstart=0
i=-1
for ustime in times:
    i=i+1
    if ustime>tgstart:
        itgstart=i
        tgstart=ustime
        print("tgstart",tgstart*ds.clk_p)
        print("itgstart",itgstart)
        break
i=-1
for ustime in times:
    i=i+1
    if ustime>tgend:
        itgend=i
        tgend=ustime
        print("tgend",tgend*ds.clk_p)
        print("itgend",itgend)
        break        
bins = np.arange(tgstart, tgend + time_bin_clk, time_bin_clk)
counts, _ = np.histogram(times[itgstart:itgend+1], bins)

sub_bursts_l = []
sub_start = tgstart
sub_istart = itgstart
for count in counts:
    sub_istop = sub_istart + count - 1
    sub_bursts_l.append(Burst(istart=sub_istart, istop=sub_istop,
                                start=sub_start, 
                                stop=times[sub_istop]))
    sub_istart += count 
    sub_start = times[sub_istart]
sub_bursts = Bursts.from_list(sub_bursts_l)

from fretbursts.phtools.burstsearch import count_ph_in_bursts
PR=[]
TIMEX=[]
gamma=0.31 

ddp=[]
adp=[]
TIMEC=[]

counts_dd = count_ph_in_bursts(sub_bursts, mask_dd)
counts_ad = count_ph_in_bursts(sub_bursts, mask_ad)
counts_aa = count_ph_in_bursts(sub_bursts, mask_aa)
pr=(counts_ad / (counts_dd*gamma + counts_ad)).tolist()
print(len(pr))
print(len(counts_dd))
binsa=np.asarray(bins)
timex=(binsa[:-1]+np.diff(bins)/2)*ds.clk_p
print(ds.clk_p)

timexl=timex.tolist()
i=-1
for tpr in pr:
    i=i+1
    if tpr!=0 and tpr!=1:
        PR.append(tpr)
        TIMEX.append(timexl[i])

timexl=timex.tolist()
i=-1
for cdd in counts_dd:
    i=i+1
    if cdd+counts_ad[i]>0:
        ddp.append(cdd)
        adp.append(counts_ad[i])
        TIMEC.append(timexl[i])


import seaborn as sns
import matplotlib.pyplot as plt

f, ax2 = plt.subplots(1)
# ax2.scatter(timex, pr)
ax2.set_xlim(TIMEX[0],TIMEX[-1])
ax2.set_ylim(0,1)
line,=ax2.plot(TIMEX[0], PR[0],label="RMSD")

# f.subplots_adjust(hspace=0)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("RMSD")
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
gifbin=3
dot,=ax2.plot(TIMEX[0],PR[0],markersize=11, marker='o', c='r',alpha=.9)
lasti=0
def update(i):
    line.set_xdata(TIMEX[0:int(i*gifbin)])
    line.set_ydata(PR[0:int(i*gifbin)])
    dot.set_xdata(TIMEX[int(i*gifbin)])
    dot.set_ydata(PR[int(i*gifbin)])
    # return paths, ax2
from matplotlib.animation import FuncAnimation    
import sys
if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(f, update, frames=np.arange(0, len(TIMEX)/gifbin-1), interval=120)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('line.gif', dpi=220, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()