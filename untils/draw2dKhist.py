import re

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
plt.style.use('classic')


import numpy as np

def readlog(filename,rlogNum):
    l=[]
    if rlogNum>0:
        for j in range(0,rlogNum):
            i=rlogNum-j
            print(filename+'.'+str(i))
            with open(filename+'.'+str(i)) as f:
                lines = f.readlines()
            l.extend(lines)
    with open(filename) as f:
        lines = f.readlines()
    l.extend(lines)
    k1=[]
    k2=[]
    chi=[]
    e1=[]
    e2=[]
    v1=[]
    v2=[]
    pChi = re.compile('(.*)chiSqr\:(\d+(\.\d+)?)')
    pK = re.compile('(.*)paramsBca\[(\d+(\.\d+)?), (\d+(\.\d+)?), (-?\d+(\.\d+)?), (-?\d+(\.\d+)?), (\d+(\.\d+)?), (\d+(\.\d+)?)\]')    
    for line in l:
        mk = pK.search(line)
        if mk!=None:
            e1.append(float(mk.group(2)) )
            e2.append(float(mk.group(4)) )
            tk1=float(mk.group(6)) 
            tk2=float(mk.group(8)) 
            v1.append(float(mk.group(10)) )
            v2.append(float(mk.group(12)) )
            k1.append(tk1)
            k2.append(tk2)
            # k1.append( 10**np.asarray(tk1)*1000)
            # k2.append( 10**np.asarray(tk2)*1000)
        mc = pChi.search(line)
        if mc!=None:
            chi.append( float(mc.group(2)) )
    rlen=min(len(k1),len(chi))
    print("Dataset len:",rlen)
    df=pd.DataFrame(dict({"k1":k1[:rlen],"k2":k2[:rlen],"chi":chi[:rlen],"e1":e1[:rlen],\
            "e2":e2[:rlen],"v1":v1[:rlen],"v2":v2[:rlen]}))
    print(df[['k1','k2','chi']].sort_values(by='chi')[:30])
    return k1[:rlen],k2[:rlen],chi[:rlen],df

def filterData(k1,k2,chi,th):
    if th<=0:
        return k1,k2,chi
    rk1=[]
    rk2=[]
    rchi=[]
    for i in range(len(k1)):
        if chi[i]<th:
            rk1.append(k1[i])
            rk2.append(k2[i])
            # rk1.append(k1[i])
            # rk2.append(k2[i])            
            rchi.append(chi[i])
    # import scipy.io as sio
    # sio.savemat('tmp.mat', {'k1':rk1,'k2':rk2,'chi':rchi})    
    return rk1,rk2,rchi

def plot3dScatter(sx,sy,sz):
    
    # fig, ax = plt.subplots(nrows=2, ncols=2)
    # # Plot the model function and the randomly selected sample points
    # ax[0,0].contourf(X, Y, T)
    # ax[0,0].scatter(sx, sy, c='k', alpha=0.2, marker='.')
    # ax[0,0].set_title('Sample points on f(X,Y)')
    # Interpolate using three different methods and plot
    # method='cubic'
    # Ti = griddata((sx, sy), sz, (X, Y), method=method)
    # plt.contourf(X, Y, Ti)
    # plt.scatter(sx, sy, alpha=0.2, marker='.')
    # plt.scatter(sx, sy, c='r', alpha=0.1, marker='.',edgecolors='none')
    minz=min(sz)
    for i in [0,minz+0.05,minz+1,minz+2,minz+3,minz+4,minz+5,minz+10]:
        fx,fy,_=filterData(sx,sy,sz,i)
        fig = plt.figure()
        plt.scatter(fx, fy, c='b', alpha=0.3, marker='o',edgecolors='none')
        plt.title(str(minz)+":"+str(i)+"#"+str(len(fx)))
    # plt.colorbar()
    plt.show()

def pickres(id,df):
    print(df.iloc[id])
    print([df.iloc[id]['e1'],df.iloc[id]['e2'],df.iloc[id]['k1'],df.iloc[id]['k2'],\
        df.iloc[id]['v1'],df.iloc[id]['v2']])

if __name__ == '__main__':
    import sys,getopt
    filename=''
    rlogNum=0
    pid=-1
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "p:l:n:", ["log=","pid=", "rlognum="])  
        for o, v in opts: 
            if o in ("-l", "--log"):
                filename = v
            if o in ("-p", "--pid"):
                pid = int(v.strip())                
            if o in ("-n", "--rlognum"):
                rlogNum = int(v.strip())                
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        sys.exit(1)
    k1,k2,chi,df=readlog(filename,rlogNum)    
    if pid>=0:
        pickres(pid,df)
    else:
        plot3dScatter(k1,k2,chi)
