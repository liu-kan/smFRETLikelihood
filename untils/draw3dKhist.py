import re

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
plt.style.use('classic')


import numpy as np

def readlog(filename,rlogNum,sNum):
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
    e=[]
    k=[]
    v=[]
    for ids in range(sNum):
        e.append([])
        for i in range(sNum-1):
            k.append([])
        v.append([])
    # print(len(k))
    chi=[]

    pChi = re.compile('(.*)chiSqr\:(\d+(\.\d+)?)')
    pKstr='(.*)paramsBca\['+'(-?\d+(\.\d+)?), '*(sNum*sNum+sNum-1)+'(-?\d+(\.\d+)?)\]'
    # pK = re.compile('(.*)paramsBca\[(\d+(\.\d+)?), (\d+(\.\d+)?), (-?\d+(\.\d+)?), (-?\d+(\.\d+)?), (\d+(\.\d+)?), (\d+(\.\d+)?)\]')
    pK = re.compile(pKstr)
    for line in l:
        mk = pK.search(line)
        if mk!=None:
            for ids in range(sNum):
                e[ids].append(float(mk.group(2*(ids+1)) ))
                for i in range(sNum-1):
                    # print(ids*(sNum-1)+i)
                    k[ids*(sNum-1)+i].append( float(mk.group(2*sNum+2*(ids*(sNum-1)+i+1))) )
                v[ids].append( float(mk.group(2*sNum*sNum+2*(ids+1))) )
        mc = pChi.search(line)
        if mc!=None:
            chi.append( float(mc.group(2)) )
    rlen=min(len(k[0]),len(chi))
    print("Dataset len:",rlen)
    df=pd.DataFrame()
    for ids in range(len(v)):
        colname="v"+str(ids+1)
        df.loc[:,colname] = v[ids][:rlen]    
    for idk in range(len(e)):
        colname="e"+str(idk+1)
        df.loc[:,colname] = e[idk][:rlen]       
    for idk in range(len(k)):
        colname="k"+str(idk+1)
        df.loc[:,colname] = k[idk][:rlen]
    df.loc[:,"chi"] = chi[:rlen]
    print(df.iloc[:,sNum:].sort_values(by='chi')[:30])
    return chi[:rlen],df

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

def pickres(id,df,sNum):
    print(df.iloc[id])
    thelist=df.iloc[id][sNum:sNum*2].tolist()
    thelist.extend(df.iloc[id][sNum*2:-1])
    thelist.extend(df.iloc[id][:sNum])
    print(thelist)
    return thelist

if __name__ == '__main__':
    import sys,getopt
    filename=''
    rlogNum=0
    pid=-1
    sn=2
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "p:l:n:s:", ["log=","pid=", "rlognum=", "state="])  
        for o, v in opts: 
            if o in ("-l", "--log"):
                filename = v
            if o in ("-p", "--pid"):
                pid = int(v.strip())                
            if o in ("-s", "--state"):
                sn = int(v.strip())                                
            if o in ("-n", "--rlognum"):
                rlogNum = int(v.strip())                
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        sys.exit(1)
    chi,df=readlog(filename,rlogNum,sn)    
    if pid>=0:
        pickres(pid,df,sn)
    else:
        plot3dScatter(df['k1'],df['k2'],chi)
