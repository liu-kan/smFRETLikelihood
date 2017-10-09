import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from pandas import DataFrame

# read .RData file as a pandas dataframe
def load_rdata_file(filename):
    r_data = robjects.r['get'](robjects.r['load'](filename))
    df = pandas2ri.ri2py(r_data)
    return df

# write pandas dataframe to an .RData file
def save_rdata_file(df, filename,name="my_df"):
    r_data = pandas2ri.py2ri(df)
    robjects.r.assign(name, r_data)
    robjects.r("save({}, file='{}')".format(name,filename))

import pickle
import numpy as np
from array import array
TauD = pickle.load( open( "../data/TauDraw.pickle", "rb" ) )
burstTauD = TauD[0]
ad=np.dot(burstTauD,1e-9)
ad=burstTauD
thebin=200
[hist,bin]=np.histogram(ad, bins=thebin)
x=array("d")
y=array("d")
sidx=np.argmax(hist) #只拟合最大值之后的数据    
lenhist=len(hist)
for vh in range(sidx,lenhist):
    if hist[vh]>0 and bin[vh]>bin[sidx]+0.0:        
        y.append(hist[vh])
        x.append(bin[vh]-bin[sidx]-0.0)     
        #print(bin[vh])           
xx=np.array(x)
yy=np.array(y)
TauDdf=DataFrame({'x':xx, 'y':yy})
save_rdata_file(TauDdf,'out.RData','TauDdf')