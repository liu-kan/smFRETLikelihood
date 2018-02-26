import pickle,sys,getopt
dbname='fret.mat'
savefn=''
try:  
    opts, args = getopt.getopt(sys.argv[1:], "i:o:b:", ["mat=", "pickle="])  
    for o, v in opts: 
        if o in ("-o", "--mat"):
            dbname=v
        if o in ("-i", "--pickle"):
            savefn = v
except getopt.GetoptError:  
    # print help information and exit:  
    print("getopt error!") 
    sys.exit(1)


F_S_Z = pickle.load( open( savefn, "rb" ) )
S=0.86

import scipy.io as sio
from array import array
n=len(F_S_Z[1])
Fsel=array("d")
for i in range(n):
    #if F_S_Z[1][i]<=S:
    Fsel.append(F_S_Z[1][i])
sio.savemat(dbname, {'fret':Fsel})