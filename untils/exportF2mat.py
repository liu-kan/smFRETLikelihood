import pickle
F_S_Z = pickle.load( open( "/home/liuk/dataZ1/smfretRes/rawRes/rsv/21c_224c_0.5_[2.3509349643187933, 1.1838989238159361, 3.3364004127243341, 3.5348338881347292].pickle", "rb" ) )
S=0.86

import scipy.io as sio
from array import array
n=len(F_S_Z[1])
Fsel=array("d")
for i in range(n):
    #if F_S_Z[1][i]<=S:
    Fsel.append(F_S_Z[1][i])
sio.savemat('fret.mat', {'fret':Fsel})