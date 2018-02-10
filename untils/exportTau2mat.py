import pickle
TauD = pickle.load( open( "/home/liuk/dataZ1/smfretRes/rawRes/rsv/21c_224c_0.5_[2.3509349643187933, 1.1838989238159361, 3.3364004127243341, 3.5348338881347292].pickle", "rb" ) )

import scipy.io as sio
from array import array
n=len(TauD[0])
Fsel=array("d")
for i in range(n):
    Fsel.append(TauD[0][i])
sio.savemat('Tau.mat', {'TauD':Fsel})