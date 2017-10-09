import pickle
TauD = pickle.load( open( "../data/TauDraw.pickle", "rb" ) )
burstTauD = TauD[0]
import scipy.io as sio
from array import array
n=len(TauD[0])
Fsel=array("d")
for i in range(n):
    Fsel.append(TauD[0][i])
sio.savemat('Tau.mat', {'TauD':Fsel})