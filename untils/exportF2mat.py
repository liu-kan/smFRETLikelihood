import pickle
F_S_Z = pickle.load( open( "objs.pickle", "rb" ) )
S=0.87

import scipy.io as sio
from array import array
n=len(F_S_Z[0])
Fsel=array("d")
for i in range(n):
    if F_S_Z[1][i]<=S:
        Fsel.append(F_S_Z[0][i])
sio.savemat('x.mat', {'x':Fsel})