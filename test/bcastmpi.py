from mpi4py import MPI
import numpy as np
import collections

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    data =dict()

    cburst = collections.namedtuple('datas', ['ntag', 'burstW','timetag','dtime','chl','e','s'])
    data["ch"]=dict({"dsf":[32,54,768,789,9834,45],"lala":2})
    data["df"]=np.ones(1,dtype=np.int)
else:
    data = dict()
    print("data",data)


data=comm.bcast(data, root=0)
print("rank,",rank,data)
r=np.zeros(1,dtype=np.int)
comm.Reduce(data['df'], r, op=MPI.SUM, root=0)
print("Reduce",r[0])
