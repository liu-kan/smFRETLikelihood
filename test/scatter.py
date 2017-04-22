#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import numpy as np
from mpi4py import MPI

def chunks(l, n):
    """Yield successive nth chunks from l."""
    lenl=len(l)
    stack=[]
    if(lenl%n==0):
        for i in range(0, lenl, lenl/n):
            stack.append(l[i:i + lenl/n])
        for i in range(0, lenl, lenl/n):
            yield stack.pop()
    else:
        for i in range(0, lenl, int(lenl/n)+1):
            stack.append(l[i:i + int(lenl/n)+1])
        for i in range(0, lenl, int(lenl/n)+1):
            yield stack.pop()


comm = MPI.COMM_WORLD
rank =comm.Get_rank()
size = comm.Get_size()
if rank==0:
    a=list(chunks(range(10), size))
    print(a)
else:
    a=list()
al=list()
# Scatter data into my_A arrays
al=comm.scatter( a, root=0)

print("After Scatter:")
for r in range(comm.size):
    if comm.rank == r:
        for i in al:
            print("~~[%d] %s [%d]" % (comm.rank, al,i))
comm.Barrier()
