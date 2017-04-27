#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import numpy as np
from mpi4py import MPI

def chunks(l, n):
    """Yield successive nth chunks from l."""
    lenl=len(l)
    stack=[]
    stackList=[]
    for i in range(0, lenl, lenl//n):
        stack.append(l[i:i + lenl//n])
    for i in range(n):
        stackList.append([stack[i]])
    if len(stack)>n:
        for i in range(n,len(stack)):
            stackList[i-n].append(stack[i])
    for i in range(0, n):
        yield stackList.pop()


comm = MPI.COMM_WORLD
rank =comm.Get_rank()
size = comm.Get_size()
if rank==0:
    a=list(chunks(range(735), size))
    print(a)
else:
    a=list()
    pass
al=list()
# Scatter data into my_A arrays
al=comm.scatter( a, root=0)

print("After Scatter:")
for r in range(comm.size):
    if comm.rank == r:
        for i in al:
            print("~~",comm.rank, al,i)
comm.Barrier()
