import numpy as np
from mpi4py import MPI

N = 100000000 # for testing

def parallel_function_caller(x,stopp,comm):
    rank = comm.Get_rank()
    size = comm.Get_size()
    step = N//size # say that N is divisible by size

    stopp[0]=comm.bcast(stopp[0], root=0)
    summ=0
    if stopp[0]==0:
        #your function here in parallel
        x=comm.bcast(x, root=0)
        array= np.arange(x[0]-N/2.+rank*step-42,x[0]-N/2.+(rank+1)*step-42,1.)
        summl=np.sum(np.square(array))
        summ=comm.reduce(summl,op=MPI.SUM, root=0)
        if rank==0:
            print ("value is "+str(summ))
    return summ
