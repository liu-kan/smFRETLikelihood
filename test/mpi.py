from scipy.optimize import fmin
from mpi4py import MPI
import numpy as np


class p():
    def __init__(self,comm,N,step):
        self.comm=comm
        self.N=N
        self.step=step
    def parallel_function_caller(self,x,stopp):
        size = self.comm.Get_size()
        rank = self.comm.Get_rank()
        N=self.N
        step=self.step
        stopp[0]=self.comm.bcast(stopp[0], root=0)
        summ=0
        if stopp[0]==0:
            #your function here in parallel

            x=self.comm.bcast(x, root=0)            
            array= np.arange(x[0]-N/2.+rank*step-42,x[0]-N/2.+(rank+1)*step-42,1.)
            summl=np.sum(np.square(array))
            summ=self.comm.reduce(summl,op=MPI.SUM, root=0)
            print(rank," calcing ",summl)
            if rank==0:
                print ("value is "+str(summ))
                return summ
            else:
                pass
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()


    N = 10 # for testing
    step = N//size # say that N is divisible by size

    pp=p(comm,N,step)
    rank = comm.Get_rank()
    if rank == 0 :
        stop=[0]
        x = np.zeros(1)
        x[0]=20
        #xs = minimize(parallel_function_caller, x, args=(stop))
        xs = fmin(pp.parallel_function_caller,x0= x, args=(stop,))
        print ("the argmin is "+str(xs))
        stop2=[1]
        pp.parallel_function_caller(x,stop2)

    else :
        stop=[0]
        x=np.zeros(1)
        idx=0
        while stop[0]==0:
            idx+=1
            print(idx,"===while rank",rank)
            pp.parallel_function_caller(x,stop)
        print(rank," done")
 