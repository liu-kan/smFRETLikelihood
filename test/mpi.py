from scipy.optimize import fmin
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

N = 10 # for testing
step = N//size # say that N is divisible by size
class p():
    def parallel_function_caller(self,x,stopp):
        stopp[0]=comm.bcast(stopp[0], root=0)
        summ=0
        if stopp[0]==0:
            #your function here in parallel

            x=comm.bcast(x, root=0)
            array= np.arange(x[0]-N/2.+rank*step-42,x[0]-N/2.+(rank+1)*step-42,1.)
            summl=np.sum(np.square(array))+rank
            summ=comm.reduce(summl,op=MPI.SUM, root=0)
            print(rank," calcing ",summl)
            if rank==0:
                print ("value is "+str(summ))
                return summ
            else:
                pass

pp=p()
if rank == 0 :
   stop=[0]
   x = np.zeros(1)
   x[0]=20
   #xs = minimize(parallel_function_caller, x, args=(stop))
   xs = fmin(pp.parallel_function_caller,x0= x, args=(stop,))
   print ("the argmin is "+str(xs))
   stop=[1]
   pp.parallel_function_caller(x,stop)

else :
   stop=[0]
   x=np.zeros(1)
   while stop[0]==0:

      pp.parallel_function_caller(x,stop)
   print(rank," done")