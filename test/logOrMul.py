import numpy as np
n=1000
x=0
import datetime
starttime = datetime.datetime.now()
for i in range(1,n):
    x+=np.log(i)
e = datetime.datetime.now()
print (e-starttime)
y=1
starttime = datetime.datetime.now()
for i in range(1,n):
    y*=i
e = datetime.datetime.now()
print (e-starttime)
print(x,y)
