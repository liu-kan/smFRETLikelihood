import datetime
starttime = datetime.datetime.now()
#long running
s=0
for i in range(60000000):
    s+=i*i
endtime = datetime.datetime.now()
print ((endtime - starttime).seconds)
