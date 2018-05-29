# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import sys
lastn=-1
fname=sys.argv[1]
change=0

with open(fname) as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content] 
lenS=len(content)
print("lenS",lenS)
lastT=0
seqS=[]
seqT=[]
for i in range(lenS):
    if lastn!=content[i]:
        change=change+1
        if lastn!=-1:
            seqS.append(lastn)
            seqT.append(lastT)
        lastn=content[i]
        lastT=1
    else:
        lastT=lastT+1
numST=len(seqS)
dicS=dict()
for i in range(numST):
    if seqS[i] not in dicS.keys():
        dicS[seqS[i]]=[]
    dicS[seqS[i]].append(seqT[i])
i=0
print("change",change)
fig, axs = plt.subplots(1, len(dicS),  tight_layout=True) #sharey=True,
for (k,v) in  dicS.items():
    i=int(k)
    axs[i].hist(v, 200,  facecolor='g', alpha=0.75)    #density=True,
plt.show()            