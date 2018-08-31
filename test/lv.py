import collections
cburst = collections.namedtuple('burst', ['ntag','e','s'])
ntag=[1,2,3]
fretE=[1,2,3]
fretS=[1,2,3]
burst=cburst(ntag,fretE,fretS)

def ch1(b):
    b.e[0]=-1

ch1(burst)
print(burst.e)
