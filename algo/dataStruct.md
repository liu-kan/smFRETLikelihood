# 各类数据结构

## calcBGrate <- BGrate

  bgRate = collections.namedtuple('bgRate', ['time', 'bgrate'])

time in sec, bgrate in cps

## ch idx
  * 1 DexAem
  * 2 DexDem
  * 3 AexAem
  * 4 AexDem

## BurstSearch

  cburst = collections.namedtuple('burst', ['ntag', 'burstW','timetag','dtime','chl','e','s'])
  burst[ch]=cburst(ntag,burstW,timetag,dtime,chl,fretE,fretS)

ntag burst中有多少光子，burstW 时间跨度，
timetag、dtime、chl 各个光子的相应属性。
