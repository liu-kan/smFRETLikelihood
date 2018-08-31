# 各类数据结构

## calcBGrate <- BGrate

  bgRate=dict({'time':timeline,'bgrate':buf})
  brcd[ch]=bgRate

time in sec, bgrate in cps

## ch idx
  * 1 DexAem
  * 2 DexDem
  * 3 AexAem
  * 4 AexDem

## BurstSearch

  cburst=dict({'ntag':ntag, 'burstW':burstW,'timetag':timetag,'dtime':dtime,\
                    'chl':chl,'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime})
  burst[ch]=cburst

ntag burst中有多少光子，burstW 时间跨度，
timetag、dtime、chl 各个光子的相应属性。

## binRawData
  binData[ch]=dict({'timetag':timetag,'dtime':dtime,'chl':chl,'binMs':binMs, \
                    'e':fretE,'s':fretS,'z':fretZ,'lifetime':lifetime})